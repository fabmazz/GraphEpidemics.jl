import StatsBase

const NO_INFECTOR = -5

struct GillespieTrans2{F<:AbstractFloat}
    idx::Int
    infector::Int
    tr_idx::StI
    rate::F
end

struct TransKey
    idx::Int
    infector::Int
    tr_idx::StI
end
struct GillespieTrans3{F<:AbstractFloat}
    trkey::TransKey
    rate::F
end
GillespieTrans3(idx, infector, tr_idx, rate) = GillespieTrans3(TransKey(idx, infector, tr_idx), rate)

GillTrans = GillespieTrans3
rate(g::GillTrans) = g.rate
key_transition(g::GillTrans) = g.trkey

draw_delays_gillespie(m::AbstractSIRModel, p, rng::AbstractRNG, nodes) = draw_delay_exp(p, rng, nodes)
function draw_truncated_exponential(p::AbstractFloat, rng::AbstractRNG, max_t::Real)
    v = rand(rng, Exponential(1/p))
    l = floor(v/max_t)
    return v - l*max_t
end


function sum_rates(trans::Vector{GillTrans})
    s = 0.0
    for t in trans
        s+= t.rate
    end
    s
end

struct TransEvent
    index::Int
    newstate::Int8
end

function check_add_infect_trans(coda::PriorityQueue,i_from::Integer, j_to::Integer,te::AbstractFloat, infector_arr::Vector{<:Real})
    
    event = TransEvent(j_to,2) # do not track infector in event

    if haskey(coda, event)
        tu = coda[event] # infection for the same j 
        if tu > te
            ## the new happens before the old one
            coda[event]=te 
            infector_arr[j_to] = i_from
        end
    else
        enqueue!(coda, event, te)
        infector_arr[j_to]=i_from
    end
    event
end

function draw_infection_delay(model, i::Integer,j::Integer,rng::AbstractRNG, infect_IorS)
    rate_infect = calc_prob_infection(model, i, j, infect_IorS)
    draw_delay_exp(rate_infect, rng, )
end

function sim_sir_gillespie(g::AbstractGraph, model::AbstractSIRModel, simdata::SIRSimData, rng::AbstractRNG, 
    patient_zeros::Vector{<:Integer}; infect_IorS::Symbol=:I, max_revive::Integer=0, debug::Bool=false, counts_func::Function =count_SIR_states)
    N = nv(g)

    infect_t = simdata.infect_time
    infect_i = simdata.infect_node
    delays = simdata.rec_delays

    pq = PriorityQueue{TransEvent,Float64}()

    t=0
    nrevive = Dict{Integer, Integer}()
    states::Vector{StI} = fill(1, N)
    states[patient_zeros] .=2
    #allevents = Tuple{Float64,TransEvent}[]
    for i in patient_zeros
        #states[i] = 2
        infect_t[i] = t ## time of infection
        infect_i[i] = -10 ## node of infection

        trec = t+delays[i]
        enqueue!(pq,TransEvent(i,3),trec)
        for j in neighbors(g, i)
            if states[j]==1
                te = draw_infection_delay(model, i,j, rng, infect_IorS)+t
                if te<=trec
                    ne = check_add_infect_trans(pq, i, j, te, infect_i)
                    #push!(allevents,(te,ne))
                end
            end
        end

    end
    ##setup complete
    #allcounts = [[N-npz, npz,0]]
    allcounts = [counts_func(states)]
    times =[0.]
    while !isempty(pq)
        ev = peek(pq)
        dequeue!(pq)
        i =ev.first.index
        s = ev.first.newstate
        t = ev.second
        ##set state
        states[i] =s
        if s == 2
            infect_t[i] = t
            trec = t+delays[i]
            enqueue!(pq, TransEvent(i,3),trec )
            # for directed graphs, neighbors(g,i) returns the outgoing neighbors
            for j in neighbors(g,i)
                if states[j]==1
                    te = draw_infection_delay(model,i,j, rng,infect_IorS)+t
                    if te <= trec
                        ne=check_add_infect_trans(pq, i, j, te, infect_i)
                        #push!(allevents,(te, ne))
                    end
                end
            end
            ## equeue event of I->R
            #println("Enqueue $i recovering at $trec")
            if i in patient_zeros
                println("Coda: $pq")
            end
 
        elseif s==3
            ## do nothing
            #println("$i has recovered: $ev")
        end
        revivepossible = max_revive>0 && !((i in keys(nrevive)) && (nrevive[i] >= max_revive))
        if isempty(pq) && (s ==3 ) && revivepossible
            ## the simulation finished with a recovery event
            ## the node has not passed the maximum number of "revive" attempt
            newdelay = draw_delay_exp(model.gamma,rng, i)
            delays[i] += newdelay
            newtrec = t + newdelay
            enqueue!(pq, ev.first, newtrec)
            ## redraw infection delays
            for j in neighbors(g, i)
                if states[j]==1
                    te = draw_infection_delay(model, i,j, rng, infect_IorS)+t
                    if te<=newtrec
                        ne = check_add_infect_trans(pq, i, j, te, infect_i)
                        #push!(allevents,(te,ne))
                    end
                end
            end
            ## mark revival
            if i in keys(nrevive)
                nrevive[i]+=1
            else
                nrevive[i] = 1
            end
            states[i] = 2
            if debug
                println("t $t: restore $i to I")
            end

        end
        ##IMPROVE: if we do not "revive", add to the count
        #push!(allcounts, StatsBase.counts(states,3))
        push!(allcounts, counts_func(states))
        push!(times, t)
        
    end

    times, allcounts
end

function run_sir_gillespie(g::AbstractGraph, model::AbstractSIRModel, rng::AbstractRNG, 
    patient_zeros::Vector{<:Integer}; dtype::DataType=Float64, kwargs...)
    ## draw recovery delays
    N= nv(g)
    nodes = collect(Int32,1:N)
    delays = draw_delays_gillespie(model, model.gamma, rng, nodes)

    data = SIRSimData(delays,N, dtype)
    times, allcounts = sim_sir_gillespie(g, model, data, rng, patient_zeros; kwargs...)

    data, times, allcounts
end


function gillespie_sir_direct(g::AbstractGraph, model::AbstractSIRModel, simdata::SIRSimData, rng::AbstractRNG, 
    patient_zeros::Vector{I}; ignore_infector::Bool=false, infect_IorS::Symbol=:I) where I<: Integer
    N = nv(g)

    infect_t = simdata.infect_time
    infect_i = simdata.infect_node
    delays = simdata.rec_delays
    #Ftype = get_type(calc_prob_infection(model,1,1, infect_IorS))
    wtree = BinaryTree{Float64,GillTrans{<:AbstractFloat}, TransKey}(3)

    t=0.0
    states::Vector{StI} = fill(1, N)
    states[patient_zeros] .=2

    infect_events=Dict{Int,Vector{GillTrans{Float64}}}()
    #allevents = Tuple{Float64,TransEvent}[]
    for i in patient_zeros
        #states[i] = 2
        infect_t[i] = t ## time of infection
        infect_i[i] = -10 ## node of infection

        #trec = t+delays[i]
        #enqueue!(pq,TransEvent(i,3),trec)
        add_event!(wtree, GillTrans(i,-1,StI(3),model.gamma))
        for j in neighbors(g, i)
            if states[j]==1
                iif = ignore_infector ? NO_INFECTOR : i
                e = GillTrans(j,iif,StI(2), calc_prob_infection(model, i,j, infect_IorS))
                add_event!(wtree,e )
                if !ignore_infector
                    if !haskey(infect_events, i)
                        infect_events[i] = [e]
                    else
                        push!(infect_events[i],e)
                    end
                end
            end
        end

    end
    ##setup complete
    npz= length(patient_zeros)
    allcounts = [[N-npz, npz,0]]
    times =[0.]
    tot_rate = total_rate(wtree)
    while (tot_rate > 0)
        
        tot_rate = total_rate(wtree)
        τ = -log( rand(rng)) / tot_rate
        if τ == 0 
            println("WARNING: τ is 0. TOTAL RATE: $tot_rate")
        elseif τ > 1e12
            println("WARNING: τ is incredibily large (order: $(round(log10(τ))) ) TOTAL RATE: $tot_rate. Occupied tree leaves: $(wtree.noccupied)")
        end
        r = rand(rng) * tot_rate
        c=0
        while (r <= 0) & (c<1000)
            r = rand(rng) * tot_rate
            c+=1
        end
        if c == 1000
            println("ERROR: the random value is 0.0, tot_rate: $tot_rate")
        end
        idx_ev = find_leaf_idx_random_draw(wtree, r)
        rate_infection = get_weight_tree_idx(wtree, idx_ev)
        event = remove_event_idx!(wtree, idx_ev)
        ekey = key_transition(event)
        i = ekey.idx
        s = ekey.tr_idx
        ##set state
        states[i] =s
        ## increase time 
        t+=τ
        if s == 2

            ## add recovery transition
            add_event!(wtree, GillTrans(i, -1, StI(3), gamma(model)))
            
            if !ignore_infector
                ##remove other infection events for this node
                remidx = Int[]
                for (k,e) in wtree.eventsByPos
                    if (e.trkey.idx == i)
                        push!(remidx,k)
                    end
                end
                for mid in remidx
                    remove_event_idx!(wtree, mid)
                end
            end

            

            infect_t[i] = t
            infect_i[i] = event.trkey.infector
            infector_r = ignore_infector ? rand(rng)*rate_infection : -2.0
            #=inf_search_sum = 0.0
            infector_found = ignore_infector ? false : true
            =#
            #=if ignore_infector
                println("$i inf, rate $rate_infection , x: $infector_r")
            end=#
            
            for j in neighbors(g,i)
                if states[j]==1
                    iif = ignore_infector ? NO_INFECTOR : i
                    e = GillTrans(j,iif,StI(2), calc_prob_infection(model, i,j, infect_IorS))
                    if ignore_infector
                        if contains_event(wtree, e)
                            increase_rate_event!(wtree, e, e.rate)
                        else
                            add_event!(wtree,e)
                        end
                        ## nothing else to do
                    else
                        add_event!(wtree,e )
                        if !haskey(infect_events, i)
                            infect_events[i] = [e]
                        else
                            push!(infect_events[i],e)
                        end
                    end
                #=elseif (!infector_found & (states[j]==2))
                    ## might have been infected by it
                    beta_ji = calc_prob_infection(model, j,i, infect_IorS)
                    inf_search_sum += beta_ji
                    #println("j $j rate $beta_ji, cumsum: $inf_search_sum")
                    if infector_r <= inf_search_sum
                        infector_found = true
                        infect_i[i] = j
                    end
                    =#
                end
            end


            #@assert i in keys(infect_events)

        elseif s==3
            #println("$i has recovered: $ev")
            delays[i] = t - infect_t[i]
            if ignore_infector
                for j in neighbors(g, i)
                    if states[j] == 1
                        e = GillTrans(j,NO_INFECTOR,StI(2), calc_prob_infection(model, i,j, infect_IorS))
                        if contains_event(wtree, e)
                            increase_rate_event!(wtree, e, -1*e.rate)
                        end
                    end
                end
            else
                if i in keys(infect_events)
                    for ev in infect_events[i]
                        if has_event(wtree,ev)
                            remove_event!(wtree, ev)
                        end
                    end
                end
            end
        end
        push!(allcounts, StatsBase.counts(states,3))
        push!(times, t)
        #reupdate total rate
        tot_rate = total_rate(wtree)
    end

    times, allcounts
end