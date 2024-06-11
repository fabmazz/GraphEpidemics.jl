import DataStructures: PriorityQueue, enqueue!, peek, delete!, dequeue!
import StatsBase


struct GillespieTrans2{F<:AbstractFloat}
    idx::Int
    infector::Int
    tr_idx::StI
    rate::F
end

GillTrans = GillespieTrans2
rate(g::GillTrans) = g.rate

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

function check_add_infect_trans(coda::PriorityQueue,i::Integer, j::Integer,te::AbstractFloat, infector_arr::Vector{<:Real})
    
    ne = TransEvent(j,2)

    if haskey(coda, ne)
        tu = coda[ne]
        if tu > te
            coda[ne]=te
            infector_arr[j] = i
        end
    else
        enqueue!(coda, ne, te)
        infector_arr[j]=i
    end
    ne
end

function draw_infection_delay(model, i::Integer,j::Integer,rng::AbstractRNG, infect_IorS)
    rate_infect = calc_prob_infection(model, i, j, infect_IorS)
    draw_delay_exp(rate_infect, rng, )
end

function sim_sir_gillespie(g::AbstractGraph, model::AbstractSIRModel, simdata::SIRSimData, rng::AbstractRNG, 
    patient_zeros::Vector{<:Integer}; infect_IorS::Symbol=:I, max_revive::Integer=0, debug::Bool=false)
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
    npz= length(patient_zeros)
    allcounts = [[N-npz, npz,0]]
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
        push!(allcounts, StatsBase.counts(states,3))
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

function gillespie_sir_direct(g::AbstractGraph, model::SIRModel, simdata::SIRSimData, rng::AbstractRNG, 
    patient_zeros::Vector{I};) where I<: Integer
    N = nv(g)

    infect_t = simdata.infect_time
    infect_i = simdata.infect_node
    delays = simdata.rec_delays

    wtree = BinaryTree{Float64,GillTrans{Float64}}(3)

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
        add_event!(wtree, GillespieTrans2(i,-1,StI(3),model.gamma))
        for j in neighbors(g, i)
            if states[j]==1
                e = GillespieTrans2(j,i,StI(2), model.beta)
                add_event!(wtree,e )
                if !haskey(infect_events, i)
                    infect_events[i] = [e]
                else
                    push!(infect_events[i],e)
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
        event = remove_event_idx!(wtree, idx_ev)

        i=event.idx
        s = event.tr_idx
        ##set state
        states[i] =s
        t+=τ
        if s == 2

            ##remove other infection events for this node
            remidx = Int[]
            for (k,e) in wtree.eventsByPos
                if (e.idx == i)
                    push!(remidx,k)
                end
            end
            for mid in remidx
                remove_event_idx!(wtree, mid)
            end

            infect_t[i] = t
            infect_i[i] = event.infector
            for j in neighbors(g,i)
                if states[j]==1
                    e = GillespieTrans2(j,i,StI(2), model.beta)
                    add_event!(wtree,e )
                    if !haskey(infect_events, i)
                        infect_events[i] = [e]
                    else
                        push!(infect_events[i],e)
                    end
                end
            end
            ## add recovery transition
            add_event!(wtree, GillespieTrans2(i, -1, StI(3), model.gamma))
            

            #@assert i in keys(infect_events)

        elseif s==3
            ## do nothing
            #println("$i has recovered: $ev")
            delays[i] = t - infect_t[i]
            if i in keys(infect_events)
                for ev in infect_events[i]
                    if has_event(wtree,ev)
                        remove_event!(wtree, ev)
                    end
                end
            end
        end
        push!(allcounts, StatsBase.counts(states,3))
        push!(times, t)

        tot_rate = total_rate(wtree)
    end

    times, allcounts
end