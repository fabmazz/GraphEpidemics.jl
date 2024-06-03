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

draw_delays_gillespie(m::SIRModel, p, rng::AbstractRNG, nodes) = draw_delay_exp(p, rng, nodes)
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

function sim_sir_gillespie(g::AbstractGraph, model::SIRModel, simdata::SIRSimData, rng::AbstractRNG, 
    patient_zeros::Vector{I};) where I<: Integer
    N = nv(g)

    infect_t = simdata.infect_time
    infect_i = simdata.infect_node
    delays = simdata.rec_delays

    pq = PriorityQueue{TransEvent,Float64}()

    t=0
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
                te = draw_delay_exp(model.beta, rng, )+t
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
            for j in neighbors(g,i)
                if states[j]==1
                    te = draw_delay_exp(model.beta, rng, )+t
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
            enqueue!(pq, TransEvent(i,3),trec )
        elseif s==3
            ## do nothing
            #println("$i has recovered: $ev")
        end
        push!(allcounts, StatsBase.counts(states,3))
        push!(times, t)
    end

    times, allcounts
end

function run_sir_gillespie(g::AbstractGraph, model::SIRModel, rng::AbstractRNG, 
    patient_zeros::Vector{<:Integer}; dtype::DataType=Float64)
    ## draw recovery delays
    N= nv(g)
    nodes = collect(Int32,1:N)
    delays = draw_delays_gillespie(model, model.gamma, rng, nodes)

    data = SIRSimData(delays,N, dtype)
    times, allcounts = sim_sir_gillespie(g, model, data, rng, patient_zeros)

    data, times, allcounts
end

function gillespie_sir_first(g::AbstractGraph, model::SIRModel, simdata::SIRSimData, rng::AbstractRNG, 
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
        idx_ev = find_leaf_idx_random_draw(wtree, r)
        event = remove_event_idx!(wtree, idx_ev)

        i=event.idx
        s = event.tr_idx
        ##set state
        states[i] =s
        t+=τ
        if s == 2
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

        elseif s==3
            ## do nothing
            #println("$i has recovered: $ev")
            delays[i] = t - infect_t[i]
            for ev in infect_events[i]
                if has_event(wtree,ev)
                    remove_event!(wtree, ev)
                end
            end
        end
        push!(allcounts, StatsBase.counts(states,3))
        push!(times, t)
    end

    times, allcounts
end