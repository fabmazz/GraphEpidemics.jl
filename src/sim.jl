#= 
 Copyright (c) 2024 Fabio Mazza
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at https://mozilla.org/MPL/2.0/.
=#

using OffsetArrays
using Statistics

beta_R0(R0::Real,g::AbstractGraph,gamma::Real) = R0*gamma/mean(degree(g))
beta_R0(R0::Real,mean_deg::Real,gamma::Real) = R0*gamma/mean_deg

struct SIRSimData{I<:Real, F<:AbstractFloat}
    rec_delays::Vector{I}
    infect_time::Vector{F}
    infect_node::Vector{F}
end
SIRSimData{I,F}(N::Integer) where {I,F}= SIRSimData(zeros(I,N), fill(F(NaN),N),fill(F(NaN),N))

"""
`explode(d::SIRSimData)`

Returns the three fields of a `SIRSimData` instance as a tuple.
"""
function explode(d::SIRSimData)
    (d.rec_delays, d.infect_time, d.infect_node)
end
"""
`explode(d::SIRSimData, t::DataType)`

Returns the fields of a `SIRSimData` object converted to a specified data type.
"""
function explode(d::SIRSimData, t::DataType)
    (convert.(t,d.rec_delays), convert.(t,d.infect_time), convert.(t,d.infect_node))
end

function SIRSimData(li::Union{Vector, Tuple})
    SIRSimData(li[1], li[2], li[3])
end
"""
`SIRSimData(rec_delays, N::Integer, typeDates::DataType)`

Constructs a `SIRSimData` object given recovery delays, number of nodes, using the specified time data type.
"""
SIRSimData(rec_delays, N::Integer, typeDates::DataType) = SIRSimData(rec_delays, 
            fill(convert(typeDates,NaN),N), 
            fill(convert(typeDates,NaN),N))

            """
`check_infection(rng::AbstractRNG, p::AbstractFloat, i::Integer, j::Integer, infect_prob_I::Symbol)`

Determines whether infection occurs between two nodes using scalar infection probability. `i`, `j` and `infect_prob_I` are ignored.
"""
function check_infection(rng::AbstractRNG, p::AbstractFloat, i::Integer, j::Integer, infect_prob_I::Symbol)
    rand(rng) < p
end

"""
`check_infection(rng::AbstractRNG, p::Vector{<:AbstractFloat}, i_I::Integer, j_S::Integer, infect_prob_IS::Symbol)`

Determines whether infection occurs using per-node probabilities, optionally depending on infecting or susceptible node.
"""
function check_infection(rng::AbstractRNG, p::Vector{<:AbstractFloat}, i_I::Integer, j_S::Integer, infect_prob_IS::Symbol)
    if infect_prob_IS == :I
        return rand(rng) < p[i_I]
    elseif infect_prob_IS == :S
        return rand(rng) < p[j_S]
    else
        return rand(rng) < sqrt(p[i_I]*p[j_S])
    end
end


"""
`sim_sir_fast(g::AbstractGraph, model::AbstractSIRModel, T::Integer, simdata::SIRSimData, rng::AbstractRNG, patient_zeros::Vector{I}; beta_IorS=:I, counter=BaseSIRStatesCounter()) where I<:Integer`

Runs a discrete-time SIR simulation over `T` time steps, with the data structures already set up.

### Optional arguments
- `beta_IorS=:I`: Whether to base infection probabilities on infectors (`:I`), susceptibles (`:S`), or their mean. Defaults to `:I`.
- `counter=BaseSIRStatesCounter()`: The state-counting struct. When using a custom struct, make sure to implement `count_states(::BaseSIRStatesCounter, states::Vector)`
- `dynstateChanger`: object used for changing the state dynamically
"""
function sim_sir_fast(G::AbstractGraph, model::AbstractSIRModel, T::Integer, simdata::SIRSimData, rng::AbstractRNG, 
    patient_zeros::Vector{I}; 
    beta_IorS::Symbol = :I, 
    counter::AbstractStatesCounter=BaseSIRStatesCounter(), 
    dynstateChanger::Union{Nothing,AbstractStateChanger}=nothing) where I<: Integer
    N = nv(G)

    infect_t = simdata.infect_time
    infect_i = simdata.infect_node
    delays = simdata.rec_delays

    if !(beta_IorS in (:I,:S,:mean))
        throw(ArgumentError("argument 'beta_IorS' must be one of :I, :S, or :mean"))
    end

    states::Vector{StI} = fill(1, N)
    for i in patient_zeros
        states[i] = 2
        infect_t[i] = -1 ## time of infection
        infect_i[i] = -10 ## node of infection
    end
    
    
    t=0
    useT = (T >= 0)
    if useT
        trace_states = Vector{SVector}(undef, T+1)
    else
        trace_states = Vector{SVector}(undef, 0)
    end
    nI = length(patient_zeros)
    mcont = true
    Iidx = Set(patient_zeros)
    new_inf = Set()
    while (mcont)
        ## set recovered state
        
        #rec_mask = @. (t >= infect_t[Iidx] + delays[Iidx]+1)
        #states[Iidx[rec_mask]] .= 3
        new_rec = Set()
        for i in Iidx
            if (t >= infect_t[i] + delays[i]+1)
                states[i] = 3
                push!(new_rec, i)
            end
        end
        #states[mask_rec] .= 3
        #nS =  (sum(states.==1))
        setdiff!(Iidx, new_rec)
        union!(Iidx, new_inf)
        nI = sum(states.==2)
        @assert nI  == length(Iidx)
        #nR= N - nI - nS 
        if useT
            # use t+1 because of indexing
            trace_states[t+1] = count_states(counter,states) #[nS, nI, nR]
            ## set here the flag
            mcont = (t+1 <= T)
        else 
            push!(trace_states, count_states(counter,states))
            mcont = nI > 0
        end

        ### State changer
        if (!isnothing(dynstateChanger))
            change_states_dyn(dynstateChanger, model, states, G, t, Iidx)
        end

        c=0
        #Iidx = findall(states.==2) # this will be used also later
        empty!(new_inf)
        for i in Iidx
            ## TODO: parallel runs need to lock the graph object and then unlock it
            for j in neighbors(G,i)
                if (states[j] == 1) ## double check to be sure
                    ## Here, multiple dispatch will help distinguish 
                    ## when it's a vector of probabilities or just a single one for everyone
                    if rand(rng) < calc_prob_infection(model, i, j, beta_IorS)
                        ## infected
                        infect_t[j] = t
                        infect_i[j] = i
                        states[j] = 2
                        push!(new_inf, j)
                        c+=1
                    end
                end
            end #for neighbors
        end #for infected
        #println("$c new infected at time $t")
        
        t+=1
    end #for time loop
    states, trace_states
end
"""
`run_sir_fast(g::AbstractGraph, model::AbstractSIRModel, T::Integer, rng::AbstractRNG, patient_zeros::Vector{<:Integer}; beta_IorS=:I, dtype=Float64, counter=BaseSIRStatesCounter())`

Sets up and runs a fast SIR and SIR-like simulation, returning simulation data and state counts.

The `beta_IorS` argument is used only with the default SIR model (::SIRModel). Other models have their own method for computation of the infection probability.

### Optional arguments
- `beta_IorS`: Infection rule — whether based on infected (`:I`), susceptible (`:S`), or their average if other value. Defaults to `:I`.
- `dtype=Float64`: Data type for output arrays.
- `counter=BaseSIRStatesCounter()`: State counting method.
- `dynstateChanger`: object used for changing the state dynamically

"""
function run_sir_fast(g::AbstractGraph, model::AbstractSIRModel, T::Integer, rng::AbstractRNG, 
    patient_zeros::Vector{<:Integer}; beta_IorS::Symbol=:I, dtype::DataType=Float64,
    counter::AbstractStatesCounter=BaseSIRStatesCounter(), dynstateChanger::Union{Nothing, AbstractStateChanger}=nothing)
    ## draw recovery delays
    N= nv(g)
    nodes = collect(Int32,1:N)
    delays = draw_delays(model, model.gamma, rng, nodes)

    data = SIRSimData(delays,N, dtype)
    endstate, cc = sim_sir_fast(g, model, T, data, rng, patient_zeros, beta_IorS=beta_IorS,counter=counter, dynstateChanger=dynstateChanger)

    data, cc
end

function calc_n_comparts(infect_times,delays,T::Integer)
    ti= infect_times .+ 1
    isinfect= collect(0:T) .>= (ti)'
    isrec = collect(0:T) .>= (ti.+delays)'
    nR = sum(isrec,dims=2)
    nI = sum(isinfect,dims=2).-nR
    nS = sum(.~isinfect,dims=2)
    hcat(nS,nI,nR)
end

function calc_n_states(infects, delays, T::Integer)
    #counts = zeros(T+1) ## from 0 to T
    counts = OffsetArray(zeros(Int,T+1,3), 0:T,1:3)

    N = length(delays)
    for i=1:N
        if isnan(infects[i])
            ## never infected
            counts[:,1].+=1
            continue
        end
        ti = convert(Int,infects[i])+1
        if ti < 0
            ## never infected
            counts[:,1] .+= 1
            continue
        end
        trec = min(ti + delays[i], T)
        counts[0:(ti-1),1] .+=1 ## times S
        counts[ti:(trec-1),2] .+= 1
        counts[trec:end,3].+=1
    end
    counts
end

function calc_nstates_all(infects,recovs,T::Integer)
    dats = [calc_n_comparts(ii,dd, T) for (ii,dd) in zip(infects,recovs)]

    counts=cat(dats..., dims=3)
end

function calc_Rarr(inf,N, T, tmax, win_size=7)
    newinf = zeros(Int,T+2)
    succinf = zeros(Int,T+2)
    
    for i=1:N
        if isnan(inf[i,1])
            continue
        end
        ti=convert(Int,inf[i,1])
        newinf[ti+2]+=1
        if !isnan(inf[i,2]) && inf[i,2]>=0
            ji = convert(Int,inf[i,2])
            
            tj = convert(Int, inf[ji,1])
            if tj >= 0
                succinf[tj+2]+=1
            end
        end
    end
    Rarr = zeros(tmax)
    ws = win_size
    for t=ws:tmax
        Rarr[t]=sum(succinf[max(t-ws,1):t])/sum(newinf[max(t-ws,1):t]) 
    end
    newinf, succinf, Rarr
end