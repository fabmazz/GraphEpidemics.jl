#= 
 Copyright (c) 2024 Fabio Mazza
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at https://mozilla.org/MPL/2.0/.
=#

import DataFrames: AbstractDataFrame

const LARGET::Int32 =10000

struct NodeHistory{I<:Integer}
    idx::Int
    i_inf::Int32
    t_inf::Int32
    trans_delays::Vector{I}
end

struct SimData{F<:AbstractFloat,I<:Integer}
    N::Int
    infect_times::Vector{F}
    infect_node::Vector{F}
    states_repr::Dict{Symbol,StI}
    last_trans_time::Vector{I}
    transition_delays::Array{I,2}
    epistate::Vector{StI}
    additional_log::Vector{NodeHistory}
end
struct IndepTrans{F<:Union{AbstractFloat,Vector{<:AbstractFloat}}, I<:Integer}
    stateto::StI
    prob::F
    matidx::I
    reverse::Bool
end
"""
`is_spreading(p::Real, rng::AbstractRNG, i::Integer, j::Integer)`

Determines whether a spreading event occurs using a scalar probability, ignoring `i` and `j`
"""
is_spreading(p::Real, rng::AbstractRNG, i::Integer, j::Integer) = rand(rng)<p

"""
`is_spreading(p::Vector{F}, rng::AbstractRNG, i::Integer, j::Integer) where F<:AbstractFloat`

Determines whether a spreading event occurs using a node-specific probability vector.
Returns true if a random draw is less than `p[i]`, indicating successful infection from node `i`.
"""
is_spreading(p::Vector{F}, rng::AbstractRNG, i::Integer, j::Integer) where F<:AbstractFloat = rand(rng)<p[i]

"""
`extract_delays(model::AbstractEpiModel, indep_transitions::Dict{StI,<:IndepTrans}, delays_trans::Matrix{II}, rng::AbstractRNG, i::Integer) where II <: Integer`

Draws and assigns delay times for all transitions associated with a specific node.

Given a model and its independent transitions, this function updates the `i`-th row of the delay matrix 
by drawing new delay values using the provided random number generator.
"""
function extract_delays(model::AbstractEpiModel,indep_transitions::Dict{StI,<:IndepTrans}, 
    delays_trans::Matrix{II}, rng::AbstractRNG, i::Integer) where II <:Integer
    for (st_check,trans) in indep_transitions
        #ns, p, c = vals
        #p = dat[3]
        delays_trans[i,trans.matidx] = draw_delays(model, trans.prob, rng, i)
    end
end

"""
`init_model_discrete(model::AbstractEpiModel, g::AbstractGraph, nodes_active::Vector{I}, init_state_inactive::Symbol, init_state_active::Symbol, rng::AbstractRNG) where I <: Integer`

Initializes the simulation state with some nodes actively infected and others inactive.

The graph `g` provides the node structure. Nodes in `nodes_active` are assigned the `init_state_active`,
while all others receive `init_state_inactive`. Delay values and infection metadata are initialized accordingly.
"""
function init_model_discrete(model::AbstractEpiModel, g::AbstractGraph, nodes_active::Vector{I}, 
    init_state_inactive::Symbol, init_state_active::Symbol, rng::AbstractRNG) where I<:Integer
    N =nv(g)
    infect_t = fill(NaN, N)
    infect_i = fill(NaN, N)
    all_states = model_states(model)
    sval::Dict{Symbol,StI} = Dict(s=>i for (i,s) in enumerate(all_states))
    init_state = sval[init_state_inactive]
    infect_state = sval[init_state_active]

    states::Vector{StI} = fill(init_state, N)
    last_trans_time::Vector{typeof(LARGET)} = fill(-1000, N)

    indep_transitions = trans_independent(model)
    n_trans = length(indep_transitions)
    delays_trans = fill(LARGET,(N,n_trans))

    
    for i in nodes_active
        states[i] = infect_state
        infect_t[i] = -1 ## time of infection
        last_trans_time[i] = 0 ##set last transition time 
        infect_i[i] = -10 ## node of infection
        ## extract delays for those in transitions
        for (c,dat) in enumerate(indep_transitions)
            p = dat[3]
            delays_trans[i,c] = draw_delays(model, p, rng, i)
        end
    end
    SimData(N, infect_t, infect_i, sval, last_trans_time, delays_trans, states, Array{Vector{Int64},1}(undef, 0))
end

"""
`process_indep_trans(model::AbstractEpiModel, sval::Dict{Symbol, StI})`

Converts model-defined independent transitions into internal transition structures.

Uses the state value dictionary `sval` to encode symbolic state names as integers.
Each independent transition is converted into an `IndepTrans` struct and stored in a dictionary.
"""
#Dict{StI,Tuple{StI,F,I}}
process_indep_trans(model::AbstractEpiModel, sval::Dict{Symbol,StI}) = Dict(sval[x[1]]=>IndepTrans(sval[x[2]], x[3], i, sval[x[2]]<sval[x[1]]) for (i,x) in enumerate(trans_independent(model)) ) 

"""
`set_state_nodes(model::AbstractEpiModel, data::SimData, rng::AbstractRNG, nodes::Vector{I}, state::Symbol, active::Bool; tset::Integer=0) where I <: Integer`

Assigns a state to specified nodes and activates them if necessary.

Updates the states of all nodes in `nodes` to `state`. If the state is marked as active or `active` is true,
it updates infection metadata and transition delays using the model and random number generator.

### Optional Arguments
- `tset`: The time to assign for infection and transition tracking (default: 0).
"""
function set_state_nodes(model::AbstractEpiModel, data::SimData, rng::AbstractRNG, 
        nodes::Vector{I}, state::Symbol, active::Bool, tset::Integer=0) where I<:Integer

    all_states = model_states(model)
    sval =data.states_repr
    for (i,st) in enumerate(all_states)
        @assert sval[st] == i
    end
    states::Vector{StI} = data.epistate

    st_int = sval[state]
    states[nodes] .= st_int
    if state in first_active_states(model)
        data.infect_times[nodes] .= tset-1
        data.infect_node[nodes] .= -10
        active = true
    end
    if active
        indep_transitions = process_indep_trans(model, sval)
        data.last_trans_time[nodes] .= tset
        for i in nodes
            extract_delays(model,indep_transitions,data.transition_delays,rng, i)
        end
    end
end

"""
`init_model_discrete(model::AbstractEpiModel, g::AbstractGraph, rng::AbstractRNG, state_base::Symbol)`

Initializes all nodes to the same base state in a simulation without any active infections.

All nodes are set to `state_base`, and default values are assigned for delays, states, and infection history.
Returns a `SimData` struct.
"""
function init_model_discrete(model::AbstractEpiModel, g::AbstractGraph, rng::AbstractRNG, state_base::Symbol) 
    N =nv(g)
    infect_t = fill(NaN, N)
    infect_i = fill(NaN, N)
    all_states = model_states(model)
    sval::Dict{Symbol,StI} = Dict(s=>i for (i,s) in enumerate(all_states))
    init_state = sval[state_base]

    states::Vector{StI} = fill(init_state, N)
    last_trans_time::Vector{typeof(LARGET)} = fill(-1000, N)

    indep_transitions = trans_independent(model)
    n_trans = length(indep_transitions)
    delays_trans = fill(LARGET,(N,n_trans))
    
    #for (i, trans) in enumerate(indep_transitions)
    #    delays_trans[:,i] = draw_delays_nodes(model, trans[3], rng, collect(1:N) )
    #end
    SimData(N, infect_t, infect_i, sval, last_trans_time, delays_trans, states, Vector{NodeHistory}(undef, 0))
end
struct StateTo{F<:Union{AbstractFloat,Vector{<:AbstractFloat}}}
    st::StI
    prob::F
end
"""
`process_spreading_states(model::AbstractEpiModel, sval::Dict{Symbol, StI})`

Builds a nested dictionary representation of spreading transitions.

Each source state maps to a dictionary of neighbor states, which map to destination states and their associated
probabilities, represented using `StateTo` structs.
"""
function process_spreading_states(model::AbstractEpiModel, sval::Dict{Symbol,StI})
    Dict(sval[k]=>Dict(sval[x[1]]=> StateTo(sval[x[2]],x[3]) for x in vals) for (k,vals) in spreading_states(model))
end

"""
`run_complex_contagion(model::AbstractEpiModel, g::AbstractGraph, T::Integer, rng::AbstractRNG, data::SimData; spreading_function::Function = is_spreading, verbose::Bool = false)`

Simulates the spread of contagion across a network over a fixed number of time steps.

Applies independent and spreading transitions from time `0` to `T`, updating the simulation state in-place.
Returns a vector of NamedTuples tracking the number of nodes in each state at each timestep.

### Optional Arguments
- `spreading_function`: A custom function to determine if infection spreads between neighbors.
- `verbose`: If true, prints additional debug information.
"""
function run_complex_contagion(model::AbstractEpiModel, g::AbstractGraph,T::Integer, rng::AbstractRNG, data::SimData;
    spreading_function::Function = is_spreading, verbose=false,
    )
    N = nv(g)
    @assert N == data.N

    #track = fill(NaN, (N,2))

    all_states = model_states(model)
    sval =data.states_repr
    for (i,st) in enumerate(all_states)
        @assert sval[st] == i
    end
  
    states::Vector{StI} = data.epistate
    infect_t = data.infect_times
    infect_i = data.infect_node

    indep_transitions = process_indep_trans(model, sval)
    delays_trans = data.transition_delays
    last_trans_time = data.last_trans_time
    @assert size(delays_trans,2) == length(keys(indep_transitions))

    if verbose
        println("Independent transitions: $indep_transitions")
    end

    spreading_trans = process_spreading_states(model, sval)
    first_act_states = Set(sval[k] for k in first_active_states(model))
    if verbose
        println("Spreading transitions: $spreading_trans")
    end
    num_states=NamedTuple[]
    for t =0:T
        ## independent transitions
        ## check 
        #c = 1
        for (st_check,trans) in indep_transitions
            for i in findall(states.==st_check)
                ##check last transition time + delay against t
                if t >= last_trans_time[i] + delays_trans[i,trans.matidx]
                    ## transitioned
                    states[i] = trans.stateto
                    last_trans_time[i] = t
                    if trans.reverse
                        ## we have to save the stats
                        #datv = [i,convert(Int, infect_i[i]), convert(Int, infect_t[i]), data.transition_delays[i,:]... ]
                        push!(data.additional_log, NodeHistory(i, convert(Int32, infect_i[i]),convert(Int32, infect_t[i]), data.transition_delays[i,:]))
                        data.transition_delays[i,:] .= LARGET
                    end
                    #println("at t=$t $i transitions from $st_check ($(all_states[st_check])) -> $ns ($(all_states[ns]))")
                end
            end
            #c+=1 ##increment counter for transitions
        end
        ##record statistics
        push!(num_states,
        (; t = t, zip(keys(sval),  (sum(states .== sval[k]) for k in keys(sval)))... ) # NamedTuple
        )

        #println("Have $nS S $(sum(states.==1)) I $nR R")
        cinf=0
        for (st_from,to_dict) in spreading_trans
            #st_to, newst, prob = ext
            outstate = Set(keys(to_dict))
            for i in findall(states.==st_from)
                for j in neighbors(g,i)
                    if states[j] in outstate #(states[j] == st_to)
                        outst = to_dict[states[j]] 
                        if spreading_function(outst.prob, rng, i, j)
                            ## infected
                            #println("$j infected, from $(states[j]) to $newst")
                            infect_t[j] = t
                            infect_i[j] = i
                            states[j] = outst.st
                            last_trans_time[j] = t+1
                            ## extract transitions
                            if outst.st in first_act_states
                                extract_delays(model,indep_transitions,delays_trans,rng, j)
                            end
                            
                            cinf+=1
                        end
                    end
                end #for neighbors
            end #for infected
        end
        #println("$c new infected at time $t")

    end #for time loop
    num_states
end

"""
`group_counts_sims(counts::Vector{<:AbstractDataFrame}, model::AbstractEpiModel)`

Aggregates state count data from multiple simulations into a single array.

Takes a list of DataFrames, each representing one simulation, and combines them into a 3D array
organized by timestep, state, and simulation index.
"""
function group_counts_sims(counts::Vector{<:AbstractDataFrame}, model::AbstractEpiModel)
    nsims = length(counts)
    T=size(counts[1],1)
    res = zeros(Int64,T,3,nsims);

    for (i,c) in enumerate(counts)
        sort!(c, [:t])
        for (n,s) in enumerate(model_states(model))
            res[:,n,i] = c[:,s]
        end
    end
    res
end