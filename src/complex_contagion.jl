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

is_spreading(p::Real, rng::AbstractRNG, i::Integer, j::Integer) = rand(rng)<p
is_spreading( p::Vector{F}, rng::AbstractRNG, i::Integer, j::Integer) where F<:AbstractFloat = rand(rng)<p[i]

function extract_delays(model::AbstractEpiModel,indep_transitions::Dict{StI,<:IndepTrans}, 
    delays_trans::Matrix{II}, rng::AbstractRNG, i::Integer) where II <:Integer
    for (st_check,trans) in indep_transitions
        #ns, p, c = vals
        #p = dat[3]
        delays_trans[i,trans.matidx] = draw_delays(model, trans.prob, rng, i)
    end
end

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
    last_trans_time::Vector{typeof(NULLT)} = fill(-1000, N)

    indep_transitions = trans_independent(model)
    n_trans = length(indep_transitions)
    delays_trans = fill(NULLT,(N,n_trans))

    
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

#Dict{StI,Tuple{StI,F,I}}
process_indep_trans(model::AbstractEpiModel, sval::Dict{Symbol,StI}) = Dict(sval[x[1]]=>IndepTrans(sval[x[2]], x[3], i, sval[x[2]]<sval[x[1]]) for (i,x) in enumerate(trans_independent(model)) ) 

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
function process_spreading_states(model::AbstractEpiModel, sval::Dict{Symbol,StI})
    Dict(sval[k]=>Dict(sval[x[1]]=> StateTo(sval[x[2]],x[3]) for x in vals) for (k,vals) in spreading_states(model))
end


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