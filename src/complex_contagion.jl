using Random
using Distributions
#=
struct EpiModel <: AbstractEpiModel
    states::Vector{Symbol}
    infection_scheme::Vector{Tuple{Symbol, Symbol, Symbol}} ## e.g. (:I,:S,:I) -> :I infects :S turning them into :I
    individual_trans::Dict{Symbol, Symbol}
end
=#
#=
Define spreading states
=#
abstract type AbstractEpiModel end

struct SIRModel{F<:AbstractFloat,EI<:Integer} <: AbstractEpiModel
    beta::F ##vector is the 
    gamma::F
    #stateType::DataType
end

const NULLT::Int32 = -100

struct SimData{F<:AbstractFloat,I<:Integer}
    N::Integer
    infect_times::Vector{F}
    infect_node::Vector{F}
    states_repr::Dict{Symbol,Int8}
    last_trans_time::Vector{I}
    transition_delays::Array{I,2}
    epistate::Vector{Int8}
    additional_log::Vector{Vector{I}} #unused, (node, node_infector, time_transitions (multiple))
end


#SIRModel(beta::F, gamma::F) = SIRModel(beta,gamma,StatesSIR)


states(x::SIRModel) = (:S,:I,:R)
function states_values(x::AbstractEpiModel)
    d::Dict{Symbol,Int8} =  Dict(s=>i for (i,s) in enumerate(states(x)))
    d
end
#spreading_state(x::SIRModel) = :I
spreading_states(x::SIRModel) = Dict(:I=>[(:S,:I, x.beta)])

trans_independent(x::SIRModel) = [(:I,:R, x.gamma)]

draw_delays(m::SIRModel, p::Real, rng::AbstractRng, i::Integer) = rand(rng, Geometric(p))+1
draw_delays(m::SIRModel, p::Vector{F}, rng::AbstractRng, i::Integer) where F<:AbstractFloat = rand(rng, Geometric(p[i]))+1
is_spreading(p::Real, rng::AbstractRng, i::Integer, j::Integer) = rand(rng)<p
is_spreading( p::Vector{F}, rng::AbstractRng, i::Integer, j::Integer) = rand(rng)<p[i]

function init_model_discrete(model::AbstractEpiModel, g::AbstractGraph, nodes_active::Vector{I}, 
    init_state_inactive::Symbol, init_state_active::Symbol, rng::AbstractRNG) where I<:Integer
    N =nv(g)
    infect_t = fill(NaN, N)
    infect_i = fill(NaN, N)
    all_states = states(model)
    sval::Dict{Symbol,Int8} = Dict(s=>i for (i,s) in enumerate(all_states))
    init_state = sval[init_state_inactive]
    infect_state = sval[init_state_active]

    states::Vector{Int8} = fill(init_state, N)
    last_trans_time = fill(-1000, N)

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
    SimData(N, infect_t, infect_i, sval, last_trans_time, delays_trans, states, [])
end

function simulate_model_discrete(model::AbstractEpiModel, g::AbstractGraph,T::Integer, rng::AbstractRNG, data::SimData,
    spreading_function::Function = is_spreading
    )
    N = nv(g)
    @assert N == data.N

    #track = fill(NaN, (N,2))

    all_states = states(model)
    sval =data.states_repr
    for (i,st) in enumerate(all_states)
        @assert sval[st] == i
    end
  
    states::Vector{Int8} = data.epistate
    infect_t = data.infect_times
    infect_i = data.infect_node

    indep_transitions = Dict(sval[x[1]]=>(sval[x[2]], x[3]) for x in trans_independent(model))
    delays_trans = data.transition_delays
    last_trans_time = data.last_trans_time
    @assert size(delays_trans,2) == length(keys(indep_transitions))

    println("Independent transitions: $indep_transitions")

    spreading_trans = Dict(sval[k]=>[(sval[x[1]], sval[x[2]],x[3]) for x in vals] for (k,vals) in spreading_states(model))
    for t =0:T
        ## independent transitions
        ## check 
        c = 1
        for (st_check,vals) in indep_transitions
            ns, p = vals
            for i in findall(states.==st_check)
                ##check last transition time + delay against t
                if t >= last_trans_time[i] + delays_trans[i,c]
                    ## transitioned
                    states[i] = ns
                end
            end
            c+=1 ##increment counter for transitions
        end

        #println("Have $nS S $(sum(states.==1)) I $nR R")
        c=0
        for (st_from,ext) in spreading_trans
            st_to, newst, prob = ext
            for i in findall(states.==st_from)
                for j in neighbors(g,i)
                    if (states[j] == st_to) 
                        if is_spreading(prob, rng, i, j)
                            ## infected
                            infect_t[j] = t
                            infect_i[j] = i
                            states[j] = newst
                            c+=1
                        end
                    end
                end #for neighbors
            end #for infected
        end
        #println("$c new infected at time $t")

    end #for time loop
    data
end