module GraphEpidemics

using Graphs
using Distributions: Geometric, Exponential
using Random

import DataStructures: DefaultDict, PriorityQueue, enqueue!, peek, delete!, dequeue!
using StaticArrays: SVector


export AbstractEpiModel, SIRModel, model_states, states_values, spreading_states, trans_independent, init_model_discrete, run_complex_contagion, draw_delays, set_state_nodes

export sim_sir_fast, calc_n_comparts, calc_nstates_all, calc_Rarr, calc_n_states, beta_R0, run_sir_fast, prob_from_rate, count_states

export sim_sir_gillespie, run_sir_gillespie, gillespie_sir_direct, BinaryTree

include("types.jl")
include("models.jl")
include("utils.jl")
include("sim.jl")
include("sim_seir.jl")
include("complex_contagion.jl")
include("recurrent.jl")
include("binary_tree.jl")
include("gillespie.jl")

export SISModel, SIRSModel, SIRModelSus, SEIRModel, SEIRHetModel
export sim_seir_fast, run_seir_fast
export AbstractStateChanger

end # module GraphEpidemics
