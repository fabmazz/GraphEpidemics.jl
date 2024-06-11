module GraphEpidemics

using Graphs
using Distributions: Geometric, Exponential
using Random

export AbstractEpiModel, SIRModel, model_states, states_values, spreading_states, trans_independent, init_model_discrete, run_complex_contagion, draw_delays, set_state_nodes

export sim_sir_fast, calc_n_comparts, calc_nstates_all, calc_Rarr, calc_n_states, beta_R0, run_sir_fast, prob_from_rate

export sim_sir_gillespie, run_sir_gillespie, gillespie_sir_direct, BinaryTree

include("models.jl")
include("sim.jl")
include("complex_contagion.jl")
include("recurrent.jl")
include("binary_tree.jl")
include("gillespie.jl")

export SISModel, SIRSModel, SIRModelSus

end # module GraphEpidemics
