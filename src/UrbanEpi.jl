module UrbanEpi

using Graphs
using Distributions
using Random

greet() = print("Hello World!")

include("sim.jl")
include("complex_contagion.jl")

export sim_sir, calc_n_comparts, calc_nstates_all, calc_Rarr, calc_n_states, beta_R0

export AbstractEpiModel, SIRModel, model_states, states_values, spreading_states, trans_independent, init_model_discrete, run_complex_contagion, draw_delays, set_state_nodes


end # module UrbanEpi
