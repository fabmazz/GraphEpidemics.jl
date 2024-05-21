module GraphEpidemics

using Graphs
using Distributions: Geometric, Xoshiro
using Random

export AbstractEpiModel, SIRModel, model_states, states_values, spreading_states, trans_independent, init_model_discrete, run_complex_contagion, draw_delays, set_state_nodes

export sim_sir_fast, calc_n_comparts, calc_nstates_all, calc_Rarr, calc_n_states, beta_R0, run_sir_fast


include("models.jl")
include("sim.jl")
include("complex_contagion.jl")
include("recurrent.jl")

export SISModel, SIRSModel

end # module GraphEpidemics
