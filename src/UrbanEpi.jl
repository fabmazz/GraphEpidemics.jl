module UrbanEpi

using Graphs
using Distributions
using Random

greet() = print("Hello World!")

include("sim.jl")

export sim_sir, calc_n_comparts, calc_nstates_all, calc_Rarr, calc_n_states


end # module UrbanEpi
