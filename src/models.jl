#= 
 Copyright (c) 2024 Fabio Mazza
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at https://mozilla.org/MPL/2.0/.
=#
"""
`AbstractEpiModel`

Abstract supertype for all epidemiological models. To be used in `run_complex_contagion`, subtypes must implement
- `model_states(::model)` which returns a tuple of symbols that describe the model states
- `spreading_states(::model)` returning a dict of transitions describing the infector state, the state of soon to be infected (from and to) and the rate, e.g. `Dict(:I=>[(:S,:I, x.beta)])`
- `trans_independent(::model)` which is a list of all the independent transitions with their rate (e.g., `[(:I,:R, x.gamma)]` for the SIR model
- `first_active_states(::model)` which describes the seed states (e.g. `(:I,)` for the SIR)
"""
abstract type AbstractEpiModel end

StI = Int8

prob_from_rate(r::Real, t::Real) = 1 - exp(-r*t)

UVFloat = Union{Vector{<:AbstractFloat}, AbstractFloat}

function states_values(x::AbstractEpiModel)
    d::Dict{Symbol,StI} =  Dict(s=>i for (i,s) in enumerate(model_states(x)))
    d
end

#=struct SIRModel{F<:AbstractFloat} <: AbstractEpiModel
    beta::Union{F,Vector{F}}##vector is the 
    gamma::Union{F,Vector{F}}
    #stateType::DataType
end
=#
"""
`AbstractSIRModel`

Abstract supertype for all SIR epidemiological models. Subtypes must implement a method for computing infection probabilities.
This model is used for `run_sir_fast` and `run_sir_gillespie`.
"""
abstract type AbstractSIRModel <: AbstractEpiModel end

"""
`SIRModel{beta, gamma} <: AbstractSIRModel`

Standard SIR model with infection (`beta`) and recovery (`gamma`) rates. Rates can be either scalars or node-specific vectors.
"""
struct SIRModel{BT<:Union{Vector{<:AbstractFloat}, AbstractFloat},
    GT<:Union{Vector{<:AbstractFloat}, AbstractFloat}} <: AbstractSIRModel
    beta::BT##vector is the 
    gamma::GT
    #stateType::DataType
end
"""
`SIRModelSus{beta, sigma, gamma} <: AbstractSIRModel`

SIR model with heterogeneous susceptibility `sigma`, where the probability of infection depends on both the infecting node and the susceptible node.

All parameters can be either vectors or scalars.
"""
struct SIRModelSus{F1<:UVFloat, F2<:UVFloat, F3<:UVFloat} <:AbstractSIRModel
    beta::F1
    sigma::F2
    gamma::F3
end

model_states(x::AbstractSIRModel) = (:S,:I,:R)
"""
`gamma(::model)`

Every model <: AbstractSIRModel has to implement this method to extract the recovery probability
"""
gamma(x::SIRModel) = x.gamma
gamma(x::SIRModelSus) = x.gamma

#spreading_state(x::SIRModel) = :I
spreading_states(x::SIRModel) = Dict(:I=>[(:S,:I, x.beta)])

trans_independent(x::SIRModel) = [(:I,:R, x.gamma)]
first_active_states(x::SIRModel) = (:I,)
AF=AbstractFloat
draw_delay_geom(p::AF, rng::AbstractRNG) = rand(rng, Geometric(p))+1
draw_delay_geom(p::AF, rng::AbstractRNG, i::Integer) = rand(rng, Geometric(p))+1
draw_delay_geom( p::AF, rng::AbstractRNG, nodes::Vector{<:Integer})= rand(rng, Geometric(p), length(nodes)) .+1


draw_delays_markov(p::AbstractFloat,rng::AbstractRNG, nodes::Union{Integer,Vector{<:Integer}}) = draw_delay_geom(p, rng, nodes)
draw_delays_markov(p::Vector{<:AbstractFloat}, rng::AbstractRNG, nodes::Vector{<:Integer}) = draw_delay_geom.(p,rng, nodes)
draw_delays_markov(p::Vector{<:AbstractFloat}, rng::AbstractRNG, i::Integer) = draw_delay_geom(p[i], rng, i)
"""
`draw_delays(m::SIRModel or SIRModelSus, p, rng::AbstractRNG, nodes::Union{Integer, Vector{<:Integer}})`

Draw delays according to the model, with `p` that can be either a vector or a single value.
"""
draw_delays(m::SIRModel, p, rng::AbstractRNG, nodes::Union{Integer, Vector{<:Integer}}) = draw_delays_markov(p, rng, nodes)
draw_delays(m::SIRModelSus, p, rng::AbstractRNG, nodes::Union{Integer, Vector{<:Integer}}) = draw_delays_markov(p, rng, nodes)


## Real time - Exponential
"""
`draw_delay_exp(p, rng::AbstractRNG)`

This function draws an exponentially distributed variable with rate p.
"""
draw_delay_exp(p::AF, rng::AbstractRNG) = rand(rng, Exponential(1/p))
"""
`draw_delay_exp(p, rng::AbstractRNG, i::Integer)`

This function draws an exponentially distributed variable with rate p, for index `i` (if `p` is a Vector)
"""
draw_delay_exp(p::AF, rng::AbstractRNG, i::Integer) = rand(rng, Exponential(1/p))
draw_delay_exp(p::Vector{<:AF}, rng::AbstractRNG, i::Integer) = rand(rng, Exponential(1/p[i]))
"""
`draw_delay_exp(p, rng::AbstractRNG, nodes::Vector)`

This function draws `m` exponentially distributed variables with rate `p`, where `m` is the length of the `nodes`
"""
draw_delay_exp(p::AF, rng::AbstractRNG, nodes)= rand(rng, Exponential(1/p), length(nodes))

#### SEIR Model ####
struct SEIRModel{F<:AbstractFloat} <: AbstractEpiModel
    eta::F
    beta::F ##vector is the 
    gamma::F
    #stateType::DataType
end

model_states(x::SEIRModel) = (:S,:E,:I,:R)
spreading_states(x::SEIRModel) = Dict(:I=>[(:S,:E, x.beta)])
trans_independent(x::SEIRModel) = [(:E,:I,x.eta),(:I,:R, x.gamma)]

first_active_states(m::SEIRModel) = (:E,)
draw_delays(m::SEIRModel, p, rng::AbstractRNG, node) = draw_delays_markov(p, rng, node )
