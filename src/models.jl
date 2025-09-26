#= 
 Copyright (c) 2024 Fabio Mazza
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at https://mozilla.org/MPL/2.0/.
=#

StI = Int8

prob_from_rate(r::Real, t::Real) = 1 - exp(-r*t)

UVFloat = Union{Vector{<:AbstractFloat}, AbstractFloat}
VFloat = Vector{<:AbstractFloat}

function states_values(x::AbstractEpiModel)
    d::Dict{Symbol,StI} =  Dict(s=>i for (i,s) in enumerate(model_states(x)))
    d
end


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


"""
`SIRModelHet{beta::Float, alpha, sigma, gamma} <: AbstractSIRModel`

SIR model where the probability of infection depends on both the infecting node and the susceptible node.

The `beta` parameter is a scalar which is the disease infectivity
Parameters `alpha` and `sigma` are the multipliers for the infectivity and susceptibility (resp)
and must be vectors

Recovery probability `gamma` can be either a scalar (for the whole population) or a vector

"""
struct SIRModelHet{F0<:AbstractFloat,F1<:VFloat, F2<:VFloat, F3<:UVFloat} <: AbstractSIRModel
    beta_d::F0
    alpha::F1
    sigma::F2
    gamma::F3
end

model_states(x::AbstractSIRModel) = (:S,:I,:R)
"""
`gamma(::model)`

Every model <: AbstractSIRModel has to implement this method to extract the recovery probability
if the γ is given in a different way (`gamma` property)
"""
gamma(x::AbstractSIRModel) = x.gamma
#gamma(x::SIRModel) = x.gamma
#gamma(x::SIRModelSus) = x.gamma



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
draw_delays(m::AbstractSIRModel, p, rng::AbstractRNG, nodes::Union{Integer, Vector{<:Integer}}) = draw_delays_markov(p, rng, nodes)
#draw_delays(m::SIRModel, p, rng::AbstractRNG, nodes::Union{Integer, Vector{<:Integer}}) = draw_delays_markov(p, rng, nodes)
#draw_delays(m::SIRModelSus, p, rng::AbstractRNG, nodes::Union{Integer, Vector{<:Integer}}) = draw_delays_markov(p, rng, nodes)


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

#=== FOR THE COMPLEX CONTAGION ROUTINES ===#

spreading_states(x::SIRModel) = Dict(:I=>[(:S,:I, x.beta)])
trans_independent(x::SIRModel) = [(:I,:R, x.gamma)]
first_active_states(x::AbstractSIRModel) = (:I,)


#======== INFECTION PROBABILITY CALCULATION ======#


"""
`get_p_infection(p::AbstractFloat, i_I::Integer, j_S::Integer, infect_prob_IS::Symbol)`

Returns the infection probability `p`. This method exists for uniformity.
"""
function get_p_infection(p::AbstractFloat, i_I::Integer, j_S::Integer, infect_prob_IS::Symbol)
    p
end

"""
`get_p_infection(p::Vector{<:AbstractFloat}, i_I::Integer, j_S::Integer, infect_prob_IS::Symbol)`

Returns the infection probability based on the infecting when `infect_prob_IS` is `:I` or susceptible node (`:S`), or their geometric mean when `infect_prob_IS` is neither value.
"""
function get_p_infection(p::Vector{<:AbstractFloat}, i_I::Integer, j_S::Integer, infect_prob_IS::Symbol)
    if infect_prob_IS == :I
        return p[i_I]
    elseif infect_prob_IS == :S
        return p[j_S]
    else
        return sqrt(p[i_I]*p[j_S])
    end
end
"""
`calc_prob_infection(model::SIRModel, i::Integer, j::Integer, infect_prob_IS::Symbol)`

Helper method to compute the probability of infection for the SIR model.
"""
function calc_prob_infection(model::SIRModel, i::Integer, j::Integer, infect_prob_IS::Symbol)
    get_p_infection(model.beta, i, j, infect_prob_IS)
end
### Average pij
function get_p_product(p1::Real, p2::Real, i1::Integer, i2::Integer)
    (p1*p2)
end
function get_p_product(p1::Vector{<:Real}, p2::Vector{<:Real}, i1::Integer, i2::Integer)
    (p1[i1]*p2[i2])
end
## TODO: conver cases when p1 is Vector, p2 is not and viceversa

"""
`calc_prob_infection(model::SIRModelSus, i::Integer, j::Integer, ignored_s::Symbol)`

Computes infection probability for the SIR model with susceptible
"""
function calc_prob_infection(model::SIRModelSus, i::Integer, j::Integer, ignored_s::Symbol)
    get_p_product(model.beta, model.sigma, i, j)
end

function calc_prob_infection(model::SIRModelHet, i::Integer, j::Integer, ignored_s::Symbol)
    model.beta_d * get_p_product(model.alpha,model.sigma, i, j)
end