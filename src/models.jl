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
abstract type AbstractSIRModel <: AbstractEpiModel end

struct SIRModel{BT<:Union{Vector{<:AbstractFloat}, AbstractFloat},
    GT<:Union{Vector{<:AbstractFloat}, AbstractFloat}} <: AbstractSIRModel
    beta::BT##vector is the 
    gamma::GT
    #stateType::DataType
end

struct SIRModelSus{F1<:UVFloat, F2<:UVFloat, F3<:UVFloat} <:AbstractSIRModel
    beta::F1
    sigma::F2
    gamma::F3
end

model_states(x::AbstractSIRModel) = (:S,:I,:R)

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

draw_delays(m::SIRModel, p, rng::AbstractRNG, nodes::Union{Integer, Vector{<:Integer}}) = draw_delays_markov(p, rng, nodes)


## Real time - Exponential
draw_delay_exp(p::AF, rng::AbstractRNG) = rand(rng, Exponential(1/p))
draw_delay_exp(p::AF, rng::AbstractRNG, i::Integer) = rand(rng, Exponential(1/p))
draw_delay_exp( p::AF, rng::AbstractRNG, nodes)= rand(rng, Exponential(1/p), length(nodes))