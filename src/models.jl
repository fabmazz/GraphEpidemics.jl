abstract type AbstractEpiModel end

function states_values(x::AbstractEpiModel)
    d::Dict{Symbol,Int8} =  Dict(s=>i for (i,s) in enumerate(model_states(x)))
    d
end

#=struct SIRModel{F<:AbstractFloat} <: AbstractEpiModel
    beta::Union{F,Vector{F}}##vector is the 
    gamma::Union{F,Vector{F}}
    #stateType::DataType
end
=#

struct SIRModel{BT<:Union{Vector{<:AbstractFloat}, AbstractFloat},
    GT<:Union{Vector{<:AbstractFloat}, AbstractFloat}} <: AbstractEpiModel
    beta::BT##vector is the 
    gamma::GT
    #stateType::DataType
end



model_states(x::SIRModel) = (:S,:I,:R)

#spreading_state(x::SIRModel) = :I
spreading_states(x::SIRModel) = Dict(:I=>[(:S,:I, x.beta)])

trans_independent(x::SIRModel) = [(:I,:R, x.gamma)]
first_active_states(x::SIRModel) = (:I,)
AF=AbstractFloat
draw_delay_mk(p::AF, rng::AbstractRNG) = rand(rng, Geometric(p))+1
draw_delay_mk(p::AF, rng::AbstractRNG, i::Integer) = rand(rng, Geometric(p))+1
draw_delay_mk( p::AF, rng::AbstractRNG, nodes::Vector{I}) where I<:Integer= rand(rng, Geometric(p), length(nodes)) .+1

draw_delay_i(m::SIRModel, p::Real, rng::AbstractRNG, i::Integer) = draw_delay_mk(p,rng,i) #rand(rng, Geometric(p))+1
draw_delay_i(m::SIRModel, p::Vector{F}, rng::AbstractRNG, i::Integer) where F<:AbstractFloat = draw_delay_mk(p[i], rng)  #rand(rng, Geometric(p[i]))+1
draw_delays_nodes(m::SIRModel, p::Real, rng::AbstractRNG, nodes::Vector{I}) where I<:Integer = rand(rng, Geometric(p), length(nodes)) .+1
draw_delays_nodes(m::SIRModel, p::Vector{<:AbstractFloat}, rng::AbstractRNG, nodes::Vector{<:Integer}) = draw_delay_mk.(p,rng, nodes)

