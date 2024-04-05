abstract type AbstractEpiModel end

function states_values(x::AbstractEpiModel)
    d::Dict{Symbol,Int8} =  Dict(s=>i for (i,s) in enumerate(model_states(x)))
    d
end

struct SIRModel{F<:AbstractFloat} <: AbstractEpiModel
    beta::Union{F,Vector{F}}##vector is the 
    gamma::Union{F,Vector{F}}
    #stateType::DataType
end




model_states(x::SIRModel) = (:S,:I,:R)

#spreading_state(x::SIRModel) = :I
spreading_states(x::SIRModel) = Dict(:I=>[(:S,:I, x.beta)])

trans_independent(x::SIRModel) = [(:I,:R, x.gamma)]
first_active_states(x::SIRModel) = (:I,)

draw_delays(m::SIRModel, p::Real, rng::AbstractRNG, i::Integer) = rand(rng, Geometric(p))+1
draw_delays(m::SIRModel, p::Vector{F}, rng::AbstractRNG, i::Integer) where F<:AbstractFloat = rand(rng, Geometric(p[i]))+1
