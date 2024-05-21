UnionFloat = Union{Vector{<:AbstractFloat}, AbstractFloat}

struct SISModel{BT<:UnionFloat, GT<:UnionFloat} <: AbstractEpiModel
    beta::BT
    gamma::GT
end

model_states(x::SISModel) = (:S,:I)
spreading_states(x::SISModel) = Dict(:I=>[(:S,:I, x.beta)])
trans_independent(x::SISModel) = [(:I,:S, x.gamma)]

first_active_states(m::SISModel) = (:I,)
draw_delays(m::SISModel, p, rng::AbstractRNG, node) = draw_delays_markov(p, rng, node )

struct SIRSModel{BT<:UnionFloat,
                GT<:UnionFloat, TT<:UnionFloat} <: AbstractEpiModel

    beta::BT
    gamma::GT
    delta::TT
end

model_states(x::SIRSModel) = (:S,:I,:R)
spreading_states(x::SIRSModel) = Dict(:I=>[(:S,:I, x.beta)])
trans_independent(x::SIRSModel) = [(:I,:R, x.gamma), (:R,:S,x.delta)]

first_active_states(m::SIRSModel) = (:I,)
draw_delays(m::SIRSModel, p, rng::AbstractRNG, node) = draw_delays_markov(p, rng, node )