get_type(r::Number) = typeof(r)
get_type(v::Vector{<:Number}) = eltype(v)

function count_SIR_states(states)
    StatsBase.counts(states,3)
end