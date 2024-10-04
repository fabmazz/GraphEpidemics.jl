get_type(r::Number) = typeof(r)
get_type(v::Vector{<:Number}) = eltype(v)

function count_SIR_states(states)
    StatsBase.counts(states,3)
end

function count_SIR_svector(states::Vector)
    c=SVector{3,Int}(
        sum(states.==i) for i=1:3
    )
end