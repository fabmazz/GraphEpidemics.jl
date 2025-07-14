#=abstract type AbstractEpiData end

struct SIREpiData{F<:AbstractFloat,I<:Integer} <: AbstractEpiData
    inf_times::Vector{F}
    infectors::Vector{F}
    rec_delays::Vector{I}
end

convert_matrix(x::Vector{<:Tuple}) = reduce(hcat, getindex.(x,i) for i in eachindex(x[1]))

abstract type InfectionDirection end
=#

abstract type AbstractStatesCounter end
abstract type AbstractStateChanger end

struct BaseSIRStatesCounter <:AbstractStatesCounter
end

function count_states(counter::BaseSIRStatesCounter,states::Vector)
    c=SVector{3,Int}(
        sum(states.==i) for i=1:3
    )
end