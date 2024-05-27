abstract type AbstractEpiData end

struct SIREpiData{F<:AbstractFloat,I<:Integer} <: AbstractEpiData
    inf_times::Vector{F}
    infectors::Vector{F}
    rec_delays::Vector{I}
end

convert_matrix(x::Vector{<:Tuple}) = reduce(hcat, getindex.(x,i) for i in eachindex(x[1]))

abstract type InfectionDirection end

struct InfectI <: InfectionDirection end
struct InfectS <: InfectionDirection end
struct InfectMeanIS <: InfectionDirection end