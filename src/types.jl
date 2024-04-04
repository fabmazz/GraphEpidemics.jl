abstract type AbstractEpiData end

struct SIREpiData{F<:AbstractFloat,I<:Integer} <: AbstractEpiData
    inf_times::Vector{F}
    infectors::Vector{F}
    rec_delays::Vector{I}
end