#= 
 Copyright (c) 2024 Fabio Mazza
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at https://mozilla.org/MPL/2.0/.
=#
import StatsBase: counts

get_type(r::Number) = typeof(r)
get_type(v::Vector{<:Number}) = eltype(v)

function count_SIR_states(states)
    counts(states,3)
end

function count_SIR_svector(states::Vector)
    c=SVector{3,Int}(
        sum(states.==i) for i=1:3
    )
end