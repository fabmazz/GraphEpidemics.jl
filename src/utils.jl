get_type(r::Number) = typeof(r)
get_type(v::Vector{<:Number}) = eltype(v)