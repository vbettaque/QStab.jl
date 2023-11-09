module Strings

using ..Binary

export weight_pauli, weight_majorana

function weight_pauli(string::AbstractVector{GF2})
    n = length(string)
    @assert iseven(n)
    weight = 0
    for i in 1:2:(n-1)
        weight += isone(string[i]) || isone(string[i+1])
    end
    return weight
end

function weight_majorana(string::AbstractVector{GF2})
    n = length(string)
    @assert iseven(n)
    return reduce(+, isone.(string); init = 0)
end

end
