module Circuits

using LinearAlgebra
using Random
using ..Binary

export rand_q_local_volume

function rand_permutation(n)
    @assert iseven(n)
    P = Matrix{GF2}(I, n, n)
    r = randperm(n)
    return P[:, r]
end

function rand_q_local_volume(n, q, D, sampler)
    @assert iseven(q) && q > 0
    @assert iszero(n % q)
    @assert D > 0
    k = n ÷ q
    C = Matrix{GF2}(I, n, n)
    for d in 1:D
        C_d = zeros(GF2, n, n)
        for i in 1:k
            range = ((i-1)*q+1):(i*q)
            C_d[range, range] = sampler(q)
        end
        C = C_d * rand_permutation(n) * C
    end
    return C
end

end
