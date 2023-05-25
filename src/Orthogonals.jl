module Orthogonals

using ..Binary
using LinearAlgebra
using Random

function group_order(n::Integer)
    @assert iseven(n) && n > 0

    k = n ÷ 2
    order = (big"2")^(k * k)
    for i in 1:(k-1)
        order *= ((big"4")^i - 1)
    end
    return order
end

function _transvect!(M::AbstractMatrix{GF2}, h::AbstractVector{GF2})
    return copyto!(M, M + h * (transpose(h) * M))
end

function _transvection(v::AbstractVector{GF2}, w::AbstractVector{GF2})
    if iszero(v ⋅ w)
        h = w - v
        return x -> _transvect!(x, h)
    else
        h = complement(w - v)
        return x -> complement!(_transvect!(x, h))
    end
end

function _transvection(h::AbstractVector{GF2})
    return x -> _transvect!(x, h)
end

function indexed_element!(M::AbstractMatrix{GF2}, i::Integer)
    n, ncols = size(M)
    @assert ncols == n
    @assert iseven(n)

    order = group_order(n)
    @assert 1 <= i <= order

    M == I || copyto!(M, I)

    p1 = (big"2")^(n - 1) # Number of odd-parity vectors v1 of length n
    p2 = (big"2")^(n - 2) - 1 # Number of odd-parity vectors v2 of length n orthogonal to v1

    if order == 2
        isone(i) || complement!(M)
        return M
    else
        i_rec = (i - 1) ÷ (p1 * p2) + 1
        indexed_element!(view(M, 3:n, 3:n), i_rec)
    end

    i1 = (i - 1) % p1 + 1
    i2 = ((i - 1) ÷ p1) % p2 + 1

    e1 = view(M, :, 1)
    e2 = view(M, :, 2)

    f1 = indexed_odd_bitvec(i1, n)
    f2 = vcat(GF2(0), indexed_odd_bitvec(i2, n - 1))

    t1! = _transvection(e1, f1)
    t2! = if iszero(f2[2])
        h = e2 + f2
        _transvection(h)
    else
        h1 = copy(e2); h2 = copy(f2)
        j = findfirst(iszero, f2[3:end]) + 2
        h1[j] = GF2(1); h2[j] = GF2(1)
        _transvection(h2) ∘ _transvection(h1)
    end

    return t1!(t2!(M))
end

function indexed_element(n::Integer, i::Integer)
    @assert iseven(n) && n > 0
    return indexed_element!(Matrix{GF2}(I, n, n), i)
end

function rand(n::Integer, rng = Random.default_rng())
    i = Random.rand(rng, 1:group_order(n))
    return indexed_element(n, i)
end

end
