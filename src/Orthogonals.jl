module Orthogonals

using ..Binary
using LinearAlgebra
using Random

function group_order(n::Integer)
    @assert n > 0

    m = n ÷ 2
    order = (big"2")^(m^2)
    for i in 1:(m-1)
        order *= ((big"4")^i - 1)
    end
    if isodd(n) && n > 1
        order *= ((big"4")^m - 1)
    end
    return order
end

function _householder!(M::AbstractMatrix{GF2}, h::AbstractVector{GF2})
    return copyto!(M, M + h * (transpose(h) * M))
end

# function _transvection(h::AbstractVector{GF2})
#     return x -> _transvect!(x, h)
# end

# function _transvection(v::AbstractVector{GF2}, w::AbstractVector{GF2})
#     @assert length(v) == length(w)
#     @assert isone(parity(v)) && isone(parity(w))

#     if iszero(v ⋅ w) || v == w
#         h = w - v
#         return _transvection(h)
#     end

#     n = length(v)
#     if iseven(n)
#         h = complement(w - v)
#         return complement! ∘ _transvection(h)
#     end

#     h = ones(GF2, n)
#     for i=1:n
#         if iszero(v[i]) && iszero(w[i])
#             h[i] = 0
#             return _transvection(w - v - h) ∘ _transvection(h)
#         end
#     end

#     h = zeros(GF2, n)
#     v_minus = v - v .* w
#     w_minus = w - v .* w
#     h[findfirst(isone, v_minus)] = 1
#     h[findfirst(isone, w_minus)] = 1

#     return _transvection(w - v - h) ∘ _transvection(h)
# end

function _householder_vector(v::AbstractVector{GF2}, w::AbstractVector{GF2})
    @assert length(v) == length(w)
    @assert isone(parity(v)) && isone(parity(w))

    if iszero(v ⋅ w) || v == w
        return (w - v, nothing)
    end

    n = length(v)
    if iseven(n)
        h = complement(w - v)
        return (h, ones(GF2, n))
    end

    h = ones(GF2, n)
    for i=1:n
        if iszero(v[i]) && iszero(w[i])
            h[i] = 0
            return (h, w - v - h)
        end
    end

    h = zeros(GF2, n)
    v_minus = v - v .* w
    w_minus = w - v .* w
    h[findfirst(isone, v_minus)] = 1
    h[findfirst(isone, w_minus)] = 1

    return (h, w - v - h)
end

# function indexed_element!(M::AbstractMatrix{GF2}, i::Integer)
#     n, ncols = size(M)
#     @assert ncols == n
#     @assert iseven(n)

#     order = group_order(n)
#     @assert 1 <= i <= order

#     M == I || copyto!(M, I)

#     p1 = (big"2")^(n - 1) # Number of odd-parity vectors v1 of length n
#     p2 = (big"2")^(n - 2) - 1 # Number of odd-parity vectors v2 of length n orthogonal to v1

#     if order == 2
#         isone(i) || complement!(M)
#         return M
#     else
#         i_rec = (i - 1) ÷ (p1 * p2) + 1
#         indexed_element!(view(M, 3:n, 3:n), i_rec)
#     end

#     i1 = (i - 1) % p1 + 1
#     i2 = ((i - 1) ÷ p1) % p2 + 1

#     e1 = view(M, :, 1)
#     e2 = view(M, :, 2)

#     f1 = indexed_odd_bitvec(i1, n)
#     f2 = vcat(GF2(0), indexed_odd_bitvec(i2, n - 1))

#     t1! = _transvection(e1, f1)
#     t2! = if iszero(f2[2])
#         h = e2 + f2
#         _transvection(h)
#     else
#         h1 = copy(e2); h2 = copy(f2)
#         j = findfirst(iszero, f2[3:end]) + 2
#         h1[j] = GF2(1); h2[j] = GF2(1)
#         _transvection(h2) ∘ _transvection(h1)
#     end

#     return t1!(t2!(M))
# end

# function indexed_element!(M::AbstractMatrix{GF2}, i::Integer)
#     n, ncols = size(M)
#     @assert ncols == n
#     @assert n > 0

#     order = group_order(n)
#     @assert 1 <= i <= order

#     M == I || copyto!(M, I)

#     p = (big"2")^(n - 1) - isodd(n) # Number of odd-parity vectors of length n

#     if order == 1
#         return M
#     else
#         i_rec = (i - 1) ÷ p + 1
#         indexed_element!(view(M, 2:n, 2:n), i_rec)
#     end

#     i_p = (i - 1) % p + 1

#     e = view(M, :, 1)

#     f = indexed_odd_bitvec(i_p, n)

#     t! = _transvection(e, f)

#     return t!(M)
# end

function indexed_element_householders(n::Integer, i::Integer)
    @assert n > 0

    hs = []
    i_k = i
    for k=2:n
        p = (big"2")^(k - 1) - isodd(k)
        i_p = (i_k - 1) % p + 1
        e = zeros(GF2, k); e[1] = 1
        f = indexed_odd_bitvec(i_p, k)
        h1, h2 = _householder_vector(e, f)
        push!(hs, [zeros(GF2, n-k); h1])
        if !isnothing(h2)
            push!(hs, [zeros(GF2, n-k); h2])
        end
        i_k = (i_k - 1) ÷ p + 1
    end
    return hs
end

function indexed_element!(M::AbstractMatrix{GF2}, i::Integer)
    n, ncols = size(M)
    @assert ncols == n
    @assert n > 0

    order = group_order(n)
    @assert 1 <= i <= order

    M == I || copyto!(M, I)

    hs = indexed_element_householders(n, i)
    for h in hs
        _householder!(M, h)
    end
    return M
end

function indexed_element(n::Integer, i::Integer)
    @assert n > 0
    return indexed_element!(Matrix{GF2}(I, n, n), i)
end

function rand(n::Integer, rng = Random.default_rng())
    i = Random.rand(rng, 1:group_order(n))
    return indexed_element(n, i)
end

end