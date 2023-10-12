module Symplectics

using ..Binary
using LinearAlgebra

function group_order(n::Integer)
    @assert iseven(n) && n > 0

    m = n ÷ 2
    order = (big"2")^(m^2)
    for i in 1:m
        order *= ((big"4")^i - 1)
    end
    return order
end

function product_form_pauli(n::Integer)
    @assert iseven(n) && n > 0
    M = Matrix{GF2}(I, n, n)
    for i in 1:2:(n-1)
        M[i, :], M[i+1, :] = M[i+1, :], M[i, :]
    end
    return M
end

function product_form_majorana(n::Integer)
    @assert iseven(n) && n > 0
    M = ones(GF2, n, n) - I
    return M
end

function inner_pauli(A::AbstractVecOrMat{T}, B::AbstractVecOrMat{T}) where T
    n = size(A, 1)
    @assert iseven(n) && n == size(B, 1)
    return view(A, 1:2:(n-1), :)' * view(B, 2:2:n, :) -
        view(A, 2:2:n, :)' * view(B, 1:2:(n-1), :)
end

∧(x::AbstractVecOrMat{T}, y::AbstractVecOrMat{T}) where T = inner_pauli(x, y)

function inner_majorana(A::AbstractVecOrMat{GF2}, B::AbstractVecOrMat{GF2})
    n = size(A, 1)
    @assert iseven(n) && n == size(B, 1)
    return transpose(A) * mapslices(parity_complement, B, dims=1)
end


function _transvect!(M::AbstractMatrix{GF2}, h::AbstractVector{GF2})
    return copyto!(M, M + h * transpose(M ∧ h))
end

function _transvection(h::AbstractVector{GF2})
    return x -> _transvect!(x, h)
end

function _transvection(v::AbstractVector{GF2}, w::AbstractVector{GF2})
    if !iszero(v ∧ w)
        h = w - v
        return _transvection(h)
    end

    n = length(v)
    u = zeros(GF2, n)

    for i in 1:2:(n-1)
        if !iszero(v[i] + v[i+1]) && !iszero(w[i] + w[i+1])
            u[i] = v[i] + w[i]
            u[i+1] = v[i+1] + w[i+1]
            if iszero(u[i] + u[i+1])
                u[i+1] = one(GF2)
                (v[i] != v[i+1]) && (u[i] = one(GF2))
            end
            h1 = u - v; h2 = w - u
            return _transvection(h2) ∘ _transvection(h1)
        end
    end

    for i in 1:2:(n-1)
        if !iszero(v[i] + v[i+1]) && iszero(w[i] + w[i+1])
            if v[i] == v[i+1]
                u[i+1] = one(GF2)
            else
                u[i+1] = v[i]
                u[i] = v[i+1]
            end
            break
        end
    end

    for i in 1:2:(n-1)
        if iszero(v[i] + v[i+1]) && !iszero(w[i] + w[i+1])
            if w[i] == w[i+1]
                u[i+1] = one(GF2)
            else
                u[i+1] = w[i]
                u[i] = w[i+1]
            end
            break
        end
    end

    h1 = u - v; h2 = w - u;
    return _transvection(h2) ∘ _transvection(h1)
end

function indexed_element!(M::AbstractMatrix{GF2}, i::Integer)
    n, ncols = size(M)
    @assert ncols == n
    @assert iseven(n)

    order = group_order(n)
    @assert 1 <= i <= order

    M == I || copyto!(M, I)

    p1 = (big"2")^n - 1 # Number of symplectic vectors v1 of length n
    p2 = (big"2")^(n - 1) # Number of symplectic vectors v2 of length n orthogonal to v1

    if n > 2
        i_rec = (i - 1) ÷ (p1 * p2) + 1
        indexed_element!(view(M, 3:n, 3:n), i_rec)
    end

    i1 = (i - 1) % p1 + 1
    i2 = ((i - 1) ÷ p1) % p2 + 1

    e1 = view(M, :, 1)
    e2 = view(M, :, 2)

    f1 = bitvec(i1, n)
    f2_coeffs = bitvec(i2, n - 1)
    f2 = vcat(f2_coeffs[1], one(GF2), f2_coeffs[2:(n-1)])

    t1! = _transvection(e1, f1)
    t2! = if !iszero(f2[1])
        h = f2 - e2
        _transvection(h)
    else
        h1 = vcat(one(GF2), zero(GF2), f2_coeffs[2:(n-1)])
        h2 = copy(e1)
        _transvection(h2) ∘ _transvection(h1)
    end

    return t1!(t2!(M))
end

function indexed_element(n::Integer, i::Integer)
    @assert iseven(n) && n > 0
    return indexed_element!(Matrix{GF2}(I, n, n), i)
end

end
