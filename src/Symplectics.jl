module Symplectics

using Random
using LinearAlgebra
using ..Binary

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

function jordan_wigner_matrix(n::Integer)
    @assert iseven(n) && n > 0
    return triu!(complement!(product_form_pauli(n)))
end

function inner_pauli(A::AbstractVecOrMat{GF2}, B::AbstractVecOrMat{GF2})
    n = size(A, 1)
    @assert iseven(n) && n == size(B, 1)
    return view(A, 1:2:(n-1), :)' * view(B, 2:2:n, :) -
        view(A, 2:2:n, :)' * view(B, 1:2:(n-1), :)
end

∧(x::AbstractVecOrMat{GF2}, y::AbstractVecOrMat{GF2}) = inner_pauli(x, y)

function inner_majorana(A::AbstractVecOrMat{GF2}, B::AbstractVecOrMat{GF2})
    n = size(A, 1)
    @assert iseven(n) && n == size(B, 1)
    return transpose(A) * parity_complement(B)
end

⊼(x::AbstractVecOrMat{GF2}, y::AbstractVecOrMat{GF2}) = inner_majorana(x, y)

function is_pauli(A::AbstractMatrix{GF2})
    n = size(A, 1)
    return A ∧ A == product_form_pauli(n)
end

function is_majorana(A::AbstractMatrix{GF2})
    n = size(A, 1)
    return A ⊼ A == product_form_majorana(n)
end

function _transvect!(v::AbstractVector{GF2}, h::AbstractVector{GF2})
    return copyto!(v, v + (v ∧ h)[1, 1] * h)
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

    if v == w
        return _transvection(v - w)
    end

    n = length(v)
    u = zeros(GF2, n)


    for i in 1:2:(n-1)
        if !iszero(v[i:i+1]) && !iszero(w[i:i+1])
            u[i] = v[i] + w[i]
            u[i+1] = v[i+1] + w[i+1]
            if iszero(u[i:i+1])
                u[i+1] = one(GF2)
                (v[i] != v[i+1]) && (u[i] = one(GF2))
            end
            h1 = u - v; h2 = w - u
            return _transvection(h2) ∘ _transvection(h1)
        end
    end

    for i in 1:2:(n-1)
        if !iszero(v[i:i+1]) && iszero(w[i:i+1])
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
        if iszero(v[i:i+1]) && !iszero(w[i:i+1])
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

function indexed_element_pauli!(M::AbstractMatrix{GF2}, i::Integer)
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
        indexed_element_pauli!(view(M, 3:n, 3:n), i_rec)
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

function indexed_element_pauli(n::Integer, i::Integer)
    @assert iseven(n) && n > 0
    return indexed_element_pauli!(Matrix{GF2}(I, n, n), i)
end

function indexed_element_majorana!(M::AbstractMatrix{GF2}, i::Integer)
    indexed_element_pauli!(M, i)
    n, _ = size(M)
    W = jordan_wigner_matrix(n)
    M = W * M * W
    return M
end

function indexed_element_majorana(n::Integer, i::Integer)
    @assert iseven(n) && n > 0
    return indexed_element_majorana!(Matrix{GF2}(I, n, n), i)
end

function rand_pauli(n::Integer, rng = Random.default_rng())
    i = Random.rand(rng, 1:group_order(n))
    return indexed_element_pauli(n, i)
end

function rand_majorana(n::Integer, rng = Random.default_rng())
    i = Random.rand(rng, 1:group_order(n))
    return indexed_element_majorana(n, i)
end

end
