module Hilbert

using Random
using LinearAlgebra
using Combinatorics

using ..Galois
using ..Symplectics
using ..Orthogonals

export pauli_x, pauli_y, pauli_z, qutrit_x, qutrit_z, α, ⊗, inout_dims, pauli_string

const pauli_x = [0 1; 1 0]
const pauli_y = [0 -im; im 0]
const pauli_z = [1 0; 0 -1]

const α = (-1/2 + sqrt(3)/2 * im)

const qutrit_x = [0 0 1; 1 0 0; 0 1 0]
const qutrit_z = [1 0 0; 0 α 0; 0 0 α^2]

⊗ = kron

function inout_dims(n)
    @assert iseven(n) && n > 0
    dims = reshape(reshape(1:n, (n ÷ 2, 2))', n)
    return collect(dims)
end

function pauli_string(v::AbstractVector{GF2})
    n = length(v)
    @assert iseven(n) && n > 0

    string = 1
    for i=1:2:(n - 1)
        pauli = if isone(v[i] * v[i+1])
            pauli_y
        elseif isone(v[i])
            pauli_z
        elseif isone(v[i+1])
            pauli_x
        else
            [1 0; 0 1]
        end
        string = string ⊗ pauli
    end

    return string
end

function pauli_string(v::AbstractVector{GF3})
    n = length(v)
    @assert iseven(n) && n > 0

    string = 1
    power = zero(GF3)
    for i=1:2:(n - 1)
        pauli = qutrit_z^Integer(v[i]) * qutrit_x^Integer(v[i+1])
        string = string ⊗ pauli
        power -= v[i] * v[i+1]
    end
    coeff = α^Integer(power / GF3(2))
    return coeff * string
end

function majorana_string(v::AbstractVector{GF2})
    n = length(v)
    @assert iseven(n) && n > 0

    v_pauli = Symplectics.jordan_wigner_matrix(n) * v
    return pauli_string(v_pauli)
end

function braid(v::AbstractVector{GF2})
    return (I + im * majorana_string(v)) / sqrt(2)
end

function indexed_pclifford(n::Integer, i_ortho::Integer, i_majo::Integer = 1)
    @assert iseven(n) && n > 0

    C = I
    hs = Orthogonals.indexed_element_householders(n, i_ortho)
    for h in hs
        C = braid(h) * C
    end
    C = C * majorana_string(indexed_even_bitvec(i_majo, n))

    return C
end

function rand_pclifford(n::Integer, rng = Random.default_rng())
    ortho_order = Orthogonals.group_order(n)
    majo_order = (big"2")^(n - 1)
    i_ortho = Random.rand(rng, 1:ortho_order)
    i_majo = Random.rand(rng, 1:majo_order)
    return indexed_pclifford(n, i_ortho, i_majo)
end

end
