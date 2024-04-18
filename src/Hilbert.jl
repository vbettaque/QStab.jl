module Hilbert

using Random
using LinearAlgebra

using ..Binary
using ..Symplectics
using ..Orthogonals

const pauli_x = [0 1; 1 0]
const pauli_y = [0 -im; im 0]
const pauli_z = [1 0; 0 -1]

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
        string = kron(string, pauli)
    end

    return string
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