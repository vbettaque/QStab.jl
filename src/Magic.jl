module Magic

using LinearAlgebra
using TensorOperations

using ..Galois, ..Symplectics, ..Hilbert

const phase_space_point_ops = [
    [1   0   0  ; 0   0   α  ; 0   α^2 0],
    [1   0   0  ; 0   0   α^2; 0   α   0],
    [0   1   0  ; 1   0   0  ; 0   0   1],
    [0   α   0  ; α^2 0   0  ; 0   0   1],
    [0   α^2 0  ; α   0   0  ; 0   0   1],
    [0   0   1  ; 0   1   0  ; 1   0   0],
    [0   0   α^2; 0   1   0  ; α   0   0],
    [0   0   α  ; 0   1   0  ; α^2 0   0],
    [1   0   0  ; 0   0   1  ; 0   1   0],
]

function phase_space_point(a::AbstractVector{GF3})
    n = length(a)
    @assert iseven(n) && n > 0
    D = 3^(n ÷ 2)
    op = zeros(ComplexF64, D, D)
    for i=1:(D^2)
        b = tritvec(i-1, n)
        op += α^Integer(a ∧ b) * pauli_string(b)'
    end
    return op' / D
end

function wigner(rho::AbstractMatrix{ComplexF64}, a::AbstractVector{GF3})
    D, D_ = size(rho)
    n = length(a)
    @assert iseven(n) && n > 0
    @assert D == D_ == 3^(n ÷ 2)
    psp = phase_space_point(a)
    wig = reduce(+, [dot(psp[i, :], rho[:, i]) for i=1:D]) / D
    @assert isreal(wig) || isapprox(imag(wig), 0, atol=1e-10)
    return real(wig)
end

function mana_old(rho::AbstractMatrix{ComplexF64})
    D, D_ = size(rho)
    n = 2 * log(3, D)
    @assert D == D_ && isinteger(n)
    n = Integer(n)
    mana = 0;
    for i=1:(D^2)
        a = tritvec(i-1, n)
        wig = wigner(rho, a)
        mana += abs(wig)
    end
    return log(3, mana[])
end

function mana(ρ::AbstractMatrix{ComplexF64})
    D, D_ = size(ρ)
    n = 2 * log(3, D)
    @assert D == D_ && isapprox(n, round(n))
    n = round(Int, n)

    ρ = reshape(ρ, repeat([3], n)...)
    ρ = permutedims(ρ, inout_dims(n))

    A = [reshape(phase_space_point_ops[i], 9) for i=1:9]
    A = reduce(vcat,transpose.(A))

    for i=1:(n÷2)
        ρ = reshape(ρ, (cat([9^(i-1), 9, 9^(n÷2-i)], dims=1)...))
        @tensor ρp[α,β,γ] := A[β,βp] * ρ[α, βp, γ]
        ρ = ρp
    end

    return log(3, ρ / D .|> abs |> sum)
end

end
