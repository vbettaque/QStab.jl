module Magic

using LinearAlgebra
using Memoize
using ..Galois, ..Symplectics, ..Hilbert

@memoize function phase_space_point(a::AbstractVector{GF3})
    n = length(a)
    @assert iseven(n) && n > 0
    D = 3^(n ÷ 2)
    op = zeros(ComplexF64, D, D)
    for i=1:(D^2)
        b = tritvec(i-1, n)
        op += qutrit_phase(a ∧ b) * pauli_string(b)'
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

function mana(rho::AbstractMatrix{ComplexF64})
    D, D_ = size(rho)
    n = 2 * log(3, D)
    @assert D == D_ && isinteger(n)
    n = Integer(n)
    mana = 0;
    for i=1:(D^2)
        print("\e[2K")
        print("\e[1G")
        print(i/D^2 * 100)
        a = tritvec(i-1, n)
        wig = wigner(rho, a)
        mana += abs(wig)
    end
    return log(3, mana[])
end

end
