module FramePotential

using LinearAlgebra

using ..Binary, ..Orthogonals

function even_parity_sector(ortho::AbstractMatrix{GF2})
    n, ncols = size(ortho)
    @assert ncols == n
    @assert ortho' * ortho == I

    even_basis = zeros(GF2, n, n-1)
    even_basis[1:(n-1), :] += I
    even_basis[2:n, :] += I

    compl_basis = zeros(GF2, n-1, n)
    compl_basis[:, 1:(n-1)] = LowerTriangular(ones(GF2, n-1, n-1))

   return compl_basis * ortho * even_basis
end

function even_parity(t::Integer, n::Integer; with_pcheck=false, max_reps=-1)
    @assert t > 0
    @assert n > 0 && iseven(n)

    reps = max_reps > 0 ? max_reps : Orthogonals.group_order(n)
    mean = 0

    for k=1:reps
        ortho = max_reps > 0 ? Orthogonals.rand(n) : Orthogonals.indexed_element(n, k)
        reduced = even_parity_sector(ortho)
        fixed_points = 2^((n - 2) - Binary.rank(reduced - I) + with_pcheck)
        mean += (fixed_points^(t-1) - mean) / k
    end

    return mean
end

end
