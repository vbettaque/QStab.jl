module FramePotential

using LinearAlgebra

using ..Binary, ..Symplectics, ..Orthogonals, ..Hilbert

function unitary(t::Integer, d::Integer)
    @assert t > 0

    d == 2 && return factorial(2 * t) / (factorial(t) * factorial(t + 1))

    d >= t && return factorial(t)

    return -1
end

function clifford_symp(t::Integer, n::Integer; max_reps=-1)
    @assert t > 0
    @assert n > 0 && iseven(n)

    symp_order = Symplectics.group_order(n)

    full_sample = max_reps <= 0 || max_reps >= symp_order

    reps = full_sample ? symp_order : max_reps

    mean = 0
    var = 0

    for k=1:reps
        symp = full_sample ? Symplectics.indexed_element_majorana(n, k) : Symplectics.rand_majorana(n)
        eigenvectors = n - Binary.rank(symp - I)
        fixed_points = 2^eigenvectors
        term = fixed_points^(t-1)

        if full_sample
            mean += term
        else
            new_mean = mean + (term - mean) / k
            var += (term - mean) * (term - new_mean)
            mean = new_mean
        end
    end

    if full_sample
        mean = mean / reps
        err = 0
    else
        var /= reps
        err = sqrt(var/reps)
    end

    return mean, err
end

function pclifford_hilbert(t::Integer, n::Integer; max_reps=-1)
    @assert t > 0
    @assert n > 0 && iseven(n)

    proj_even = (I + Hilbert.majorana_string(ones(GF2, n))) / 2

    mean = 0
    var = 0

    ortho_order = Orthogonals.group_order(n)
    majo_order = (big"2")^(n - 1)

    if (max_reps <= 0 || max_reps >= (majo_order * ortho_order))
        for i=1:ortho_order
            mean_majo = 0
            for j=1:majo_order
                U = Hilbert.indexed_pclifford(n, i, j)
                trace = round(abs(tr(proj_even * U))^(2 * t))
                mean_majo += (trace - mean_majo) / j
            end
            mean += (mean_majo - mean) / i
        end
    else
        for i=1:max_reps
            U = Hilbert.rand_pclifford(n)
            trace = round(abs(tr(proj_even * U))^(2 * t))

            new_mean = mean + (trace - mean) / i
            var += (trace - mean) * (trace - new_mean)
            mean = new_mean
        end
    end

    var /= max_reps
    err = sqrt(var/max_reps)

    return mean, err
end

# function orthogonal(t::Integer, n::Integer; max_reps=-1)
#     @assert t > 0
#     @assert n > 0 && iseven(n)

#     reps = max_reps > 0 ? max_reps : Orthogonals.group_order(n)
#     mean = 0

#     for k=1:reps
#         ortho = max_reps > 0 ? Orthogonals.rand(n) : Orthogonals.indexed_element(n, k)
#         eigenvectors = n - Binary.rank(ortho - I)
#         fixed_points = 2^eigenvectors
#         mean += (fixed_points^(t-1) - mean) / k
#     end

#     return mean
# end

# function even_parity_sector_matrix(ortho::AbstractMatrix{GF2})
#     n, ncols = size(ortho)
#     @assert ncols == n
#     @assert ortho' * ortho == I

#     even_basis = zeros(GF2, n, n-1)
#     even_basis[1:(n-1), :] += I
#     even_basis[2:n, :] += I

#     compl_basis = zeros(GF2, n-1, n)
#     compl_basis[:, 1:(n-1)] = LowerTriangular(ones(GF2, n-1, n-1))

#    return compl_basis * ortho * even_basis
# end

# function orthogonal_even(t::Integer, n::Integer; with_pcheck=false, max_reps=-1)
#     @assert t > 0
#     @assert n > 0 && iseven(n)

#     reps = max_reps > 0 ? max_reps : Orthogonals.group_order(n)
#     mean = 0

#     for k=1:reps
#         ortho = max_reps > 0 ? Orthogonals.rand(n) : Orthogonals.indexed_element(n, k)
#         reduced = even_parity_sector_matrix(ortho)
#         eigenvectors = (n - 1) - Binary.rank(reduced - I)
#         fixed_points = 2^eigenvectors - (1 - with_pcheck)
#         mean += (fixed_points^(t-1) - mean) / k
#     end

#     return mean
# end

end
