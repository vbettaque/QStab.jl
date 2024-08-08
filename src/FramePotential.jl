module FramePotential

using LinearAlgebra, Statistics

using ..Galois, ..Symplectics, ..Orthogonals, ..Hilbert, ..Utils

function unitary(t::Integer, d::Integer)
    @assert t > 0
    d == 2 && return factorial(2 * t) / (factorial(t) * factorial(t + 1))
    d >= t && return factorial(t)
    return -1
end

function fixed_points_symp(symp::AbstractMatrix{GF2})
    n, m = size(symp)
    @assert iseven(n) && n == m

    kernel_dim = n - Binary.rank(symp - I)
    return 2^kernel_dim
end

function fixed_points_symp(n::Integer; max_reps=-1)
    @assert n > 0 && iseven(n)

    symp_order = Symplectics.group_order(n)
    full_sample = max_reps <= 0 || max_reps >= symp_order
    reps = full_sample ? symp_order : max_reps

    fixed_points = zeros(Int, reps)

    done = Threads.Atomic{Int}(0)
    Threads.@threads for i=1:reps
        symp = full_sample ?
            Symplectics.indexed_element_majorana(n, i) : Symplectics.rand_majorana(n)
        fixed_points[i] = fixed_points_symp(symp)
        Threads.atomic_add!(done, 1)
        print(Int(round(done[]/reps * 10000))/100, "% \u001b[1000D")
    end

    return fixed_points
end

function fixed_points_ortho(ortho::AbstractMatrix{GF2})
    n, m = size(ortho)
    @assert iseven(n) && n == m

    ortho_reduced = even_parity_sector_matrix(ortho)
    j_reduced = GF2.([isodd(i) for i=1:(n-1)])

    eqs = rref(hcat(ortho_reduced - I, j_reduced))
    r = Binary.rank(eqs)
    has_complements = !iszero(eqs[r, 1:(n-1)])
    r -= !has_complements

    kernel_dim = n - r - 1

    return 2^(kernel_dim - 1 + has_complements)
end

function fixed_points_ortho(n::Integer; max_reps=-1)
    @assert n > 0 && iseven(n)

    ortho_order = Orthogonals.group_order(n)
    full_sample = max_reps <= 0 || max_reps >= ortho_order
    reps = full_sample ? ortho_order : max_reps

    fixed_points = zeros(Int, reps)

    done = Threads.Atomic{Int}(0)
    Threads.@threads for i=1:reps
        ortho = full_sample ?
            Orthogonals.indexed_element(n, i) : Orthogonals.rand(n)
        fixed_points[i] = fixed_points_ortho(ortho)
        Threads.atomic_add!(done, 1)
        print(Int(round(done[]/reps * 10000))/100, "% \u001b[1000D")
    end

    return fixed_points
end

function from_fixed_points(fixed_points, t::Integer; exact = false)
    powers = fixed_points .^ (t-1)
    pot_mean = mean(powers)
    pot_std = exact ? 0 : std(powers; mean=pot_mean)
    return pot_mean, pot_std / sqrt(length(fixed_points))
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

function even_parity_sector_matrix(ortho::AbstractMatrix{GF2})
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

end
