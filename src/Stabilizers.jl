module Stabilizers

using Random, Combinatorics, LinearAlgebra
using ..Binary, ..Orthogonals, ..Symplectics, ..Utils

export canon_stab_matrix, rand_stab_matrix, entangled_stab_matrix, weight_pauli, weight_majorana, entropy, mutual_info, stab_subgroup, code_distance

function canon_stab_matrix(n; odd=false)
    @assert iseven(n) && n > 0
    k = n ÷ 2
    !odd && return kron(Matrix{GF2}(I, k, k), [1 1])
    return triu!(ones(GF2, n, n), 1)[1:2:(n-1), :]
end

function rand_stab_matrix(n; orthogonal=false)
    @assert iseven(n) && n > 0
    A = orthogonal ? Orthogonals.rand(n) : Symplectics.rand_majorana(n)
    return canon_stab_matrix(n) * A'
end

function entangled_stab_matrix(n_party)
    @assert iseven(n_party) && n_party > 0
    return kron([1 1], Matrix{GF2}(I, n_party, n_party))
end

function weight_pauli(stab::AbstractMatrix{GF2})
    k, n = size(stab)
    @assert iseven(n) && n > 0
    weights = Vector{Int}(undef, k)
    for i in 1:2:(n-1)
        weights += isone.(stab[:, i]).||(isone.(stab[:, i+1]))
    end
    return weights
end

function weight_pauli(string::AbstractVector{GF2})
    return weight_pauli(string')[1]
end

function weight_majorana(stab::AbstractMatrix{GF2})
    k, n = size(stab)
    @assert iseven(n) && n > 0
    return reduce(+, isone.(stab); dims = 2, init = 0)
end

function weight_majorana(string::AbstractVector{GF2})
    return weight_majorana(string')[1]
end

function entropy(stab::AbstractMatrix{GF2})
    k, n = size(stab)
    @assert iseven(n)
    n_half = n ÷ 2
    @assert k ≤ n_half
    return n_half - k
end

function mutual_info(stab::AbstractMatrix{GF2}, ids_A, ids_B)
    ids_AB = [ids_A; ids_B]
    stab_A = stab_subgroup(stab, ids_A)
    stab_B = stab_subgroup(stab, ids_B)
    stab_AB = stab_subgroup(stab, ids_AB)
    return entropy(stab_A) + entropy(stab_B) - entropy(stab_AB)
end

function stab_subgroup(stab::AbstractMatrix{GF2}, ids::AbstractVector{<:Integer})
    k, n = size(stab)
    @assert iseven(n)
    n_half = n ÷ 2
    @assert k ≤ n_half
    m = length(ids)
    @assert iseven(m)

    stab_perm = permute_to_back(stab, ids)
    rref!(stab_perm)
    sub_elems = [iszero(stab_perm[i, 1:n-m]) for i=1:k]
    return stab_perm[sub_elems, (n-m+1):n]
end

function code_distance(code_cliff, k; mc_max= 0)
    @assert iseven(k) && k > 0
    n, _  = size(code_cliff)
    @assert iseven(n) && n > 0

    init_stab = canon_stab_matrix(n-k) ⊕ entangled_stab_matrix(k)
    stab = init_stab * (code_cliff' ⊕ Matrix{GF2}(I, k, k))
    R = (n+1):(n+k)

    d_max = n - k + 2

    for d in (d_max-2):-2:2
        non_zero_info = false
        if mc_max <= 0 || mc_max >= try binomial(n, d) catch; binomial(big(n), k) end
            subsets = Combinatorics.combinations(1:n, d)
            for A in subsets
                info_AR = mutual_info(stab, A, R)
                if info_AR > 0
                    non_zero_info = true
                    break  
                end
            end
        else
            for _ in 1:mc_max
                A = randperm(n)[1:d]
                info_AR = mutual_info(stab, A, R)
                if info_AR > 0
                    non_zero_info = true
                    break  
                end
            end
        end
        !non_zero_info && return d + 2
    end
    return 2
end

end
