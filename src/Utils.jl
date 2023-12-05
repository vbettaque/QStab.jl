module Utils

using Base: swaprows!

export directsum, ⊕, setcomplement, permute_to_back!, permute_to_back, rref!, rref

function directsum(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
    return [A zeros(T, size(A, 1), size(B, 2)); zeros(T, size(B, 1), size(A, 2)) B]
end

⊕(A, B) = directsum(A, B)

function setcomplement(U, A)
    A_compl = Vector{eltype(U)}(undef, length(U) - length(A))

    i = 1
    @inbounds for el in U
        el ∈ A && continue
        A_compl[i] = el
        i += 1
    end

    return A_compl
end

function permute_to_back(M::AbstractMatrix, ids::AbstractVector{<:Integer})
    _, n = size(M)

    @assert allunique(ids)
    @assert all(i -> (1 ≤ i ≤ n), ids)

    ids_compl = setcomplement(collect(1:n), ids)

    perm = [ids_compl; ids]

    return M[:, perm]
end

function permute_to_back!(M::AbstractMatrix, ids::AbstractVector{<:Integer})
    M = permute_to_front(M, ids)
    return M
end

function rref!(M::AbstractMatrix)
    lead = 1
    n_rows, n_columns = size(M)
	for r in 1:n_rows
        lead > n_columns && return M

        i = findnext(!iszero, M[:, lead], r)
        while isnothing(i)
            lead += 1
            lead > n_columns && return M
            i = findnext(!iszero, M[:, lead], r)
        end

        i != r && swaprows!(M, i, r)

        M[r, :] /= M[r, lead]

        js = 1:n_rows .!= r
        M[js, :] -= M[js, lead] * M[r, :]'

        lead += 1
	end
    return M
end

function rref(M::AbstractMatrix)
    return rref!(copy(M))
end

end
