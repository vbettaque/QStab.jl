module Utils

using Base: swaprows!

export directsum, ⊕, rref

function directsum(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
    return [A zeros(T, size(A, 1), size(B, 2)); zeros(T, size(B, 1), size(A, 2)) B]
end

⊕(A, B) = directsum(A, B)

function rref(M::AbstractMatrix)
    lead = 1
    n_rows, n_columns = size(M)
    M_rref = copy(M)
	for r in 1:n_rows
        lead > n_columns && return M_rref

        i = findnext(!iszero, M_rref[:, lead], r)
        while isnothing(i)
            lead += 1
            lead > n_columns && return M_rref
            i = findnext(!iszero, M_rref[:, lead], r)
        end

        i != r && swaprows!(M_rref, i, r)

        M_rref[r, :] /= M_rref[r, lead]

        js = 1:n_rows .!= r
        M_rref[js, :] -= M_rref[js, lead] * M_rref[r, :]'

        lead += 1
	end
    return M_rref
end

end
