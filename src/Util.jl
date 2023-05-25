module Util

export directsum, ⊕

function directsum(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T <: Number
    return [A zeros(T, size(A, 1), size(B, 2)); zeros(T, size(B, 1), size(A, 2)) B]
end

⊕(A, B) = directsum(A, B)

function rref(M)
    lead = 1
    n_rows, n_columns = size(M)
    M_rref = copy(M)
	for r in 1:n_rows
        if n_columns < lead
            return M_rref
		end
        i = r
        while iszero(M_rref[i, lead])
            i += 1
            if i > n_rows
                i = r
                lead += 1
                if lead > n_columns
                    return M_rref
				end
			end
		end
        if i != r
            temp = copy(M_rref[i, :])
            M_rref[i, :] = M_rref[r, :]
            M_rref[r, :] = temp
		end
        M_rref[r] /= M_rref[r, lead]
        for j in 1:n_rows
            if j != r
                M_rref[j, :] -= M_rref[j, lead] * M_rref[r, :]
			end
		end
        lead += 1
	end
    return M_rref
end

end
