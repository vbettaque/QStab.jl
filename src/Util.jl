module Util

export directsum, ⊕

function directsum(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T <: Number
    return [A zeros(T, size(A, 1), size(B, 2)); zeros(T, size(B, 1), size(A, 2)) B]
end

⊕(A, B) = directsum(A, B)

end
