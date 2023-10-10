module Symplectics

using ..Binary
using LinearAlgebra

export group_order, symp, ∧

function group_order(n::Integer)
    @assert iseven(n) && n > 0

    m = n ÷ 2
    order = (big"2")^(m^2)
    for i in 1:m
        order *= ((big"4")^i - 1)
    end
    return order
end

function symp(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
    n = size(A, 1)
    @assert iseven(n) && n == size(B, 1)
    return view(A, 1:2:(n-1), :)' * view(B, 2:2:n, :) -
        view(A, 2:2:n, :)' * view(B, 1:2:(n-1), :)
end

function symp(v::AbstractVector{T}, w::AbstractVector{T}) where T
    return symp(reshape(v, :, 1), reshape(w, :, 1))[1]
end

∧(x::AbstractVecOrMat{T}, y::AbstractVecOrMat{T}) where T = symp(x, y)

end
