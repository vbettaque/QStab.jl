module Cliffords

using Base
using ..Binary
using LinearAlgebra
using ..Util

export Clifford, I, tensor

struct Clifford{N}
    _M::Matrix{GF2}

    function Clifford{N}(M::AbstractMatrix{GF2}) where N
        @assert N isa Integer
        @assert iseven(N)
        nrows, ncols = size(M)
        @assert nrows == N == ncols
        @assert transpose(M) * M == I # TODO: Check Symplecity
        return new(M)
    end
end

function Clifford{N}(U::UniformScaling{<:Integer}) where N
    return Clifford{N}(GF2.(U(N)))
end

function Clifford(M::AbstractMatrix{<:Integer})
    N, _ = size(M)
    return Clifford{N}(GF2.(M))
end

Base.convert(::Type{Clifford}, M::AbstractMatrix{<:Integer}) = Clifford(M)
Base.convert(::Type{Clifford{N}}, U::UniformScaling{<:Integer}) where N = Clifford{N}(U)

Base.:*(c1::Clifford{N}, c2::Clifford{N}) where N = Clifford{N}(c1._M * c2._M)

tensor(c1::Clifford{N}, c2::Clifford{M}) where {N, M} = Clifford{N + M}(c1._M ⊕ c2._M)

end
