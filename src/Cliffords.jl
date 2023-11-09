module Cliffords

using Base
using ..Binary
using LinearAlgebra
using ..Utils

export Clifford, Pauli, Majorana, I, tensor, ⊗

IntOrGF2 = Union{Integer, GF2}

@enum CliffordType Pauli Majorana

struct Clifford{T,N}
    _M::Matrix{GF2}
    function Clifford{T,N}(M::AbstractMatrix{GF2}) where {T,N}
        @assert T isa CliffordType
        @assert N isa Integer
        @assert iseven(N)
        nrows, ncols = size(M)
        @assert nrows == N == ncols
        @assert transpose(M) * M == I # TODO: Check Symplecity
        return new(M)
    end
end

function Clifford{T,N}(U::UniformScaling{<:IntOrGF2}) where {T,N}
    return Clifford{T,N}(GF2.(U(N)))
end

function Clifford(M::AbstractMatrix{<:IntOrGF2})
    N, _ = size(M)
    return Clifford{N}(GF2.(M))
end


Base.convert(::Type{Clifford}, M::AbstractMatrix{<:IntOrGF2}) = Clifford(M)
Base.convert(::Type{Clifford{N}}, U::UniformScaling{<:IntOrGF2}) where N = Clifford{N}(U)

Base.:*(c1::Clifford{N}, c2::Clifford{N}) where N = Clifford{N}(c1._M * c2._M)

tensor(c1::Clifford{N}, c2::Clifford{M}) where {N, M} = Clifford{N + M}(c1._M ⊕ c2._M)
⊗(c1::Clifford{N}, c2::Clifford{M}) where {N, M} = tensor(c1, c2)

end
