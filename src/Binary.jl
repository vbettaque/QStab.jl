module Binary

using GaloisFields

export GF2, parity, complement, complement!, bitvec, indexed_odd_bitvec

const GF2 = @GaloisField 2

Base.conj(x::GF2) = x

parity(vec::AbstractVector{GF2}) = reduce(+, vec; init = GF2(0))

complement(arr::AbstractArray{GF2}) = arr .+ GF2(1)

function complement!(arr::AbstractArray{GF2})
    return copyto!(arr, complement(arr))
end

function bitvec(n::Integer, len::Integer)
    @assert n >= 0
    @assert len > 0

    vec = zeros(Integer, len)
    digits!(vec, n, base=2)
    return GF2.(vec)
end

function indexed_odd_bitvec(i::Integer, len::Integer)
    @assert len > 1
    @assert 1 <= i <= (big"2")^(len - 1)

    n = i - 1
    vec = zeros(GF2, len)
    vec[1] = GF2(count_ones(n) + 1)
    vec[2:len] = bitvec(n, len - 1)
    return vec
end

end
