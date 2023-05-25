module QStab

include("Util.jl")
include("Binary.jl")
include("Cliffords.jl")
include("Orthogonals.jl")

using BenchmarkTools

using .Cliffords

using .Binary

using .Orthogonals


function main()
    @benchmark begin
    c1::Clifford{10} = I
    c2::Clifford{10} = I
    c1 * c2
    end
end

end
