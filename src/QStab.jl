module QStab

include("Utils.jl")
include("Binary.jl")
include("Cliffords.jl")
include("Orthogonals.jl")
include("Symplectics.jl")

using .Cliffords

using .Binary

using .Orthogonals
using .Symplectics

using .Utils


function main()
    v1 = bitvec(3, 6)
    v2 = bitvec(2, 6)
    display(v1)
    display(v2)
    display(v1 ∧ v2)
end

end
