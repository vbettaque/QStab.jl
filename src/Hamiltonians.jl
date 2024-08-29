module Hamiltonians

using Random, Distributions
using ..Galois, ..Hilbert

function tensor_qutrits(n)
    @assert n > 0
    zs = zeros(ComplexF64, (n, 3^n, 3^n))
    xs = zeros(ComplexF64, (n, 3^n, 3^n))
    for i=1:n
        z_string = zeros(GF3, 2*n)
        z_string[2*i-1] = 1
        x_string = zeros(GF3, 2*n)
        x_string[2*i] = 1
        zs[i, :, :] = pauli_string(z_string)
        xs[i, :, :] = pauli_string(x_string)
    end
    return zs, xs
end

function qutrit_sherrington_kirkpatrick(N, J, g, λ=0)
    @assert N > 1
    zs, xs = tensor_qutrits(N)
    gauss = Normal(0, sqrt(J/N))
    H = zeros(ComplexF64, 3^N, 3^N)
    for i=1:(N-1)
        for j=(i+1):N
            H -= rand(gauss) * (zs[i, :, :]' * zs[j, :, :] + zs[j, :, :]' * zs[i, :, :])
        end
    end
    for i=1:N
        H -= g * (xs[i, :, :] + xs[i, :, :]') + λ * (zs[i, :, :] + zs[i, :, :]')
    end
    return H
end

end
