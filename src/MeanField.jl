module MeanField

using Random, Distributions, JuMP, Ipopt

using ..Magic

function product_sherrington_kirkpatrick(N, J, g)
    gauss = Normal(0, sqrt(J/N))
    J_sum = reduce(+, rand(gauss, N^2))
    Z_mean(α) = (3 * α^2 - 1) / 2
    X_mean(α) = α * sqrt(2 * (1 - α^2)) + (1 - α^2) / 2
    E_mean(α) = - 2 * J_sum * Z_mean(α)^2 - 2 * g * N * X_mean(α)
    return E_mean
end

product_state(α) = [α; sqrt((1 - α^2)/2); sqrt((1 - α^2)/2)]

function minimize_energy(E)
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, 0 <= α <= 1)
    @objective(model, Min, E(α))
    optimize!(model)
    α_min = value.(α)
    if E(0) < E(α_min)
        α_min = 0
    end
    if E(1) < E(α_min)
        α_min = 1
    end
    return α_min
end

function mean_field_sk_mana(N, J, g_max, iters, steps)
    gs = []
    manas = []
    for g in LinRange(0, g_max, steps)
        println("g = ", g)
        append!(gs, g)
        avg_mana = 0
        for i=1:iters
            E = product_sherrington_kirkpatrick(N, J, g)
            α_min = minimize_energy(E)
            min_state = product_state(α_min)
            rho = ComplexF64.(min_state * min_state')
            avg_mana += Magic.mana(rho) * N
        end
        avg_mana /= iters
        append!(manas, avg_mana)
    end
    return gs, manas
end

end
