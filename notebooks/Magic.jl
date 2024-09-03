### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ e402a58a-ff73-4558-856e-47875c861c69
# ╠═╡ show_logs = false
begin
	using Pkg
	using LinearAlgebra
	using Plots
	using DataFrames
	using CSV
	using LaTeXStrings
	

	Pkg.activate(temp=true)
 	Pkg.develop(path="/home/vbettaque/Development/QStab.jl")
 	using QStab
	using QStab.Hamiltonians
	using QStab.Magic
end

# ╔═╡ 3c918d16-6a06-11ef-3eca-e104c47e8a00
function spin_glass_mana(N, J, g_max, iters, steps)
    gs = []
    manas = []
    for g in LinRange(0, g_max, steps)
        println("g = ", g)
        append!(gs, g)
        avg_mana = 0
        for i = 1:iters
            H = Hamiltonians.qutrit_sherrington_kirkpatrick(N, J, g, 2.0^(-17))
            # display(H)
            e = eigen(H)
            # display(e.values)
            ground_state = normalize(e.vectors[:,1])
            # println(round.(ground_state; digits=2))
            # display(e.vectors)
            rho = ground_state * ground_state'
            avg_mana += Magic.mana(rho)
        end
        avg_mana /= iters
        # @assert avg_mana < 0.1
        append!(manas, avg_mana)
    end
    return gs, manas
end

# ╔═╡ d827144b-0e8b-46db-9446-977ffdd7493a
begin
Ns = 2:6
J = 1
g_max = 10
iters = 100
steps = 10000
path = "../data/magic/transverse_sk/"
!ispath(path) && mkpath(path)
for N = Ns
    gs, manas = spin_glass_mana(N, J, g_max, iters, steps)
    filename = "sk_N" * string(N) * "J" * string(J) * "g" * string(g_max) * "s" * string(steps) * "i" * string(iters) * ".csv"
    labels = ["g", "M"]
    frame = DataFrame(hcat(gs, manas), labels)
    CSV.write(path * filename, frame)
end
end

# ╔═╡ 11b16d8e-8f8a-43b3-9b53-2858ada36434
begin
	data_2 = CSV.read("../data/magic/transverse_sk/sk_N2J1g100s10000i1.csv", DataFrame)
	data_3 = CSV.read("../data/magic/transverse_sk/sk_N3J1g100s10000i1.csv", DataFrame)
	data_4 = CSV.read("../data/magic/transverse_sk/sk_N4J1g100s10000i1.csv", DataFrame)
	data_5 = CSV.read("../data/magic/transverse_sk/sk_N5J1g100s10000i1.csv", DataFrame)
	data_6 = CSV.read("../data/magic/transverse_sk/sk_N6J1g100s10000i1.csv", DataFrame)
end

# ╔═╡ b84a20b1-1dd3-41f9-b718-8459d03807a8
begin
	g_range = 1:1000
	scatter(data_2[!, 1][g_range], data_2[!, 2][g_range], ms=3, label=L"N=2")
	scatter!(data_3[!, 1][g_range], data_3[!, 2][g_range], ms=3, label=L"N=3")
	scatter!(data_4[!, 1][g_range], data_4[!, 2][g_range], ms=3, label=L"N=4")
	scatter!(data_5[!, 1][g_range], data_5[!, 2][g_range], ms=3, label=L"N=5")
	scatter!(data_6[!, 1][g_range], data_6[!, 2][g_range], ms=3, label=L"N=6")
	xlabel!(L"g")
	ylabel!(L"\mathcal{M}")
end

# ╔═╡ Cell order:
# ╠═e402a58a-ff73-4558-856e-47875c861c69
# ╠═3c918d16-6a06-11ef-3eca-e104c47e8a00
# ╠═d827144b-0e8b-46db-9446-977ffdd7493a
# ╠═11b16d8e-8f8a-43b3-9b53-2858ada36434
# ╠═b84a20b1-1dd3-41f9-b718-8459d03807a8
