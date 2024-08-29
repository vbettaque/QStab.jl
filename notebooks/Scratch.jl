using DataFrames
using CSV

using Random
using LinearAlgebra
using QStab
using QStab.Symplectics, QStab.Orthogonals, QStab.Stabilizers, QStab.Galois, QStab.Utils, QStab.Circuits, QStab.FramePotential, QStab.Hilbert, QStab.Magic, QStab.Hamiltonians

function single_clifford_single_string_data()
    reps = 100000
    ns = [22 24 26 28 30 32 34 36 38 40]

    for n in ns
        println("n = ", n)

        init_1 = zeros(GF2, n)
        init_1[1] = 1

        init_2 = zeros(GF2, n)
        init_2[1] = 1
        init_2[2] = 1

        data = zeros(Int, 4, reps)

        Threads.@threads for i in 1:reps
            O = Orthogonals.rand(n)
            S = Symplectics.rand_majorana(n)

            string_symp_1 = S * init_1
            string_symp_2 = S * init_2
            string_ortho_1 = O * init_1
            string_ortho_2 = O * init_2

            data[1, i] = weight_majorana(string_symp_1)
            data[2, i] = weight_majorana(string_symp_2)
            data[3, i] = weight_majorana(string_ortho_1)
            data[4, i] = weight_majorana(string_ortho_2)

            print("\e[2K")
            print("\e[1G")
            print(Int(floor(i / reps * 100)), "%")
        end

        file_name = "weights_n" * string(n) * "r" * string(reps) * ".csv"
        path = "data/single_clifford/single_string/"
        !ispath(path) && mkpath(path)
        frame = DataFrame(data', ["Symp 1", "Symp 2", "Ortho 1", "Ortho 2"])
        CSV.write(path * file_name, frame)

        println("")
    end

end

function single_clifford_stabilizer_data()
    reps = 100000
    ns = [22 24 26 28 30 32 34 36 38 40]

    for n in ns
        println("n = ", n)
        k = n ÷ 2

        init_even = canon_stab_matrix(n)
        init_odd = canon_stab_matrix(n, odd=true)

        data = zeros(Int, 8 * k, reps)

        Threads.@threads for i in 1:reps
            O = Orthogonals.rand(n)
            S = Symplectics.rand_majorana(n)

            stab_symp_even = init_even * S'
            stab_symp_odd = init_odd * S'
            stab_ortho_even = init_even * O'
            stab_ortho_odd = init_odd * O'

            rref_symp_even = rref(stab_symp_even)
            rref_symp_odd = rref(stab_symp_odd)
            rref_ortho_even = rref(stab_ortho_even)
            rref_ortho_odd = rref(stab_ortho_odd)

            data[1:k, i] = weight_majorana(stab_symp_even)
            data[(k+1):(2*k), i] = weight_majorana(rref_symp_even)

            data[(2*k+1):(3*k), i] = weight_majorana(stab_symp_odd)
            data[(3*k+1):(4*k), i] = weight_majorana(rref_symp_odd)

            data[(4*k+1):(5*k), i] = weight_majorana(stab_ortho_even)
            data[(5*k+1):(6*k), i] = weight_majorana(rref_ortho_even)

            data[(6*k+1):(7*k), i] = weight_majorana(stab_ortho_odd)
            data[(7*k+1):(8*k), i] = weight_majorana(rref_ortho_odd)

            print("\e[2K")
            print("\e[1G")
            print(Int(floor(i / reps * 100)), "%")
        end

        filename = "weights_n" * string(n) * "r" * string(reps) * ".csv"
        path = "data/single_clifford/stabilizer/"
        !ispath(path) && mkpath(path)

        labels_symp_even = map(i -> "Symp Even " * string(i), 1:k)
        labels_symp_even_rref = map(i -> "Symp Even RREF " * string(i), 1:k)
        labels_symp_odd = map(i -> "Symp Odd " * string(i), 1:k)
        labels_symp_odd_rref = map(i -> "Symp Odd RREF " * string(i), 1:k)
        labels_ortho_even = map(i -> "Ortho Even " * string(i), 1:k)
        labels_ortho_even_rref = map(i -> "Ortho Even RREF " * string(i), 1:k)
        labels_ortho_odd = map(i -> "Ortho Odd " * string(i), 1:k)
        labels_ortho_odd_rref = map(i -> "Ortho Odd RREF " * string(i), 1:k)

        labels = vcat(labels_symp_even, labels_symp_even_rref, labels_symp_odd, labels_symp_odd_rref,
            labels_ortho_even, labels_ortho_even_rref, labels_ortho_odd, labels_ortho_odd_rref)

        frame = DataFrame(data', labels)
        CSV.write(path * filename, frame)

        println("")
    end

end

function q_local_single_string_data()
    reps = 100000

    qs = 2:2:10

    Ds = 1:10

    for q in qs
        ns = (1:10).*q

        println("q = ", q)

        for n in ns

            println("n = ", n)

            init_1 = zeros(GF2, n)
            init_1[1] = 1

            init_2 = zeros(GF2, n)
            init_2[1] = 1
            init_2[2] = 1

            for D in Ds

                println("D = ", D)

                data = zeros(Int, 2, reps)

                for i in 1:reps
                    O = rand_q_local_volume(n, q, D, Orthogonals.rand)

                    string_ortho_1 = O * init_1
                    string_ortho_2 = O * init_2

                    data[1, i] = weight_majorana(string_ortho_1)
                    data[2, i] = weight_majorana(string_ortho_2)

                    print("\e[2K")
                    print("\e[1G")
                    print(Int(floor(i / reps * 100)), "%")
                end

                filename = "weights_n" * string(n) * "r" * string(reps) * ".csv"
                path = "/Users/vbettaque/Development/Julia/QStab.jl/data/q_local/single_string/q" * string(q) * "/D" * string(D) * "/" * filename
                frame = DataFrame(data', ["Ortho 1", "Ortho 2"])
                CSV.write(path, frame)

                println("")

            end

        end

    end

end

function q_local_stabilizer_data()
    reps = 100000

    qs = 2:2:10

    Ds = 1:10

    for q in qs
        ns = (1:10).*q

        println("q = ", q)

        for n in ns

            println("n = ", n)
            k = n ÷ 2

            init_even = canon_stab_matrix(n)
            init_odd = canon_stab_matrix(n, odd=true)

            for D in Ds

                println("D = ", D)

                data = zeros(Int, 4 * k, reps)


                Threads.@threads for i in 1:reps
                    O = rand_q_local_volume(n, q, D, Orthogonals.rand)

                    stab_ortho_even = init_even * O'
                    stab_ortho_odd = init_odd * O'

                    rref_ortho_even = rref(stab_ortho_even)
                    rref_ortho_odd = rref(stab_ortho_odd)

                    data[1:k, i] = weight_majorana(stab_ortho_even)
                    data[(k+1):(2*k), i] = weight_majorana(rref_ortho_even)

                    data[(2*k+1):(3*k), i] = weight_majorana(stab_ortho_odd)
                    data[(3*k+1):(4*k), i] = weight_majorana(rref_ortho_odd)

                    #print("\e[2K")
                    #print("\e[1G")
                    #print(Int(floor(i / reps * 100)), "%")
                end

                filename = "weights_n" * string(n) * "r" * string(reps) * ".csv"
                path = "/Users/vbettaque/Development/Julia/QStab.jl/data/q_local/stabilizer/q" * string(q) * "/D" * string(D) * "/"
                !ispath(path) && mkpath(path)

                labels_ortho_even = map(i -> "Ortho Even " * string(i), 1:k)
                labels_ortho_even_rref = map(i -> "Ortho Even RREF " * string(i), 1:k)
                labels_ortho_odd = map(i -> "Ortho Odd " * string(i), 1:k)
                labels_ortho_odd_rref = map(i -> "Ortho Odd RREF " * string(i), 1:k)

                labels = vcat(labels_ortho_even, labels_ortho_even_rref, labels_ortho_odd, labels_ortho_odd_rref)


                frame = DataFrame(data', labels)
                CSV.write(path * filename, frame)


                println("")

            end

        end

    end

end

function single_clifford_distance()
    reps = 100000
    n_mc = 1000
    ns = 52:2:128
    k = 2

    for n in ns
        println("n = ", n)

        data = zeros(Int, 2, reps)

        Threads.@threads for i in 1:reps
            O = Orthogonals.rand(n)
            S = Symplectics.rand_majorana(n)

            data[1, i] = code_distance(O, k, mc_max = n_mc)
            data[2, i] = code_distance(S, k, mc_max = n_mc)
        end

        filename = "distance_n" * string(n) * "k" * string(k) * "mc" * string(n_mc) * "r" * string(reps) * ".csv"
        path = "/home/vbettaque/Development/QStab.jl/data/single_clifford/distance/k" * string(k) * "/"
        !ispath(path) && mkpath(path)

        labels = ["Ortho", "Symp"]

        frame = DataFrame(data', labels)
        CSV.write(path * filename, frame)
    end

end

function q_local_distance()
    reps = 100000
    n_mc = 1000
    k = 2

    qs = 4:2:4

    Ds = 1:10

    for q in qs
        ns = (22:30).*q

        println("q = ", q)

        for n in ns

            if n <= k continue end

            println("n = ", n)

            for D in Ds

                println("D = ", D)

                data = zeros(Int, 1, reps)

                Threads.@threads for i in 1:reps
                    O = rand_q_local_volume(n, q, D, Orthogonals.rand)

                    data[1, i] = code_distance(O, k, mc_max = n_mc)
                end

                filename = "distance_n" * string(n) * "r" * string(reps) * ".csv"
                path = "/home/vbettaque/Development/QStab.jl/data/q_local/stabilizer/q" * string(q) * "/D" * string(D) * "/"
                !ispath(path) && mkpath(path)

                frame = DataFrame(data', ["Ortho"])
                CSV.write(path * filename, frame)
            end

            println("")

        end

    end

end

function even_frame_potentials(n_max, t_max; max_reps = -1)
    @assert iseven(n_max) && n_max > 2

    unit_data = zeros(t_max, n_max÷2-1)
    cliff_data = zeros(t_max, n_max÷2-1)
    pcliff_data = zeros(t_max, n_max÷2-1)


    for n=4:2:n_max
        println("n = ", n)
        for t=1:t_max
            println("t = ", t)
            unit_data[t, n÷2-1] = FramePotential.unitary(t, 2^(n÷2-1))
            cliff_data[t, n÷2-1], _ = FramePotential.clifford_symp(t, n-2; max_reps=max_reps)
            pcliff_data[t, n÷2-1], _ = FramePotential.pclifford_hilbert(t, n; max_reps=max_reps)
        end
        println("")
    end

    labels_unit = map(n -> "d = " * string(2^(n÷2)), 4:2:n_max)
    labels_cliff = map(n -> "n = " * string(n), 4:2:n_max)

    unit_frame = DataFrame(unit_data, labels_unit)
    cliff_frame = DataFrame(cliff_data, labels_cliff)
    pcliff_frame = DataFrame(pcliff_data, labels_cliff)

    filename_unit = "unitary_d" * string(2^(n_max÷2)) * "t" * string(t_max) * ".csv"
    filename_cliff = "clifford_n" * string(n_max) * "t" * string(t_max) * (max_reps > 0 ? "r"*string(max_reps) : "")* ".csv"
    filename_pcliff = "pclifford_n" * string(n_max) * "t" * string(t_max) * (max_reps > 0 ? "r"*string(max_reps) : "")* ".csv"

    path = "data/frame_potential/even/"
    !ispath(path) && mkpath(path)

    CSV.write(path * filename_unit, unit_frame)
    CSV.write(path * filename_cliff, cliff_frame)
    CSV.write(path * filename_pcliff, pcliff_frame)
end

# single_clifford_stabilizer_data()

# ns = 14:2:20
# reps = 1000000

# for n in ns
#     fixed_points = FramePotential.fixed_points_ortho(n; max_reps=reps)
#     data_frame = DataFrame(x = fixed_points)
#     file_name = "orthogonal_n" * string(n) * "r" * string(reps) * ".csv"
#     path = "data/fixed_points/orthogonal/"
#     !ispath(path) && mkpath(path)
#     CSV.write(path * file_name, data_frame, header=false)

#     fixed_points = FramePotential.fixed_points_symp(n-2; max_reps=reps)
#     data_frame = DataFrame(x = fixed_points)
#     file_name = "symplectic_n" * string(n-2) * "r" * string(reps) * ".csv"
#     path = "data/fixed_points/symplectic/"
#     !ispath(path) && mkpath(path)
#     CSV.write(path * file_name, data_frame, header=false)
# end

# Orthogonals.group_order(8)

# FramePotential.from_fixed_points(fixed_points, 5; exact=true)
# x = 1

# n = 8
# println("start")
# for t=1:10
#     println("symp: ", FramePotential.clifford_symp(t, n-2; max_reps = 1000000))
#     println("ortho: ", FramePotential.pclifford_symp(t, n; max_reps = 1000000))
#     println("")
# end


# n = 6
# majo_order = (big"2")^(n - 1)
# ortho = Orthogonals.rand(n)

# ortho_reduced = FramePotential.even_parity_sector_matrix(ortho)

# j_reduced = [GF2(isodd(i)) for i=1:(n-1)]

# eqs = hcat(ortho_reduced-I, j_reduced)
# display(rref(eqs))


# for i=1:majo_order
#     v = indexed_even_bitvec(i, n)
#     if ortho * v == complement(v)
#         println("has complement!")
#         break
#     end
# end
# println("done")


# n = 4
# order = Orthogonals.group_order(n)
# proj_even = (I + Hilbert.majorana_string(ones(GF2, n))) / 2
# for i=1:order
#     cliff = Hilbert.indexed_pclifford(n, i)
#     ortho = Orthogonals.indexed_element(n, i)
#     trace = Int(round(abs(tr(proj_even * cliff)^2)))

#     fixed = 0
#     complements = 0

#     for j=1:(big"2")^(n - 1)
#         v = indexed_even_bitvec(j, n)
#         if ortho * v == v
#             fixed += 1
#             # println("fixed: ", v)
#         end
#         if ortho * v == complement(v)
#             complements += 1
#             # println("complement: ", v)
#             println(v, " has ", (-1)^isone(parity(v[2:2:n])))
#         end
#     end
#     (fixed == complements == trace) && display(ortho)
#     println("trace = ", trace, " fixed = ", fixed, " complements = ", complements)
# end


# n = 6
# order = Orthogonals.group_order(n)
# orthos = Set([])
# for i=1:order
#     ortho =  Orthogonals.indexed_element(n, i)
#     if ortho ∈ orthos
#         display(ortho)
#     else
#         push!(orthos, ortho)
#     end
# end
# println(order, " ", length(orthos))

# even_frame_potentials(10, 10; max_reps = 100000)

#q_local_distance()


#FramePotential.even_parity(1, 10; with_pcheck=false, max_reps = 10000) * 2

# generate_data()

# N = 4
# for i=1:Orthogonals.group_order(N)
#     O = Orthogonals.indexed_element(N, i)
#     display(O)
#     @assert(O' * O == I)
# end

# Orthogonals.group_order(4)

# U = Hilbert.indexed_pclifford(4, 10)
# O = Orthogonals.indexed_element(4, 10)

# N = 4
# proj_even = (I + Hilbert.majorana_string(ones(GF2, N))) / 2
# println("start")
# for i=1:Orthogonals.group_order(N)
#     O = Orthogonals.indexed_element(N, i)
#     O_even = FramePotential.even_parity_sector_matrix(O)
#     U = Hilbert.indexed_pclifford(N, i)
#     trace = Int(round(abs(tr(proj_even * U)^2)))
#     fixed = 2^(N - 1 - Binary.rank(O_even - I)) ÷ 2
#     println(trace, " ", fixed)
#     if trace != fixed
#         display(O)
#     end
# end

# n = 8
# for t=1:10
#     println("t = ", t)
#     println("cliff: ", tdesign_even(t, n, 1000000))
#     # println("haar: ", factorial(2*t) / ( factorial(t) * factorial(t+1)))
#     println("haar: ", factorial(t))
# end

# n = 20
# size = (big"2")^(n - 1)
# for r = 1:10000
#     v = Galois.indexed_even_bitvec(rand(1:size), n)
#     w = Galois.indexed_even_bitvec(rand(1:size), n)
#     if all(isone.(v)) || all(isone.(w)) || all(iszero.(v)) || all(iszero.(w))
#         continue
#     end

#     h1, h2 = Orthogonals._householder_vector(v, w)
#     z = v + (v ⋅ h1) * h1
#     if !isnothing(h2)
#         z = z + (z ⋅ h2) * h2
#     end
#     @assert z == w
# end

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

using Plots
gs, manas = spin_glass_mana(4, 1, 20, 400, 400)
plot(gs, manas)

# for i=1:9
#     a = tritvec(i, 2)
#     @assert isapprox(Magic.phase_space_point(a), Magic.phase_space_point_ops[i], atol=1e-8)
# end
# println("done")

# H = Hamiltonians.qutrit_sherrington_kirkpatrick(4, 1, 0)
# display(H)

#  display(map(x -> isapprox(x, 0, atol=1e-8) ? 0 : x, Magic.phase_space_point(GF3.([1, 1, 2, 1]))))

#  Magic.phase_space_point(GF3.([1, 1, 2, 1]))