using DataFrames
using CSV

using LinearAlgebra
using QStab
using QStab.Symplectics, QStab.Orthogonals, QStab.Stabilizers, QStab.Binary, QStab.Utils, QStab.Circuits, QStab.FramePotential

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

        for i in 1:reps
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

        filename = "weights_n" * string(n) * "r" * string(reps) * ".csv"
        path = "/Users/vbettaque/Development/Julia/QStab.jl/data/" * filename
        frame = DataFrame(data', ["Symp 1", "Symp 2", "Ortho 1", "Ortho 2"])
        CSV.write(path, frame)

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
        path = "/Users/vbettaque/Development/Julia/QStab.jl/data/single_clifford/stabilizer/" * filename

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
        CSV.write(path, frame)

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

function frame_potentials(n, t_max; max_reps = -1)
    data = zeros(t_max, 5)

    for t=1:t_max
        unitary = factorial(t)
        symp = FramePotential.symplectic(t, n; max_reps = max_reps)
        ortho = FramePotential.orthogonal(t, n; max_reps = max_reps)
        ortho_even = FramePotential.orthogonal_even(t, n; with_pcheck = false, max_reps = max_reps)
        ortho_even_pcheck = FramePotential.orthogonal_even(t, n; with_pcheck = true, max_reps = max_reps)

        data[t, :] = [unitary symp ortho ortho_even ortho_even_pcheck]
    end

    frame = DataFrame(data, ["Unitary", "Symplectic", "Orthogonal", "Orthogonal (Even)", "Orthogonal (Even, p-Check)"])

    filename = "frame_n" * string(n) * "t" * string(t_max) * (max_reps > 0 ? "r"*string(max_reps) : "")* ".csv"

    path = "/home/vbettaque/Development/QStab.jl/data/frame_potential/"
    !ispath(path) && mkpath(path)

    CSV.write(path * filename, frame)
end

# frame_potentials(6, 10; max_reps = -1)

#q_local_distance()


#FramePotential.even_parity(1, 10; with_pcheck=false, max_reps = 10000) * 2

# generate_data()

N = 2
for i=1:Orthogonals.group_order(N)
    O = Orthogonals.indexed_element(N, i)
    display(O)
    @assert(O' * O == I)
end