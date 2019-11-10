include("backTracker.jl")
include("dsets.jl")


function Base.show(io::IO, ds::DSet)
    print(io, "<1.1:", size(ds))
    if dim(ds) != 2
        print(io, " ", dim(ds))
    end
    print(io, ":")

    for i in 0 : dim(ds)
        if i > 0
            print(",")
        end
        for D in 1 : size(ds)
            E = get(ds, i, D)
            if E == 0 || E >= D
                if D > 1
                    print(io, " ")
                end
                print(io, E)
            end
        end
    end
    print(io, ":")

    for i in 0 : dim(ds) - 1
        if i > 0
            print(",")
        end
        j = i + 1
        seen = falses(size(ds))

        for D in 1 : size(ds)
            if !seen[D]
                if D > 1
                    print(io, " ")
                end
                E, k = D, 0
                while true
                    Ei = get(ds, i, E)
                    E = Ei == 0 ? E : Ei
                    seen[E] = true
                    Ej = get(ds, j, E)
                    E = Ej == 0 ? E : Ej
                    seen[E] = true
                    k += 1

                    if E == D
                        break
                    end
                end

                print(io, k * Int(ceil(3 / k)))
            end
        end
    end

    print(io, ">")
end


for ds in DSetGenerator(2, parse(Int, ARGS[1]))
    println(ds)
end
