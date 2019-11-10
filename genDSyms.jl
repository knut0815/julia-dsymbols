include("backTracker.jl")
include("dsets.jl")


struct Orbit
    index::Int
    elements::Vector{Int}
    isChain::Bool
end

Base.length(orb::Orbit) = length(orb.elements)

r(orb::Orbit) = orb.isChain ? length(orb) : div(length(orb) + 1, 2)

minV(orb::Orbit) = r(orb) * Int(ceil(3 / r(orb)))


function orbits(ds::DSet, i::Int)
    seen = falses(size(ds))
    result::Vector{Orbit} = []

    for D in 1 : size(ds)
        if !seen[D]
            orb = [D]
            seen[D] = true
            isChain::Bool = false

            E = D
            k = i

            while true
                Ek = get(ds, k, E)
                isChain = isChain || Ek == E
                E = Ek == 0 ? E : Ek
                k = i + i + 1 - k

                if !seen[E]
                    seen[E] = true
                    push!(orb, E)
                end

                if E == D && k == i
                    break
                end
            end

            push!(result, Orbit(i, orb, isChain))
        end
    end

    return result
end


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

        for orb in orbits(ds, i)
            if first(orb.elements) > 1
                print(io, " ")
            end
            print(io, minV(orb))
        end
    end

    print(io, ">")
end


for ds in DSetGenerator(2, parse(Int, ARGS[1]))
    println(ds)
end
