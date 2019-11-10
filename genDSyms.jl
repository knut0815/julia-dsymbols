include("backTracker.jl")
include("dsets.jl")


struct Orbit
    index::Int
    elements::Vector{Int}
    isChain::Bool
end

Base.length(orb::Orbit) = length(orb.elements)

r(orb::Orbit) = orb.isChain ? length(orb) : div(length(orb) + 1, 2)

minV(orb::Orbit) = Int(ceil(3 / r(orb)))


function curvature(ds::DSet, orbs::Vector{Orbit}, vs::Vector{Int})
    result = -size(ds)//2

    for i in 1 : length(orbs)
        result += (orbs[i].isChain ? 1 : 2) // vs[i]
    end

    return result
end


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


function morphism(ds::DSet, D0::Int)
    m = zeros(Int, size(ds))
    m[1] = D0
    q = [(1, D0)]

    while length(q) > 0
        D, E = popfirst!(q)

        for i in 0 : dim(ds)
            Di = get(ds, i, D)
            Ei = get(ds, i, E)

            if Di > 0 || Ei > 0
                if m[Di] == 0
                    m[Di] = Ei
                    push!(q, (Di, Ei))
                elseif m[Di] != Ei
                    return nothing
                end
            end
        end
    end

    return m
end


function automorphisms(ds::DSet)
    result::Vector{Vector{Int}} = []

    for D in 1 : size(ds)
        map = morphism(ds, D)
        if map != nothing
            push!(result, map)
        end
    end

    return result
end


function onOrbits(map::Vector{Int}, orbs::Vector{Orbit}, ds::DSet)
    inOrb = zeros(Int, dim(ds), size(ds))

    for i in 1 : length(orbs)
        for D in orbs[i].elements
            inOrb[orbs[i].index + 1, D] = i
        end
    end

    orbMap = zeros(Int, length(orbs))

    for D in 1 : size(ds)
        for i in 0 : dim(ds) - 1
            orbMap[inOrb[i + 1, D]] = inOrb[i + 1, map[D]]
        end
    end

    return orbMap
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
            print(io, r(orb) * minV(orb))
        end
    end

    print(io, ">")
end


for ds in DSetGenerator(2, parse(Int, ARGS[1]))
    println(ds)

    orbs = vcat(orbits(ds, 0), orbits(ds, 1))
    vs = map(minV, orbs)
    println("# $(vs)")

    curv = curvature(ds, orbs, vs)
    println("# $(curv)")

    elmMaps = automorphisms(ds)

    orbMaps = Set{Vector{Int}}()
    for m in elmMaps
        push!(orbMaps, onOrbits(m, orbs, ds))
    end
    println("# $(orbMaps)")

    println()
end
