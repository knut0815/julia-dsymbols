struct DSet
    op::Array{Int64,2}
end


Base.size(ds::DSet) = first(size(ds.op))

dim(ds::DSet) = last(size(ds.op)) - 1

get(ds::DSet, i::Int64, D::Int64) = 1 <= D <= size(ds) ? ds.op[D, i + 1] : 0


function set!(ds::DSet, i::Int64, D::Int64, E::Int64)
    ds.op[D, i + 1] = E
    ds.op[E, i + 1] = D
end


function partialOrientation(ds::DSet)
    ori = zeros(Int64, size(ds))
    ori[1] = 1
    queue::Vector{Int64} = [1]

    while length(queue) > 0
        D = pop!(queue)

        for i in 0 : dim(ds)
            Di = get(ds, i, D)
            if ori[Di] == 0
                ori[Di] = -ori[D]
                push!(queue, Di)
            end
        end
    end

    return ori
end


function isLoopless(ds::DSet)
    for i in 0 : dim(ds)
        for D in 1 : size(ds)
            if get(ds, i, D) == D
                return false
            end
        end
    end

    return true
end


function isWeaklyOriented(ds::DSet)
    ori = partialOrientation(ds)

    for i in 0 : dim(ds)
        for D in 1 : size(ds)
            Di = get(ds, i, D)
            if Di != D && ori[Di] == ori[D]
                return false
            end
        end
    end

    return true
end


struct Orbit
    index::Int64
    elements::Vector{Int64}
    isChain::Bool
end

Base.length(orb::Orbit) = length(orb.elements)

r(orb::Orbit) = orb.isChain ? length(orb) : div(length(orb) + 1, 2)

minV(orb::Orbit) = Int64(ceil(3 / r(orb)))


function orbits(ds::DSet, i::Int64, j::Int64)
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
                k = i + j - k

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

orbits(ds::DSet) = vcat(orbits(ds, 0, 1), orbits(ds, 1, 2))


function morphism(ds::DSet, D0::Int64)
    m = zeros(Int64, size(ds))
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
    result::Vector{Vector{Int64}} = []

    for D in 1 : size(ds)
        map = morphism(ds, D)
        if map != nothing
            push!(result, map)
        end
    end

    return result
end
