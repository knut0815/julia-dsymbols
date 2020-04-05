abstract type AbstractDelaneySet end



struct DelaneySetUnderConstruction <: AbstractDelaneySet
    op::Array{Int64,2}
end

DelaneySetUnderConstruction(size, dim) =
    DelaneySetUnderConstruction(zeros(Int64, size, dim + 1))

Base.size(ds::DelaneySetUnderConstruction) = first(size(ds.op))

dim(ds::DelaneySetUnderConstruction) = last(size(ds.op)) - 1

get(ds::DelaneySetUnderConstruction, i::Int64, D::Int64) =
    1 <= D <= size(ds) ? ds.op[D, i + 1] : 0

function set!(ds::DelaneySetUnderConstruction, i::Int64, D::Int64, E::Int64)
    ds.op[D, i + 1] = E
    ds.op[E, i + 1] = D
end



struct DelaneySet <: AbstractDelaneySet
    op::Array{Int64,2}

    function DelaneySet(op::Array{Int64,2})
        # TODO add consistency checks here
        new(op)
    end
end

DelaneySet(ds::DelaneySetUnderConstruction) = DelaneySet(ds.op)

Base.size(ds::DelaneySet) = first(size(ds.op))

dim(ds::DelaneySet) = last(size(ds.op)) - 1

get(ds::DelaneySet, i::Int64, D::Int64) =
    1 <= D <= size(ds) ? ds.op[D, i + 1] : 0



struct Orbit
    indices::Vector{Int64}
    elements::Vector{Int64}
    isChain::Bool
end

Base.length(orb::Orbit) = length(orb.elements)

r(orb::Orbit) = orb.isChain ? length(orb) : div(length(orb) + 1, 2)

minV(orb::Orbit) = Int64(ceil(3 / r(orb)))



function partialOrientation(ds::AbstractDelaneySet)
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


function isLoopless(ds::AbstractDelaneySet)
    for i in 0 : dim(ds)
        for D in 1 : size(ds)
            if get(ds, i, D) == D
                return false
            end
        end
    end

    return true
end


function isWeaklyOriented(ds::AbstractDelaneySet)
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


isOriented(ds::AbstractDelaneySet) = isLoopless(ds) && isWeaklyOriented(ds)


function orbits(ds::AbstractDelaneySet, i::Int64, j::Int64)
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

            push!(result, Orbit([i, j], orb, isChain))
        end
    end

    return result
end


function morphism(ds::AbstractDelaneySet, D0::Int64)
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


function automorphisms(ds::AbstractDelaneySet)
    result::Vector{Vector{Int64}} = []

    for D in 1 : size(ds)
        map = morphism(ds, D)
        if map != nothing
            push!(result, map)
        end
    end

    return result
end


function orientedCover(ds::AbstractDelaneySet)
    if isOriented(ds)
        return ds
    else
        sz = size(ds)
        ori = partialOrientation(ds)
        cov = DelaneySetUnderConstruction(2 * sz, dim(ds))

        for i in 0 : dim(ds)
            for D in 1 : sz
                Di = get(ds, i, D)
                if ori[Di] != ori[D]
                    set!(cov, i, D, Di)
                    set!(cov, i, D + sz, Di + sz)
                else
                    set!(cov, i, D, Di + sz)
                    set!(cov, i, D + sz, Di)
                end
            end
        end

        return DelaneySet(cov)
    end
end
