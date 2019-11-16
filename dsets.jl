include("backTracker.jl")


struct DSet
    op::Array{Int,2}
end


Base.size(ds::DSet) = first(size(ds.op))

dim(ds::DSet) = last(size(ds.op)) - 1

get(ds::DSet, i::Int, D::Int) = 1 <= D <= size(ds) ? ds.op[D, i + 1] : 0


function set!(ds::DSet, i::Int, D::Int, E::Int)
    ds.op[D, i + 1] = E
    ds.op[E, i + 1] = D
end


function partialOrientation(ds::DSet)
    ori = zeros(Int, size(ds))
    ori[1] = 1

    for D in 1 : size(ds)
        for i in 0 : dim(ds)
            Di = get(ds, i, D)
            if Di != D
                ori[Di] = -ori[D]
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


function orbits(ds::DSet, i::Int, j::Int)
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


struct DSetGenerator <: BackTracker{DSet, DSet}
    dim::Int
    maxSize::Int
end


root(g::DSetGenerator) =
    DSet(zeros(Int, 1, g.dim + 1))

extract(g::DSetGenerator, ds::DSet) =
    firstUndefined(ds) == nothing ? ds : nothing


function children(g::DSetGenerator, ds::DSet)
    result = []
    undef = firstUndefined(ds)

    if undef != nothing
        D, i = undef

        for E in D : min(size(ds) + 1, g.maxSize)
            if get(ds, i, E) == 0
                if E > size(ds)
                    out = DSet(vcat(ds.op, zeros(Int, 1, dim(ds) + 1)))
                else
                    out = DSet(copy(ds.op))
                end

                set!(out, i, D, E)

                head, tail, gap, k = scan02Orbit(out, D)

                if gap == 1
                    set!(out, k, head, tail)
                elseif gap == 0 && head != tail
                    continue;
                end

                if isCanonical(out)
                    push!(result, out)
                end
            end
        end
    end

    return result
end


function firstUndefined(ds::DSet)
    for D in 1 : size(ds)
        for i in 0 : dim(ds)
            if get(ds, i, D) == 0
                return D, i
            end
        end
    end

    return nothing
end


function scan02Orbit(ds::DSet, D::Int)
    head, i = scan(ds, [0, 2, 0, 2], D, 4)
    tail, j = scan(ds, [2, 0, 2, 0], D, 4 - i)

    return head, tail, 4 - i - j, 2 * (i % 2)
end


function scan(ds::DSet, w::Vector{Int}, D::Int, limit::Int)
    E, k = D, 1

    while k <= limit && get(ds, w[k], E) != 0
        E = get(ds, w[k], E)
        k += 1
    end

    return E, k - 1
end


function isCanonical(ds::DSet)
    n2o = zeros(Int, size(ds))
    o2n = zeros(Int, size(ds))

    for D in 1 : size(ds)
        if compareRenumberedFrom(ds, D, n2o, o2n) < 0
            return false
        end
    end

    return true
end


function compareRenumberedFrom(
    ds::DSet, D0::Int, n2o::Vector{Int}, o2n::Vector{Int}
)
    for D in 1 : size(ds)
        n2o[D] = 0
        o2n[D] = 0
    end

    n2o[1] = D0
    o2n[D0] = 1
    next = 2

    for D in 1 : size(ds)
        for i in 0 : dim(ds)
            oval = get(ds, i, D)
            E = get(ds, i, n2o[D])

            if E == 0
                if oval != 0
                    return 1
                end
            else
                if o2n[E] == 0
                    o2n[E] = next
                    n2o[next] = E
                    next += 1
                end

                if oval == 0
                    return -1
                else
                    d = o2n[E] - oval
                    if d != 0
                        return d
                    end
                end
            end
        end
    end

    return 0
end
