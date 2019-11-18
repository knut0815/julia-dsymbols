include("backTracker.jl")


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


struct DSetState
    dset::DSet
    isRemapStart::BitVector
end


struct DSetGenerator <: BackTracker{DSet, DSetState}
    dim::Int64
    maxSize::Int64
end


root(g::DSetGenerator) =
    DSetState(DSet(zeros(Int64, 1, g.dim + 1)), falses(g.maxSize))

extract(g::DSetGenerator, st::DSetState) =
    firstUndefined(st.dset) == nothing ? st.dset : nothing


function children(g::DSetGenerator, st::DSetState)
    result = []
    ds = st.dset
    undef = firstUndefined(ds)

    if undef != nothing
        D, i = undef

        for E in D : min(size(ds) + 1, g.maxSize)
            if get(ds, i, E) == 0
                isRemapStart = copy(st.isRemapStart)

                if E > size(ds)
                    dset = DSet(vcat(ds.op, zeros(Int64, 1, dim(ds) + 1)))
                    isRemapStart[E] = true
                else
                    dset = DSet(copy(ds.op))
                end

                set!(dset, i, D, E)

                head, tail, gap, k = scan02Orbit(dset, D)

                if gap == 1
                    set!(dset, k, head, tail)
                elseif gap == 0 && head != tail
                    continue;
                end

                if checkCanonicity!(dset, isRemapStart)
                    push!(result, DSetState(dset, isRemapStart))
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


function scan02Orbit(ds::DSet, D::Int64)
    head, i = scan(ds, [0, 2, 0, 2], D, 4)
    tail, j = scan(ds, [2, 0, 2, 0], D, 4 - i)

    return head, tail, 4 - i - j, 2 * (i % 2)
end


function scan(ds::DSet, w::Vector{Int64}, D::Int64, limit::Int64)
    E, k = D, 1

    while k <= limit && get(ds, w[k], E) != 0
        E = get(ds, w[k], E)
        k += 1
    end

    return E, k - 1
end


function checkCanonicity!(ds::DSet, isRemapStart::BitVector)
    n2o = zeros(Int64, size(ds))
    o2n = zeros(Int64, size(ds))

    for D in 1 : size(ds)
        if isRemapStart[D]
            d = compareRenumberedFrom(ds, D, n2o, o2n)
            if d < 0
                return false
            elseif d > 0
                isRemapStart[D] = false
            end
        end
    end

    return true
end


function compareRenumberedFrom(
    ds::DSet, D0::Int64, n2o::Vector{Int64}, o2n::Vector{Int64}
)
    fill!(n2o, zero(Int64))
    fill!(o2n, zero(Int64))

    n2o[1] = D0
    o2n[D0] = 1
    next = 2

    for D in 1 : size(ds)
        for i in 0 : dim(ds)
            Ei = get(ds, i, n2o[D])

            if Ei == 0
                return 0
            else
                if o2n[Ei] == 0
                    o2n[Ei] = next
                    n2o[next] = Ei
                    next += 1
                end

                Di = get(ds, i, D)

                if Di == 0
                    return 0
                elseif o2n[Ei] != Di
                    return o2n[Ei] - Di
                end
            end
        end
    end

    return 0
end
