include("backTracker.jl")
include("dsyms.jl")


struct DSymState
    vs::Vector{Int64}
    curv::Rational{Int64}
    next::Int64
end


struct DSymGenerator <: BackTracker{DSym, DSymState}
    dset::DSet
    orbs::Vector{Orbit}
    orbMaps::Set{Vector{Int64}}

    function DSymGenerator(dset::DSet)
        orbs = orbits(dset)

        orbMaps = Set{Vector{Int64}}()
        for m in automorphisms(dset)
            push!(orbMaps, onOrbits(m, orbs, dset))
        end

        return new(dset, orbs, orbMaps)
    end
end


function root(g::DSymGenerator)
    vs = map(minV, g.orbs)
    curv = curvature(g.dset, g.orbs, vs)
    return DSymState(vs, curv, 1)
end


function extract(g::DSymGenerator, st::DSymState)
    if st.next > length(g.orbs) && goodResult(g, st) && isCanonical(g, st)
        return DSym(g.dset, st.vs)
    end

    return nothing
end


function children(g::DSymGenerator, st::DSymState)
    result = []

    if st.next <= length(g.orbs)
        if st.curv < 0
            push!(result, DSymState(st.vs, st.curv, length(g.orbs) + 1))
        else
            orb = g.orbs[st.next]

            for v in st.vs[st.next] : 7
                vs = copy(st.vs)
                vs[st.next] = v
                curv = curvature(g.dset, g.orbs, vs)

                if curv >= 0 || isMinimallyHyperbolic(g.dset, g.orbs, vs)
                    push!(result, DSymState(vs, curv, st.next + 1))
                end

                if curv < 0
                    break;
                end
            end
        end
    end

    return result
end


function goodResult(g::DSymGenerator, st::DSymState)
    if st.curv <= 0
        return true
    else
        cones::Vector{Int64} = []
        corners::Vector{Int64} = []

        for orb in orbits(g.dset, 0, 2)
            if orb.isChain
                if length(orb.elements) == 1
                    push!(corners, 2)
                end
            else
                if length(orb.elements) == 2
                    push!(cones, 2)
                end
            end
        end

        for i in 1 : length(g.orbs)
            if st.vs[i] > 1
                if g.orbs[i].isChain
                    push!(corners, st.vs[i])
                else
                    push!(cones, st.vs[i])
                end
            end
        end
    end

    front = join(reverse(sort(cones)), "")
    middle = isLoopless(g.dset) ? "" : "*"
    back = join(reverse(sort(corners)), "")
    cross = isWeaklyOriented(g.dset) ? "" : "x"
    key = front * middle * back * cross

    goodKeys = [
        "", "*", "x",
        "532", "432", "332",
        "422", "322", "222",
        "44", "33", "22",
        "*532", "*432", "*332", "3*2",
        "*422", "*322", "*222", "2*4", "2*3", "2*2",
        "*44", "*33", "*22", "4*", "3*", "2*", "4x", "3x", "2x"
    ]

    return key in goodKeys
end


function isCanonical(g::DSymGenerator, st::DSymState)
    vs = st.vs

    for m in g.orbMaps
        if map(i -> vs[m[i]], 1 : length(vs)) > vs
            return false
        end
    end

    return true
end


function isMinimallyHyperbolic(ds::DSet, orbs::Vector{Orbit}, vs::Vector{Int64})
    curv = curvature(ds, orbs, vs)

    if curv >= 0
        return false
    else
        for i in 1 : length(orbs)
            k = orbs[i].isChain ? 1 : 2
            v = vs[i]
            if v > minV(orbs[i]) && curv - k // v + k // (v - 1) < 0
                return false
            end
        end
    end

    return true
end


function onOrbits(map::Vector{Int64}, orbs::Vector{Orbit}, ds::DSet)
    inOrb = zeros(Int64, dim(ds), size(ds))

    for i in 1 : length(orbs)
        for D in orbs[i].elements
            inOrb[orbs[i].index + 1, D] = i
        end
    end

    orbMap = zeros(Int64, length(orbs))

    for D in 1 : size(ds)
        for i in 0 : dim(ds) - 1
            orbMap[inOrb[i + 1, D]] = inOrb[i + 1, map[D]]
        end
    end

    return orbMap
end
