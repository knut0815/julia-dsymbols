include("dsyms.jl")


struct ExplicitDelaneySymbol <: AbstractDelaneySymbol
    op::Array{Int64, 2}
    v::Array{Int64, 2}
end


Base.size(ds::ExplicitDelaneySymbol) = first(size(ds.op))

dim(ds::ExplicitDelaneySymbol) = last(size(ds.op)) - 1

get(ds::ExplicitDelaneySymbol, i::Int64, D::Int64) =
    1 <= D <= size(ds) ? ds.op[D, i + 1] : 0


function v(ds::ExplicitDelaneySymbol, i::Int64, j::Int64, D::Int64)
    if !(1 <= D <= size(ds))
        return 0
    elseif j == i + 1
        return ds.v[D, j]
    elseif i == j + 1
        return ds.v[D, i]
    elseif j != i && get(ds, i, D) == get(ds, j, D)
        return 2
    else
        return 1
    end
end



function isPseudoConvex(dsRaw::AbstractDelaneySymbol)
    ds = makeOriented(dsRaw)
    ori = partialOrientation(ds)
    hasHandles = eulerCharacteristic(ds) <= 0

    for A1 in filter(D -> ori[D] > 0, 1 : size(ds))
        if findDiskBoundingTwoCut(ds, hasHandles, A1)
            return false
        end
        if findDiskBoundingFourCut(ds, hasHandles, A1)
            return false
        end
    end

    return true
end


function makeOriented(ds::AbstractDelaneySymbol)
    s = size(ds)
    d = dim(ds)

    if isOriented(ds)
        op = zeros(Int64, s, d + 1)
        vs = zeros(Int64, s, d)

        for D in 1 : s
            for i in 0 : d
                op[D, i + 1] = get(ds, i, D)
            end
            for i in 1 : d
                vs[D, i] = v(ds, i - 1, i, D)
            end
        end

        return ExplicitDelaneySymbol(op, vs)
    else
        op = zeros(Int64, 2 * s, d + 1)
        vs = zeros(Int64, 2 * s, d)

        for D in 1 : s
            for i in 0 : d
                Di = get(ds, i, D)
                op[D, i + 1] = Di + s
                op[D + s, i + 1] = Di
            end
            for i in 1 : d
                vs[D, i] = vs[D + s, i] = v(ds, i - 1, i, D)
            end
        end

        return ExplicitDelaneySymbol(op, vs)
    end
end


function findDiskBoundingTwoCut(
    ds::AbstractDelaneySymbol, hasHandles::Bool, A::Int64
)
    seen = falses(size(ds))
    seen[A] = seen[get(ds, 2, A)] = true

    B = get(ds, 1, A)
    while true
        B = get(ds, 0, get(ds, 1, B))
        if seen[B]
            break
        elseif B == get(ds, 0, A)
            seen[B] = seen[get(ds, 1, B)] = true
            continue
        end

        seen[B] = seen[get(ds, 1, B)] = true

        T = get(ds, 1, B)
        while true
            T = get(ds, 2, get(ds, 1, T))
            if T == A
                return boundsDisk(ds, hasHandles, A, B)
            elseif seen[T]
                break
            end
        end
    end

    return false
end


function findDiskBoundingFourCut(
    ds::AbstractDelaneySymbol, hasHandles::Bool, A1::Int64
)
    seen1 = falses(size(ds))
    seen1[A1] = seen1[get(ds, 2, A1)] = true

    A2 = get(ds, 1, A1)
    while true
        A2 = get(ds, 0, get(ds, 1, A2))
        if seen1[A2]
            break
        elseif A2 == get(ds, 0, A1)
            seen1[A2] = seen1[get(ds, 1, A2)] = true
            continue
        end

        seen2 = copy(seen1)
        seen2[A2] = seen1[A2] = seen1[get(ds, 1, A2)] = true

        B2 = get(ds, 1, A2)
        while true
            B2 = get(ds, 2, get(ds, 1, B2))
            if seen2[B2]
                break
            elseif B2 < A1 || B2 == get(ds, 2, A2)
                seen2[B2] = seen2[get(ds, 1, B2)] = true
                continue
            end

            seen3 = copy(seen2)
            seen3[B2] = seen2[B2] = seen2[get(ds, 1, B2)] = true

            B1 = get(ds, 1, B2)
            while true
                B1 = get(ds, 0, get(ds, 1, B1))
                if seen3[B1]
                    break
                end

                seen3[B1] = seen3[get(ds, 1, B1)] = true

                if B1 == get(ds, 0, B2)
                    continue
                end

                T = get(ds, 1, B1)
                while true
                    T = get(ds, 2, get(ds, 1, T))
                    if T == A1 && boundsDisk(ds, hasHandles, A1, A2, B2, B1)
                        return true
                    elseif seen3[T]
                        break
                    end
                end
            end
        end
    end

    return false
end


function boundsDisk(
    ds::AbstractDelaneySymbol, hasHandles::Bool,
    A::Int64, B::Int64
)
    if !hasHandles
        A1 = get(ds, 1, A)
        B1 = get(ds, 1, B)

        if A1 == B
            v01 = v(ds, 0, 1, A)
            v12 = v(ds, 1, 2, A)

            if (v01 == 1 && v12 <= 2) || (v01 <= 2 && v12 == 1)
                return false
            end
        end

        if (
            get(ds, 0, A1) == B1 && get(ds, 2, A1) == B1 &&
            v(ds, 0, 1, A) == 1 && v(ds, 1, 2, A) == 1
        )
            return false
        end
    end

    return checkPatch(ds, [A, B], true)
end


function boundsDisk(
    ds::AbstractDelaneySymbol, hasHandles::Bool,
    A::Int64, B::Int64, C::Int64, D::Int64
)
    if !hasHandles
        A1 = get(ds, 1, A)
        B1 = get(ds, 1, B)
        C1 = get(ds, 1, C)
        D1 = get(ds, 1, D)

        if (
            ((A1 == B && C1 == D) || (A1 == D && B1 == C)) &&
            v(ds, 0, 1, A) == 1 &&
            v(ds, 1, 2, A) == 1 &&
            v(ds, 1, 2, B) == 1
        )
            return false
        end

        if (
            get(ds, 0, A1) == B1 &&
            get(ds, 2, B1) == C1 &&
            get(ds, 0, C1) == D1 &&
            v(ds, 0, 1, A) == 1 && v(ds, 1, 2, B) == 1 &&
            v(ds, 0, 1, C) == 1 && v(ds, 1, 2, D) == 1
        )
            return false
        end
    end

    return checkPatch(ds, [A, B, C, D], false)
end


function checkPatch(
    ds::AbstractDelaneySymbol, cut::Vector{Int64}, allow2Cone::Bool
)
    seed = cut[1]

    inCut = falses(size(ds))
    outside = falses(size(ds))
    for D in cut
        D1 = get(ds, 1, D)
        inCut[D] = inCut[D1] = true
        if !(D1 in cut)
            outside[D1] = true
        end
    end

    inPatch = falses(size(ds))
    inPatch[seed] = true
    elements = [seed]
    next = 1
    nrLoops = 0

    while next <= length(elements)
        D = elements[next]
        next += 1

        for i in 0 : dim(ds)
            Di = (i == 1 && inCut[D]) ? D : get(ds, i, D)
            if Di == D
                nrLoops += 1
            elseif outside[Di]
                return false
            elseif !inPatch[Di]
                inPatch[Di] = true
                push!(elements, Di)
            end
        end
    end

    seen2Cone = false
    eulerChar = -div(length(elements) + nrLoops, 2)

    for i in 0 : dim(ds) - 1
        for j in i + 1 : dim(ds)
            seen = falses(size(ds))
            for D in elements
                if !seen[D]
                    eulerChar += 1
                    isChain = false

                    E = D
                    k = i

                    while true
                        Ek = (k == 1 && inCut[E]) ? E : get(ds, k, E)
                        seen[Ek] = true
                        if Ek == E
                            isChain = true
                        end

                        E = Ek
                        k = i + j - k
                        if E == D && k == i
                            break
                        end
                    end

                    if !isChain
                        vD = v(ds, i, j, D)
                        if vD > 2
                            return false
                        elseif vD == 2
                            if !allow2Cone || seen2Cone
                                return false
                            else
                                seen2Cone = true
                            end
                        end
                    end
                end
            end
        end
    end

    return eulerChar == 1
end


function eulerCharacteristic(ds::AbstractDelaneySet)
    nrLoops(i) = count(D -> get(ds, i, D) == D, 1 : size(ds))
    nrOrbits(i, j) = length(orbits(ds, i, j))

    nf = size(ds)
    ne = div(3 * nf + nrLoops(0) + nrLoops(1) + nrLoops(2), 2)
    nv = nrOrbits(0, 1) + nrOrbits(0, 2) + nrOrbits(1, 2)

    return nf - ne + nv
end
