include("dsyms.jl")


function isPseudoConvex(dsRaw::AbstractDelaneySymbol)
    ds = orientedCover(dsRaw)
    ori = partialOrientation(ds)
    hasHandles = eulerCharacteristic(ds) <= 0

    for A1 in filter(D -> ori[D] > 0, 1 : size(ds))
        seen1 = falses(size(ds))
        seen1[A1] = true

        A2 = get(ds, 1, A1)
        while true
            A2 = get(ds, 0, get(ds, 1, A2))
            if seen1[A2]
                break
            elseif A2 == get(ds, 0, A1) || A2 == get(ds, 2, A1)
                seen1[A2] = seen1[get(ds, 1, A2)] = true
                continue
            end

            seen2 = copy(seen1)
            seen2[A2] = seen1[A2] = seen1[get(ds, 1, A2)] = true

            B2 = get(ds, 1, A2)
            while true
                B2 = get(ds, 2, get(ds, 1, B2))
                if seen2[B2]
                    if B2 == A1 && cutsOffDisk(ds, hasHandles, A1, A2)
                        return false
                    end
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

                    if B1 == get(ds, 0, B2) || B1 == get(ds, 2, A1)
                        continue
                    end

                    T = get(ds, 1, B1)
                    while true
                        T = get(ds, 2, get(ds, 1, T))
                        if seen3[T]
                            if (
                                T == A1 &&
                                cutsOffDisk(ds, hasHandles, A1, A2, B2, B1)
                            )
                                return false
                            end
                            break
                        end
                    end
                end
            end
        end
    end

    return true
end


function cutsOffDisk(
    ds::AbstractDelaneySymbol, hasHandles::Bool,
    A::Int64, B::Int64
)
    goodCones::Vector{Vector{Int64}} = [[], [2]]

    if !hasHandles
        A1 = get(ds, 1, A)
        B1 = get(ds, 1, B)

        if A1 == B
            vs = [v(ds, 0, 1, A), v(ds, 1, 2, A)]
            if filter(v -> v > 1, vs) in goodCones
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

    return checkPatch(ds, [A, B], goodCones)
end


function cutsOffDisk(
    ds::AbstractDelaneySymbol, hasHandles::Bool,
    A::Int64, B::Int64, C::Int64, D::Int64
)
    goodCones::Vector{Vector{Int64}} = [[]]

    if !hasHandles
        A1 = get(ds, 1, A)
        B1 = get(ds, 1, B)
        C1 = get(ds, 1, C)
        D1 = get(ds, 1, D)

        if (A1 == B && C1 == D) || (A1 == D && B1 == C)
            vs = [v(ds, 0, 1, A), v(ds, 1, 2, A), v(ds, 1, 2, B)]
            if filter(v -> v > 1, vs) in goodCones
                return false
            end
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

    return checkPatch(ds, [A, B, C, D], goodCones)
end


function checkPatch(
    ds::AbstractDelaneySymbol, cut::Vector{Int64},
    goodCones::Vector{Vector{Int64}}
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

    cones = []
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
                        push!(cones, v(ds, i, j, D))
                    end
                end
            end
        end
    end

    return eulerChar == 1 && filter(v -> v > 1, cones) in goodCones
end


function eulerCharacteristic(ds::AbstractDelaneySet)
    nrLoops(i) = count(D -> get(ds, i, D) == D, 1 : size(ds))
    nrOrbits(i, j) = length(orbits(ds, i, j))

    nf = size(ds)
    ne = div(3 * nf + nrLoops(0) + nrLoops(1) + nrLoops(2), 2)
    nv = nrOrbits(0, 1) + nrOrbits(0, 2) + nrOrbits(1, 2)

    return nf - ne + nv
end
