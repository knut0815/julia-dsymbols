include("dsyms.jl")


struct ExplicitDelaneySymbol <: AbstractDelaneySymbol
    op::Array{Int64, 2}
    v::Array{Int64, 2}
end


Base.size(ds::ExplicitDelaneySymbol) = first(size(ds.op))

dim(ds::ExplicitDelaneySymbol) = last(size(ds.op)) - 1

get(ds::ExplicitDelaneySymbol, i::Int64, D::Int64) = ds.op[D, i + 1]


function v(ds::ExplicitDelaneySymbol, i::Int64, j::Int64, D::Int64)
    if j == i + 1
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

    nv = countOrbits(ds, 0, 1) + countOrbits(ds, 0, 2) + countOrbits(ds, 1, 2)
    hasHandles = 2 * nv <= size(ds)

    for D in 1 : size(ds)
        if ori[D] > 0
            if findDiskBoundingTwoCut(ds, hasHandles, D)
                return false
            end
            if findDiskBoundingFourCut(ds, hasHandles, D)
                return false
            end
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


function countOrbits(ds::AbstractDelaneySymbol, i::Int64, j::Int64)
    n = 0
    seen = fill(false, size(ds))

    for D in 1 : size(ds)
        if !seen[D]
            n += 1
            E = D
            while !seen[E]
                seen[E] = seen[get(ds, i, E)] = true
                E = get(ds, j, get(ds, i, E))
            end
        end
    end

    return n
end


# The remaining functions in this file are specialized helpers for
# isPseudoConvex() which assume that the D-symbol they work on is oriented.

function findDiskBoundingTwoCut(
    ds::AbstractDelaneySymbol, hasHandles::Bool, A::Int64
)
    B = get(ds, 1, get(ds, 0, A))

    while B != A
        B = get(ds, 1, get(ds, 0, B))

        T = get(ds, 2, get(ds, 1, B))
        while T != B
            T = get(ds, 2, get(ds, 1, T))
            if T == A
                return boundsDisk(ds, hasHandles, A, get(ds, 1, B))
            end
        end
    end

    return false
end


function findDiskBoundingFourCut(
    ds::AbstractDelaneySymbol, hasHandles::Bool, A1::Int64
)
    seen1 = fill(false, size(ds))
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
            if v(ds, 0, 1, A) == 1 || v(ds, 1, 2, A) == 1
                return false
            end
        elseif get(ds, 0, A1) == B1 && get(ds, 2, A1) == B1
            if v(ds, 0, 1, A) == 1 && v(ds, 1, 2, A) == 1
                return false
            end
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
    elements = patchElements(ds, cut)
    if elements == nothing
        return false
    end

    seen2Cone = false
    eulerChar = div(length(cut) - length(elements), 2)

    for i in 0 : dim(ds) - 1
        for j in i + 1 : dim(ds)
            seen = fill(false, size(ds))

            if i == 1 || j == 1
                for D in cut
                    if !seen[D]
                        E = D
                        k = i + j - 1
                        seen[E] = true

                        while !(k == 1 && E in cut)
                            E = get(ds, k, E)
                            k = i + j - k
                            seen[E] = true
                        end
                    end
                end
            end

            for D in elements
                if !seen[D]
                    eulerChar += 1
                    E = D
                    while !seen[E]
                        seen[E] = seen[get(ds, i, E)] = true
                        E = get(ds, j, get(ds, i, E))
                    end

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

    return eulerChar == 1
end


function patchElements(ds::AbstractDelaneySymbol, cut::Vector{Int64})
    inPatch = fill(false, size(ds))
    elements = fill(0, size(ds))

    seed = cut[1]
    inPatch[seed] = true
    elements[1] = seed
    next = nrElements = 1

    while next <= nrElements
        D = elements[next]
        next += 1

        for i in 0 : dim(ds)
            Di = get(ds, i, D)
            if (i == 1 && D in cut)
                continue
            elseif (i == 1 && Di in cut)
                return nothing
            elseif !inPatch[Di]
                inPatch[Di] = true
                nrElements += 1
                elements[nrElements] = Di
            end
        end
    end

    return elements[1 : nrElements]
end
