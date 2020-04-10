include("dsyms.jl")


struct DoubleStep
    ds::AbstractDelaneySymbol
    i::Int64
    j::Int64
    start::Int64
end


function Base.iterate(spec::DoubleStep, state=nothing)
    if state == nothing
        return (spec.start, spec.start)
    else
        next = get(spec.ds, spec.i, get(spec.ds, spec.j, state))
        if next == spec.start
            return nothing
        else
            return (next, next)
        end
    end
end


function isPseudoConvex(dsRaw::AbstractDelaneySymbol)
    ds = orientedCover(dsRaw)
    ori = partialOrientation(ds)

    for A1 in filter(D -> ori[D] > 0, 1 : size(ds))
        seen1 = falses(size(ds))
        seen1[A1] = true

        for A2 in DoubleStep(ds, 0, 1, get(ds, 0, A1))
            if seen1[A2]
                break
            end

            seen2 = copy(seen1)
            seen2[A2] = seen1[A2] = seen1[get(ds, 1, A2)] = true

            for B2 in DoubleStep(ds, 2, 1, get(ds, 2, A2))
                if seen2[B2]
                    if B2 == A1 && cutsOffDisk(ds, [A1, A2], true)
                        return false
                    else
                        break
                    end
                end

                seen3 = copy(seen2)
                seen3[B2] = seen2[B2] = seen2[get(ds, 1, B2)] = true

                for B1 in DoubleStep(ds, 0, 1, get(ds, 0, B2))
                    if seen3[B1]
                        break
                    end

                    seen3[B1] = seen3[get(ds, 1, B1)] = true
                    seen4 = copy(seen3)

                    for T in DoubleStep(ds, 2, 1, get(ds, 2, B1))
                        if seen4[T]
                            if T == A1 && cutsOffDisk(ds, [A1, A2, B2, B1])
                                return false
                            else
                                break
                            end
                        end

                        seen4[T] = seen4[get(ds, 1, T)] = true
                    end
                end
            end
        end
    end

    return true
end


function cutsOffDisk(
    ds::AbstractDelaneySymbol, cut::Vector{Int64}, allow2Cone::Bool=false
)
    checkCones(cones) = cones == [] || (allow2Cone && cones == [2])

    patch, rest = splitAlong(ds, cut)

    if size(patch) == length(cut)
        return false
    end

    if eulerCharacteristic(ds) > 0
        if size(patch) == size(ds)
            vs = [v(ds, 0, 1, cut[1]), v(ds, 1, 2, cut[1])]
            if length(cut) > 2
                push!(vs, v(ds, 1, 2, cut[2]))
            end

            if checkCones(filter(v -> v > 1, vs))
                return false
            end
        end

        if (
            size(patch) == size(ds) - length(cut) &&
            all(D -> v(ds, 0, 1, D) == 1 && v(ds, 1, 2, D) == 1, cut) &&
            checkCones(coneDegrees(rest))
        )
            return false
        end
    end

    return (
        isWeaklyOriented(patch) &&
        eulerCharacteristic(patch) == 1 &&
        checkCones(coneDegrees(patch))
    )
end


function splitAlong(ds::AbstractDelaneySymbol, cut::Vector{Int64})
    inCut = falses(size(ds))
    for D in cut
        inCut[D] = inCut[get(ds, 1, D)] = true
    end

    result = []
    for seed in [first(cut), get(ds, 1, first(cut))]
        src2img = zeros(Int64, size(ds))
        img2src = zeros(Int64, size(ds))
        count = 0
        queue = [seed]

        while length(queue) > 0
            D = popfirst!(queue)
            if src2img[D] == 0
                count += 1
                src2img[D] = count
                img2src[count] = D

                for i in 0 : dim(ds)
                    if i != 1 || !inCut[D]
                        push!(queue, get(ds, i, D))
                    end
                end
            end
        end

        dset = DelaneySetUnderConstruction(count, dim(ds))
        for D in 1 : count
            for i in 0 : dim(ds)
                E = img2src[D]
                if i == 1 && inCut[E]
                    set!(dset, i, D, D)
                else
                    set!(dset, i, D, src2img[get(ds, i, E)])
                end
            end
        end

        part = DelaneySymbolUnderConstruction(dset)
        for D in 1 : count
            for i in 1 : dim(ds)
                setV!(part, i - 1, i, D, v(ds, i - 1, i, img2src[D]))
            end
        end

        push!(result, DelaneySymbol(part))
    end

    return result
end


function eulerCharacteristic(ds::AbstractDelaneySet)
    nrLoops(i) = count(D -> get(ds, i, D) == D, 1 : size(ds))
    nrOrbits(i, j) = length(orbits(ds, i, j))

    nf = size(ds)
    ne = div(3 * nf + nrLoops(0) + nrLoops(1) + nrLoops(2), 2)
    nv = nrOrbits(0, 1) + nrOrbits(0, 2) + nrOrbits(1, 2)

    return nf - ne + nv
end


function coneDegrees(ds::AbstractDelaneySymbol)
    result = []

    for i in 0 : dim(ds) - 1
        for j in i + 1 : dim(ds)
            for orb in orbits(ds, i, j)
                vOrb = v(ds, i, j, first(orb.elements))
                if vOrb > 1 && !orb.isChain
                    push!(result, vOrb)
                end
            end
        end
    end

    return result
end
