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
    return false
end
