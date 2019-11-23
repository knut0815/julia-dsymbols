include("backTracker.jl")
include("dsets.jl")


struct DSetGenState
    dset::DSet
    isRemapStart::BitVector
end


struct DSetGenerator <: BackTracker{DSet, DSetGenState}
    dim::Int64
    maxSize::Int64
end


extract(g::DSetGenerator, st::DSetGenState) =
    firstUndefined(st.dset) == nothing ? st.dset : nothing

root(g::DSetGenerator) =
    DSetGenState(DSet(zeros(Int64, 1, g.dim + 1)), falses(g.maxSize))


function children(g::DSetGenerator, st::DSetGenState)
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
                    push!(result, DSetGenState(dset, isRemapStart))
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
