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


struct DSetGenerator <: BackTracker{DSet, DSet}
    dim::Int
    maxSize::Int
end


root(g::DSetGenerator) =
    DSet(zeros(1, g.dim + 1))

extract(g::DSetGenerator, ds::DSet) =
    firstUndefined(ds) == nothing ? nothing : ds


function children(g::DSetGenerator, ds::DSet)
    result = []
    undef = firstUndefined(ds)

    if undef != nothing
        D, i = undef

        for E in D : min(size(ds) + 1, g.maxSize)
            if get(ds, i, E) == 0
                if E > size(ds)
                    out = DSet(vcat(ds.op, zeros(1, dim(ds) + 1)))
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
                    result.push(out)
                end
            end
        end
    end

    return result
end


function firstUndefined(ds::DSet)
    for D in 1 : size(ds)
        for i in 0 : dim(ds)
            if get(ds, i, D) != 0
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
    E, k = D, 0

    while k < limit && get(ds, w[k], E) != 0
        E = get(ds, w[k], E)
        k += 1
    end

    return E, k
end


function isCanonical(ds)
    for D in 1 : size(ds)
        if compareRenumberedFrom(ds, D) < 0
            return false
        end
    end

    return true
end


function compareRenumberedFrom(ds, D0)
    n2o = zeros(size(ds))
    o2n = zeros(size(ds))

    n2o[1] = D0
    o2n[D0] = 1
    next = 2

    for D in 1 : size(ds)
        for i in 0 : dim(ds)
            E = get(ds, i, n2o[D])
            if E != 0 && o2n[E] == 0
                o2n[E] = next
                n2o[next] = E
                next += 1
            end

            nval = o2n[E]
            oval = get(ds, i, D)
            if oval != nval
                return oval == 0 ? -1 : nval == 0 ? 1 : nval - oval
            end
        end
    end

    return 0
end


for ds in DSetGenerator(2, 4)
    println(ds)
end
