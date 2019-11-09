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
    for D in 1 : size(ds)
        if compareRenumberedFrom(ds, D) < 0
            return false
        end
    end

    return true
end


function compareRenumberedFrom(ds::DSet, D0::Int)
    n2o = zeros(Int, size(ds))
    o2n = zeros(Int, size(ds))

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

            nval = E == 0 ? 0 : o2n[E]
            oval = get(ds, i, D)
            if oval != nval
                return oval == 0 ? -1 : nval == 0 ? 1 : nval - oval
            end
        end
    end

    return 0
end


function Base.show(io::IO, ds::DSet)
    print(io, "<1.1:", size(ds))
    if dim(ds) != 2
        print(io, " ", dim(ds))
    end
    print(io, ":")

    for i in 0 : dim(ds)
        if i > 0
            print(",")
        end
        for D in 1 : size(ds)
            E = get(ds, i, D)
            if E == 0 || E >= D
                if D > 1
                    print(io, " ")
                end
                print(io, E)
            end
        end
    end
    print(io, ":")

    for i in 0 : dim(ds) - 1
        if i > 0
            print(",")
        end
        j = i + 1
        seen = falses(size(ds))

        for D in 1 : size(ds)
            if !seen[D]
                if D > 1
                    print(io, " ")
                end
                E, k = D, 0
                while true
                    Ei = get(ds, i, E)
                    E = Ei == 0 ? E : Ei
                    seen[E] = true
                    Ej = get(ds, j, E)
                    E = Ej == 0 ? E : Ej
                    seen[E] = true
                    k += 1

                    if E == D
                        break
                    end
                end

                print(io, k * Int(ceil(3 / k)))
            end
        end
    end

    print(io, ">")
end


for ds in DSetGenerator(2, parse(Int, ARGS[1]))
    println(ds)
end
