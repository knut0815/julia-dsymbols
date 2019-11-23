include("dsets.jl")


struct DSym
    dset::DSet
    vs::Vector{Int64}
end


Base.size(ds::DSym) = size(ds.dset)

dim(ds::DSym) = dim(ds.dset)

get(ds::DSym, i::Int64, D::Int64) = get(ds.dset, i, D)

orbits(ds::DSym) = orbits(ds.dset)


struct NumberedDSym
    dsym::DSym
    count1::Int64
    count2::Int64
end


NumberedDSym(dset::DSet, vs::Vector{Int64}, count1::Int64, count2::Int64) =
    NumberedDSym(DSym(dset, vs), count1, count2)

Base.size(ds::NumberedDSym) = size(ds.dsym)

dim(ds::NumberedDSym) = dim(ds.dsym)

get(ds::NumberedDSym, i::Int64, D::Int64) = get(ds.dsym, i, D)

orbits(ds::NumberedDSym) = orbits(ds.dsym)


function curvature(ds::DSet, orbs::Vector{Orbit}, vs::Vector{Int64})
    result = -size(ds)//2

    for i in 1 : length(orbs)
        result += (orbs[i].isChain ? 1 : 2) // vs[i]
    end

    return result
end


function Base.show(io::IO, ds::NumberedDSym)
    print(io, "<$(ds.count1).$(ds.count2):$(size(ds))")
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

    orbs = orbits(ds)

    for i in 0 : dim(ds) - 1
        if i > 0
            print(",")
        end

        for k in 1 : length(orbs)
            if orbs[k].index == i
                if first(orbs[k].elements) > 1
                    print(io, " ")
                end
                print(io, r(orbs[k]) * ds.dsym.vs[k])
            end
        end
    end

    print(io, ">")
end
