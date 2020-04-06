include("dsets.jl")


abstract type AbstractDelaneySymbol <: AbstractDelaneySet end

setCount(ds::AbstractDelaneySymbol) = 1

symbolCount(ds::AbstractDelaneySymbol) = 1



function collectOrbits(ds::DelaneySet)
    allOrbits::Vector{Orbit} = []
    orbitIndex = zeros(Int64, dim(ds), size(ds))

    for i in 1 : dim(ds)
        for orb in orbits(ds, i - 1, i)
            push!(allOrbits, orb)
            for D in orb.elements
                orbitIndex[i, D] = length(allOrbits)
            end
        end
    end

    return (allOrbits, orbitIndex)
end



struct DelaneySymbolUnderConstruction <: AbstractDelaneySymbol
    dset::DelaneySet
    orbits::Vector{Orbit}
    orbitIndex::Array{Int64, 2}
    vs::Vector{Int64}

    function DelaneySymbolUnderConstruction(dset::DelaneySet)
        (allOrbits, orbitIndex) = collectOrbits(dset)
        new(dset, allOrbits, orbitIndex, zeros(Int64, length(allOrbits)))
    end
end


Base.size(ds::DelaneySymbolUnderConstruction) = size(ds.dset)

dim(ds::DelaneySymbolUnderConstruction) = dim(ds.dset)

get(ds::DelaneySymbolUnderConstruction, i::Int64, D::Int64) =
    get(ds.dset, i, D)



struct DelaneySymbol <: AbstractDelaneySymbol
    dset::DelaneySet
    vs::Vector{Int64}
end


Base.size(ds::DelaneySymbol) = size(ds.dset)

dim(ds::DelaneySymbol) = dim(ds.dset)

get(ds::DelaneySymbol, i::Int64, D::Int64) = get(ds.dset, i, D)



struct NumberedDelaneySymbol <: AbstractDelaneySymbol
    dsym::DelaneySymbol
    setCount::Int64
    symbolCount::Int64
end


Base.size(ds::NumberedDelaneySymbol) = size(ds.dsym)

dim(ds::NumberedDelaneySymbol) = dim(ds.dsym)

get(ds::NumberedDelaneySymbol, i::Int64, D::Int64) = get(ds.dsym, i, D)

setCount(ds::NumberedDelaneySymbol) = ds.setCount

symbolCount(ds::NumberedDelaneySymbol) = ds.symbolCount



function curvature(
    ds::AbstractDelaneySet, orbs::Vector{Orbit}, vs::Vector{Int64}
)
    result = -size(ds)//2

    for i in 1 : length(orbs)
        result += (orbs[i].isChain ? 1 : 2) // vs[i]
    end

    return result
end


function Base.show(io::IO, ds::AbstractDelaneySymbol)
    print(io, "<$(setCount(ds)).$(symbolCount(ds)):$(size(ds))")
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

    orbs = vcat(orbits(ds, 0, 1), orbits(ds, 1, 2))

    for i in 0 : dim(ds) - 1
        if i > 0
            print(",")
        end

        for k in 1 : length(orbs)
            if orbs[k].indices == [i, i + 1]
                if first(orbs[k].elements) > 1
                    print(io, " ")
                end
                print(io, r(orbs[k]) * ds.dsym.vs[k])
            end
        end
    end

    print(io, ">")
end
