include("dsets.jl")


abstract type AbstractDelaneySymbol <: AbstractDelaneySet end

setCount(ds::AbstractDelaneySymbol) = 1

symbolCount(ds::AbstractDelaneySymbol) = 1

m(ds::AbstractDelaneySymbol, i::Int64, j::Int64, D::Int64) =
    r(ds, i, j, D) * v(ds, i, j, D)



function collectOrbits(ds::AbstractDelaneySet)
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



struct DelaneySymbol <: AbstractDelaneySymbol
    dset::DelaneySet
    orbits::Vector{Orbit}
    orbitIndex::Array{Int64, 2}
    vs::Vector{Int64}

    function DelaneySymbol(dset::DelaneySet)
        (allOrbits, orbitIndex) = collectOrbits(dset)
        new(dset, allOrbits, orbitIndex, zeros(Int64, length(allOrbits)))
    end

    function DelaneySymbol(dset::DelaneySet, vs::Vector{Int64})
        (allOrbits, orbitIndex) = collectOrbits(dset)
        new(dset, allOrbits, orbitIndex, copy(vs))
    end
end


Base.size(ds::DelaneySymbol) = size(ds.dset)

dim(ds::DelaneySymbol) = dim(ds.dset)

get(ds::DelaneySymbol, i::Int64, D::Int64) = get(ds.dset, i, D)


function v(ds::DelaneySymbol, i::Int64, j::Int64, D::Int64)
    if !(1 <= D <= size(ds))
        return 0
    elseif j == i + 1
        return ds.vs[ds.orbitIndex[j, D]]
    elseif i == j + 1
        return ds.vs[ds.orbitIndex[i, D]]
    elseif j != i && get(ds, i, D) == get(ds, j, D)
        return 2
    else
        return 1
    end
end


function r(ds::DelaneySymbol, i::Int64, j::Int64, D::Int64)
    if !(1 <= D <= size(ds))
        return 0
    elseif j == i + 1
        return r(ds.orbits[ds.orbitIndex[j, D]])
    elseif i == j + 1
        return r(ds.orbits[ds.orbitIndex[i, D]])
    elseif j != i && get(ds, i, D) == get(ds, j, D)
        return 1
    else
        return 2
    end
end



struct DelaneySymbolUnderConstruction <: AbstractDelaneySymbol
    ds::DelaneySymbol

    function DelaneySymbolUnderConstruction(ds::DelaneySymbol)
        new(DelaneySymbol(ds.dset, copy(ds.vs)))
    end
end

DelaneySymbolUnderConstruction(dset::DelaneySet) =
    DelaneySymbolUnderConstruction(DelaneySymbol(dset))

DelaneySymbolUnderConstruction(dset::DelaneySet, vs::Vector{Int64}) =
    DelaneySymbolUnderConstruction(DelaneySymbol(dset, vs))

DelaneySymbol(ds::DelaneySymbolUnderConstruction) = ds.ds


Base.size(ds::DelaneySymbolUnderConstruction) = size(ds.ds)

dim(ds::DelaneySymbolUnderConstruction) = dim(ds.ds)

get(ds::DelaneySymbolUnderConstruction, i::Int64, D::Int64) = get(ds.ds, i, D)

v(ds::DelaneySymbolUnderConstruction, i::Int64, j::Int64, D::Int64) =
    v(ds.ds, i, j, D)

r(ds::DelaneySymbolUnderConstruction, i::Int64, j::Int64, D::Int64) =
    r(ds.ds, i, j, D)


function setV!(
    ds::DelaneySymbolUnderConstruction, i::Int64, j::Int64, D::Int64, v::Int64
)
    if j == i + 1
        ds.ds.vs[ds.ds.orbitIndex[j, D]] = v
    elseif i == j + 1
        ds.ds.vs[ds.ds.orbitIndex[i, D]] = v
    end
end



struct NumberedDelaneySymbol <: AbstractDelaneySymbol
    ds::AbstractDelaneySymbol
    setCount::Int64
    symbolCount::Int64
end


Base.size(ds::NumberedDelaneySymbol) = size(ds.ds)

dim(ds::NumberedDelaneySymbol) = dim(ds.ds)

get(ds::NumberedDelaneySymbol, i::Int64, D::Int64) = get(ds.ds, i, D)

v(ds::NumberedDelaneySymbol, i::Int64, j::Int64, D::Int64) =
    v(ds.ds, i, j, D)

r(ds::NumberedDelaneySymbol, i::Int64, j::Int64, D::Int64) =
    r(ds.ds, i, j, D)

setCount(ds::NumberedDelaneySymbol) = ds.setCount

symbolCount(ds::NumberedDelaneySymbol) = ds.symbolCount



function orientedCover(ds::AbstractDelaneySymbol)
    if isOriented(ds)
        return ds
    else
        dset = invoke(orientedCover, Tuple{AbstractDelaneySet}, ds)
        cov = DelaneySymbolUnderConstruction(dset)

        for i in 1 : dim(ds)
            for orb in orbits(cov, i - 1, i)
                D = first(orb.elements)
                E = (D - 1) % size(ds) + 1
                vD = div(m(ds, i - 1, i, E), r(cov, i - 1, i, D))
                setV!(cov, i - 1, i, D, vD)
            end
        end

        return DelaneySymbol(cov)
    end
end
