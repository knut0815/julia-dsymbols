include("dsetGenerator.jl")
include("dsymGenerator.jl")
include("properties.jl")


filterPseudoConvex = false
maxSize = 6

for s in ARGS
    if s == "-c"
        global filterPseudoConvex
        filterPseudoConvex = true
    else
        n = tryparse(Int64, s)
        if n != nothing
            global maxSize
            maxSize = n
        end
    end
end


printIfGood(ds) = (!filterPseudoConvex || isPseudoConvex(ds)) && println(ds)

for (count1, dset) in enumerate(DSetGenerator(2, maxSize))
    orbs = vcat(orbits(dset, 0, 1), orbits(dset, 1, 2))
    vs = map(minV, orbs)
    curv = curvature(dset, orbs, vs)

    if curv < 0
        printIfGood(NumberedDelaneySymbol(DelaneySymbol(dset, vs), count1, 1))
    else
        for (count2, ds) in enumerate(DSymGenerator(dset))
            printIfGood(NumberedDelaneySymbol(ds, count1, count2))
        end
    end
end
