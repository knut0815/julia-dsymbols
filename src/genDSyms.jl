include("dsetGenerator.jl")
include("dsymGenerator.jl")


for (count1, dset) in enumerate(DSetGenerator(2, parse(Int64, ARGS[1])))
    orbs = vcat(orbits(dset, 0, 1), orbits(dset, 1, 2))
    vs = map(minV, orbs)
    curv = curvature(dset, orbs, vs)

    if curv < 0
        println(NumberedDelaneySymbol(DelaneySymbol(dset, vs), count1, 1))
    else
        for (count2, ds) in enumerate(DSymGenerator(dset))
            println(NumberedDelaneySymbol(ds, count1, count2))
        end
    end
end
