include("dsyms.jl")


for (count1, dset) in enumerate(DSetGenerator(2, parse(Int64, ARGS[1])))
    orbs = orbits(dset)
    vs = map(minV, orbs)
    curv = curvature(dset, orbs, vs)

    if curv < 0
        println(NumberedDSym(dset, vs, count1, 1))
    else
        for (count2, ds) in enumerate(DSymGenerator(dset))
            println(NumberedDSym(ds, count1, count2))
        end
    end
end
