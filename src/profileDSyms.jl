using Profile

include("dsetGenerator.jl")
include("dsymGenerator.jl")
include("properties.jl")


function gen(n)
    for (count1, dset) in enumerate(DSetGenerator(2, n))
        orbs = vcat(orbits(dset, 0, 1), orbits(dset, 1, 2))
        vs = map(minV, orbs)
        curv = curvature(dset, orbs, vs)

        if curv < 0
            isPseudoConvex(DelaneySymbol(dset, vs))
            #println(NumberedDelaneySymbol(DelaneySymbol(dset, vs), count1, 1))
        else
            for (count2, ds) in enumerate(DSymGenerator(dset))
                isPseudoConvex(ds)
                #println(NumberedDelaneySymbol(ds, count1, count2))
            end
        end
    end
end


println("Warmup run...")
gen(8)

Profile.init(n = 10000000)
Profile.clear()

println("Profile run...")
@profile gen(12)

println("Profiling result:")
Profile.print(format=:tree, combine=true, mincount=200)
#Profile.print(format=:flat, combine=true, sortedby=:count)
