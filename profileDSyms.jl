using Profile

include("dsetGenerator.jl")
include("dsymGenerator.jl")


function gen(n)
    for (count1, dset) in enumerate(DSetGenerator(2, n))
        orbs = orbits(dset)
        vs = map(minV, orbs)
        curv = curvature(dset, orbs, vs)

        if curv < 0
            #println(NumberedDSym(dset, vs, count1, 1))
        else
            for (count2, ds) in enumerate(DSymGenerator(dset))
                #println(NumberedDSym(ds, count1, count2))
            end
        end
    end
end


println("Warmup run...")
gen(12)

Profile.init(n = 10000000)
Profile.clear()

println("Profile run...")
@profile gen(16)

println("Profiling result:")
Profile.print(format=:tree, combine=true, mincount=100)
#Profile.print(format=:flat, combine=true, sortedby=:count)
