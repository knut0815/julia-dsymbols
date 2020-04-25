using Profile

include("dsetGenerator.jl")
include("dsymGenerator.jl")
include("properties.jl")


function gen(n)
    for (count1, dset) in enumerate(DSetGenerator(2, n))
        for (count2, ds) in enumerate(DSymGenerator(dset))
            isPseudoConvex(ds)
            #println(NumberedDelaneySymbol(ds, count1, count2))
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
Profile.print(format=:tree, combine=true, mincount=750)
#Profile.print(format=:flat, combine=true, sortedby=:count)
