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


for (count1, dset) in enumerate(DSetGenerator(2, maxSize))
    for (count2, ds) in enumerate(DSymGenerator(dset))
        if !filterPseudoConvex || !isPseudoConvex(ds)
            println(NumberedDelaneySymbol(ds, count1, count2))
        end
    end
end
