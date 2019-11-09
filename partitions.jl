include("backTracker.jl")


struct PartitionState
    xs::Vector{Int}
    sz::Int
    mx::Int
end


struct Partitions <: BackTracker{Vector{Int}, PartitionState}
    n::Int
end


root(p::Partitions) = PartitionState([], 0, 1)

extract(p::Partitions, st::PartitionState) = st.sz == p.n ? st.xs : nothing

function children(p::Partitions, st::PartitionState)
    ch::Vector{PartitionState} = []

    for i in st.mx : p.n - st.sz
        push!(ch, PartitionState(vcat(st.xs, [i]), st.sz + i, max(st.mx, i)))
    end

    return ch
end


for p in Partitions(70)
    println(p)
end
