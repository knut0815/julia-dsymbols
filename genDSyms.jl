import Base: iterate

abstract type BackTracker{T} end


function iterate(bt::BackTracker{T}, stack=nothing) where T
    if stack == nothing
        stack = [[root(bt)]]
    end

    while length(stack) > 0
        current = last(last(stack))
        value = extract(bt, current)

        next = children(bt, current)
        if length(next) > 0
            push!(stack, reverse(next))
        else
            while length(stack) > 0 && length(last(stack)) < 2
                pop!(stack)
            end

            if length(stack) > 0
                pop!(last(stack))
            end
        end

        if value != nothing
            return (value, stack)
        end
    end

    return nothing
end


struct PartitionState
    xs::Vector{Int}
    sz::Int
    mx::Int
end


struct Partitions <: BackTracker{PartitionState}
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


for p in Partitions(100)
    println(p)
end
