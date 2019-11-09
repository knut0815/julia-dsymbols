abstract type BackTracker{R, S} end


function Base.iterate(
    bt::BackTracker{R, S},
    stack::Union{Nothing, Vector{Vector{S}}}=nothing
)::Union{Nothing, Tuple{R, Vector{Vector{S}}}} where {R, S}

    if stack == nothing
        stack::Vector{Vector{S}} = [[root(bt)]]
    end

    while length(stack) > 0
        current::S = last(last(stack))
        value::Union{Nothing, R} = extract(bt, current)

        next::Vector{S} = children(bt, current)
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


Base.eltype(::Type{BackTracker{R, S}}) where {R, S} = R
