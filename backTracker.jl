abstract type BackTracker{R, S} end


function root(bt::BackTracker{R, S})::S where {R, S}
    return S()
end


function extract(bt::BackTracker{R, S}, state::S)::R where {R, S}
    return R()
end


function children(bt::BackTracker{R, S}, state::S)::Vector{S} where {R, S}
    return []
end


function Base.iterate(
    bt::BackTracker{R, S},
    stack::Vector{Vector{S}}=[[root(bt)]]
)::Union{Nothing, Tuple{R, Vector{Vector{S}}}} where {R, S}

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
