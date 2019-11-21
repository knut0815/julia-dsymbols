"""
    abstract type BackTracker{R, S}

Defines a generic interface for enumeration algorithms.

The approach is to traverse an implicit tree with an instance of type `S`
(the node state) assigned to each node, such that the current node's
children and their states only depend on the current node's state.  The
result type `R`, which can be identical to `S` or different from it,
determines the type of objects to be enumerated.

The file `partitions.jl` contains a simple example in the form of a
backtracker that enumerates integer partitions.
"""
abstract type BackTracker{R, S} end


"""
    extract(bt::BackTracker{R, S}, st::S)::R

Return the result for the state `st`, if available, otherwise nothing.

The state `st` in the context of the backtracker `bt` may correspond to a
complete enumeration result, which is then returned, or a partial result, in
which case `nothing` is returned.
"""
function extract(bt::BackTracker{R, S}, st::S)::R where {R, S}
    return R()
end


"""
The root() function returns the root of the enumeration tree for the
backtracker instance defined by the argument bt.
"""
function root(bt::BackTracker{R, S})::S where {R, S}
    return S()
end


"""
The children() function returns a list of states defining the children of
the node with state st in the context of the backtracker bt.
"""
function children(bt::BackTracker{R, S}, st::S)::Vector{S} where {R, S}
    return []
end


"""
The following is a basic implementation of Julia's Iteration interface for
backtrackers.  This enables one to write code like

    for result in bt
        println(result)
    end

As by Julia's conventions, the iterate() function performs a single
iteration step, taking as its first argument an object representing the
collection to iterate over, and as its optional second argument a state
object representing the current state of the iteration process.  It returns
either a pair consisting of a result value and a new state, or nothing if
the end of the iteration has been reached.
"""
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
