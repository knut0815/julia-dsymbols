"""
    abstract type BackTracker{R, S}

Defines a generic interface for tree-based enumeration algorithms.

An implicit tree with nodes labelled by instances of the parameter type `S`
(the node state) is defined by specifying the state of the root and, for any
node state appearing in the tree, the list of states of its children.  The
enumeration is performed by traversing this tree, producing extracted
elements of type `R` for nodes that correspond to results.

See also: [`extract`](@ref) [`root`](@ref) [`children`](@ref)

# Example

    # A backtracking enumerator for integer partitions.

    struct PState
        xs::Vector{Int}
        left::Int
        top::Int
    end

    struct Partitions <: BackTracker{Vector{Int}, PState}
        n::Int
    end

    extract(p::Partitions, st::PState) = st.left == 0 ? st.xs : nothing

    root(p::Partitions) = PState([], p.n, 1)

    children(p::Partitions, st::PState) = map(
        i -> PState(vcat(st.xs, [i]), st.left - i, max(st.top, i)),
        st.top : st.left
    )

    # Print the partitions of the number 10.

    for p in Partitions(10)
        println(p)
    end
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

    for result in backtracker
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
