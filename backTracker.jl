"""
    BackTracker{R, S}

Supertype for enumeration algorithms based on tree traversal.

An implicit enumeration tree with nodes labelled by instances of the
parameter type `S` (the node state) is defined by the state of the root and
a function that computes a list of child states for a given node state.  The
enumeration is performed by traversing this tree, producing extracted
elements of type `R` for only those nodes that correspond to finished
results.

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
    extract(bt::BackTracker{R, S}, st::S) -> R

Return the result for the state `st`, if available, otherwise `nothing`.

The state `st` in the context of the backtracker `bt` may correspond to a
finished enumeration result, which is then returned, or a partial result, in
which case `nothing` is returned.

See also: [`BackTracker`](@ref) [`root`](@ref) [`children`](@ref)
"""
function extract(bt::BackTracker{R, S}, st::S)::R where {R, S}
    return R()
end


"""
    root(bt::BackTracker{R, S}) -> S

Return the root state of the enumeration tree for the backtracker `bt`.

See also: [`BackTracker`](@ref) [`extract`](@ref) [`children`](@ref)
"""
function root(bt::BackTracker{R, S})::S where {R, S}
    return S()
end


"""
    children(bt::BackTracker{R, S}, st::S) -> Vector{S}

Return the list of child states for the node state `st`.

See also: [`BackTracker`](@ref) [`extract`](@ref) [`root`](@ref)
"""
function children(bt::BackTracker{R, S}, st::S)::Vector{S} where {R, S}
    return []
end


"""
    iterate(bt::BackTracker{R, S} [, stack::Vector{Vector{S}}]) ->
        Union{Nothing, Tuple{R, Vector{Vector{S}}}}

Provides a basic backtracker implementation of Julia's iteration protocol.

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
