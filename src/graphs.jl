"""Ordered trivalent plumbing graph. Slot 1/2/3 is infinity/one/zero.

Positive slot entries are internal edges; zero is an external puncture.
Half-edge IDs are 3(vertex-1)+slot. Edge endpoints retain their sewing order.
"""
struct PlumbingGraph
    slots::Vector{NTuple{3,Int}}
    ends::Vector{NTuple{2,Int}}
    inversions::Vector{NTuple{2,Int}}
end

"""Validate graph incidence and precompute the Koszul inversion pairs."""
function PlumbingGraph(slots::Vector{NTuple{3,Int}};
                       ends::Union{Nothing,Vector{NTuple{2,Int}}}=nothing)
    isempty(slots) && throw(ArgumentError("graph has no vertices"))
    all(e->e>=0, Iterators.flatten(slots)) || throw(ArgumentError("negative edge ID"))
    ne = maximum(maximum, slots)
    occurrences = [Int[] for _ in 1:ne]
    for (v, vertex) in enumerate(slots), j in 1:3
        vertex[j] == 0 || push!(occurrences[vertex[j]], 3(v-1)+j)
    end
    all(x->length(x)==2, occurrences) || throw(ArgumentError("each internal edge needs two ends"))
    pairs = isnothing(ends) ? [Tuple(x) for x in occurrences] : copy(ends)
    length(pairs) == ne || throw(ArgumentError("wrong number of endpoint pairs"))
    for e in 1:ne
        sort(collect(pairs[e])) == occurrences[e] || throw(ArgumentError("endpoint mismatch"))
    end
    order = collect(Iterators.flatten(pairs))
    edgeof(h) = slots[div(h-1,3)+1][mod(h-1,3)+1]
    inversions = NTuple{2,Int}[]
    for i in eachindex(order), j in i+1:length(order)
        order[i] > order[j] && push!(inversions, (edgeof(order[i]), edgeof(order[j])))
    end
    PlumbingGraph(copy(slots), pairs, inversions)
end

theta_graph() = PlumbingGraph([(1,2,3), (1,2,3)])
glasses_graph() = PlumbingGraph([(1,2,2), (1,3,3)])
tetrahedron_graph() = PlumbingGraph([(1,2,3), (1,4,5), (2,4,6), (3,5,6)])

"""Koszul sign from edge-pair order to ordered vertex slots."""
function sewing_sign(graph::PlumbingGraph, mask::Integer)
    exponent = 0
    for (e,f) in graph.inversions
        exponent ⊻= Int(((mask >> (e-1)) & 1) & ((mask >> (f-1)) & 1))
    end
    parity_sign(exponent)
end

"""Split marked edges with a new (right,external,left) insertion vertex.

Returns the expanded graph and original-edge => (left,right) edge IDs.
"""
function split_graph(graph::PlumbingGraph, marked::AbstractVector{<:Integer})
    length(unique(marked)) == length(marked) || throw(ArgumentError("edge split twice"))
    slots, ends = copy(graph.slots), copy(graph.ends)
    splits = Dict{Int,NTuple{2,Int}}()
    for e in marked
        1 <= e <= length(graph.ends) || throw(ArgumentError("invalid split edge"))
        right, v = length(ends)+1, length(slots)+1
        lh, rh = ends[e]
        rv, rs = div(rh-1,3)+1, mod(rh-1,3)+1
        slots[rv] = Base.setindex(slots[rv], right, rs)
        push!(slots, (right, 0, e))
        ends[e] = (lh, 3(v-1)+3)
        push!(ends, (3(v-1)+1, rh))
        splits[e] = (e, right)
    end
    PlumbingGraph(slots; ends=ends), splits
end

"""Return Ramond cycles as (edge IDs, vertex IDs); validate NSRR slot order."""
function ramond_cycles(graph::PlumbingGraph, ramond::AbstractVector{Bool})
    length(ramond) == length(graph.ends) || throw(DimensionMismatch("sector count"))
    byvertex = [Int[] for _ in graph.slots]
    for (v, slots) in enumerate(graph.slots)
        for j in 1:3
            e = slots[j]
            e > 0 && ramond[e] && push!(byvertex[v], e)
        end
        length(byvertex[v]) in (0,2) || throw(ArgumentError("vertex is not NSNSNS/NSRR"))
        isempty(byvertex[v]) || (slots[1] == 0 || !ramond[slots[1]]) ||
            throw(ArgumentError("NSRR vertex requires NS first slot"))
    end
    cycles = Tuple{Vector{Int},Vector{Int}}[]
    unseen = Set(findall(ramond))
    while !isempty(unseen)
        edges, vertices = Int[], Set{Int}()
        pending = [minimum(unseen)]
        while !isempty(pending)
            e = pop!(pending)
            e in unseen || continue
            delete!(unseen, e)
            push!(edges, e)
            for h in graph.ends[e]
                v = div(h-1,3)+1
                push!(vertices, v)
                append!(pending, byvertex[v])
            end
        end
        push!(cycles, (sort(edges), sort(collect(vertices))))
    end
    cycles
end

"""Mixed Koszul sign for closed ordinary auxiliary/physical convolution.

This includes full tensor-Gram and both ordered-vertex signs. Split insertion
graphs require their external operator parity and must be assembled separately.
"""
function mixed_convolution_sign(graph::PlumbingGraph,physical::Integer,auxiliary::Integer)
    all(>(0),Iterators.flatten(graph.slots)) ||
        throw(ArgumentError("external insertions require the split-graph parity rule"))
    exponent=0
    for e in eachindex(graph.ends)
        exponent+=Int(((physical>>(e-1))&1)*((auxiliary>>(e-1))&1))
    end
    for slots in graph.slots
        exponent+=Int(((physical>>(slots[1]-1))&1)*((auxiliary>>(slots[1]-1))&1))
        exponent+=Int(((physical>>(slots[2]-1))&1)*((auxiliary>>(slots[3]-1))&1))
    end
    sewing_sign(graph,xor(physical,auxiliary))*sewing_sign(graph,physical)*
        sewing_sign(graph,auxiliary)*parity_sign(exponent)
end

"""Choose one deterministic edge on each cycle with product of vertex eta=-1."""
function insertion_edges(graph::PlumbingGraph, ramond::AbstractVector{Bool},
                         eta::AbstractVector{<:Integer})
    length(eta) == length(graph.slots) || throw(DimensionMismatch("vertex signs"))
    all(x->x in (-1,1), eta) || throw(ArgumentError("eta must be ±1"))
    [first(edges) for (edges,vertices) in ramond_cycles(graph,ramond)
                  if prod(eta[v] for v in vertices) == -1]
end
