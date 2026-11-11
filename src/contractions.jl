"""Dense network tensor with explicit integer half-edge labels; contractions are bilinear."""
struct NetworkTensor{T}
    data::Array{T}
    labels::Vector{Int}
    function NetworkTensor(data::Array{T},labels::Vector{Int}) where T
        ndims(data)==length(labels) || throw(DimensionMismatch("tensor labels"))
        length(unique(labels))==length(labels) || throw(ArgumentError("repeated tensor label"))
        new{T}(data,labels)
    end
end

"""Contract every shared index using ordinary matrix multiplication, never adjoints."""
function contract_pair(a::NetworkTensor{T},b::NetworkTensor{T}) where T
    shared=intersect(a.labels,b.labels)
    afree=findall(x->x ∉ shared,a.labels)
    bfree=findall(x->x ∉ shared,b.labels)
    ai=[something(findfirst(==(x),a.labels)) for x in shared]
    bi=[something(findfirst(==(x),b.labels)) for x in shared]
    all(i->size(a.data,ai[i])==size(b.data,bi[i]),eachindex(shared)) ||
        throw(DimensionMismatch("contracted edge dimension"))
    adims=Tuple(size(a.data,i) for i in afree)
    bdims=Tuple(size(b.data,i) for i in bfree)
    inner=prod(size(a.data,i) for i in ai;init=1)
    am=reshape(permutedims(a.data,(afree...,ai...)),prod(adims),inner)
    bm=reshape(permutedims(b.data,(bi...,bfree...)),inner,prod(bdims))
    data=reshape(am*bm,(adims...,bdims...))
    NetworkTensor(data,[a.labels[afree];b.labels[bfree]])
end

"""Greedy contraction cost in log units, avoiding integer overflow for large networks."""
function contraction_cost(a::NetworkTensor,b::NetworkTensor)
    isempty(intersect(a.labels,b.labels)) && !isempty(a.labels) && !isempty(b.labels) && return Inf
    sum(log(max(1,size(a.data,i))) for i in eachindex(a.labels);init=0.0)+
    sum(log(max(1,size(b.data,i))) for i in eachindex(b.labels) if b.labels[i] ∉ a.labels;init=0.0)
end

"""Contract a closed network, selecting small shared-index contractions first."""
function contract_network(tensors::Vector{NetworkTensor{T}}) where T
    isempty(tensors) && return one(T)
    work=copy(tensors)
    while length(work)>1
        best=(1,2)
        cost=Inf
        for i in eachindex(work),j in i+1:length(work)
            candidate=contraction_cost(work[i],work[j])
            if candidate<cost
                cost,best=candidate,(i,j)
            end
        end
        i,j=best
        product=contract_pair(work[i],work[j])
        deleteat!(work,j)
        work[i]=product
    end
    isempty(only(work).labels) || throw(ArgumentError("network contains open indices"))
    only(only(work).data)
end

"""Sew ordered three-point tensors with full supplied inverse Grams.

Selections restrict both ends of each edge to the same parity sector. External
slots have dimension one. This generic contraction is also used by independent
PBW checks; it contains no recursion or branching data.
"""
function sew_tensors(graph::PlumbingGraph,tensors::Vector{Array{T,3}},
                     inverses::AbstractVector,selections::AbstractVector) where T
    length(tensors)==length(graph.slots) || throw(DimensionMismatch("vertex count"))
    length(inverses)==length(graph.ends)==length(selections) || throw(DimensionMismatch("edge count"))
    any(isempty,selections) && return zero(T)
    network=NetworkTensor{T}[]
    for (v,slots) in enumerate(graph.slots)
        indices=ntuple(j->slots[j]==0 ? [1] : selections[slots[j]],3)
        localdata=tensors[v][indices...]
        active=findall(>(0),slots)
        data=reshape(localdata,Tuple(size(localdata,j) for j in active))
        push!(network,NetworkTensor(data,[3(v-1)+j for j in active]))
    end
    for (e,ends) in enumerate(graph.ends)
        push!(network,NetworkTensor(Matrix{T}(inverses[e][selections[e],selections[e]]),collect(ends)))
    end
    contract_network(network)
end

"""Physical PBW coefficient in all edge-parity sectors, using independent local Ward forms."""
function pbw_coefficient(graph::PlumbingGraph,vertices::AbstractVector{<:WardVertex{T}},
                         levels::NTuple{N,Int},forms::AbstractVector{<:Integer},
                         eta::AbstractVector{<:Integer}=ones(Int,length(vertices))) where {T,N}
    N==length(graph.ends) || throw(DimensionMismatch("edge levels"))
    length(vertices)==length(graph.slots)==length(forms)==length(eta) || throw(DimensionMismatch("vertex data"))
    all(>(0),Iterators.flatten(graph.slots)) || throw(ArgumentError("physical closed-graph PBW requires internal slots"))
    modules=[vertices[div(graph.ends[e][1]-1,3)+1].modules[mod(graph.ends[e][1]-1,3)+1] for e in 1:N]
    for (e,ends) in enumerate(graph.ends)
        other=vertices[div(ends[2]-1,3)+1].modules[mod(ends[2]-1,3)+1]
        m=modules[e]
        (m.ramond,m.c,m.h,m.beta,m.primary_parity)==
            (other.ramond,other.c,other.h,other.beta,other.primary_parity) ||
            throw(ArgumentError("module mismatch across edge $e"))
    end
    bases=[basis(modules[e],levels[e]) for e in 1:N]
    out=zeros(Complex{T},2^N)
    any(isempty,bases) && return out
    tensors=[vertex_tensor(vertices[v],ntuple(j->levels[graph.slots[v][j]],3),forms[v],eta[v]) for v in eachindex(vertices)]
    inverses=[inverse_gram(modules[e],levels[e]) for e in 1:N]
    selections=[ntuple(p->findall(s->Int(state_parity(modules[e],s))==p-1,bases[e]),2) for e in 1:N]
    for mask in 0:2^N-1
        selected=[selections[e][((mask>>(e-1))&1)+1] for e in 1:N]
        out[mask+1]=sewing_sign(graph,mask)*sew_tensors(graph,tensors,inverses,selected)
    end
    out
end

"""Evaluate parity coefficients at fixed edge spin signs, including intrinsic NS parity."""
function evaluate_spins(coefficients::AbstractVector{T},signs::AbstractVector{<:Integer}) where T
    length(coefficients)==2^length(signs) || throw(DimensionMismatch("parity coefficients"))
    all(s->s in (-1,1),signs) || throw(ArgumentError("spin signs must be ±1"))
    out=zero(T)
    for mask in 0:length(coefficients)-1
        sign=prod(signs[e] for e in eachindex(signs) if isodd(mask>>(e-1));init=1)
        out+=sign*coefficients[mask+1]
    end
    out
end
