"""Independent finite-c Virasoro reference, using only L modes and direct Gram sewing.

No CCY pole, fusion, global-seed, Schottky or production coefficients enter this
calculation. It shares the elementary SCA L-mode algebra and generic contraction.
"""
module VirasoroPBWReference
using LinearAlgebra
using SuperconformalBlocks

"""Reference Virasoro graph with separate local descendant and Gram caches."""
mutable struct Reference{T<:AbstractFloat}
    graph::PlumbingGraph
    modules::Vector{SCAModule{T}}
    vertices::Vector{WardVertex{T,1}}
    bases::Dict{Int,Vector{PBWState}}
    inverses::Dict{Tuple{Int,Int},Matrix{Complex{T}}}
    tensors::Dict{NTuple{4,Int},Array{Complex{T},3}}
end

"""Construct direct finite-c data, with per-puncture external weights."""
function Reference(::Type{T},graph::PlumbingGraph,c,weights;external=Dict()) where T<:AbstractFloat
    plan=Combinatorics()
    modules=[SCAModule(T,:NS,c;h=h,plan=plan) for h in weights]
    vertices=WardVertex{T,1}[]
    for (v,slots) in enumerate(graph.slots)
        localmodules=ntuple(j->slots[j]>0 ? modules[slots[j]] :
            SCAModule(T,:NS,c;h=external[3(v-1)+j],plan=plan),3)
        push!(vertices,ns_vertex(localmodules))
    end
    Reference{T}(graph,modules,vertices,Dict(),Dict(),Dict())
end

"""Pure L-descendant basis at an integer Virasoro level."""
function vir_basis(reference::Reference,n::Int)
    get!(reference.bases,n) do
        [PBWState(word,()) for word in partitions(reference.modules[1].plan,n)]
    end
end

"""Restricted Virasoro Gram inverse, with no SCA fermionic descendants."""
function vir_inverse(reference::Reference{T},edge::Int,n::Int) where T
    get!(reference.inverses,(edge,n)) do
        states=vir_basis(reference,n)
        matrix=Complex{T}[pairing(reference.modules[edge],a,b) for a in states,b in states]
        lu(matrix)\Matrix{Complex{T}}(I,length(states),length(states))
    end
end

"""Reference descendant vertex obtained entirely through Virasoro Ward reductions."""
function vir_tensor(reference::Reference{T},v::Int,levels::NTuple{3,Int}) where T
    get!(reference.tensors,(v,levels...)) do
        bases=map(n->vir_basis(reference,n),levels)
        out=zeros(Complex{T},map(length,bases))
        for ix in CartesianIndices(out)
            out[ix]=rho(reference.vertices[v],0,ntuple(j->bases[j][ix[j]],3))
        end
        out
    end
end

"""One independent direct coefficient at the specified integer edge levels."""
function coefficient(reference::Reference{T},levels::NTuple{N,Int}) where {T,N}
    graph=reference.graph
    tensors=[vir_tensor(reference,v,ntuple(j->graph.slots[v][j]==0 ? 0 : levels[graph.slots[v][j]],3))
             for v in eachindex(graph.slots)]
    inverses=[vir_inverse(reference,e,levels[e]) for e in 1:N]
    selections=[collect(1:length(vir_basis(reference,levels[e]))) for e in 1:N]
    sew_tensors(graph,tensors,inverses,selections)
end
end
