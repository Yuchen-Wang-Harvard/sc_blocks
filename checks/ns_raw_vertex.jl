"""Independent raw all-NS chi-vector three-point reference, using physical G Ward forms."""
module NSRawReference
using SuperconformalBlocks

mutable struct Reference{T<:AbstractFloat}
    modules::NTuple{3,FreeField{T}}
    physical::WardVertex{T,1}
    auxiliary::FermionVertex{T}
    ordered::Dict{Tuple{Int,Int},Dict{Tuple{FermionState,PBWState},Complex{T}}}
end

function Reference(::Type{T},b,momenta::NTuple{3};parities=(0,0,0)) where T<:AbstractFloat
    modules=ntuple(j->FreeField(T,:NS,b,momenta[j];parity=parities[j]),3)
    Reference{T}(modules,ns_vertex(map(m->m.physical,modules)),FermionVertex(map(m->m.auxiliary,modules)),Dict())
end

"""Raw local value without the closed Appendix-B formula or a norm-root convention."""
function coefficient(reference::Reference{T},labels::NTuple{3,Int},a::Int) where T
    vectors=ntuple(j->get!(reference.ordered,(j,labels[j])) do
        ordered_pbw(reference.modules[j],ff_primary(reference.modules[j],2labels[j]))
    end,3)
    value=zero(Complex{T})
    for (s1,c1) in vectors[1],(s2,c2) in vectors[2],(s3,c3) in vectors[3]
        states=(s1,s2,s3)
        value+=c1*c2*c3*enlarged_rho(reference.auxiliary,reference.physical,map(first,states),map(last,states),a)
    end
    value
end
end
