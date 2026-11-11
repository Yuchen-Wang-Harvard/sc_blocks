"""Independent enlarged-module reference using only SCA PBW and auxiliary Clifford actions."""
module TensorPBWReference
using SuperconformalBlocks
const SCB=SuperconformalBlocks
const State=Tuple{FermionState,PBWState}

"""Physical and auxiliary modules and independent tensor-action caches."""
mutable struct Reference{T<:AbstractFloat}
    physical::SCAModule{T}
    auxiliary::FermionModule{T}
    b::Complex{T}
    bases::Dict{Int,Vector{State}}
    actions::Dict{Tuple{Int,Int,State},Dict{State,Complex{T}}}
    auxiliary_actions::Dict{Tuple{Int,FermionState},Dict{FermionState,Complex{T}}}
end

"""Construct the reference independently of free-field oscillator data."""
function Reference(::Type{T},sector::Symbol,b,p;parity=0) where T<:AbstractFloat
    C=Complex{T}
    bb,pp=C(b),C(p)
    q=bb+inv(bb)
    c=C(3)/2+3q*q
    physical=sector==:NS ? SCAModule(T,:NS,c;h=q*q/8-pp*pp/2,parity=parity) :
                          SCAModule(T,:R,c;beta=pp/sqrt(T(2)))
    Reference{T}(physical,FermionModule(T,sector),bb,Dict(),Dict(),Dict())
end

"""Auxiliary-first tensor PBW basis at doubled grade n."""
function tensor_basis(reference::Reference,n::Int)
    get!(reference.bases,n) do
        [(a,x) for j in 0:n for a in basis(reference.auxiliary,j) for x in basis(reference.physical,n-j)]
    end
end

"""Auxiliary Virasoro action evaluated by a bounded full mode sum."""
function auxiliary_virasoro(reference::Reference{T},n::Int,s::FermionState) where T
    get!(reference.auxiliary_actions,(n,s)) do
        C=Complex{T}
        out=Dict{FermionState,C}()
        2n>level2(s) && return out
        if n==0
            out[s]=C(level2(s))/2+(reference.auxiliary.ramond ? C(1)/16 : zero(C))
            return out
        end
        for r2 in 2n-level2(s):level2(s)
            iseven(r2)==reference.auxiliary.ramond || continue
            for (x,a) in act_psi(reference.auxiliary,r2,s),(y,b) in act_psi(reference.auxiliary,2n-r2,x)
                out[y]=get(out,y,zero(C))+r2*a*b/4
            end
        end
        SCB.prune_zeros!(out)
    end
end

"""Embedded Virasoro mode evaluated directly in the physical PBW tensor basis."""
function embedded_action(reference::Reference{T},copy::Int,n::Int,state::State) where T
    get!(reference.actions,(copy,n,state)) do
        C=Complex{T}
        auxiliary,physical=state
        out=Dict{State,C}()
        2n>level2(auxiliary)+level2(physical) && return out
        b=reference.b
        scale,fs,denominator=copy==1 ? (inv(b),-(inv(b)+2b),inv(b)-b) :
                                      (b,-(b+2inv(b)),b-inv(b))
        for (y,value) in act(reference.physical,L(n),physical)
            out[(auxiliary,y)]=get(out,(auxiliary,y),zero(C))+scale*value/denominator
        end
        for (a,value) in auxiliary_virasoro(reference,n,auxiliary)
            out[(a,physical)]=get(out,(a,physical),zero(C))+fs*value/denominator
        end
        for r2 in 2n-level2(auxiliary):level2(physical)
            iseven(r2)==reference.physical.ramond || continue
            for (y,cg) in act(reference.physical,G(r2),physical),(a,cp) in act_psi(reference.auxiliary,2n-r2,auxiliary)
                out[(a,y)]=get(out,(a,y),zero(C))+
                    SCB.parity_sign(Int(relative_parity(auxiliary)))*cg*cp/denominator
            end
        end
        SCB.prune_zeros!(out)
    end
end

"""Reference embedded mode acting on an ordered tensor vector."""
function embedded_apply(reference::Reference{T},copy::Int,n::Int,vector::Dict{State,Complex{T}}) where T
    out=Dict{State,Complex{T}}()
    for (s,value) in vector
        SCB.sparse_add!(out,embedded_action(reference,copy,n,s),value)
    end
    SCB.prune_zeros!(out)
end
end
