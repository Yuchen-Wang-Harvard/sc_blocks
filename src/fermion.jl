"""Auxiliary Fock state: distinct descending doubled modes and Ramond ground bit."""
struct FermionState
    modes::Tuple{Vararg{Int16}}
    ground::Bool
end
FermionState(modes::Tuple=(),ground::Integer=0)=FermionState(map(Int16,modes),Bool(ground))
level2(s::FermionState)=sum(Int,s.modes;init=0)
relative_parity(s::FermionState)=xor(isodd(length(s.modes)),s.ground)

"""Auxiliary Majorana module; shares exact combinatorics but not physical caches."""
mutable struct FermionModule{T<:AbstractFloat}
    ramond::Bool
    plan::Combinatorics
    bases::Dict{Int,Vector{FermionState}}
    actions::Dict{Tuple{Int,FermionState},Dict{FermionState,Complex{T}}}
end
function FermionModule(::Type{T},sector::Symbol;plan=Combinatorics()) where T<:AbstractFloat
    sector in (:NS,:R) || throw(ArgumentError("sector must be :NS or :R"))
    FermionModule{T}(sector==:R,plan,Dict(),Dict())
end
state_parity(::FermionModule,s::FermionState)=relative_parity(s)

"""Auxiliary Fock basis at doubled descendant grade."""
function basis(m::FermionModule,n::Integer)
    get!(m.bases,Int(n)) do
        [FermionState(word,a) for word in fermion_partitions(m.plan,n,m.ramond)
                             for a in 0:Int(m.ramond)]
    end
end

"""Clifford mode action; mode is doubled, including the normalized Ramond zero mode."""
function act_psi(m::FermionModule{T},mode::Integer,s::FermionState) where T
    iseven(mode)==m.ramond || throw(ArgumentError("psi mode has wrong sector"))
    get!(m.actions,(Int(mode),s)) do
        C=Complex{T}
        out=Dict{FermionState,C}()
        if mode==0
            out[FermionState(s.modes,!s.ground)]=C(parity_sign(length(s.modes))/sqrt(T(2)))
        elseif mode>0
            j=findfirst(==(mode),s.modes)
            isnothing(j) || (out[FermionState((s.modes[1:j-1]...,s.modes[j+1:end]...),s.ground)]=C(parity_sign(j-1)))
        elseif -mode ∉ s.modes
            j=count(>(-mode),s.modes)
            out[FermionState((s.modes[1:j]...,Int16(-mode),s.modes[j+1:end]...),s.ground)]=C(parity_sign(j))
        end
        out
    end
end

"""Diagonal auxiliary BPZ norm, distinguished from its ordered three-point value."""
fermion_norm(::FermionModule{T},s::FermionState) where T =
    relative_parity(s) ? Complex{T}(0,-1) : one(Complex{T})
pairing(m::FermionModule{T},a::FermionState,b::FermionState) where T =
    a==b ? fermion_norm(m,a) : zero(Complex{T})
gram(m::FermionModule,n::Integer)=Diagonal([fermion_norm(m,s) for s in basis(m,n)])
inverse_gram(m::FermionModule,n::Integer)=Diagonal([inv(fermion_norm(m,s)) for s in basis(m,n)])

"""Even auxiliary three-point form in NSNSNS or NSRR order."""
mutable struct FermionVertex{T<:AbstractFloat}
    modules::NTuple{3,FermionModule{T}}
    values::Dict{NTuple{3,FermionState},Complex{T}}
end
function FermionVertex(modules::NTuple{3,FermionModule{T}}) where T
    map(m->m.ramond,modules) in ((false,false,false),(false,true,true)) ||
        throw(ArgumentError("auxiliary slots must be NSNSNS or NSRR"))
    FermionVertex{T}(modules,Dict())
end

"""Auxiliary mode inserted in one ordered slot."""
function psi_value(v::FermionVertex{T},slot::Integer,mode::Integer,states::NTuple{3,FermionState}) where T
    sum((coef*rho(v,Base.setindex(states,s,slot)) for (s,coef) in act_psi(v.modules[slot],mode,states[slot]));init=zero(Complex{T}))
end

"""Finite auxiliary binomial Ward sum in doubled-mode units."""
function psi_sum(v::FermionVertex{T},slot::Integer,start::Integer,states::NTuple{3,FermionState},a::Rational;
                 alternating::Bool=false,first::Integer=0) where T
    out=zero(Complex{T})
    for j in first:fld(level2(states[slot])-start,2)
        coefficient=exact_binomial(v.modules[1].plan,a,j)
        iszero(coefficient) && continue
        out+=coefficient*(alternating ? parity_sign(j) : 1)*psi_value(v,slot,start+2j,states)
    end
    out
end
psi_sum(v,slot,start,states,a::Integer;kwargs...)=psi_sum(v,slot,start,states,a//1;kwargs...)

"""Normalized auxiliary vertex from its own psi Ward identity, without physical forms."""
function rho(v::FermionVertex{T},states::NTuple{3,FermionState}) where T
    C=Complex{T}
    isodd(sum(Int(relative_parity(s)) for s in states)) && return zero(C)
    get!(v.values,states) do
        ii=C(0,1)
        if !isempty(states[1].modes)
            r2=Int(first(states[1].modes))
            base=Base.setindex(states,FermionState(Base.tail(states[1].modes)),1)
            e1,e2=Int(relative_parity(base[1])),Int(relative_parity(base[2]))
            if !v.modules[2].ramond
                k=div(r2-1,2)
                return -ii*parity_sign(e1)*psi_sum(v,2,1,base,k)-
                        ii*parity_sign(e1+e2)*psi_value(v,3,2k+1,base)
            end
            k=div(r2+1,2)
            return -psi_sum(v,1,-r2,base,-1//2;alternating=true,first=1)-
                    ii*parity_sign(e1)*psi_sum(v,2,0,base,(2k-1)//2)-
                    parity_sign(e1+e2)*psi_sum(v,3,2k,base,-1//2;alternating=true)
        elseif !isempty(states[2].modes)
            r2=Int(first(states[2].modes))
            base=Base.setindex(states,FermionState(Base.tail(states[2].modes),states[2].ground),2)
            e2=Int(relative_parity(base[2]))
            if !v.modules[2].ramond
                k=div(r2-1,2)
                return -parity_sign(e2-k-1)*psi_sum(v,3,1,base,-k-1;alternating=true)
            end
            k=div(r2,2)
            return -psi_sum(v,2,-2k,base,-1//2;first=1)+
                    ii*parity_sign(e2-k)*psi_sum(v,3,0,base,(-2k-1)//2;alternating=true)
        elseif !isempty(states[3].modes)
            v.modules[3].ramond || return zero(C)
            k=div(Int(first(states[3].modes)),2)
            base=Base.setindex(states,FermionState(Base.tail(states[3].modes),states[3].ground),3)
            return -psi_sum(v,3,-2k,base,-1//2;alternating=true,first=1)-
                    ii*parity_sign(Int(relative_parity(states[2])))*psi_value(v,2,0,base)
        end
        !v.modules[2].ramond && return one(C)
        states[2].ground==states[3].ground || return zero(C)
        states[2].ground ? ii : one(C)
    end
end

"""Independent master psi identity residual for arbitrary contour integers m,n."""
function master_ward_residual(v::FermionVertex{T},states::NTuple{3,FermionState},m::Integer,n::Integer) where T
    e1,e2=Int(relative_parity(states[1])),Int(relative_parity(states[2]))
    ii=Complex{T}(0,1)
    if v.modules[2].ramond
        lhs=psi_sum(v,1,1-2m-2n,states,(2n-1)//2;alternating=true)
        rhs=-ii*parity_sign(e1)*psi_sum(v,2,2n,states,(2m-1)//2)-
             parity_sign(e1+e2+n)*psi_sum(v,3,2m,states,(2n-1)//2;alternating=true)
    else
        lhs=psi_sum(v,1,-1-2m-2n,states,n;alternating=true)
        rhs=-ii*parity_sign(e1)*psi_sum(v,2,2n+1,states,m)-
             ii*parity_sign(e1+e2+n)*psi_sum(v,3,2m+1,states,n;alternating=true)
    end
    residual_pair((lhs,),(rhs,))
end

"""Full auxiliary-first graded tensor BPZ pairing."""
function tensor_pairing(aux::FermionModule,physical::SCAModule,a::FermionState,x::PBWState,b::FermionState,y::PBWState)
    parity_sign(Int(relative_parity(a))*Int(state_parity(physical,x)))*pairing(aux,a,b)*pairing(physical,x,y)
end

"""Ordered enlarged form including infinity dualization and middle/third reordering."""
function enlarged_rho(aux::FermionVertex,physical::WardVertex{T,K},
                      astates::NTuple{3,FermionState},states::NTuple{3,PBWState},f::Integer,eta::Integer=1) where {T,K}
    exponent=Int(relative_parity(astates[1]))*Int(state_parity(physical.modules[1],states[1]))+
             Int(relative_parity(astates[3]))*Int(state_parity(physical.modules[2],states[2]))
    value=K==1 ? rho(physical,f,states) : rho(physical,f,eta,states)
    parity_sign(exponent)*rho(aux,astates)*value
end
