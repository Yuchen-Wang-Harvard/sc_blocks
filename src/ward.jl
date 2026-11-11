"""Ordered physical three-point form. K=1 for NSNSNS, K=2 for the two NSRR eta signs."""
mutable struct WardVertex{T<:AbstractFloat,K}
    modules::NTuple{3,SCAModule{T}}
    values::Dict{NTuple{3,PBWState},NTuple{K,Complex{T}}}
    tensors::Dict{NTuple{5,Int},Array{Complex{T},3}}
end

"""All-NS form with relative-descendant labels a and explicit primary parities."""
function ns_vertex(modules::NTuple{3,SCAModule{T}}) where T
    any(m->m.ramond,modules) && throw(ArgumentError("all three modules must be NS"))
    all(m->m.c==modules[1].c,modules) || throw(ArgumentError("central charges differ"))
    WardVertex{T,1}(modules,Dict(),Dict())
end

"""NSRR form, computing eta=+1 and eta=-1 together from the current Ward identity."""
function nsrr_vertex(modules::NTuple{3,SCAModule{T}}) where T
    map(m->m.ramond,modules)==(false,true,true) ||
        throw(ArgumentError("slots must be ordered NS,R,R"))
    all(m->m.c==modules[1].c,modules) || throw(ArgumentError("central charges differ"))
    WardVertex{T,2}(modules,Dict(),Dict())
end

@inline form_zero(::WardVertex{T,K}) where {T,K} = ntuple(_->zero(Complex{T}),Val(K))
@inline form_add(a::NTuple{K,C}, b::NTuple{K,C}) where {K,C} = map(+,a,b)
@inline form_scale(a::NTuple{K,C}, b) where {K,C} = map(x->C(b*x),a)
@inline form_sub(a,b) = form_add(a,form_scale(b,-1))

"""Insert one mode into a form, using translation before descendant expansion."""
function mode_value(v::WardVertex{T,K}, slot::Integer, mode::Mode,
                    states::NTuple{3,PBWState}) where {T,K}
    if mode==L(-1) && slot==2
        weights=ntuple(j->v.modules[j].h+T(level2(states[j]))/2,3)
        return form_scale(ward_values(v,states),weights[1]-weights[2]-weights[3])
    elseif mode==L(-1) && slot==3
        return form_sub(mode_value(v,1,L(1),states),mode_value(v,2,L(-1),states))
    end
    out=form_zero(v)
    for (state,coefficient) in act(v.modules[slot],mode,states[slot])
        out=form_add(out,form_scale(ward_values(v,Base.setindex(states,state,slot)),coefficient))
    end
    out
end

"""Finite binomial mode sum; positive modes above the state grade annihilate it."""
function ward_sum(v::WardVertex, slot::Integer, odd::Bool, start::Integer,
                  states::NTuple{3,PBWState}, a::Rational;
                  alternating::Bool=false, first::Integer=0)
    finish=odd ? fld(level2(states[slot])-start,2) : fld(level2(states[slot]),2)-start
    out=form_zero(v)
    for j in first:finish
        coefficient=exact_binomial(v.modules[1].plan,a,j)
        iszero(coefficient) && continue
        alternating && (coefficient*=parity_sign(j))
        mode=odd ? G(start+2j) : L(start+j)
        out=form_add(out,form_scale(mode_value(v,slot,mode,states),coefficient))
    end
    out
end
ward_sum(v,slot,odd,start,states,a::Integer;kwargs...) =
    ward_sum(v,slot,odd,start,states,a//1;kwargs...)

"""Shared Virasoro reduction of the first remaining L lowering mode."""
function reduce_l(v::WardVertex,states::NTuple{3,PBWState},slot::Int)
    old=states[slot]
    k=Int(first(old.l))
    base=Base.setindex(states,PBWState(Base.tail(old.l),old.g,old.ground),slot)
    if slot==1
        return form_add(ward_sum(v,2,false,-1,base,k+1),mode_value(v,3,L(k),base))
    elseif slot==3
        return form_sub(mode_value(v,1,L(k),base),ward_sum(v,2,false,-1,base,1-k))
    elseif k==1
        return mode_value(v,2,L(-1),base)
    end
    form_sub(ward_sum(v,1,false,k,base,1-k;alternating=true),
             form_scale(ward_sum(v,3,false,-1,base,1-k;alternating=true),parity_sign(1-k)))
end

"""No-f NSRR fermionic reduction; signs use TOTAL state parity."""
function reduce_g(v::WardVertex{T,2},states::NTuple{3,PBWState},slot::Int) where T
    old=states[slot]
    r2=Int(first(old.g))
    base=Base.setindex(states,PBWState(old.l,Base.tail(old.g),old.ground),slot)
    e1,e2=Int(state_parity(v.modules[1],base[1])),Int(state_parity(v.modules[2],base[2]))
    ii=Complex{T}(0,1)
    if slot==1
        k=div(r2-1,2)
        return form_sub(form_sub(
            form_scale(ward_sum(v,2,true,0,base,(2k+1)//2),ii*parity_sign(e1)),
            form_scale(ward_sum(v,3,true,2k,base,1//2;alternating=true),parity_sign(e1+e2))),
            ward_sum(v,1,true,-r2,base,1//2;alternating=true,first=1))
    end
    k=div(r2,2)
    if slot==2
        return form_sub(form_scale(form_add(
            ward_sum(v,1,true,2k-1,base,(1-2k)//2;alternating=true),
            form_scale(ward_sum(v,3,true,0,base,(1-2k)//2;alternating=true),parity_sign(e1+e2-k))),
            inv(ii*parity_sign(e1))),ward_sum(v,2,true,-2k,base,1//2;first=1))
    end
    form_sub(form_scale(form_sub(
        form_scale(ward_sum(v,2,true,0,base,(1-2k)//2),ii*parity_sign(e1)),
        ward_sum(v,1,true,2k-1,base,1//2;alternating=true)),parity_sign(e1+e2)),
        ward_sum(v,3,true,-2k,base,1//2;alternating=true,first=1))
end

"""All-NS fermionic reduction, retaining the normalized middle G_-1/2 anchor."""
function reduce_g(v::WardVertex{T,1},states::NTuple{3,PBWState},slot::Int) where T
    old=states[slot]
    k=div(Int(first(old.g))-1,2)
    slot==2 && k==0 && return (one(Complex{T}),)
    base=Base.setindex(states,PBWState(old.l,Base.tail(old.g),old.ground),slot)
    e1,e2=Int(state_parity(v.modules[1],base[1])),Int(state_parity(v.modules[2],base[2]))
    ii=Complex{T}(0,1)
    if slot==1
        return form_add(form_scale(ward_sum(v,2,true,-1,base,k+1),ii*parity_sign(e1)),
                        form_scale(mode_value(v,3,G(2k+1),base),ii*parity_sign(e1+e2)))
    elseif slot==3
        return form_scale(form_sub(mode_value(v,1,G(2k+1),base),
            form_scale(ward_sum(v,2,true,-1,base,-k),ii*parity_sign(e1))),inv(ii*parity_sign(e1+e2)))
    end
    form_scale(form_sub(ward_sum(v,1,true,2k+1,base,-k;alternating=true),
        form_scale(ward_sum(v,3,true,-1,base,-k;alternating=true),ii*parity_sign(e1+e2-k))),
        inv(ii*parity_sign(e1)))
end

"""Memoized normalized form values in the relative-parity channel of the states."""
function ward_values(v::WardVertex{T,K},states::NTuple{3,PBWState}) where {T,K}
    get!(v.values,states) do
        for slot in (1,3,2)
            isempty(states[slot].l) || return reduce_l(v,states,slot)
        end
        for slot in (K==1 ? (1,3,2) : (1,2,3))
            isempty(states[slot].g) || return reduce_g(v,states,slot)
        end
        C=Complex{T}
        K==1 && return (one(C),)
        a,b=states[2].ground,states[3].ground
        !a && return (one(C),one(C))
        b ? (one(C),-one(C)) : (C(0,1),C(0,-1))
    end
end
ward_values(v::WardVertex,s1::PBWState,s2::PBWState,s3::PBWState) = ward_values(v,(s1,s2,s3))

"""Select the normalized NSRR rho_f^(eta), with f a RELATIVE parity label."""
function rho(v::WardVertex{T,2},f::Integer,eta::Integer,states::NTuple{3,PBWState}) where T
    f in (0,1) && eta in (-1,1) || throw(ArgumentError("invalid form label"))
    mod(sum(Int(relative_parity(s)) for s in states),2)==f || return zero(Complex{T})
    ward_values(v,states)[eta==1 ? 1 : 2]
end

"""Select normalized all-NS rho_a, with a the parity of lowering operators."""
function rho(v::WardVertex{T,1},a::Integer,states::NTuple{3,PBWState}) where T
    a in (0,1) || throw(ArgumentError("invalid form label"))
    mod(sum(Int(relative_parity(s)) for s in states),2)==a || return zero(Complex{T})
    only(ward_values(v,states))
end

"""Cache a complete descendant vertex once for reuse in graph contractions."""
function vertex_tensor(v::WardVertex{T,K},levels::NTuple{3,Int},f::Integer,eta::Integer=1) where {T,K}
    get!(v.tensors,(levels...,Int(f),Int(eta))) do
        bases=ntuple(j->basis(v.modules[j],levels[j]),3)
        out=zeros(Complex{T},map(length,bases))
        for ix in CartesianIndices(out)
            states=ntuple(j->bases[j][ix[j]],3)
            out[ix]=K==1 ? rho(v,f,states) : rho(v,f,eta,states)
        end
        out
    end
end

"""Return absolute/scaled residual of an arbitrary master G Ward identity."""
function master_ward_residual(v::WardVertex{T,K},states::NTuple{3,PBWState},m::Integer,n::Integer) where {T,K}
    e1,e2=Int(state_parity(v.modules[1],states[1])),Int(state_parity(v.modules[2],states[2]))
    ii=Complex{T}(0,1)
    if K==2
        lhs=ward_sum(v,1,true,-2m-2n-1,states,(2n+1)//2;alternating=true)
        rhs=form_sub(form_scale(ward_sum(v,2,true,2n,states,(2m+1)//2),ii*parity_sign(e1)),
            form_scale(ward_sum(v,3,true,2m,states,(2n+1)//2;alternating=true),parity_sign(e1+e2+n)))
    else
        lhs=ward_sum(v,1,true,-2m-2n+1,states,n;alternating=true)
        rhs=form_add(form_scale(ward_sum(v,2,true,2n-1,states,m),ii*parity_sign(e1)),
            form_scale(ward_sum(v,3,true,2m-1,states,n;alternating=true),ii*parity_sign(e1+e2+n)))
    end
    residual_pair(lhs,rhs)
end

"""Return absolute/scaled residual of an arbitrary master Virasoro Ward identity."""
function virasoro_ward_residual(v::WardVertex,states::NTuple{3,PBWState},m::Integer,n::Integer)
    lhs=ward_sum(v,1,false,1-m-n,states,n;alternating=true)
    rhs=form_add(ward_sum(v,2,false,n-1,states,m),
        form_scale(ward_sum(v,3,false,m-1,states,n;alternating=true),parity_sign(n)))
    residual_pair(lhs,rhs)
end

"""Residual scaled by max(1, component magnitudes), without precision conversion."""
function residual_pair(lhs,rhs)
    absolute=maximum(abs,form_sub(lhs,rhs))
    absolute, absolute/max(one(absolute),maximum(abs,lhs),maximum(abs,rhs))
end
