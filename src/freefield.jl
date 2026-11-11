"""Interned free-field state, ordered auxiliary first, then bosons and physical fermions."""
struct OscillatorState
    bosons::Tuple{Vararg{Int16}}
    fermions::FermionState
    auxiliary::FermionState
end
OscillatorState(bosons::Tuple=(),fermions::FermionState=FermionState(),auxiliary::FermionState=FermionState())=
    OscillatorState(map(Int16,bosons),fermions,auxiliary)
physical_level2(s::OscillatorState)=2sum(Int,s.bosons;init=0)+level2(s.fermions)
level2(s::OscillatorState)=physical_level2(s)+level2(s.auxiliary)

"""One physical PBW-to-oscillator change of basis, with reusable factorization."""
struct FockMapping{T<:AbstractFloat}
    states::Vector{UInt32}
    index::Dict{UInt32,Int}
    matrix::Matrix{Complex{T}}
    factor::LU{Complex{T},Matrix{Complex{T}},Vector{Int}}
end

"""Free-field realization with current Ramond zero-mode phases and interned sparse actions.

Physical action caches omit auxiliary spectators, and auxiliary Virasoro caches
omit physical spectators. Embedded actions and common descendant suffixes are
shared. PBW conversion is local basis conversion, never a final block conversion.
"""
mutable struct FreeField{T<:AbstractFloat}
    ramond::Bool
    b::Complex{T}
    p::Complex{T}
    q::Complex{T}
    physical::SCAModule{T}
    auxiliary::FermionModule{T}
    states::Vector{OscillatorState}
    ids::Dict{OscillatorState,UInt32}
    physical_actions::Dict{Tuple{Mode,UInt32},Dict{UInt32,Complex{T}}}
    auxiliary_actions::Dict{Tuple{Int,FermionState},Dict{FermionState,Complex{T}}}
    embedded_actions::Dict{NTuple{3,Int},Dict{UInt32,Complex{T}}}
    pbw_images::Dict{PBWState,Dict{UInt32,Complex{T}}}
    mappings::Dict{Int,FockMapping{T}}
    primaries::Dict{Tuple{Int,Int},Dict{UInt32,Complex{T}}}
    descendants::Dict{Tuple{Int,Int,Tuple{Vararg{Int16}},Tuple{Vararg{Int16}}},Dict{UInt32,Complex{T}}}
    reflected::Union{Nothing,FreeField{T}}
end

"""Construct one chart at momentum P, with beta=P/sqrt(2) in the R sector."""
function FreeField(::Type{T},sector::Symbol,b,p;parity::Integer=0,plan=Combinatorics()) where T<:AbstractFloat
    C=Complex{T}
    bb,pp=C(b),C(p)
    bb in (zero(C),one(C),-one(C)) && throw(ArgumentError("b must differ from 0,±1"))
    q=bb+inv(bb)
    c=C(3)/2+3q*q
    physical=sector==:NS ? SCAModule(T,sector,c;h=q*q/8-pp*pp/2,parity=parity,plan=plan) :
                          SCAModule(T,sector,c;beta=pp/sqrt(T(2)),parity=parity,plan=plan)
    module_=FreeField{T}(physical.ramond,bb,pp,q,physical,FermionModule(T,sector;plan=plan),
        OscillatorState[],Dict(),Dict(),Dict(),Dict(),Dict(),Dict(),Dict(),Dict(),nothing)
    intern!(module_,OscillatorState())
    module_
end

"""Assign a compact stable ID to a canonical oscillator state."""
function intern!(module_::FreeField,s::OscillatorState)
    get!(module_.ids,s) do
        length(module_.states)<typemax(UInt32) || throw(OverflowError("oscillator IDs exhausted"))
        push!(module_.states,s)
        UInt32(length(module_.states))
    end
end

"""Accumulate a single oscillator state into an interned sparse vector."""
function ff_add!(module_::FreeField{T},out::Dict{UInt32,Complex{T}},s::OscillatorState,value) where T
    iszero(value) && return out
    id=intern!(module_,s)
    out[id]=get(out,id,zero(Complex{T}))+value
    out
end

"""Bosonic oscillator action; repeated creators are stored as a descending tuple."""
function boson_action(n::Int,s::OscillatorState)
    n==0 && throw(ArgumentError("c_0 is not an oscillator"))
    if n>0
        multiplicity=count(==(n),s.bosons)
        multiplicity==0 && return nothing,0
        ix=something(findfirst(==(n),s.bosons))
        bosons=(s.bosons[1:ix-1]...,s.bosons[ix+1:end]...)
        return OscillatorState(bosons,s.fermions,s.auxiliary),n*multiplicity
    end
    ix=count(>(-n),s.bosons)
    bosons=(s.bosons[1:ix]...,Int16(-n),s.bosons[ix+1:end]...)
    OscillatorState(bosons,s.fermions,s.auxiliary),1
end

"""Physical free fermion eta action, including eta_0=-exp(∓iπ/4)/sqrt(2)."""
function eta_action(module_::FreeField{T},r2::Int,s::FermionState) where T
    if r2!=0
        return act_psi(module_.auxiliary,r2,s)
    end
    module_.ramond || throw(ArgumentError("eta zero mode in NS sector"))
    # exp(∓iπ/4)/sqrt(2)=(1∓i)/2, avoiding numerical pi.
    coefficient=-Complex{T}(1,s.ground ? 1 : -1)*parity_sign(length(s.modes))/2
    Dict(FermionState(s.modes,!s.ground)=>coefficient)
end

"""Physical L/G action with auxiliary spectators removed from its cache key."""
function ff_physical_action(module_::FreeField{T},mode::Mode,id::UInt32) where T
    state=module_.states[id]
    clean=OscillatorState(state.bosons,state.fermions)
    cleanid=intern!(module_,clean)
    get!(module_.physical_actions,(mode,cleanid)) do
        C=Complex{T}
        out=Dict{UInt32,C}()
        shift=mode.odd ? Int(mode.index) : 2Int(mode.index)
        shift>physical_level2(clean) && return out
        if !mode.odd && mode.index==0
            out[cleanid]=module_.physical.h+T(physical_level2(clean))/2
            return out
        end
        ii=C(0,1)
        if !mode.odd
            n=Int(mode.index)
            candidates=Set{Int}()
            for k in clean.bosons
                push!(candidates,Int(k),n-Int(k))
            end
            n<0 && union!(candidates,n+1:-1)
            for k in candidates
                (k==0 || k==n) && continue
                first,a=boson_action(k,clean)
                iszero(a) && continue
                second,b=boson_action(n-k,first)
                iszero(b) || ff_add!(module_,out,second,C(a*b)/2)
            end
            empty!(candidates)
            for r2 in clean.fermions.modes
                push!(candidates,Int(r2),2n-Int(r2))
            end
            if n<0
                union!(candidates,module_.ramond ? (2n:2:0) : (2n+1:2:-1))
            end
            for r2 in candidates,(first,a) in eta_action(module_,r2,clean.fermions),
                (second,b) in eta_action(module_,2n-r2,first)
                ff_add!(module_,out,OscillatorState(clean.bosons,second),C(r2)*a*b/4)
            end
            changed,coefficient=boson_action(n,clean)
            iszero(coefficient) || ff_add!(module_,out,changed,ii*(module_.q*n-2module_.p)*coefficient/2)
        else
            r2=Int(mode.index)
            iseven(r2)==module_.ramond || throw(ArgumentError("G mode has wrong sector"))
            candidates=Set{Int}(Int.(clean.bosons))
            for s2 in clean.fermions.modes
                push!(candidates,div(r2-Int(s2),2))
            end
            module_.ramond && push!(candidates,div(r2,2))
            r2<0 && union!(candidates,fld(r2,2):-1)
            for k in candidates
                k==0 && continue
                for (fermion,a) in eta_action(module_,r2-2k,clean.fermions)
                    changed,b=boson_action(k,OscillatorState(clean.bosons,fermion))
                    iszero(b) || ff_add!(module_,out,changed,a*b)
                end
            end
            for (fermion,coefficient) in eta_action(module_,r2,clean.fermions)
                ff_add!(module_,out,OscillatorState(clean.bosons,fermion),ii*(module_.q*r2/2-module_.p)*coefficient)
            end
        end
        prune_zeros!(out)
    end
end

"""Apply a physical-factor mode to an oscillator vector, reattaching auxiliary spectators."""
function ff_physical_apply(module_::FreeField{T},mode::Mode,vector::Dict{UInt32,Complex{T}}) where T
    out=Dict{UInt32,Complex{T}}()
    for (id,value) in vector
        auxiliary=module_.states[id].auxiliary
        for (dest,coefficient) in ff_physical_action(module_,mode,id)
            state=module_.states[dest]
            ff_add!(module_,out,OscillatorState(state.bosons,state.fermions,auxiliary),value*coefficient)
        end
    end
    prune_zeros!(out)
end

"""Auxiliary Virasoro action, cached independently of all physical spectators."""
function ff_auxiliary_action(module_::FreeField{T},n::Int,state::FermionState) where T
    get!(module_.auxiliary_actions,(n,state)) do
        C=Complex{T}
        out=Dict{FermionState,C}()
        2n>level2(state) && return out
        if n==0
            out[state]=T(level2(state))/2+(module_.ramond ? C(1)/16 : zero(C))
            return out
        end
        candidates=Set{Int}()
        for r2 in state.modes
            push!(candidates,Int(r2),2n-Int(r2))
        end
        n<0 && union!(candidates,module_.ramond ? (2n:2:0) : (2n+1:2:-1))
        for r2 in candidates,(first,a) in act_psi(module_.auxiliary,r2,state),
            (second,b) in act_psi(module_.auxiliary,2n-r2,first)
            out[second]=get(out,second,zero(C))+r2*a*b/4
        end
        prune_zeros!(out)
    end
end

"""One embedded Virasoro mode in the current auxiliary-first tensor convention."""
function ff_embedded_action(module_::FreeField{T},copy::Int,n::Int,id::UInt32) where T
    copy in (1,2) || throw(ArgumentError("Virasoro copy must be 1 or 2"))
    get!(module_.embedded_actions,(copy,n,Int(id))) do
        C=Complex{T}
        out=Dict{UInt32,C}()
        state=module_.states[id]
        2n>level2(state) && return out
        t=copy==1 ? inv(module_.b) : module_.b
        opposite=inv(t)
        denominator=t-opposite
        for (dest,value) in ff_physical_action(module_,L(n),id)
            target=module_.states[dest]
            ff_add!(module_,out,OscillatorState(target.bosons,target.fermions,state.auxiliary),t*value/denominator)
        end
        for (auxiliary,value) in ff_auxiliary_action(module_,n,state.auxiliary)
            ff_add!(module_,out,OscillatorState(state.bosons,state.fermions,auxiliary),-(t+2opposite)*value/denominator)
        end
        for r2 in 2n-level2(state.auxiliary):physical_level2(state)
            iseven(r2)==module_.ramond || continue
            for (auxiliary,psi) in act_psi(module_.auxiliary,2n-r2,state.auxiliary),
                (dest,value) in ff_physical_action(module_,G(r2),id)
                target=module_.states[dest]
                ff_add!(module_,out,OscillatorState(target.bosons,target.fermions,auxiliary),
                    parity_sign(Int(relative_parity(state.auxiliary)))*psi*value/denominator)
            end
        end
        prune_zeros!(out)
    end
end

"""Embedded mode on a sparse vector, reusing every basis-state action."""
function ff_embedded_apply(module_::FreeField{T},copy::Int,n::Int,vector::Dict{UInt32,Complex{T}}) where T
    out=Dict{UInt32,Complex{T}}()
    for (id,value) in vector
        sparse_add!(out,ff_embedded_action(module_,copy,n,id),value)
    end
    prune_zeros!(out)
end

"""Chi=psi-i eta action; reflected_zero supplies chi_0(-P) in the same ground basis."""
function chi_apply(module_::FreeField{T},r2::Int,vector::Dict{UInt32,Complex{T}};reflected_zero::Bool=false) where T
    reflected_zero && r2!=0 && throw(ArgumentError("only the reflected zero mode shares this chart"))
    out=Dict{UInt32,Complex{T}}()
    ii=Complex{T}(0,1)
    for (id,value) in vector
        state=module_.states[id]
        for (auxiliary,coefficient) in act_psi(module_.auxiliary,r2,state.auxiliary)
            ff_add!(module_,out,OscillatorState(state.bosons,state.fermions,auxiliary),value*coefficient)
        end
        for (fermions,coefficient) in eta_action(module_,r2,state.fermions)
            ff_add!(module_,out,OscillatorState(state.bosons,fermions,state.auxiliary),
                    value*(reflected_zero ? ii : -ii)*parity_sign(Int(relative_parity(state.auxiliary)))*coefficient)
        end
    end
    prune_zeros!(out)
end

"""Physical PBW word realized in the oscillator chart, caching common word tails."""
function ff_pbw_image(module_::FreeField{T},state::PBWState) where T
    get!(module_.pbw_images,state) do
        if !isempty(state.l)
            tail=PBWState(Base.tail(state.l),state.g,state.ground)
            return ff_physical_apply(module_,L(-first(state.l)),ff_pbw_image(module_,tail))
        elseif !isempty(state.g)
            tail=PBWState((),Base.tail(state.g),state.ground)
            return ff_physical_apply(module_,G(-first(state.g)),ff_pbw_image(module_,tail))
        end
        id=intern!(module_,OscillatorState((),FermionState((),state.ground)))
        Dict(id=>one(Complex{T}))
    end
end

"""Factor one physical PBW-to-Fock map for use by every auxiliary spectator group."""
function fock_mapping(module_::FreeField{T},n::Int) where T
    get!(module_.mappings,n) do
        basis_=basis(module_.physical,n)
        states=[intern!(module_,OscillatorState(s.l,FermionState(s.g,s.ground))) for s in basis_]
        index=Dict(id=>j for (j,id) in enumerate(states))
        matrix=zeros(Complex{T},length(states),length(states))
        for (j,state) in enumerate(basis_),(id,value) in ff_pbw_image(module_,state)
            matrix[index[id],j]=value
        end
        FockMapping{T}(states,index,matrix,lu(matrix))
    end
end

"""Convert a local oscillator vector to the ordered auxiliary-times-PBW basis."""
function ordered_pbw(module_::FreeField{T},vector::Dict{UInt32,Complex{T}}) where T
    C=Complex{T}
    groups=Dict{Tuple{FermionState,Int},Dict{UInt32,C}}()
    for (id,value) in vector
        state=module_.states[id]
        key=(state.auxiliary,physical_level2(state))
        physical_id=intern!(module_,OscillatorState(state.bosons,state.fermions))
        group=get!(() -> Dict{UInt32,C}(),groups,key)
        group[physical_id]=get(group,physical_id,zero(C))+value
    end
    out=Dict{Tuple{FermionState,PBWState},C}()
    for ((auxiliary,n),group) in groups
        mapping=fock_mapping(module_,n)
        rhs=zeros(C,length(mapping.states))
        for (id,value) in group
            rhs[mapping.index[id]]=value
        end
        coefficients=mapping.factor\rhs
        for (j,state) in enumerate(basis(module_.physical,n))
            iszero(coefficients[j]) || (out[(auxiliary,state)]=coefficients[j])
        end
    end
    out
end

"""Evaluate an ordered auxiliary-times-PBW vector in the chosen oscillator chart."""
function from_ordered_pbw(module_::FreeField{T},vector::Dict{Tuple{FermionState,PBWState},Complex{T}}) where T
    out=Dict{UInt32,Complex{T}}()
    for ((auxiliary,physical),value) in vector,(id,coefficient) in ff_pbw_image(module_,physical)
        state=module_.states[id]
        ff_add!(module_,out,OscillatorState(state.bosons,state.fermions,auxiliary),value*coefficient)
    end
    prune_zeros!(out)
end

"""The manuscript chi-vector primary, n4=4n; negative labels use its explicit reflection."""
function ff_primary(module_::FreeField{T},n4::Int,alpha::Int=0) where T
    (module_.ramond ? isodd(n4) && alpha in (0,1) : iseven(n4) && alpha==0) ||
        throw(ArgumentError("invalid free-field branch label"))
    get!(module_.primaries,(n4,alpha)) do
        C=Complex{T}
        if n4<0
            if isnothing(module_.reflected)
                module_.reflected=FreeField(T,module_.ramond ? :R : :NS,module_.b,-module_.p;
                    parity=Int(module_.physical.primary_parity),plan=module_.physical.plan)
            end
            reflected=module_.reflected
            ordered=ordered_pbw(reflected,ff_primary(reflected,-n4,alpha))
            if module_.ramond
                for (key,value) in ordered
                    ordered[key]=parity_sign(Int(last(key).ground))*value
                end
            end
            return from_ordered_pbw(module_,ordered)
        end
        vector=Dict(UInt32(1)=>one(C))
        if !module_.ramond
            for r2 in 1:2:n4-1
                vector=chi_apply(module_,-r2,vector)
            end
        else
            m=div(n4-1,2)
            mod(m+1,2)!=alpha && (vector=chi_apply(module_,0,vector;reflected_zero=true))
            for r in m:-1:1
                vector=chi_apply(module_,-2r,vector)
            end
            vector=chi_apply(module_,0,vector)
        end
        vector
    end
end

"""Descendant of a branching primary with shared suffixes in both Virasoro copies."""
function ff_descendant(module_::FreeField,n4::Int,alpha::Int,first::Tuple,second::Tuple)
    a,b=map(Int16,first),map(Int16,second)
    get!(module_.descendants,(n4,alpha,a,b)) do
        isempty(a) && isempty(b) && return ff_primary(module_,n4,alpha)
        if !isempty(a) && (isempty(b) || a[1]>=b[1])
            ff_embedded_apply(module_,1,-Int(a[1]),ff_descendant(module_,n4,alpha,Base.tail(a),b))
        else
            ff_embedded_apply(module_,2,-Int(b[1]),ff_descendant(module_,n4,alpha,a,Base.tail(b)))
        end
    end
end

"""Current full tensor BPZ pairing of oscillator vectors via independent physical Grams."""
function ff_pairing(module_::FreeField{T},left::Dict{UInt32,Complex{T}},right::Dict{UInt32,Complex{T}}) where T
    ordered_left,ordered_right=ordered_pbw(module_,left),ordered_pbw(module_,right)
    value=zero(Complex{T})
    for ((a,x),vx) in ordered_left,((b,y),vy) in ordered_right
        a==b || continue
        value+=vx*vy*tensor_pairing(module_.auxiliary,module_.physical,a,x,b,y)
    end
    value
end
"""Current full tensor BPZ norm of a local oscillator vector."""
ff_norm(module_::FreeField,vector)=ff_pairing(module_,vector,vector)
