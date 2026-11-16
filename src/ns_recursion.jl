"""NS superconformal central-charge pole and its native residue prefactor."""
struct NSPole{T<:AbstractFloat}
    central::Complex{T}
    residue::Complex{T}
    x::Complex{T}
    id::Int
    fusion_id::Int
end

"""Independent all-NS c-recursion, using superglobal seeds and native BPZ null transport."""
mutable struct NSRecursionPlan{T<:AbstractFloat,N}
    graph::PlumbingGraph
    weights::NTuple{N,Complex{T}}
    parities::NTuple{N,Int}
    primary_mask::Int
    phases::Vector{Int}
    centers::Vector{Complex{T}}
    fusion_ids::Dict{Complex{T},Int}
    poles::Dict{NTuple{4,Int},NSPole{T}}
    fusions::Dict{NTuple{8,Int},Complex{T}}
    denominators::Dict{NTuple{2,Int},Complex{T}}
    vertices::Dict{NTuple{7,Int},Complex{T}}
    norms::Dict{NTuple{3,Int},Complex{T}}
    inverse_norms::Dict{NTuple{3,Int},Complex{T}}
    compact_theta::Bool
    transitions::Int
    seed_terms::Int
end

"""Exponent of the local BPZ phase of an all-NS superglobal state, without its Koszul sign."""
function ns_mode_phase(graph::PlumbingGraph,parities,mask::Int)
    exponent=0
    for e in eachindex(graph.ends)
        bit=(mask>>(e-1))&1
        exponent-=bit*(1+2parities[e])
    end
    for slots in graph.slots
        firstbit=(mask>>(slots[1]-1))&1
        thirdbit=(mask>>(slots[3]-1))&1
        exponent+=firstbit*(1+2parities[slots[1]])+2parities[slots[2]]*thirdbit
    end
    mod(exponent,4)
end

function NSRecursionPlan(::Type{T},graph::PlumbingGraph,central,weights::NTuple{N};
                         primary_parities::NTuple{N,Int}=ntuple(_->0,N),compact_theta::Bool=true) where {T<:AbstractFloat,N}
    N==length(graph.ends) || throw(DimensionMismatch("edge weights"))
    all(>(0),Iterators.flatten(graph.slots)) || throw(ArgumentError("NS reference currently requires closed graphs"))
    all(p->p in (0,1),primary_parities) || throw(ArgumentError("invalid intrinsic parity"))
    all(ends->div(ends[1]-1,3)!=div(ends[2]-1,3),graph.ends) ||
        throw(ArgumentError("NS c-recursion self-loop residues are not implemented"))
    C=Complex{T}
    mask=sum(primary_parities[e]<<(e-1) for e in 1:N)
    NSRecursionPlan{T,N}(graph,map(C,weights),primary_parities,mask,
        [ns_mode_phase(graph,primary_parities,m) for m in 0:2^N-1],[C(central)],
        Dict(),Dict(),Dict(),Dict(),Dict(),Dict(),Dict(),
        compact_theta && length(graph.slots)==2 && graph.slots[1]==graph.slots[2],0,0)
end

"""Stable NS Kac pole with removable x±1 factors canceled algebraically."""
function ns_c_pole(plan::NSRecursionPlan{T},edge::Int,shift::Int,r::Int,s::Int) where T
    r>=2 && s>=1 && iseven(r+s) || throw(ArgumentError("invalid NS Kac labels"))
    get!(plan.poles,(edge,shift,r,s)) do
        C=Complex{T}
        h=plan.weights[edge]+C(shift)/2
        affine=4h+r*s-1
        rad=sqrt(16h*h+8(r*s-1)*h+(r-s)^2)
        x1,x2=(affine+rad)/(1-r*r),(affine-rad)/(1-r*r)
        x=abs(x1)>=abs(x2) ? x1 : x2
        iszero(x) && throw(DomainError((edge,shift,r,s),"NS pole at infinity"))
        central=C(15//2)+3(x+inv(x))
        numerator=-12C(2)^(r*s-1)*x^(r*s-1)
        denominator=(1-r*r)*x*x-(1-s*s)
        factors=[(a,t) for a in 1-r:r for t in 1-s:s if iseven(a+t) && (a,t)!=(0,0) && (a,t)!=(r,s)]
        for target in (-1,1)
            index=findfirst(pair->pair[1]!=0 && pair[2]==pair[1]*target,factors)
            if isnothing(index)
                numerator*=x+target
            else
                denominator*=factors[index][1]
                deleteat!(factors,index)
            end
        end
        for (a,t) in factors
            denominator*=a*x+t
        end
        iszero(denominator) && throw(DomainError((edge,shift,r,s),"confluent NS pole requires a limiting prescription"))
        id=findfirst(==(central),plan.centers)
        if isnothing(id)
            push!(plan.centers,central)
            id=length(plan.centers)
        end
        # The central charge is invariant under x -> 1/x; fusion data are not.
        fusion_id=get!(plan.fusion_ids,x) do
            length(plan.fusion_ids)+1
        end
        NSPole{T}(central,numerator/denominator,x,id,fusion_id)
    end
end
"""NS fusion polynomial, pairing opposite momentum factors before evaluation."""
function ns_c_fusion(plan::NSRecursionPlan{T},pole::NSPole{T},r::Int,s::Int,a::Int,da::Int,b::Int,db::Int,f::Int) where T
    get!(plan.fusions,(pole.fusion_id,r,s,a,da,b,db,f)) do
        ha,hb=plan.weights[a]+Complex{T}(da)/2,plan.weights[b]+Complex{T}(db)/2
        x=pole.x
        lambda=x+2+inv(x)-8ha
        difference=8(hb-ha)
        value=one(Complex{T})
        for u in 1-r:2:r-1,v in 1-s:2:s-1
            mod(u+v-r-s,4)==(f==1 ? 0 : 2) || continue
            if u==0 && v==0
                value*=hb-ha
            elseif u>0 || (u==0 && v>0)
                shift=u*u*x+2u*v+v*v/x
                value*=((difference+shift)^2-4lambda*shift)/64
            end
        end
        value
    end
end

"""Cached reciprocal central-charge difference for NS null transport."""
function ns_c_denominator(plan::NSRecursionPlan{T},incoming::Int,outgoing::Int) where T
    get!(plan.denominators,(incoming,outgoing)) do
        a,b=plan.centers[incoming],plan.centers[outgoing]
        abs(a-b)>64eps(T)*max(one(T),abs(a),abs(b)) || throw(DomainError((a,b),"unresolved NS poles"))
        inv(a-b)
    end
end

"""Independent osp(1|2) descendant polynomial before the local BPZ phase."""
function ns_superglobal_core(weights::NTuple{3,C},levels::NTuple{3,Int}) where C<:Number
    h1,h2,h3=weights
    a,b,c=map(x->mod(x,2),levels)
    bits=4a+2b+c
    bottom=bits==6 ? h1+h2-h3 : bits==5 ? h1-h2+h3 : bits==3 ? h1-h2-h3 :
        bits==1 ? -one(C) : bits==7 ? -(h1+h2+h3-C(1//2)) : one(C)
    h1+=C(a)/2
    h2+=C(b)/2
    h3+=C(c)/2
    i,j,k=map(x->div(x,2),levels)
    value=zero(C)
    for p in 0:min(i,k)
        term=C(binomial(big(i),big(p)))
        for t in 0:p-1
            term*=(k-t)*(2h3+k-1-t)
        end
        value+=term*rising(h2+h3-h1,k-p)*rising(h1+h2-h3+p-k,i-p)
    end
    bottom*value*rising(h1-h2-h3+i-j+1-k,j)
end

"""Current ordered superglobal vertex, using real arithmetic only for exactly real weights."""
function ns_superglobal_vertex(weights::NTuple{3,C},levels::NTuple{3,Int},parities::NTuple{3,Int}) where C<:Number
    value=all(isreal,weights) ? ns_superglobal_core(map(real,weights),levels) : ns_superglobal_core(weights,levels)
    a,c=mod(levels[1],2),mod(levels[3],2)
    typeof(complex(zero(C)))(0,1)^a*parity_sign(parities[1]*a+parities[2]*c)*value
end

"""Cached global vertex with shifted weights and fixed intrinsic-primary convention."""
function ns_c_vertex(plan::NSRecursionPlan{T,N},v::Int,shift::NTuple{N,Int},level::NTuple{N,Int}) where {T,N}
    slots=plan.graph.slots[v]
    localshift=ntuple(j->shift[slots[j]],3)
    locallevel=ntuple(j->level[slots[j]],3)
    # On theta every (shift,level) is visited once and both vertices coincide.
    # Retaining these entries consumes memory without any later cache hits.
    plan.compact_theta && return ns_superglobal_vertex(
        ntuple(j->plan.weights[slots[j]]+Complex{T}(localshift[j])/2,3),
        locallevel,ntuple(j->plan.parities[slots[j]],3))
    get!(plan.vertices,(v,localshift...,locallevel...)) do
        ns_superglobal_vertex(ntuple(j->plan.weights[slots[j]]+Complex{T}(localshift[j])/2,3),
            locallevel,ntuple(j->plan.parities[slots[j]],3))
    end
end

"""Cached inverse native BPZ superglobal norm; each division is performed only once."""
function ns_c_inverse_norm(plan::NSRecursionPlan{T},edge::Int,shift::Int,level::Int) where T
    get!(plan.inverse_norms,(edge,shift,level)) do
        norm=get!(plan.norms,(edge,shift,level)) do
            n,bit=div(level,2),mod(level,2)
            weight=2plan.weights[edge]+shift
            value=isreal(weight) ? rising(real(weight),n+bit) : rising(weight,n+bit)
            for j in 2:n
                value*=j
            end
            value*=Complex{T}(0,1)^bit*parity_sign(plan.parities[edge]*bit)
            iszero(value) && throw(DomainError((edge,shift,level),"degenerate superglobal norm"))
            value
        end
        inv(norm)
    end
end

"""Full native-BPZ superglobal graph seed; no finite-c Gram block enters."""
function ns_c_global(plan::NSRecursionPlan{T,N},shift::NTuple{N,Int},level::NTuple{N,Int}) where {T,N}
    relative_mask=sum(mod(level[e],2)<<(e-1) for e in 1:N)
    value=Complex{T}(sewing_sign(plan.graph,xor(relative_mask,plan.primary_mask)))
    if plan.compact_theta
        vertex=ns_c_vertex(plan,1,shift,level)
        value*=vertex*vertex
    else
        for v in eachindex(plan.graph.slots)
            value*=ns_c_vertex(plan,v,shift,level)
        end
    end
    for e in 1:N
        inverse=ns_c_inverse_norm(plan,e,shift[e],level[e])
        if isreal(inverse)
            value*=real(inverse)
        elseif iszero(real(inverse))
            value=Complex{T}(-imag(value)*imag(inverse),real(value)*imag(inverse))
        else
            value*=inverse
        end
    end
    value
end

"""NS null residue with current vertex, inverse-Gram and Koszul phases included locally."""
function ns_c_residue(plan::NSRecursionPlan{T,N},pole::NSPole{T},shift::NTuple{N,Int},
                      edge::Int,r::Int,s::Int,remaining_mask::Int) where {T,N}
    value=pole.residue
    exponent=0
    odd=isodd(r*s)
    for half in plan.graph.ends[edge]
        vertex,slot=div(half-1,3)+1,mod(half-1,3)+1
        slots=plan.graph.slots[vertex]
        a,b=slot==1 ? (slots[3],slots[2]) : slot==2 ? (slots[3],slots[1]) : (slots[1],slots[2])
        form=mod(sum((remaining_mask>>(e-1))&1 for e in slots),2)
        value*=ns_c_fusion(plan,pole,r,s,a,shift[a],b,shift[b],form)
        if odd
            slot==2 && (exponent+=Int((remaining_mask>>(slots[3]-1))&1))
            slot==3 && (exponent+=1)
        end
    end
    child_mask=xor(remaining_mask,Int(odd)<<(edge-1))
    value*parity_sign(exponent)*sewing_sign(plan.graph,xor(remaining_mask,plan.primary_mask))*
        sewing_sign(plan.graph,xor(child_mask,plan.primary_mask))*
        Complex{T}(0,1)^(plan.phases[remaining_mask+1]-plan.phases[child_mask+1])
end

"""Batched independent SCA c-recursion, merged by output parity, null shift and incoming pole."""
function ns_c_reduced(plan::NSRecursionPlan{T,N},domain::CoefficientDomain{N};progress=(state)->nothing) where {T,N}
    C=Complex{T}
    origin=ntuple(_->0,N)
    out=Series{N,C}()
    maskof(key)=sum(mod(key[e],2)<<(e-1) for e in 1:N)
    for output_mask in unique(maskof.(domain.indices))
        progress((stage="poles",output_mask=output_mask,transitions=plan.transitions,
            seed_terms=plan.seed_terms,cached_vertices=length(plan.vertices)))
        amplitudes=Dict{NTuple{N,Int},Dict{Int,C}}(origin=>Dict(1=>one(C)))
        totals=Series{N,C}()
        for shift in domain.indices
            incoming=pop!(amplitudes,shift,nothing)
            isnothing(incoming) && continue
            totals[shift]=sum(values(incoming);init=zero(C))
            remaining_mask=xor(output_mask,maskof(shift))
            for edge in 1:N,r in 2:domain.limits[edge]-shift[edge],s in 1:div(domain.limits[edge]-shift[edge],r)
                iseven(r+s) || continue
                changed=Base.setindex(shift,shift[edge]+r*s,edge)
                contains(domain,changed) || continue
                pole=ns_c_pole(plan,edge,shift[edge],r,s)
                residue=ns_c_residue(plan,pole,shift,edge,r,s,remaining_mask)
                iszero(residue) && continue
                value=zero(C)
                for (id,amplitude) in incoming
                    iszero(amplitude) && continue
                    value+=amplitude*ns_c_denominator(plan,id,pole.id)
                    plan.transitions+=1
                end
                target=get!(() -> Dict{Int,C}(),amplitudes,changed)
                target[pole.id]=get(target,pole.id,zero(C))+residue*value
            end
        end
        progress((stage="seeds",output_mask=output_mask,transitions=plan.transitions,
            seed_terms=plan.seed_terms,cached_vertices=length(plan.vertices)))
        for level in domain.indices
            maskof(level)==output_mask || continue
            value=zero(C)
            for (shift,amplitude) in totals
                all(e->shift[e]<=level[e],1:N) || continue
                iszero(amplitude) && continue
                value+=amplitude*ns_c_global(plan,shift,ntuple(e->level[e]-shift[e],N))
                plan.seed_terms+=1
            end
            out[level]=value
        end
        progress((stage="completed_parity",output_mask=output_mask,transitions=plan.transitions,
            seed_terms=plan.seed_terms,cached_vertices=length(plan.vertices)))
    end
    out
end
