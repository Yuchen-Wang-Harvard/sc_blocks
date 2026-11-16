"""Central-charge pole geometry and its universal derivative/residue prefactor."""
struct CCYPole{T<:AbstractFloat}
    central::Complex{T}
    residue::Complex{T}
    x::Complex{T}
    id::Int
    fusion_id::Int
end

"""Momentum-local graph CCY plan with shared pole, fusion, denominator and global caches."""
mutable struct CCYPlan{T<:AbstractFloat,N}
    graph::PlumbingGraph
    weights::Vector{Complex{T}}
    slots::Vector{NTuple{3,Int}}
    vertex_representatives::Vector{Int}
    norm_owners::Vector{Vector{Int}}
    centers::Vector{Complex{T}}
    fusion_ids::Dict{Complex{T},Int}
    geometry::Dict{NTuple{4,Int},CCYPole{T}}
    fusions::Dict{NTuple{7,Int},Complex{T}}
    denominators::Dict{NTuple{2,Int},Complex{T}}
    norms::Dict{NTuple{3,Int},Complex{T}}
    vertices::Dict{NTuple{7,Int},Complex{T}}
    vertex_cores::Dict{NTuple{6,Int},Complex{T}}
    normalized_vertices::Dict{NTuple{7,Int},Complex{T}}
    transitions::Int
    seed_terms::Int
    pole_seconds::Float64
    seed_seconds::Float64
end

"""Build a standalone CCY plan; external weights are keyed by half-edge ID."""
function CCYPlan(::Type{T},graph::PlumbingGraph,central,weights::AbstractVector;
                 external=Dict{Int,Int}()) where T<:AbstractFloat
    N=length(graph.ends)
    length(weights)==N || throw(DimensionMismatch("internal weights"))
    C=Complex{T}
    allweights=C.(weights)
    slots=copy(graph.slots)
    for v in eachindex(slots),j in 1:3
        slots[v][j]==0 || continue
        half=3(v-1)+j
        haskey(external,half) || throw(ArgumentError("missing external weight at half-edge $half"))
        push!(allweights,C(external[half]))
        slots[v]=Base.setindex(slots[v],length(allweights),j)
    end
    representatives=[something(findfirst(==(slots[v]),slots)) for v in eachindex(slots)]
    owners=[Int[] for _ in slots]
    for e in 1:N
        push!(owners[div(last(graph.ends[e])-1,3)+1],e)
    end
    CCYPlan{T,N}(graph,allweights,slots,representatives,owners,[C(central)],
        Dict(),Dict(),Dict(),Dict(),Dict(),Dict(),Dict(),Dict(),0,0,0.0,0.0)
end

"""Cached c_(r,s)(h), b_(r,s)^2 and regularized universal Kac residue."""
function ccy_pole(plan::CCYPlan{T},edge::Int,shift::Int,r::Int,s::Int) where T
    r>=2 && s>=1 || throw(ArgumentError("CCY poles require r>=2,s>=1"))
    get!(plan.geometry,(edge,shift,r,s)) do
        C=Complex{T}
        h=plan.weights[edge]+shift
        rad=C((r-s)^2)+4(r*s-1)*h+4h^2
        base=C(r*s-1)+2h
        x1,x2=(base+sqrt(rad))/(1-r*r),(base-sqrt(rad))/(1-r*r)
        x=abs(x1)>=abs(x2) ? x1 : x2
        iszero(x) && throw(DomainError((edge,shift,r,s),"CCY pole at infinity"))
        central=C(13)+6(x+inv(x))
        numerator=-12x^(2r*s-1)
        denominator=(1-r*r)*x*x-(1-s*s)
        factors=[(a,b) for a in 1-r:r for b in 1-s:s if (a,b)!=(0,0) && (a,b)!=(r,s)]
        # Cancel x±1 algebraically before numerical evaluation, as in the
        # optimized engine. This resolves removable zeros at these factors.
        for target in (-1,1)
            ix=findfirst(pair->first(pair)!=0 && last(pair)==first(pair)*target,factors)
            if isnothing(ix)
                numerator*=x+target
            else
                denominator*=first(factors[ix])
                deleteat!(factors,ix)
            end
        end
        for (a,b) in factors
            denominator*=a*x+b
        end
        iszero(denominator) && throw(DomainError((edge,shift,r,s),"confluent Kac residue requires a limiting prescription"))
        id=findfirst(==(central),plan.centers)
        if isnothing(id)
            push!(plan.centers,central)
            id=length(plan.centers)
        end
        # Reciprocal x values share a central charge but not their fusion data.
        # Keep integer fusion keys and reuse them only for identical x.
        fusion_id=get!(plan.fusion_ids,x) do
            length(plan.fusion_ids)+1
        end
        CCYPole{T}(central,numerator/denominator,x,id,fusion_id)
    end
end

"""Cached fusion polynomial, using paired factors to avoid auxiliary square-root choices."""
function ccy_fusion(plan::CCYPlan{T},pole::CCYPole{T},r::Int,s::Int,
                    top::Int,ts::Int,bottom::Int,bs::Int) where T
    get!(plan.fusions,(pole.fusion_id,r,s,top,ts,bottom,bs)) do
        a,b=plan.weights[top]+ts,plan.weights[bottom]+bs
        x=pole.x
        lambda=x+2+inv(x)-4a
        difference=4(b-a)
        value=one(Complex{T})
        for u in 1-r:2:r-1,t in 1-s:2:s-1
            if u==0 && t==0
                value*=b-a
            elseif u>0 || (u==0 && t>0)
                shift=u*u*x+2u*t+t*t/x
                value*=((difference+shift)^2-4lambda*shift)/16
            end
        end
        value
    end
end

"""Cache reciprocal pole differences, rejecting numerically unresolved collisions."""
function ccy_denominator(plan::CCYPlan{T},incoming::Int,outgoing::Int) where T
    get!(plan.denominators,(incoming,outgoing)) do
        a,b=plan.centers[incoming],plan.centers[outgoing]
        difference=a-b
        scale=max(one(T),abs(a),abs(b))
        abs(difference)>64eps(T)*scale ||
            throw(DomainError((a,b),"unresolved CCY poles; increase precision or specify a generic limit"))
        inv(difference)
    end
end

"""Global sl(2) descendant norm n!(2h)_n, cached across all graph seed terms."""
function global_norm(plan::CCYPlan{T},edge::Int,shift::Int,n::Int) where T
    get!(plan.norms,(edge,shift,n)) do
        n==0 && return one(Complex{T})
        h=plan.weights[edge]+shift
        value=isreal(h) ? rising(2real(h),n) : rising(2h,n)
        for j in 2:n
            value*=j
        end
        iszero(value) && throw(DomainError((edge,shift,n),"degenerate global Gram norm"))
        Complex{T}(value)
    end
end

"""Middle-level-independent core of the global three-point form."""
function global_vertex_core(weights::NTuple{3,C},i::Int,k::Int) where C<:Number
    if C<:Complex && all(isreal,weights)
        return C(global_vertex_core(map(real,weights),i,k))
    end
    h1,h2,h3=weights
    i==0 && k==0 && return one(C)
    i==0 && return rising(h3+h2-h1,k)
    k==0 && return rising(h1+h2-h3,i)
    value=zero(C)
    rise=Vector{C}(undef,k+1)
    suffix=Vector{C}(undef,i+1)
    rise[1],suffix[end]=one(C),one(C)
    rise_base,suffix_base,fall_base=h3+h2-h1,h1+h2-h3-k,2h3+k-1
    for t in 0:k-1
        rise[t+2]=rise[t+1]*(rise_base+t)
    end
    for t in i-1:-1:0
        suffix[t+1]=suffix[t+2]*(suffix_base+t)
    end
    falling_h,falling_k,choose=one(C),big(1),big(1)
    for p in 0:min(i,k)
        value+=(choose*falling_k)*falling_h*rise[k-p+1]*suffix[p+1]
        p==min(i,k) && break
        falling_h*=fall_base-p
        falling_k*=k-p
        choose=div(choose*(i-p),p+1)
    end
    value
end

"""Closed primary-normalized three-point value of L_-1 descendants."""
function global_vertex(weights::NTuple{3,C},levels::NTuple{3,Int}) where C<:Number
    C<:Complex && all(isreal,weights) && return C(global_vertex(map(real,weights),levels))
    h1,h2,h3=weights
    i,j,k=levels
    core=global_vertex_core(weights,i,k)
    j==0 ? core : core*rising(h1-h2-h3+i-j+1-k,j)
end

"""Complete local global factor, shared by all spectator-edge choices."""
function ccy_vertex(plan::CCYPlan{T,N},v::Int,shift::NTuple{N,Int},level::NTuple{N,Int}) where {T,N}
    edges=plan.slots[v]
    shifts=ntuple(j->edges[j]<=N ? shift[edges[j]] : 0,3)
    levels=ntuple(j->edges[j]<=N ? level[edges[j]] : 0,3)
    get!(plan.vertices,(plan.vertex_representatives[v],shifts...,levels...)) do
        weights=ntuple(j->plan.weights[edges[j]]+shifts[j],3)
        i,j,k=levels
        core=get!(plan.vertex_cores,(plan.vertex_representatives[v],shifts...,i,k)) do
            global_vertex_core(weights,i,k)
        end
        j==0 && return core
        if all(isreal,weights)
            return Complex{T}(real(core)*rising(real(weights[1])-real(weights[2])-real(weights[3])+i-j+1-k,j))
        end
        core*rising(weights[1]-weights[2]-weights[3]+i-j+1-k,j)
    end
end

"""Cache each local vertex together with its assigned inverse edge norms.

Every internal edge is assigned once, including self-loops. Normalization is
independent of spectator edges and is never repeated in the global seed sum.
"""
function ccy_normalized_vertex(plan::CCYPlan{T,N},v::Int,shift::NTuple{N,Int},level::NTuple{N,Int}) where {T,N}
    edges=plan.slots[v]
    shifts=ntuple(j->edges[j]<=N ? shift[edges[j]] : 0,3)
    levels=ntuple(j->edges[j]<=N ? level[edges[j]] : 0,3)
    get!(plan.normalized_vertices,(v,shifts...,levels...)) do
        value=ccy_vertex(plan,v,shift,level)
        for e in plan.norm_owners[v]
            level[e]==0 && continue
            norm=global_norm(plan,e,shift[e],level[e])
            value=isreal(norm) ? value/real(norm) : value/norm
        end
        value
    end
end

"""Weight-dependent global seed; the Schottky vacuum is restored separately once."""
function ccy_global(plan::CCYPlan{T,N},shift::NTuple{N,Int},level::NTuple{N,Int}) where {T,N}
    value=ccy_normalized_vertex(plan,1,shift,level)
    for v in 2:length(plan.slots)
        value*=ccy_normalized_vertex(plan,v,shift,level)
    end
    value
end

"""Graph residue, including the paired fusion factors for a Ramond-style self-loop ordering."""
function ccy_residue(plan::CCYPlan{T,N},pole::CCYPole{T},shift::NTuple{N,Int},edge::Int,r::Int,s::Int) where {T,N}
    endpoints=plan.graph.ends[edge]
    v1,v2=map(h->div(h-1,3)+1,endpoints)
    value=pole.residue
    if v1==v2
        sort([mod(h-1,3)+1 for h in endpoints])==[2,3] ||
            throw(ArgumentError("self-loop CCY chart requires loop slots 2 and 3"))
        external=plan.slots[v1][1]
        es=external<=N ? shift[external] : 0
        return parity_sign(r*s)*value*
            ccy_fusion(plan,pole,r,s,external,es,edge,shift[edge]+r*s)*
            ccy_fusion(plan,pole,r,s,external,es,edge,shift[edge])
    end
    for half in endpoints
        v,slot=div(half-1,3)+1,mod(half-1,3)+1
        edges=plan.slots[v]
        a,b=slot==1 ? (edges[3],edges[2]) : slot==2 ? (edges[3],edges[1]) : (edges[1],edges[2])
        as,bs=a<=N ? shift[a] : 0,b<=N ? shift[b] : 0
        value*=ccy_fusion(plan,pole,r,s,a,as,b,bs)
    end
    value
end

"""Attach global seeds with ordinary arithmetic for non-BigFloat plans."""
function ccy_seed_coefficient(plan::CCYPlan{T,N},totals,n::NTuple{N,Int},::Nothing) where {T,N}
    value=zero(Complex{T})
    for (shift,amplitude) in totals
        all(e->shift[e]<=n[e],1:N) || continue
        iszero(amplitude) && continue
        value+=amplitude*ccy_global(plan,shift,ntuple(e->n[e]-shift[e],N))
        plan.seed_terms+=1
    end
    value
end

"""Attach cached normalized seeds while reusing MPFR mantissas in the inner sum."""
function ccy_seed_coefficient(plan::CCYPlan{BigFloat,N},totals,n::NTuple{N,Int},work::BigComplexScratch) where N
    value=Complex(BigFloat(0),BigFloat(0))
    for (shift,amplitude) in totals
        all(e->shift[e]<=n[e],1:N) || continue
        iszero(amplitude) && continue
        level=ntuple(e->n[e]-shift[e],N)
        scratch_product!(work.product,amplitude,ccy_normalized_vertex(plan,1,shift,level),work)
        for v in 2:length(plan.slots)
            scratch_product!(work.product,work.product,ccy_normalized_vertex(plan,v,shift,level),work)
        end
        scratch_add!(value,work.product)
        plan.seed_terms+=1
    end
    value
end

"""Batched forward CCY recursion, merging every path with the same shift and incoming pole."""
function ccy_reduced(plan::CCYPlan{T,N},domain::CoefficientDomain{N}) where {T,N}
    started=time_ns()
    C=Complex{T}
    origin=ntuple(_->0,N)
    amplitudes=Dict{NTuple{N,Int},Dict{Int,C}}(origin=>Dict(1=>one(C)))
    totals=Series{N,C}()
    work=complex_scratch(C)
    for shift in domain.indices
        incoming=pop!(amplitudes,shift,nothing)
        isnothing(incoming) && continue
        totals[shift]=sum(values(incoming);init=zero(C))
        for edge in 1:N,r in 2:domain.limits[edge]-shift[edge],
            s in 1:div(domain.limits[edge]-shift[edge],r)
            changed=Base.setindex(shift,shift[edge]+r*s,edge)
            contains(domain,changed) || continue
            pole=ccy_pole(plan,edge,shift[edge],r,s)
            residue=ccy_residue(plan,pole,shift,edge,r,s)
            iszero(residue) && continue
            value=zero(C)
            for (id,amplitude) in incoming
                iszero(amplitude) && continue
                value=sum_product(value,amplitude,ccy_denominator(plan,id,pole.id),work)
                plan.transitions+=1
            end
            target=get!(() -> Dict{Int,C}(),amplitudes,changed)
            previous=get!(() -> zero(C),target,pole.id)
            target[pole.id]=sum_product(previous,residue,value,work)
        end
    end
    plan.pole_seconds+=(time_ns()-started)/1e9
    started=time_ns()
    out=Series{N,C}()
    for n in domain.indices
        out[n]=ccy_seed_coefficient(plan,totals,n,work)
    end
    plan.seed_seconds+=(time_ns()-started)/1e9
    out
end

"""Full Virasoro graph block: forward CCY followed by one exact Schottky convolution."""
function virasoro_block(plan::CCYPlan{T,N},domain::CoefficientDomain{N};vacuum=nothing) where {T,N}
    reduced=ccy_reduced(plan,domain)
    qvac=isnothing(vacuum) ? schottky_vacuum(plan.graph,domain) : vacuum
    cvac=Series{N,Complex{T}}(key=>Complex{T}(value) for (key,value) in qvac)
    series_multiply(reduced,cvac,domain)
end
