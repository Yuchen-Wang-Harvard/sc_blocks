"""Graph plan with NSNSNS/NSRR vertices, sharing the ordinary and inserted paths."""
mutable struct RamondBlockPlan{T<:AbstractFloat,N}
    graph::PlumbingGraph
    b::Complex{T}
    momenta::NTuple{N,Complex{T}}
    ramond::NTuple{N,Bool}
    primary_parities::NTuple{N,Int}
    forms::Vector{Int}
    eta::Vector{Int}
    marked::Vector{Int}
    actions::NTuple{N,BranchActions{T}}
    vertices::Vector{Union{Nothing,NSRRBranching{T}}}
    ns_vertices::Vector{Dict{NTuple{3,Int},Complex{T}}}
    ell_values::Dict{Tuple{Complex{T},Int},Complex{T}}
    auxiliary::FermionVertex{T}
    auxiliary_ns::FermionVertex{T}
    auxiliary_tensors::Dict{NTuple{4,Int},Array{Complex{T},3}}
    norms::Dict{NTuple{3,Int},Complex{T}}
    weights::Dict{NTuple{3,Int},Complex{T}}
    cache_products::Bool
    product_domain::Vector{NTuple{N,Int}}
    products::Dict{Tuple{Tuple,Tuple},Series{N,Complex{T}}}
    product_cache_hits::Int
    branches::Int
    ccy_transitions::Int
    forward_residual::T
    timings::Dict{String,Float64}
    recovery_tolerance::T
    recovery_policy::Symbol
    recovery_diagnostics::Dict{String,Any}
end

"""Accumulate wall time for a named numerical stage without changing its arithmetic."""
function record_stage!(plan::RamondBlockPlan,key::String,started::UInt64)
    plan.timings[key]=get(plan.timings,key,0.0)+(time_ns()-started)/1e9
    nothing
end

"""Accumulate CCY pole/seed substage measurements (included in the overall CCY time)."""
function record_ccy_stages!(plan::RamondBlockPlan,copies)
    plan.timings["ccy_poles"]=get(plan.timings,"ccy_poles",0.0)+sum(c.pole_seconds for c in copies)
    plan.timings["ccy_seeds"]=get(plan.timings,"ccy_seeds",0.0)+sum(c.seed_seconds for c in copies)
end

"""Change vertex labels while retaining momentum-local actions and all local caches."""
function set_vertex_labels!(plan::RamondBlockPlan,forms::AbstractVector{<:Integer},eta::AbstractVector{<:Integer})
    length(forms)==length(eta)==length(plan.graph.slots) || throw(DimensionMismatch("vertex labels"))
    all(x->x in (0,1),forms) && all(x->x in (-1,1),eta) || throw(ArgumentError("invalid vertex labels"))
    plan.forms.=forms
    plan.eta.=eta
    plan.marked=insertion_edges(plan.graph,collect(plan.ramond),eta)
    plan
end

"""Invalidate optional Virasoro-product reuse when the requested coefficient domain changes.

Products depend on momentum, branch labels, insertion topology and cutoff, but
not on physical vertex-form labels. The owning plan already fixes momentum and
the ordered graph. Copy the domain indices so later caller mutations cannot
silently reuse a product computed at another cutoff.
"""
function prepare_product_cache!(plan::RamondBlockPlan,domain::CoefficientDomain)
    plan.cache_products || return nothing
    if plan.product_domain != domain.indices
        empty!(plan.products)
        empty!(plan.product_domain)
        append!(plan.product_domain,domain.indices)
    end
    nothing
end

"""Construct a mixed graph and select one insertion edge on each odd Ramond cycle."""
function RamondBlockPlan(::Type{T},graph::PlumbingGraph,b,momenta::NTuple{N},ramond::NTuple{N,Bool},
                         forms::AbstractVector{<:Integer},eta::AbstractVector{<:Integer};
                         primary_parities::NTuple{N,Int}=ntuple(_->0,N),
                         recovery_tolerance=inv(T(10)^8),recovery_policy::Symbol=:error,
                         cache_products::Bool=false) where {T<:AbstractFloat,N}
    N==length(graph.ends) && length(forms)==length(eta)==length(graph.slots) || throw(DimensionMismatch("graph data"))
    all(x->x in (0,1),forms) && all(x->x in (0,1),primary_parities) || throw(ArgumentError("invalid parity"))
    all(x->x in (-1,1),eta) || throw(ArgumentError("eta must be ±1"))
    recovery_policy in (:error,:record) || throw(ArgumentError("recovery policy must be error or record"))
    isfinite(recovery_tolerance) && recovery_tolerance>0 || throw(ArgumentError("positive finite recovery tolerance required"))
    all(>(0),Iterators.flatten(graph.slots)) || throw(ArgumentError("ordinary plan requires a closed graph"))
    for slots in graph.slots
        map(e->ramond[e],slots) in ((false,false,false),(false,true,true)) ||
            throw(ArgumentError("each vertex must be NSNSNS or NSRR, with NS in slot one"))
    end
    marked=insertion_edges(graph,collect(ramond),eta)
    all(e->!ramond[e] || primary_parities[e]==0,1:N) || throw(ArgumentError("Ramond parity belongs to its ground bit"))
    C=Complex{T}
    coupling=C(b)
    momenta_=map(C,momenta)
    combinatorics=Combinatorics()
    actions=ntuple(e->BranchActions(T,ramond[e] ? :R : :NS,coupling,momenta_[e];
        parity=primary_parities[e],plan=combinatorics),N)
    local_vertices=Dict{NTuple{3,Int},NSRRBranching{T}}()
    vertices=Union{Nothing,NSRRBranching{T}}[ramond[slots[2]] ?
        get!(() -> NSRRBranching(ntuple(j->actions[slots[j]],3)),local_vertices,slots) : nothing
        for slots in graph.slots]
    auxns,auxr=FermionModule(T,:NS;plan=combinatorics),FermionModule(T,:R;plan=combinatorics)
    RamondBlockPlan{T,N}(graph,coupling,momenta_,ramond,primary_parities,Int.(forms),Int.(eta),marked,actions,vertices,
        [Dict{NTuple{3,Int},C}() for _ in graph.slots],Dict(),
        FermionVertex((auxns,auxr,auxr)),FermionVertex((auxns,auxns,auxns)),Dict(),Dict(),Dict(),
        cache_products,NTuple{N,Int}[],Dict(),0,0,0,zero(T),Dict(),
        T(recovery_tolerance),recovery_policy,Dict())
end

"""Doubled descendant level of a branching primary with label n4=4n."""
branch_level2(ramond::Bool,n4::Int)=ramond ? div(n4*n4-1,4) : div(n4*n4,4)

"""Enumerate only branch label tuples whose primary levels fit the requested domain."""
function ramond_branch_labels(plan::RamondBlockPlan{T,N},domain::CoefficientDomain{N}) where {T,N}
    choices=ntuple(e->begin
        if plan.ramond[e]
            bound=isqrt(4domain.limits[e]+1)
            iseven(bound) && (bound-=1)
            collect(-bound:2:bound)
        else
            bound=2isqrt(domain.limits[e])
            collect(-bound:2:bound)
        end
    end,N)
    out=NTuple{N,Int}[]
    for labels in Iterators.product(choices...)
        shift=ntuple(e->branch_level2(plan.ramond[e],labels[e]),N)
        contains(domain,shift) && push!(out,labels)
    end
    out
end

"""Cached unnormalized chi-vector norm for graph sewing."""
function ramond_plan_norm(plan::RamondBlockPlan,e::Int,label::Int,alpha::Int)
    get!(plan.norms,(e,label,alpha)) do
        plan.ramond[e] ? ramond_branch_norm(plan.b,plan.momenta[e],label,alpha) :
            ns_branch_norm(plan.b,plan.momenta[e],div(label,2),plan.primary_parities[e])
    end
end

"""Cached embedded branch weight independent of total Ramond parity."""
function ramond_plan_weight(plan::RamondBlockPlan,e::Int,label::Int,copy::Int)
    get!(plan.weights,(e,label,copy)) do
        branch_weight(plan.b,plan.momenta[e],label,copy)
    end
end

"""Reuse the universal ordered auxiliary NSRR tensor at each local grade triple."""
function ramond_auxiliary_tensor(plan::RamondBlockPlan{T},levels::NTuple{3,Int};ramond::Bool=true) where T
    get!(plan.auxiliary_tensors,(Int(ramond),levels...)) do
        vertex=ramond ? plan.auxiliary : plan.auxiliary_ns
        bases=ntuple(j->basis(vertex.modules[j],levels[j]),3)
        tensor=zeros(Complex{T},map(length,bases))
        for ix in CartesianIndices(tensor)
            tensor[ix]=rho(vertex,ntuple(j->bases[j][ix[j]],3))
        end
        tensor
    end
end

"""Direct auxiliary Fock block, retaining Ramond ground parity components and spin support."""
function ramond_auxiliary_block(plan::RamondBlockPlan{T,N},domain::CoefficientDomain{N}) where {T,N}
    out=ParitySeries{N,T}()
    modules=ntuple(e->plan.auxiliary.modules[plan.ramond[e] ? 2 : 1],N)
    for level in domain.indices
        any(e->plan.ramond[e] && isodd(level[e]),1:N) && continue
        bases=ntuple(e->basis(modules[e],level[e]),N)
        any(isempty,bases) && continue
        inverses=[inverse_gram(modules[e],level[e]) for e in 1:N]
        tensors=[ramond_auxiliary_tensor(plan,ntuple(j->level[slots[j]],3);ramond=plan.ramond[slots[2]])
            for slots in plan.graph.slots]
        selection=[ntuple(p->findall(s->Int(relative_parity(s))==p-1,bases[e]),2) for e in 1:N]
        values=zeros(Complex{T},2^N)
        for mask in 0:2^N-1
            all(slots->iseven(sum((mask>>(e-1))&1 for e in slots)),plan.graph.slots) || continue
            selected=[selection[e][((mask>>(e-1))&1)+1] for e in 1:N]
            values[mask+1]=sewing_sign(plan.graph,mask)*sew_tensors(plan.graph,tensors,inverses,selected)
        end
        any(!iszero,values) && (out[level]=values)
    end
    out
end

"""Raw local vertex in a common chi basis, with closed NS and recursive NSRR evaluation."""
function graph_branch_value(plan::RamondBlockPlan,v::Int,labels::NTuple{3,Int},parities::NTuple{3,Int})
    vertex=plan.vertices[v]
    if isnothing(vertex)
        relative_form=mod(sum(div(label,2) for label in labels),2)
        relative_form==plan.forms[v] || return zero(plan.b)
        return get!(plan.ns_vertices[v],labels) do
            slots=plan.graph.slots[v]
            ns_branching_raw(plan.b,ntuple(j->plan.momenta[slots[j]],3),map(n->div(n,2),labels),
                ntuple(j->plan.primary_parities[slots[j]],3),relative_form;ell_values=plan.ell_values)
        end
    end
    branching_raw(vertex,labels,parities[2],parities[3],plan.forms[v],plan.eta[v])
end

"""Raw branching vertex product divided by one full branch norm per edge.

No square-root gauge enters: each edge carries its unnormalized primary at both
ends, and its BPZ norm is divided out exactly once.
"""
function ramond_branch_coefficients(plan::RamondBlockPlan{T,N},labels::NTuple{N,Int}) where {T,N}
    out=zeros(Complex{T},2^N)
    for mask in 0:2^N-1
        all(e->plan.ramond[e] || ((mask>>(e-1))&1)==mod(div(labels[e],2)+plan.primary_parities[e],2),1:N) || continue
        coefficient=one(Complex{T})
        for (v,slots) in enumerate(plan.graph.slots)
            alpha,gamma=Int((mask>>(slots[2]-1))&1),Int((mask>>(slots[3]-1))&1)
            coefficient*=graph_branch_value(plan,v,ntuple(j->labels[slots[j]],3),(0,alpha,gamma))
            iszero(coefficient) && break
        end
        iszero(coefficient) && continue
        for e in 1:N
            alpha=plan.ramond[e] ? Int((mask>>(e-1))&1) : 0
            coefficient/=ramond_plan_norm(plan,e,labels[e],alpha)
        end
        out[mask+1]=sewing_sign(plan.graph,mask)*coefficient
    end
    out
end

"""Ordinary NSRR graph block from branching recurrence and two forward CCY factors.

Both Ramond parities reuse each Virasoro computation. The algorithm returns the
physical parity coefficients directly after allowed-sector auxiliary recovery;
there is no post-processing convention conversion.
"""
function ramond_block(plan::RamondBlockPlan{T,N},domain::CoefficientDomain{N}) where {T,N}
    empty!(plan.timings)
    isempty(plan.marked) || return inserted_ramond_block(SplitBlockPlan(plan),domain)
    prepare_product_cache!(plan,domain)
    C=Complex{T}
    started=time_ns()
    labels_=ramond_branch_labels(plan,domain)
    for (v,slots) in enumerate(plan.graph.slots)
        isnothing(plan.vertices[v]) && continue
        requests=unique([ntuple(j->labels[slots[j]],3) for labels in labels_])
        prepare_branching!(plan.vertices[v],requests)
    end
    record_stage!(plan,"branching",started)
    numerator=ParitySeries{N,T}()
    for labels in labels_
        started=time_ns()
        coefficients=ramond_branch_coefficients(plan,labels)
        record_stage!(plan,"branch_coefficients",started)
        any(!iszero,coefficients) || continue
        shift=ntuple(e->branch_level2(plan.ramond[e],labels[e]),N)
        key=((),labels)
        product=plan.cache_products ? get(plan.products,key,nothing) : nothing
        if isnothing(product)
            started=time_ns()
            remaining=residual_domain(domain,shift;scale=2)
            copies=ntuple(copy->CCYPlan(T,plan.graph,branch_central(plan.b,copy),
                [ramond_plan_weight(plan,e,labels[e],copy) for e in 1:N]),2)
            first,second=ccy_reduced(copies[1],remaining),ccy_reduced(copies[2],remaining)
            record_stage!(plan,"ccy",started)
            record_ccy_stages!(plan,copies)
            plan.ccy_transitions+=copies[1].transitions+copies[2].transitions
            started=time_ns()
            product=series_multiply(first,second,remaining)
            plan.cache_products && (plan.products[key]=product)
            record_stage!(plan,"products",started)
        else
            plan.product_cache_hits+=1
        end
        started=time_ns()
        for (degree,value) in product
            target=ntuple(e->shift[e]+2degree[e],N)
            values=get!(() -> zeros(C,2^N),numerator,target)
            values .+= value.*coefficients
        end
        record_stage!(plan,"assembly",started)
        plan.branches+=1
    end
    started=time_ns()
    auxiliary=ramond_auxiliary_block(plan,domain)
    record_stage!(plan,"auxiliary",started)
    started=time_ns()
    reduced,residual=recover_allowed_sector(numerator,auxiliary,domain;
        tolerance=plan.recovery_tolerance,policy=plan.recovery_policy,diagnostics=plan.recovery_diagnostics)
    record_stage!(plan,"recovery_and_forward_check",started)
    plan.forward_residual=max(plan.forward_residual,residual)
    started=time_ns()
    integer_domain=residual_domain(domain,ntuple(_->0,N);scale=2)
    vacuum=schottky_vacuum(plan.graph,integer_domain)
    squared=series_multiply(vacuum,vacuum,integer_domain)
    lifted=Series{N,C}(map(n->2n,key)=>C(value) for (key,value) in squared)
    result=parity_scalar_product(reduced,lifted,domain)
    record_stage!(plan,"vacuum_and_restoration",started)
    result
end
