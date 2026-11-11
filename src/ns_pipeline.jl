"""All-NS production plan with closed branching, graph CCY and direct auxiliary recovery."""
mutable struct NSBlockPlan{T<:AbstractFloat,N}
    graph::PlumbingGraph
    b::Complex{T}
    momenta::NTuple{N,Complex{T}}
    primary_parities::NTuple{N,Int}
    forms::Vector{Int}
    vertices::Vector{Dict{NTuple{3,Int},Complex{T}}}
    ell_values::Dict{Tuple{Complex{T},Int},Complex{T}}
    branch_weights::Dict{NTuple{3,Int},Complex{T}}
    auxiliary::FermionVertex{T}
    auxiliary_tensors::Dict{NTuple{3,Int},Array{Complex{T},3}}
    branches::Int
    ccy_transitions::Int
end

"""Construct an all-NS closed-graph plan with explicit intrinsic edge parities."""
function NSBlockPlan(::Type{T},graph::PlumbingGraph,b,momenta::NTuple{N},forms::AbstractVector{<:Integer};
                     primary_parities::NTuple{N,Int}=ntuple(_->0,N)) where {T<:AbstractFloat,N}
    N==length(graph.ends) && length(forms)==length(graph.slots) || throw(DimensionMismatch("graph data"))
    all(>(0),Iterators.flatten(graph.slots)) || throw(ArgumentError("NS block plan currently requires a closed graph"))
    all(x->x in (0,1),forms) && all(x->x in (0,1),primary_parities) || throw(ArgumentError("invalid NS parity"))
    C=Complex{T}
    coupling=C(b)
    coupling in (zero(C),one(C),-one(C)) && throw(ArgumentError("b must differ from 0,±1"))
    aux=FermionModule(T,:NS)
    NSBlockPlan{T,N}(graph,coupling,map(C,momenta),primary_parities,Int.(forms),
        [Dict{NTuple{3,Int},C}() for _ in graph.slots],Dict(),Dict(),FermionVertex((aux,aux,aux)),Dict(),0,0)
end

"""Cached closed all-NS branching coefficient for one vertex and its incident labels."""
function ns_local_branch(plan::NSBlockPlan,v::Int,n2::NTuple{3,Int})
    get!(plan.vertices[v],n2) do
        edges=plan.graph.slots[v]
        ns_branching(plan.b,ntuple(j->plan.momenta[edges[j]],3),n2,
                      ntuple(j->plan.primary_parities[edges[j]],3),plan.forms[v];ell_values=plan.ell_values)
    end
end

"""Cache embedded edge weights across every branch with the same local label."""
function ns_branch_weight(plan::NSBlockPlan,e::Int,n2::Int,copy::Int)
    get!(plan.branch_weights,(e,n2,copy)) do
        branch_weight(plan.b,plan.momenta[e],2n2,copy)
    end
end

"""Reuse auxiliary three-point tensors at each triple of doubled grades."""
function ns_auxiliary_tensor(plan::NSBlockPlan{T},levels::NTuple{3,Int}) where T
    get!(plan.auxiliary_tensors,levels) do
        bases=ntuple(j->basis(plan.auxiliary.modules[j],levels[j]),3)
        out=zeros(Complex{T},map(length,bases))
        for ix in CartesianIndices(out)
            out[ix]=rho(plan.auxiliary,ntuple(j->bases[j][ix[j]],3))
        end
        out
    end
end

"""Direct auxiliary all-NS block in doubled plumbing powers, with its own Gram and spin parity."""
function ns_auxiliary_block(plan::NSBlockPlan{T,N},domain::CoefficientDomain{N}) where {T,N}
    C=Complex{T}
    out=Series{N,C}()
    module_=plan.auxiliary.modules[1]
    for level in domain.indices
        all(slots->iseven(sum(level[e] for e in slots)),plan.graph.slots) || continue
        bases=[basis(module_,level[e]) for e in 1:N]
        any(isempty,bases) && continue
        tensors=[ns_auxiliary_tensor(plan,ntuple(j->level[slots[j]],3)) for slots in plan.graph.slots]
        inverses=[inverse_gram(module_,level[e]) for e in 1:N]
        selections=[collect(eachindex(bases[e])) for e in 1:N]
        mask=sum(mod(level[e],2)<<(e-1) for e in 1:N)
        value=sewing_sign(plan.graph,mask)*sew_tensors(plan.graph,tensors,inverses,selections)
        iszero(value) || (out[level]=value)
    end
    out
end

"""Triangular coefficientwise division by a unit auxiliary series."""
function recover_unit(numerator::Series{N,T},auxiliary::Series{N,T},domain::CoefficientDomain{N}) where {N,T}
    origin=ntuple(_->0,N)
    constant=get(auxiliary,origin,zero(T))
    iszero(constant) && throw(ArgumentError("auxiliary constant is zero; allowed-sector or split recovery is required"))
    out=Series{N,T}()
    for level in domain.indices
        value=get(numerator,level,zero(T))
        for (shift,coefficient) in auxiliary
            shift==origin && continue
            all(e->shift[e]<=level[e],1:N) || continue
            value-=coefficient*get(out,ntuple(e->level[e]-shift[e],N),zero(T))
        end
        out[level]=value/constant
    end
    out
end

"""Production all-NS graph block; indices are doubled descendant levels.

Closed Appendix-B vertices feed two forward CCY recursions per allowed branch.
Each receives only its remaining domain. The two universal vacuum factors are
restored once after branch summation and auxiliary division. Returned scalar
coefficients include the graph Koszul sign; edge spin evaluation is separate.
"""
function ns_block(plan::NSBlockPlan{T,N},domain::CoefficientDomain{N}) where {T,N}
    C=Complex{T}
    reduced=Series{N,C}()
    bounds=ntuple(e->isqrt(domain.limits[e]),N)
    ranges=ntuple(e->-bounds[e]:bounds[e],N)
    for labels in Iterators.product(ranges...)
        shift=ntuple(e->labels[e]^2,N)
        contains(domain,shift) || continue
        coefficient=one(C)
        for v in eachindex(plan.graph.slots)
            edges=plan.graph.slots[v]
            coefficient*=ns_local_branch(plan,v,ntuple(j->labels[edges[j]],3))
            iszero(coefficient) && break
        end
        iszero(coefficient) && continue
        mask=sum(mod(labels[e]+plan.primary_parities[e],2)<<(e-1) for e in 1:N)
        coefficient*=sewing_sign(plan.graph,mask)
        remaining=residual_domain(domain,shift;scale=2)
        weights=ntuple(copy->[ns_branch_weight(plan,e,labels[e],copy) for e in 1:N],2)
        copies=ntuple(copy->CCYPlan(T,plan.graph,branch_central(plan.b,copy),weights[copy]),2)
        product=series_multiply(ccy_reduced(copies[1],remaining),ccy_reduced(copies[2],remaining),remaining)
        for (degree,value) in product
            target=ntuple(e->shift[e]+2degree[e],N)
            reduced[target]=get(reduced,target,zero(C))+coefficient*value
        end
        plan.branches+=1
        plan.ccy_transitions+=copies[1].transitions+copies[2].transitions
    end
    auxiliary=ns_auxiliary_block(plan,domain)
    physical_reduced=recover_unit(reduced,auxiliary,domain)
    integer_domain=residual_domain(domain,ntuple(_->0,N);scale=2)
    vacuum=schottky_vacuum(plan.graph,integer_domain)
    squared=series_multiply(vacuum,vacuum,integer_domain)
    lifted=Series{N,C}(map(n->2n,key)=>C(value) for (key,value) in squared)
    series_multiply(physical_reduced,lifted,domain)
end

"""Evaluate one all-NS coefficient at edge spin signs, including primary parity."""
function ns_spin_coefficient(plan::NSBlockPlan{T,N},coefficient,level::NTuple{N,Int},signs::NTuple{N,Int}) where {T,N}
    all(s->s in (-1,1),signs) || throw(ArgumentError("spin signs must be ±1"))
    coefficient*prod(signs[e]^mod(level[e]+plan.primary_parities[e],2) for e in 1:N)
end
