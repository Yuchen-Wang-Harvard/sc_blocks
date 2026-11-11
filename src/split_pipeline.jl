"""Expanded mixed-sector graph with one Theta psi vertex per odd Ramond cycle."""
mutable struct SplitBlockPlan{T<:AbstractFloat,N,M}
    original::RamondBlockPlan{T,N}
    graph::PlumbingGraph
    splits::Dict{Int,NTuple{2,Int}}
    owners::NTuple{M,Int}
    splittings::Dict{Int,SplittingPlan{T}}
    auxiliary_insertions::Dict{NTuple{2,Int},Array{Complex{T},3}}
    domains::Dict{NTuple{M,Int},CoefficientDomain{M}}
    unequal_terms::Int
    transpose_reuse::Int
end

function SplitBlockPlan(original::RamondBlockPlan{T,N}) where {T,N}
    isempty(original.marked) && throw(ArgumentError("split plan needs an odd Ramond cycle"))
    graph,splits=split_graph(original.graph,original.marked)
    M=length(graph.ends)
    owners=collect(1:N)
    append!(owners,original.marked)
    SplitBlockPlan{T,N,M}(original,graph,splits,Tuple(owners),
        Dict(e=>SplittingPlan(original.actions[e]) for e in original.marked),Dict(),Dict(),0,0)
end

"""Downward closure of equal-power output targets, retaining every unequal intermediate grade."""
function split_domain(domain::CoefficientDomain{N},splits::Dict{Int,NTuple{2,Int}}) where N
    M=N+length(splits)
    sort([right for (_,right) in values(splits)])==collect(N+1:M) ||
        throw(ArgumentError("split edge IDs must append to original edges"))
    indices=NTuple{M,Int}[]
    marked=sort(collect(keys(splits)))
    work=zeros(Int,M)
    function visit(key,position)
        if position>length(marked)
            push!(indices,Tuple(work))
            return
        end
        edge=marked[position]
        left,right=splits[edge]
        left==edge || throw(ArgumentError("left split edge must retain original ID"))
        highest=key[edge]
        # max(left,right)=highest makes the enumeration disjoint.
        for lower in 0:highest
            work[left],work[right]=highest,lower
            visit(key,position+1)
            if lower<highest
                work[left],work[right]=lower,highest
                visit(key,position+1)
            end
        end
    end
    for key in domain.indices
        work[1:N].=key
        visit(key,1)
    end
    CoefficientDomain(indices)
end

"""Copy each original edge parity onto both halves of its split edge."""
function split_mask(plan::SplitBlockPlan,mask::Int)
    result=mask
    for (edge,(_,right)) in plan.splits
        result|=((mask>>(edge-1))&1)<<(right-1)
    end
    result
end

"""Koszul sign with the inserted bra/ket evaluated as an operator kernel.

An insertion is a matrix element, not a third ordered trilinear form. Converting
its dual input to the (right,external,left) tensor axes reverses the two equal-
parity internal legs once. The extra (-1)^parity is forced by the identity
kernel: two inverse BPZ pairings with an identity between them must reproduce
one original propagator, without inserting a parity operator. This is applied
to both auxiliary and enlarged sewing, before coefficient recovery.
"""
function split_sewing_sign(plan::SplitBlockPlan,mask::Int)
    sewing_sign(plan.graph,mask)*parity_sign(sum(Int((mask>>(edge-1))&1) for edge in plan.original.marked))
end

"""Extract equal split powers/parities after complete expanded-graph assembly."""
function split_diagonal(plan::SplitBlockPlan{T,N,M},expanded::ParitySeries{M,T}) where {T,N,M}
    out=ParitySeries{N,T}()
    for (degree,coefficients) in expanded
        all(pair->degree[pair[1]]==degree[pair[2]],values(plan.splits)) || continue
        out[ntuple(e->degree[e],N)]=[coefficients[split_mask(plan,mask)+1] for mask in 0:2^N-1]
    end
    out
end

"""Branch labels on the split graph, imposing the exact primary splitting support."""
function split_branch_labels(plan::SplitBlockPlan{T,N,M},domain::CoefficientDomain{M}) where {T,N,M}
    base=plan.original
    choices=ntuple(e->begin
        original=plan.owners[e]
        if base.ramond[original]
            bound=isqrt(4domain.limits[e]+1)
            iseven(bound) && (bound-=1)
            collect(-bound:2:bound)
        else
            bound=2isqrt(domain.limits[e])
            collect(-bound:2:bound)
        end
    end,M)
    labels_=NTuple{M,Int}[]
    for labels in Iterators.product(choices...)
        all(pair->abs(labels[pair[1]]-labels[pair[2]])==2,values(plan.splits)) || continue
        shift=ntuple(e->branch_level2(base.ramond[plan.owners[e]],labels[e]),M)
        contains(domain,shift) && push!(labels_,labels)
    end
    labels_
end

"""Expanded raw vertex product, insertion matrix elements and full inverse branch norms."""
function split_branch_coefficients(plan::SplitBlockPlan{T,N,M},labels::NTuple{M,Int}) where {T,N,M}
    base=plan.original
    out=zeros(Complex{T},2^M)
    for original_mask in 0:2^N-1
        mask=split_mask(plan,original_mask)
        all(e->base.ramond[e] || ((mask>>(e-1))&1)==mod(div(labels[e],2)+base.primary_parities[e],2),1:N) || continue
        coefficient=one(Complex{T})
        for v in eachindex(base.graph.slots)
            slots=plan.graph.slots[v]
            alpha,gamma=Int((mask>>(slots[2]-1))&1),Int((mask>>(slots[3]-1))&1)
            coefficient*=graph_branch_value(base,v,ntuple(j->labels[slots[j]],3),(0,alpha,gamma))
            iszero(coefficient) && break
        end
        iszero(coefficient) && continue
        for (edge,(left,right)) in plan.splits
            coefficient*=splitting_raw(plan.splittings[edge],labels[right],labels[left],Int((mask>>(edge-1))&1))
        end
        for e in 1:M
            original=plan.owners[e]
            alpha=base.ramond[original] ? Int((mask>>(e-1))&1) : 0
            coefficient/=ramond_plan_norm(base,original,labels[e],alpha)
        end
        out[mask+1]=split_sewing_sign(plan,mask)*coefficient
    end
    out
end

"""Auxiliary insertion tensor <right|Theta psi(1)|left> at independent split grades."""
function auxiliary_insertion_tensor(plan::SplitBlockPlan{T},rightlevel::Int,leftlevel::Int) where T
    get!(plan.auxiliary_insertions,(rightlevel,leftlevel)) do
        module_=plan.original.auxiliary.modules[2]
        right,left=basis(module_,rightlevel),basis(module_,leftlevel)
        out=zeros(Complex{T},length(right),1,length(left))
        for (j,ket) in enumerate(left),(target,value) in theta_psi_action(module_,leftlevel-rightlevel,ket)
            i=findfirst(==(target),right)
            isnothing(i) || (out[i,1,j]+=fermion_norm(module_,target)*value)
        end
        out
    end
end

"""Direct auxiliary sewing on the split graph; the domain may contain unequal grades."""
function split_auxiliary_block(plan::SplitBlockPlan{T,N,M},domain::CoefficientDomain{M};diagonal_only::Bool=false) where {T,N,M}
    base=plan.original
    out=ParitySeries{M,T}()
    modules=ntuple(e->base.auxiliary.modules[base.ramond[plan.owners[e]] ? 2 : 1],M)
    for level in domain.indices
        diagonal_only && any(pair->level[pair[1]]!=level[pair[2]],values(plan.splits)) && continue
        any(e->modules[e].ramond && isodd(level[e]),1:M) && continue
        bases=ntuple(e->basis(modules[e],level[e]),M)
        any(isempty,bases) && continue
        tensors=[ramond_auxiliary_tensor(base,ntuple(j->level[slots[j]],3);ramond=base.ramond[plan.owners[slots[2]]])
            for slots in plan.graph.slots[1:length(base.graph.slots)]]
        for edge in base.marked
            left,right=plan.splits[edge]
            push!(tensors,auxiliary_insertion_tensor(plan,level[right],level[left]))
        end
        inverses=[inverse_gram(modules[e],level[e]) for e in 1:M]
        selections=[ntuple(p->findall(s->Int(relative_parity(s))==p-1,bases[e]),2) for e in 1:M]
        coefficients=zeros(Complex{T},2^M)
        for original_mask in 0:2^N-1
            mask=split_mask(plan,original_mask)
            all(slots->iseven(sum((mask>>(e-1))&1 for e in slots)),plan.graph.slots[1:length(base.graph.slots)]) || continue
            selected=[selections[e][((mask>>(e-1))&1)+1] for e in 1:M]
            coefficients[mask+1]=split_sewing_sign(plan,mask)*sew_tensors(plan.graph,tensors,inverses,selected)
        end
        any(!iszero,coefficients) && (out[level]=coefficients)
    end
    out
end

"""Compute only final diagonal product targets, using every unequal intermediate power.

The two CCY factors already contain the full downward closure on independent
split edges. This restricts output coefficients of their product, never the
intermediate degrees being convolved. Returned keys include the branch shift.
"""
function split_diagonal_product(plan::SplitBlockPlan{T,N,M},first::Series{M,Complex{T}},
                                second::Series{M,Complex{T}},shift::NTuple{M,Int},
                                domain::CoefficientDomain{N}) where {T,N,M}
    minimal=ntuple(e->maximum(shift[j] for j in 1:M if plan.owners[j]==e),N)
    out=Series{N,Complex{T}}()
    contains(domain,minimal) || return out
    remaining=residual_domain(domain,minimal;scale=2)
    z=zero(Complex{T})
    work=complex_scratch(Complex{T})
    for key in remaining.indices
        target=ntuple(e->minimal[e]+2key[e],N)
        all(j->iseven(target[plan.owners[j]]-shift[j]),1:M) || continue
        degree=ntuple(j->div(target[plan.owners[j]]-shift[j],2),M)
        value=zero(Complex{T})
        for left in Iterators.product(ntuple(j->0:degree[j],M)...)
            right=ntuple(j->degree[j]-left[j],M)
            a,b=get(first,left,z),get(second,right,z)
            (iszero(a) || iszero(b)) && continue
            value=sum_product(value,a,b,work)
            any(pair->left[pair[1]]!=left[pair[2]],values(plan.splits)) && (plan.unequal_terms+=1)
        end
        iszero(value) || (out[target]=value)
    end
    out
end

"""Theta symmetry under exchanging the two halves of the inserted Ramond edge."""
function split_transpose_supported(plan::SplitBlockPlan)
    plan.original.graph.slots==[(1,2,3),(1,2,3)] && plan.original.marked==[2]
end

"""Assemble diagonal enlarged coefficients with local and transposed-factor reuse."""
function split_enlarged_diagonal(plan::SplitBlockPlan{T,N,M},domain::CoefficientDomain{N},
                                 expanded_domain::CoefficientDomain{M};reuse_transpose::Bool=true) where {T,N,M}
    base=plan.original
    prepare_product_cache!(base,domain)
    started=time_ns()
    labels_=split_branch_labels(plan,expanded_domain)
    empty!(plan.domains)
    requests=IdDict{NSRRBranching{T},Vector{NTuple{3,Int}}}()
    for v in eachindex(base.graph.slots)
        isnothing(base.vertices[v]) && continue
        slots=plan.graph.slots[v]
        append!(get!(() -> NTuple{3,Int}[],requests,base.vertices[v]),
            [ntuple(j->labels[slots[j]],3) for labels in labels_])
    end
    for (vertex,labels) in requests
        prepare_branching!(vertex,unique!(labels))
    end
    record_stage!(base,"branching",started)
    numerator=ParitySeries{N,T}()
    products=Dict{NTuple{M,Int},Series{N,Complex{T}}}()
    transpose=reuse_transpose && split_transpose_supported(plan)
    external_weights=psi_weights(base.b)
    for labels in labels_
        started=time_ns()
        expanded_coefficients=split_branch_coefficients(plan,labels)
        coefficients=[expanded_coefficients[split_mask(plan,mask)+1] for mask in 0:2^N-1]
        record_stage!(base,"branch_coefficients",started)
        any(!iszero,coefficients) || continue
        shift=ntuple(e->branch_level2(base.ramond[plan.owners[e]],labels[e]),M)
        key=transpose ? ntuple(e->e==2 ? min(labels[2],labels[M]) : e==M ? max(labels[2],labels[M]) : labels[e],M) : labels
        cachekey=(Tuple(base.marked),key)
        stored=base.cache_products ? get(base.products,cachekey,nothing) : transpose ? pop!(products,key,nothing) : nothing
        if isnothing(stored)
            started=time_ns()
            remaining=get!(plan.domains,shift) do
                residual_domain(expanded_domain,shift;scale=2)
            end
            copies=ntuple(copy->CCYPlan(T,plan.graph,branch_central(base.b,copy),
                [ramond_plan_weight(base,plan.owners[e],labels[e],copy) for e in 1:M];
                external=Dict(3(v-1)+2=>external_weights[copy] for v in length(base.graph.slots)+1:length(plan.graph.slots))),2)
            first,second=ccy_reduced(copies[1],remaining),ccy_reduced(copies[2],remaining)
            record_stage!(base,"ccy",started)
            record_ccy_stages!(base,copies)
            base.ccy_transitions+=copies[1].transitions+copies[2].transitions
            started=time_ns()
            product=split_diagonal_product(plan,first,second,shift,domain)
            if base.cache_products
                base.products[cachekey]=product
            elseif transpose
                products[key]=product
            end
            record_stage!(base,"products",started)
        else
            product=stored
            base.cache_products && (base.product_cache_hits+=1)
            transpose && (plan.transpose_reuse+=1)
        end
        started=time_ns()
        for (target,value) in product
            destination=get!(() -> zeros(Complex{T},2^N),numerator,target)
            for i in eachindex(coefficients)
                iszero(coefficients[i]) || (destination[i]+=value*coefficients[i])
            end
        end
        record_stage!(base,"assembly",started)
        base.branches+=1
    end
    numerator
end

"""Two CCY factors on independent split edges, retaining unequal grades until their product is formed."""
function split_enlarged_reduced(plan::SplitBlockPlan{T,N,M},domain::CoefficientDomain{M}) where {T,N,M}
    base=plan.original
    started=time_ns()
    labels_=split_branch_labels(plan,domain)
    empty!(plan.domains)
    for v in eachindex(base.graph.slots)
        isnothing(base.vertices[v]) && continue
        slots=plan.graph.slots[v]
        requests=unique([ntuple(j->labels[slots[j]],3) for labels in labels_])
        prepare_branching!(base.vertices[v],requests)
    end
    record_stage!(base,"branching",started)
    numerator=ParitySeries{M,T}()
    external_weights=psi_weights(base.b)
    for labels in labels_
        started=time_ns()
        coefficients=split_branch_coefficients(plan,labels)
        record_stage!(base,"branch_coefficients",started)
        any(!iszero,coefficients) || continue
        started=time_ns()
        shift=ntuple(e->branch_level2(base.ramond[plan.owners[e]],labels[e]),M)
        remaining=get!(plan.domains,shift) do
            residual_domain(domain,shift;scale=2)
        end
        copies=ntuple(copy->CCYPlan(T,plan.graph,branch_central(base.b,copy),
            [ramond_plan_weight(base,plan.owners[e],labels[e],copy) for e in 1:M];
            external=Dict(3(v-1)+2=>external_weights[copy] for v in length(base.graph.slots)+1:length(plan.graph.slots))),2)
        first,second=ccy_reduced(copies[1],remaining),ccy_reduced(copies[2],remaining)
        record_stage!(base,"ccy",started)
        record_ccy_stages!(base,copies)
        started=time_ns()
        product=series_multiply(first,second,remaining)
        record_stage!(base,"products",started)
        started=time_ns()
        for (degree,value) in product
            target=ntuple(e->shift[e]+2degree[e],M)
            any(pair->target[pair[1]]!=target[pair[2]],values(plan.splits)) && (plan.unequal_terms+=1)
            destination=get!(() -> zeros(Complex{T},2^M),numerator,target)
            destination .+= value.*coefficients
        end
        record_stage!(base,"assembly",started)
        base.branches+=1
        base.ccy_transitions+=copies[1].transitions+copies[2].transitions
    end
    numerator
end

"""Inserted physical block with equal-power extraction followed by allowed-sector recovery.

The insertion sphere's two sewing maps compose to the original map with
q=q_left*q_right. The Schottky vacuum therefore depends only on that product;
its square is restored once after diagonal extraction and auxiliary division.
"""
function inserted_ramond_block(plan::SplitBlockPlan{T,N,M},domain::CoefficientDomain{N}) where {T,N,M}
    empty!(plan.original.timings)
    started=time_ns()
    expanded_domain=split_domain(domain,plan.splits)
    record_stage!(plan.original,"split_domain",started)
    numerator=split_enlarged_diagonal(plan,domain,expanded_domain)
    started=time_ns()
    auxiliary=split_diagonal(plan,split_auxiliary_block(plan,expanded_domain;diagonal_only=true))
    record_stage!(plan.original,"auxiliary",started)
    started=time_ns()
    reduced,residual=recover_allowed_sector(numerator,auxiliary,domain;
        tolerance=plan.original.recovery_tolerance,policy=plan.original.recovery_policy,
        diagnostics=plan.original.recovery_diagnostics)
    record_stage!(plan.original,"recovery_and_forward_check",started)
    plan.original.forward_residual=max(plan.original.forward_residual,residual)
    started=time_ns()
    integer_domain=residual_domain(domain,ntuple(_->0,N);scale=2)
    vacuum=schottky_vacuum(plan.original.graph,integer_domain)
    squared=series_multiply(vacuum,vacuum,integer_domain)
    lifted=Series{N,Complex{T}}(map(n->2n,key)=>Complex{T}(value) for (key,value) in squared)
    result=parity_scalar_product(reduced,lifted,domain)
    record_stage!(plan.original,"vacuum_and_restoration",started)
    result
end
