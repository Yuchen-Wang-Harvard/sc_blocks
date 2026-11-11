"""One physical-mode action term in the two-Virasoro branching basis."""
struct BranchActionTerm{T<:AbstractFloat}
    label::Int
    first::Tuple{Vararg{Int16}}
    second::Tuple{Vararg{Int16}}
    coefficient::Complex{T}
end

"""Cached free-field actions in positive/reflected charts, including alternative Ward modes."""
mutable struct BranchActions{T<:AbstractFloat}
    positive::FreeField{T}
    reflected::Union{Nothing,FreeField{T}}
    actions::Dict{NTuple{3,Int},Vector{BranchActionTerm{T}}}
    evidence::Dict{NTuple{3,Int},NamedTuple}
    parity_reuse::Bool
end

"""Build a per-edge action plan; parity reuse is separately checked against direct alpha=1 fits."""
function BranchActions(::Type{T},sector::Symbol,b,p;parity=0,parity_reuse=true,plan=Combinatorics()) where T<:AbstractFloat
    BranchActions{T}(FreeField(T,sector,b,p;parity=parity,plan=plan),nothing,Dict(),Dict(),parity_reuse)
end

"""All pairs of Virasoro partitions of a fixed combined integer level."""
function partition_pairs(plan::Combinatorics,n::Int)
    [(a,b) for j in 0:n for a in partitions(plan,j) for b in partitions(plan,n-j)]
end

"""List the local branching span prescribed by the physical L_±1/L_0 support rules."""
function action_descriptors(module_::FreeField,n4::Int,mode::Int)
    n4>=0 || throw(ArgumentError("construct action spans in a positive chart"))
    mode in (-1,0,1) || throw(ArgumentError("branching action supports L_-1,L_0,L_1"))
    out=Tuple{Int,Tuple{Vararg{Int16}},Tuple{Vararg{Int16}}}[]
    if mode==1
        (module_.ramond ? n4==1 : n4<=2) && return out
        for (a,b) in partition_pairs(module_.physical.plan,n4-3)
            push!(out,(n4-4,a,b))
        end
    elseif mode==-1
        push!(out,(n4,(Int16(1),),()),(n4,(),(Int16(1),)))
        n4==0 && return out
        for (a,b) in partition_pairs(module_.physical.plan,n4-1)
            push!(out,(n4-4,a,b))
        end
    else
        push!(out,(n4,(),()))
        if module_.ramond && n4>1
            for (a,b) in partition_pairs(module_.physical.plan,n4-2)
                push!(out,(n4-4,a,b))
            end
        elseif !module_.ramond && n4>2
            for (a,b) in partition_pairs(module_.physical.plan,n4-2)
                push!(out,(n4-4,a,b))
            end
        elseif !module_.ramond && n4==2
            push!(out,(-2,(),()))
        end
    end
    out
end

"""Solve one physical-mode action from oscillator columns, checking every Fock row."""
function direct_branch_action(module_::FreeField{T},n4::Int,alpha::Int,mode::Int) where T
    descriptors=action_descriptors(module_,n4,mode)
    if isempty(descriptors)
        return BranchActionTerm{T}[],(rows=0,columns=0,iterations=0,residual=zero(T),
            smallest_machine_singular=1.0,precise_fallback=false)
    end
    primary=ff_primary(module_,n4,alpha)
    target=ff_physical_apply(module_,L(mode),primary)
    columns=[ff_descendant(module_,label,alpha,a,b) for (label,a,b) in descriptors]
    solver=SpanSolver(columns;extra=collect(keys(target)))
    coefficients,evidence=solve_span(solver,target)
    terms=[BranchActionTerm{T}(label,a,b,coefficients[j]) for (j,(label,a,b)) in enumerate(descriptors)]
    terms,evidence
end

"""Return cached current-convention physical L-mode branching coefficients."""
function branch_action(plan::BranchActions{T},n4::Int,alpha::Int,mode::Int) where T
    key=(n4,alpha,mode)
    haskey(plan.actions,key) && return plan.actions[key]
    if plan.positive.ramond && alpha==1 && plan.parity_reuse
        base=branch_action(plan,n4,0,mode)
        source=parity_sign(div(abs(n4)-1,2))
        out=[BranchActionTerm{T}(term.label,term.first,term.second,
              Complex{T}(2)^div(parity_sign(div(abs(term.label)-1,2))-source,2)*term.coefficient) for term in base]
        plan.actions[key]=out
        plan.evidence[key]=plan.evidence[(n4,0,mode)]
        return out
    end
    module_=plan.positive
    if n4<0
        if isnothing(plan.reflected)
            plan.reflected=FreeField(T,module_.ramond ? :R : :NS,module_.b,-module_.p;
                parity=Int(module_.physical.primary_parity),plan=module_.physical.plan)
        end
        module_=plan.reflected
    end
    terms,evidence=direct_branch_action(module_,abs(n4),alpha,mode)
    if n4<0
        terms=[BranchActionTerm{T}(-term.label,term.first,term.second,term.coefficient) for term in terms]
    end
    plan.actions[key],plan.evidence[key]=terms,evidence
    terms
end

"""Reconstruct an action in a requested chart for directed independent residual checks."""
function reconstruct_action(module_::FreeField{T},terms::Vector{BranchActionTerm{T}},alpha::Int) where T
    out=Dict{UInt32,Complex{T}}()
    for term in terms
        sparse_add!(out,ff_descendant(module_,term.label,alpha,term.first,term.second),term.coefficient)
    end
    prune_zeros!(out)
end

"""Closed one-leg Virasoro descendant form needed by local branching recurrences."""
function single_leg_vertex(slot::Int,word::Tuple,weights::NTuple{3,C}) where C<:Number
    slot in (1,2,3) || throw(ArgumentError("invalid ordered slot"))
    value=one(C)
    level=0
    h1,h2,h3=weights
    for n in reverse(word)
        factor=slot==1 ? h1+level+n*h2-h3 : slot==2 ?
            parity_sign(n)*(h2+level+n*h3-h1) : h3+level+n*h2-h1
        value*=factor
        level+=n
    end
    value
end
