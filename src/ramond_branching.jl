"""One local physical Virasoro Ward row, with cancellation bounds for its coefficients."""
struct BranchWardRow{T<:AbstractFloat}
    coefficients::Dict{NTuple{3,Int},Complex{T}}
    bounds::Dict{NTuple{3,Int},T}
end

"""NSRR branching recurrence: direct boundary data, cached actions and small local systems."""
mutable struct NSRRBranching{T<:AbstractFloat}
    actions::NTuple{3,BranchActions{T}}
    physical::WardVertex{T,2}
    auxiliary::FermionVertex{T}
    ordered::NTuple{3,Dict{Tuple{Int,Int},Dict{Tuple{FermionState,PBWState},Complex{T}}}}
    tables::Dict{NTuple{2,Int},Dict{NTuple{3,Int},NTuple{2,Complex{T}}}}
    rows::Dict{NTuple{6,Int},BranchWardRow{T}}
    factors::Dict{Tuple{NTuple{3,Int},Int,Int,Tuple{Vararg{Int16}}},Tuple{Complex{T},T}}
    weights::Dict{NTuple{3,Int},Complex{T}}
    direct_boundaries::Int
    recursive_values::Int
    recursive_groups::Int
    second_ward_rows::Int
    maximum_local_unknowns::Int
    maximum_ward_residual::T
end

"""Construct a local NSRR plan, sharing existing per-edge action caches when supplied."""
function NSRRBranching(actions::NTuple{3,BranchActions{T}}) where T
    modules=map(a->a.positive,actions)
    map(m->m.ramond,modules)==(false,true,true) || throw(ArgumentError("branching slots must be NS,R,R"))
    all(m->m.b==modules[1].b,modules) || throw(ArgumentError("couplings differ"))
    C=Complex{T}
    ordered=ntuple(_->Dict{Tuple{Int,Int},Dict{Tuple{FermionState,PBWState},C}}(),3)
    NSRRBranching{T}(actions,nsrr_vertex(map(m->m.physical,modules)),
        FermionVertex(map(m->m.auxiliary,modules)),ordered,
        Dict((a,g)=>Dict{NTuple{3,Int},NTuple{2,C}}() for a in 0:1 for g in 0:1),
        Dict(),Dict(),Dict(),0,0,0,0,0,zero(T))
end
function NSRRBranching(::Type{T},b,momenta::NTuple{3};parity::Int=0) where T<:AbstractFloat
    NSRRBranching(ntuple(j->BranchActions(T,j==1 ? :NS : :R,b,momenta[j];parity=j==1 ? parity : 0),3))
end

"""Height decreases under every bulk dependency; boundary couplings stay within one group."""
branch_height(labels::NTuple{3,Int})=div(abs(labels[1]),4)+div(abs(labels[2])-1,4)+div(abs(labels[3])-1,4)
branch_boundary(labels::NTuple{3,Int})=abs(labels[1])<=2 && abs(labels[2])<=3 && abs(labels[3])<=3

"""Group the explicit NS ±1/2 and Ramond 1/4↔-3/4 boundary couplings."""
function branch_group(labels::NTuple{3,Int})
    (abs(labels[1])==2 ? 2 : labels[1],
     abs(labels[2])<=3 ? (labels[2] in (1,-3) ? 1 : -1) : labels[2],
     abs(labels[3])<=3 ? (labels[3] in (1,-3) ? 1 : -1) : labels[3])
end

"""All labels in one boundary-coupled local group, independently of requested output support."""
function branch_group_members(labels::NTuple{3,Int})
    key=branch_group(labels)
    choices=ntuple(j->j==1 ? (abs(key[j])==2 ? (-2,2) : (key[j],)) :
        abs(key[j])<=3 ? (key[j],-3key[j]) : (key[j],),3)
    collect(Iterators.product(choices...))[:]
end

"""Cached ordered PBW expression of a local chi vector, used only for boundaries and checks."""
function branching_ordered(plan::NSRRBranching,slot::Int,label::Int,alpha::Int)
    get!(plan.ordered[slot],(label,alpha)) do
        module_=plan.actions[slot].positive
        ordered_pbw(module_,ff_primary(module_,label,alpha))
    end
end

"""Direct local NSRR chi-vector value from the current physical/auxiliary Ward forms.

Production calls this only in the fixed boundary box. High-label calls are
reserved for independent recurrence checks. Both eta signs share every term.
"""
function direct_branch_values(plan::NSRRBranching{T},labels::NTuple{3,Int},alpha::Int,gamma::Int) where T
    vectors=(branching_ordered(plan,1,labels[1],0),branching_ordered(plan,2,labels[2],alpha),
             branching_ordered(plan,3,labels[3],gamma))
    out=(zero(Complex{T}),zero(Complex{T}))
    for ((a1,x1),v1) in vectors[1],((a2,x2),v2) in vectors[2],((a3,x3),v3) in vectors[3]
        auxiliary=rho(plan.auxiliary,(a1,a2,a3))
        iszero(auxiliary) && continue
        exponent=Int(relative_parity(a1))*Int(state_parity(plan.physical.modules[1],x1))+
                 Int(relative_parity(a3))*Int(state_parity(plan.physical.modules[2],x2))
        coefficient=v1*v2*v3*auxiliary*parity_sign(exponent)
        out=form_add(out,form_scale(ward_values(plan.physical,(x1,x2,x3)),coefficient))
    end
    out
end

"""Read one established raw value, evaluating only genuine boundary data on demand."""
function known_branch_value(plan::NSRRBranching,labels::NTuple{3,Int},alpha::Int,gamma::Int)
    table=plan.tables[(alpha,gamma)]
    haskey(table,labels) && return table[labels]
    branch_boundary(labels) || throw(ArgumentError("missing lower branching coefficient $labels"))
    value=direct_branch_values(plan,labels,alpha,gamma)
    table[labels]=value
    plan.direct_boundaries+=1
    value
end

"""Cached embedded weight per slot, label and Virasoro copy."""
function branching_weight(plan::NSRRBranching,slot::Int,label::Int,copy::Int)
    get!(plan.weights,(slot,label,copy)) do
        module_=plan.actions[slot].positive
        branch_weight(module_.b,module_.p,label,copy)
    end
end

"""Absolute product bound for a one-leg Virasoro factor, before cancellations."""
function single_leg_bound(slot::Int,word::Tuple,weights)
    value=one(real(first(weights)))
    level=0
    for n in reverse(word)
        value*=abs(weights[slot])+level+n*abs(weights[slot==2 ? 3 : 2])+abs(weights[slot==1 ? 3 : 1])
        level+=n
    end
    value
end

"""Reduce an action term to a raw branching coefficient times cached descendant factors."""
function branch_term_factor(plan::NSRRBranching{T},labels::NTuple{3,Int},slot::Int,term::BranchActionTerm{T}) where T
    changed=Base.setindex(labels,term.label,slot)
    value,bound=term.coefficient,abs(term.coefficient)
    for copy in 1:2
        word=copy==1 ? term.first : term.second
        isempty(word) && continue
        factor,estimate=get!(plan.factors,(changed,slot,copy,word)) do
            weights=ntuple(j->branching_weight(plan,j,changed[j],copy),3)
            single_leg_vertex(slot,word,weights),single_leg_bound(slot,word,weights)
        end
        value*=factor
        bound*=estimate
    end
    changed,value,bound
end

"""Build either physical L_1 Ward identity; the second supplies an alternative pivot."""
function branch_ward_row(plan::NSRRBranching{T},labels::NTuple{3,Int},alpha::Int,gamma::Int;second::Bool=false) where T
    get!(plan.rows,(labels...,alpha,gamma,Int(second))) do
        coefficients=Dict{NTuple{3,Int},Complex{T}}()
        bounds=Dict{NTuple{3,Int},T}()
        function add_terms(slot,mode,scale)
            parity=slot==1 ? 0 : slot==2 ? alpha : gamma
            for term in branch_action(plan.actions[slot],labels[slot],parity,mode)
                changed,value,bound=branch_term_factor(plan,labels,slot,term)
                coefficients[changed]=get(coefficients,changed,zero(Complex{T}))+scale*value
                bounds[changed]=get(bounds,changed,zero(T))+abs(scale)*bound
            end
        end
        if !second
            add_terms(1,1,1)
            add_terms(2,-1,-1)
            add_terms(3,-1,-1)
        else
            add_terms(1,-1,1)
            add_terms(2,-1,-1)
            add_terms(2,0,-2)
            add_terms(2,1,-1)
            add_terms(3,1,-1)
        end
        BranchWardRow{T}(coefficients,bounds)
    end
end

"""Relative pivot threshold uses working precision and local cancellation bounds."""
ward_pivot_tolerance(::Type{T}) where T<:AbstractFloat=max(T(256)*eps(T),eps(T)^(T(2)/3))

"""Solve at most four local unknowns and both eta right sides with full pivoting."""
function solve_branch_group(plan::NSRRBranching{T},nodes::Vector{NTuple{3,Int}},rows::Vector{BranchWardRow{T}},
                            alpha::Int,gamma::Int) where T
    C=Complex{T}
    n,m=length(nodes),length(rows)
    index=Dict(node=>j for (j,node) in enumerate(nodes))
    matrix=zeros(C,m,n+2)
    tolerance=ward_pivot_tolerance(T)
    for (i,row) in enumerate(rows)
        for (changed,value) in row.coefficients
            if haskey(index,changed)
                abs(value)>tolerance*row.bounds[changed] && (matrix[i,index[changed]]+=value)
            else
                branch_height(changed)<branch_height(first(nodes)) ||
                    throw(ArgumentError("nondecreasing branching dependency $(first(nodes)) -> $changed"))
                known=known_branch_value(plan,changed,alpha,gamma)
                matrix[i,n+1]-=value*known[1]
                matrix[i,n+2]-=value*known[2]
            end
        end
        scale=sum(abs,view(matrix,i,1:n))
        iszero(scale) || (matrix[i,:]./=scale)
    end
    columns=collect(1:n)
    for k in 1:n
        block=view(matrix,k:m,k:n)
        pivot=argmax(abs.(block))
        pi,pj=k+pivot[1]-1,k+pivot[2]-1
        abs(matrix[pi,pj])>tolerance || return nothing
        matrix[k,:],matrix[pi,:]=copy(matrix[pi,:]),copy(matrix[k,:])
        matrix[:,k],matrix[:,pj]=copy(matrix[:,pj]),copy(matrix[:,k])
        columns[k],columns[pj]=columns[pj],columns[k]
        matrix[k,k:n+2]./=matrix[k,k]
        for i in k+1:m
            factor=matrix[i,k]
            matrix[i,k:n+2].-=factor.*matrix[k,k:n+2]
        end
    end
    values=zeros(C,n,2)
    for eta in 1:2,i in n:-1:1
        values[i,eta]=matrix[i,n+eta]-sum((matrix[i,j]*values[j,eta] for j in i+1:n);init=zero(C))
    end
    out=Dict{NTuple{3,Int},NTuple{2,C}}()
    for j in 1:n
        out[nodes[columns[j]]]=(values[j,1],values[j,2])
    end
    out
end

"""Scaled residual of a physical Ward row after all local values are established."""
function branch_row_residual(plan::NSRRBranching{T},row::BranchWardRow{T},alpha::Int,gamma::Int) where T
    residual=(zero(Complex{T}),zero(Complex{T}))
    scales=(one(T),one(T))
    for (changed,coefficient) in row.coefficients
        term=form_scale(known_branch_value(plan,changed,alpha,gamma),coefficient)
        residual=form_add(residual,term)
        scales=map(+,scales,map(abs,term))
    end
    maximum(map((r,s)->abs(r)/s,residual,scales))
end

"""Close requested labels under first-Ward dependencies and complete boundary groups."""
function branch_support(plan::NSRRBranching,requests::AbstractVector{NTuple{3,Int}})
    support=Set(requests)
    pending=collect(requests)
    while !isempty(pending)
        labels=pop!(pending)
        iseven(labels[1]) && isodd(labels[2]) && isodd(labels[3]) || throw(ArgumentError("invalid NSRR labels"))
        dependencies=branch_group_members(labels)
        for slot in 1:3,term in branch_action(plan.actions[slot],labels[slot],0,slot==1 ? 1 : -1)
            push!(dependencies,Base.setindex(labels,term.label,slot))
        end
        for changed in dependencies
            if changed ∉ support
                push!(support,changed)
                push!(pending,changed)
            end
        end
    end
    support
end

"""Prepare requested NSRR branching coefficients by local L_1 recurrence, sharing both eta signs."""
function prepare_branching!(plan::NSRRBranching{T},requests::AbstractVector{NTuple{3,Int}}) where T
    groups=Dict{NTuple{3,Int},Vector{NTuple{3,Int}}}()
    for labels in branch_support(plan,requests)
        push!(get!(() -> NTuple{3,Int}[],groups,branch_group(labels)),labels)
    end
    ordered=sort!(collect(keys(groups));by=k->(branch_height(k),k))
    for key in ordered
        nodes=sort!(groups[key])
        for alpha in 0:1,gamma in 0:1
            table=plan.tables[(alpha,gamma)]
            all(node->haskey(table,node),nodes) && continue
            if branch_boundary(first(nodes))
                for node in nodes
                    known_branch_value(plan,node,alpha,gamma)
                end
                continue
            end
            length(nodes)<=4 || throw(ArgumentError("unexpected branching group of $(length(nodes)) unknowns"))
            plan.maximum_local_unknowns=max(plan.maximum_local_unknowns,length(nodes))
            rows=[branch_ward_row(plan,node,alpha,gamma) for node in nodes]
            values=solve_branch_group(plan,nodes,rows,alpha,gamma)
            if isnothing(values)
                append!(rows,[branch_ward_row(plan,node,alpha,gamma;second=true) for node in nodes])
                plan.second_ward_rows+=length(nodes)
                values=solve_branch_group(plan,nodes,rows,alpha,gamma)
            end
            isnothing(values) && throw(ErrorException("both branching Ward pivots failed at $nodes, alpha=$alpha, gamma=$gamma"))
            merge!(table,values)
            plan.recursive_values+=length(nodes)
            plan.recursive_groups+=1
            for row in rows
                residual=branch_row_residual(plan,row,alpha,gamma)
                plan.maximum_ward_residual=max(plan.maximum_ward_residual,residual)
                residual<=64ward_pivot_tolerance(T) || throw(ErrorException("branching Ward residual $residual at $nodes"))
            end
        end
    end
    plan
end

"""Read an unnormalized NSRR branching value, preserving the manuscript's relative form label."""
function branching_raw(plan::NSRRBranching{T},labels::NTuple{3,Int},alpha::Int,gamma::Int,f::Int,eta::Int) where T
    f in (0,1) && eta in (-1,1) && alpha in (0,1) && gamma in (0,1) || throw(ArgumentError("invalid branching label"))
    mod(div(labels[1],2)+alpha+gamma,2)==f || return zero(Complex{T})
    haskey(plan.tables[(alpha,gamma)],labels) || prepare_branching!(plan,[labels])
    plan.tables[(alpha,gamma)][labels][eta==1 ? 1 : 2]
end
