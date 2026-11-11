"""Precision-dependent accuracy target for verified local span solves."""
span_tolerance(::Type{T}) where T<:AbstractFloat=max(T(256)*eps(T),eps(T)^(T(7)/8))

"""Reusable normalized sparse span, with a machine LU preconditioner and exact-precision fallback."""
mutable struct SpanSolver{T<:AbstractFloat,K}
    keys::Vector{K}
    index::Dict{K,Int}
    matrix::Matrix{Complex{T}}
    norms::Vector{T}
    selected::Vector{Int}
    machine_factor::Union{Nothing,LU{ComplexF64,Matrix{ComplexF64},Vector{Int}}}
    precise_factor::Union{Nothing,LU{Complex{T},Matrix{Complex{T}},Vector{Int}}}
    smallest_machine_singular::Float64
end

"""Construct a span once; use machine arithmetic only to select rows/precondition.

All matrix entries, residuals and final coefficients remain at working precision.
A numerically deficient machine restriction falls back to working-precision LU
row selection. Every solve checks all rows, including those not selected.
"""
function SpanSolver(columns::Vector{Dict{K,Complex{T}}};extra::Vector{K}=K[]) where {K,T<:AbstractFloat}
    ids=Set(extra)
    for column in columns
        union!(ids,keys(column))
    end
    rows=sort!(collect(ids))
    m,n=length(rows),length(columns)
    m>=n>0 || throw(ArgumentError("invalid span dimensions: $m rows, $n columns"))
    index=Dict(key=>i for (i,key) in enumerate(rows))
    matrix=zeros(Complex{T},m,n)
    norms=T[sqrt(sum(abs2,values(column);init=zero(T))) for column in columns]
    any(iszero,norms) && throw(ArgumentError("zero span column"))
    for (j,column) in enumerate(columns),(key,value) in column
        matrix[index[key],j]=value/norms[j]
    end
    shadow=ComplexF64.(matrix)
    pivoted=qr(Matrix(transpose(shadow)),ColumnNorm())
    selected=Vector{Int}(pivoted.p[1:n])
    restriction=shadow[selected,:]
    smallest=minimum(svdvals(restriction))
    machine_factor=smallest>1e-11 ? lu(restriction) : nothing
    precise_factor=nothing
    if isnothing(machine_factor)
        precise_pivots=lu(matrix;check=false)
        selected=Vector{Int}(precise_pivots.p[1:n])
        precise_factor=lu(matrix[selected,:])
    end
    SpanSolver{T,K}(rows,index,matrix,norms,selected,machine_factor,precise_factor,smallest)
end

"""Solve a sparse span and return verified coefficients plus numerical evidence."""
function solve_span(solver::SpanSolver{T,K},target::Dict{K,Complex{T}}) where {T,K}
    C=Complex{T}
    rhs=zeros(C,length(solver.keys))
    for (key,value) in target
        haskey(solver.index,key) || iszero(value) || throw(ArgumentError("target row absent from span"))
        haskey(solver.index,key) && (rhs[solver.index[key]]=value)
    end
    n=length(solver.norms)
    selected_rhs=rhs[solver.selected]
    restriction=solver.matrix[solver.selected,:]
    target_norm=sqrt(sum(abs2,rhs))
    tolerance=span_tolerance(T)
    values=zeros(C,n)
    iterations=0
    used_precise=isnothing(solver.machine_factor)
    if !used_precise
        for iteration in 1:32
            residual=selected_rhs-restriction*values
            scale=maximum(abs,residual;init=zero(T))
            scale<=tolerance*max(target_norm,eps(T))/64 && break
            correction=solver.machine_factor\ComplexF64.(residual/scale)
            values .+= scale.*C.(correction)
            iterations+=1
        end
    else
        values=solver.precise_factor\selected_rhs
    end
    error=sqrt(sum(abs2,solver.matrix*values-rhs))/max(target_norm,eps(T))
    if error>tolerance && !used_precise
        isnothing(solver.precise_factor) && (solver.precise_factor=lu(restriction))
        values=solver.precise_factor\selected_rhs
        error=sqrt(sum(abs2,solver.matrix*values-rhs))/max(target_norm,eps(T))
        used_precise=true
    end
    isfinite(error) && error<=64tolerance ||
        throw(ErrorException("full-row span residual $error exceeds tolerance $(64tolerance); increase precision or check the proposed span"))
    coefficients=values./solver.norms
    evidence=(rows=size(solver.matrix,1),columns=n,iterations=iterations,residual=error,
              smallest_machine_singular=solver.smallest_machine_singular,precise_fallback=used_precise)
    coefficients,evidence
end
