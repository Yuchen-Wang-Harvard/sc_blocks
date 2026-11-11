"""Momentum-independent caches. Populate per worker or before concurrent reads."""
mutable struct Combinatorics
    partitions::Dict{Tuple{Int,Int},Vector{Tuple{Vararg{Int16}}}}
    fermions::Dict{Tuple{Int,Bool,Int},Vector{Tuple{Vararg{Int16}}}}
    binomials::Dict{Tuple{Rational{BigInt},Int},Rational{BigInt}}
end
Combinatorics() = Combinatorics(Dict(), Dict(), Dict())

"""Descending integer partitions; cached independently of all module weights."""
function partitions(plan::Combinatorics, n::Integer, maximum::Integer=n)
    key = (Int(n), Int(min(n, maximum)))
    get!(plan.partitions, key) do
        out = Tuple{Vararg{Int16}}[]
        n == 0 && return push!(out, ())
        n < 0 && return out
        maximum > typemax(Int16) && throw(ArgumentError("mode exceeds Int16 storage"))
        for k in min(n, maximum):-1:1, rest in partitions(plan, n-k, k)
            push!(out, (Int16(k), rest...))
        end
        out
    end
end

"""Distinct descending doubled fermion modes: odd in NS, positive even in R."""
function fermion_partitions(plan::Combinatorics, n::Integer, ramond::Bool,
                            maximum::Integer=n)
    key = (Int(n), ramond, Int(min(n, maximum)))
    get!(plan.fermions, key) do
        out = Tuple{Vararg{Int16}}[]
        n == 0 && return push!(out, ())
        n < 0 && return out
        maximum > typemax(Int16) && throw(ArgumentError("mode exceeds Int16 storage"))
        for k in min(n, maximum):-1:1
            iseven(k) == ramond || continue
            for rest in fermion_partitions(plan, n-k, ramond, k-1)
                push!(out, (Int16(k), rest...))
            end
        end
        out
    end
end

"""Exact generalized binomial, including negative and half-integer arguments."""
function exact_binomial(plan::Combinatorics, a::Rational, j::Integer)
    j < 0 && return zero(Rational{BigInt})
    key = (Rational{BigInt}(a), Int(j))
    get!(plan.binomials, key) do
        value = one(Rational{BigInt})
        for k in 0:j-1
            value *= (a-k) // (k+1)
        end
        value
    end
end
exact_binomial(plan::Combinatorics, a::Integer, j::Integer) =
    exact_binomial(plan, a//1, j)

"""A downward-closed domain of nonnegative exponents in caller-chosen units."""
struct CoefficientDomain{N}
    limits::NTuple{N,Int}
    total::Int
    indices::Vector{NTuple{N,Int}}
    index::Dict{NTuple{N,Int},Int}
    regular::Bool
end

"""Enumerate a box intersected with a total-degree bound, once per cutoff."""
function CoefficientDomain(limits::NTuple{N,Int}, total::Integer=sum(limits)) where N
    all(>=(0), limits) || throw(ArgumentError("negative edge cutoff"))
    total >= 0 || throw(ArgumentError("negative total cutoff"))
    indices = NTuple{N,Int}[]
    work = zeros(Int, N)
    function visit(e, budget)
        if e > N
            push!(indices, Tuple(work))
        else
            for k in 0:min(limits[e], budget)
                work[e] = k
                visit(e+1, budget-k)
            end
        end
    end
    visit(1, Int(total))
    sort!(indices; by=x->(sum(x), x))
    CoefficientDomain(limits, Int(total), indices,
                      Dict(x=>i for (i,x) in enumerate(indices)),true)
end

"""Construct an explicit downward-closed domain, including unequal split-edge budgets."""
function CoefficientDomain(indices::AbstractVector{NTuple{N,Int}}) where N
    isempty(indices) && throw(ArgumentError("empty coefficient domain"))
    all(x->all(>=(0),x),indices) || throw(ArgumentError("negative exponent"))
    ordered=sort!(unique(indices);by=x->(sum(x),x))
    index=Dict(x=>i for (i,x) in enumerate(ordered))
    for key in ordered,e in 1:N
        key[e]==0 || haskey(index,Base.setindex(key,key[e]-1,e)) ||
            throw(ArgumentError("coefficient domain is not downward closed"))
    end
    limits=ntuple(e->maximum(x[e] for x in ordered),N)
    CoefficientDomain(limits,maximum(sum,ordered),ordered,index,false)
end
total_domain(edges::Integer, total::Integer) =
    CoefficientDomain(ntuple(_->Int(total), edges), total)
box_domain(limits::NTuple{N,Int}) where N = CoefficientDomain(limits)
contains(d::CoefficientDomain{N}, x::NTuple{N,Int}) where N =
    d.regular ? all(e->0 <= x[e] <= d.limits[e], 1:N) && sum(x) <= d.total : haskey(d.index,x)

"""Return (-1)^n for signed integers without negative integer exponent errors."""
@inline parity_sign(n::Integer) = isodd(n) ? -1 : 1

"""Rising factorial in the scalar's own arithmetic, without Float64 constants."""
function rising(x, n::Integer)
    n >= 0 || throw(ArgumentError("negative rising-factorial length"))
    out = one(x)
    for j in 0:n-1
        out *= x+j
    end
    out
end
