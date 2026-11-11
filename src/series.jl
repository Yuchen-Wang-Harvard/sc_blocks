"""Sparse multivariate coefficient series; exponents use the caller's chosen units."""
const Series{N,T}=Dict{NTuple{N,Int},T}

"""Constant series in a fixed coefficient domain."""
function series_constant(d::CoefficientDomain{N},value::T) where {N,T}
    out=Series{N,T}()
    iszero(value) || (out[ntuple(_->0,N)]=value)
    out
end

"""Add sparse series without retaining exact zero entries."""
function series_add(a::Series{N,T},b::Series{N,T}) where {N,T}
    out=copy(a)
    sparse_add!(out,b)
    prune_zeros!(out)
end

"""Scale sparse coefficients in their existing numeric type."""
function series_scale(a::Series{N,T},scale) where {N,T}
    iszero(scale) && return Series{N,T}()
    out=Series{N,T}(key=>T(value*scale) for (key,value) in a)
    prune_zeros!(out)
end

"""Truncated convolution with independent-edge and total cutoffs checked together."""
function series_multiply(a::Series{N,T},b::Series{N,T},d::CoefficientDomain{N}) where {N,T}
    out=Series{N,T}()
    for (ka,va) in a,(kb,vb) in b
        key=ntuple(e->ka[e]+kb[e],N)
        contains(d,key) || continue
        out[key]=get(out,key,zero(T))+va*vb
    end
    prune_zeros!(out)
end

"""Lowest surviving total degree, or total cutoff plus one for the zero series."""
series_order(a::Series,d::CoefficientDomain)=minimum((sum(k) for (k,v) in a if !iszero(v));init=d.total+1)

"""Shift a series and discard terms outside the original domain."""
function series_shift(a::Series{N,T},shift::NTuple{N,Int},d::CoefficientDomain{N}) where {N,T}
    out=Series{N,T}()
    for (k,v) in a
        target=ntuple(e->k[e]+shift[e],N)
        contains(d,target) && (out[target]=v)
    end
    out
end

"""Remaining downward-closed domain after a fixed leading monomial."""
function residual_domain(d::CoefficientDomain{N},shift::NTuple{N,Int};scale::Int=1) where N
    contains(d,shift) || throw(ArgumentError("leading monomial outside domain"))
    scale>0 || throw(ArgumentError("nonpositive exponent scale"))
    d.regular && return CoefficientDomain(ntuple(e->div(d.limits[e]-shift[e],scale),N),div(d.total-sum(shift),scale))
    indices=NTuple{N,Int}[]
    for key in d.indices
        all(e->key[e]>=shift[e] && (key[e]-shift[e])%scale==0,1:N) || continue
        push!(indices,ntuple(e->div(key[e]-shift[e],scale),N))
    end
    CoefficientDomain(indices)
end

"""Formal inverse of a series with a nonzero constant term."""
function series_inverse(a::Series{N,T},d::CoefficientDomain{N}) where {N,T}
    constant=get(a,ntuple(_->0,N),zero(T))
    iszero(constant) && throw(ArgumentError("series has no invertible constant term"))
    unit=series_constant(d,one(T))
    rem=series_add(unit,series_scale(a,-inv(constant)))
    out,power=copy(unit),unit
    for j in 1:div(d.total,series_order(rem,d))
        power=series_multiply(power,rem,d)
        out=series_add(out,power)
    end
    series_scale(out,inv(constant))
end

"""Formal exponential for a zero-constant logarithm."""
function series_exp(a::Series{N,T},d::CoefficientDomain{N}) where {N,T}
    iszero(get(a,ntuple(_->0,N),zero(T))) || throw(ArgumentError("logarithm must have zero constant term"))
    out=series_constant(d,one(T))
    power=copy(out)
    for j in 1:div(d.total,series_order(a,d))
        power=series_scale(series_multiply(power,a,d),one(T)/j)
        out=series_add(out,power)
    end
    out
end

"""Row-major 2x2 multiplication over truncated sparse coefficient series."""
function series_matrix_multiply(a::NTuple{4,Series{N,T}},b::NTuple{4,Series{N,T}},d::CoefficientDomain{N}) where {N,T}
    (series_add(series_multiply(a[1],b[1],d),series_multiply(a[2],b[3],d)),
     series_add(series_multiply(a[1],b[2],d),series_multiply(a[2],b[4],d)),
     series_add(series_multiply(a[3],b[1],d),series_multiply(a[4],b[3],d)),
     series_add(series_multiply(a[3],b[2],d),series_multiply(a[4],b[4],d)))
end
