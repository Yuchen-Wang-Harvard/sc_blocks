"""Sparse multigrade series whose coefficients retain every edge-parity component."""
const ParitySeries{N,T}=Dict{NTuple{N,Int},Vector{Complex{T}}}

"""In-place Walsh transform between edge parities and fixed spin signs."""
function walsh!(values::AbstractVector)
    n=length(values)
    ispow2(n) || throw(ArgumentError("Walsh vector length must be a power of two"))
    half=1
    while half<n
        for start in 1:2half:n,j in 0:half-1
            a,b=values[start+j],values[start+j+half]
            values[start+j],values[start+j+half]=a+b,a-b
        end
        half*=2
    end
    values
end

"""Ordinary parity convolution: plumbing powers add and parity masks combine by XOR."""
function parity_product(a::ParitySeries{N,T},b::ParitySeries{N,T},domain::CoefficientDomain{N}) where {N,T}
    out=ParitySeries{N,T}()
    length_=2^N
    for (ka,va) in a,(kb,vb) in b
        key=ntuple(e->ka[e]+kb[e],N)
        contains(domain,key) || continue
        dest=get!(() -> zeros(Complex{T},length_),out,key)
        for x in 0:length_-1,y in 0:length_-1
            iszero(va[x+1]) || iszero(vb[y+1]) || (dest[xor(x,y)+1]+=va[x+1]*vb[y+1])
        end
    end
    out
end

"""Convolve a scalar even-descendant vacuum factor into each parity component."""
function parity_scalar_product(a::ParitySeries{N,T},b::Series{N,S},domain::CoefficientDomain{N}) where {N,T,S}
    out=ParitySeries{N,T}()
    for (ka,va) in a,(kb,vb) in b
        key=ntuple(e->ka[e]+kb[e],N)
        contains(domain,key) || continue
        dest=get!(() -> zeros(Complex{T},2^N),out,key)
        for j in eachindex(va)
            dest[j]+=va[j]*vb
        end
    end
    out
end

"""Coefficientwise division in the recoverable spin sector after any split diagonal extraction.

The caller must insert Theta psi on every odd Ramond cycle before invoking this method. It
rejects a nonzero numerator in an annihilated spin sector and verifies the full
forward convolution after recovering the allowed sector.
"""
function recover_allowed_sector(numerator::ParitySeries{N,T},auxiliary::ParitySeries{N,T},domain::CoefficientDomain{N};
                                tolerance::T=inv(T(10)^8),policy::Symbol=:error,
                                diagnostics=Dict{String,Any}()) where {N,T}
    policy in (:error,:record) || throw(ArgumentError("recovery policy must be error or record"))
    isfinite(tolerance) && tolerance>0 || throw(ArgumentError("positive finite recovery tolerance required"))
    C=Complex{T}
    count=2^N
    origin=ntuple(_->0,N)
    haskey(auxiliary,origin) || throw(ArgumentError("missing auxiliary ground coefficient"))
    ground=walsh!(copy(auxiliary[origin]))
    numerator_spin=Dict(key=>walsh!(copy(value)) for (key,value) in numerator)
    auxiliary_spin=Dict(key=>walsh!(copy(value)) for (key,value) in auxiliary)
    result_spin=ParitySeries{N,T}()
    empty!(diagnostics)
    diagnostics["tolerance"]=string(tolerance)
    diagnostics["policy"]=string(policy)
    maximum_sector_error=zero(T)
    for level in domain.indices
        input=get(numerator_spin,level,zeros(C,count))
        value=copy(input)
        for (shift,coefficient) in auxiliary_spin
            shift==origin && continue
            all(e->shift[e]<=level[e],1:N) || continue
            lower=get(result_spin,ntuple(e->level[e]-shift[e],N),nothing)
            isnothing(lower) || (value .-= coefficient.*lower)
        end
        scale=max(one(T),maximum(abs,input))
        for j in 1:count
            if iszero(ground[j])
                error=abs(value[j])/scale
                if error>maximum_sector_error
                    maximum_sector_error=error
                    diagnostics["maximum_sector_residual"]=string(error)
                    diagnostics["worst_level2"]=collect(level)
                    diagnostics["worst_spin_mask"]=j-1
                end
                policy==:record || error<=tolerance ||
                    throw(ErrorException("auxiliary sector residual $error exceeds tolerance $tolerance at grade $level, spin mask $(j-1); check precision, insertion and conventions"))
                value[j]=zero(C)
            else
                value[j]/=ground[j]
            end
        end
        result_spin[level]=value
    end
    out=ParitySeries{N,T}(key=>walsh!(value)/count for (key,value) in result_spin)
    forward=parity_product(out,auxiliary,domain)
    maximum_error=zero(T)
    for level in domain.indices
        expected=get(numerator,level,zeros(C,count))
        actual=get(forward,level,zeros(C,count))
        error=maximum(abs,actual-expected)/max(one(T),maximum(abs,expected))
        maximum_error=max(maximum_error,error)
    end
    diagnostics["maximum_sector_residual"]=string(maximum_sector_error)
    diagnostics["maximum_forward_residual"]=string(maximum_error)
    diagnostics["passed"]=maximum_sector_error<=tolerance && maximum_error<=16tolerance
    policy==:record || maximum_error<=16tolerance || throw(ErrorException("forward convolution residual $maximum_error exceeds $(16tolerance)"))
    out,maximum_error
end
