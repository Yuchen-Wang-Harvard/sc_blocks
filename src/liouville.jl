"""Exact even Bernoulli numbers for the independently evaluated complex log-Gamma."""
function bernoulli_even(count::Int)
    count>0 || throw(ArgumentError("positive Bernoulli count required"))
    work=Rational{BigInt}[]
    out=Rational{BigInt}[]
    for m in 0:2count
        push!(work,big(1)//(m+1))
        for j in m:-1:1
            work[j]=j*(work[j]-work[j+1])
        end
        m>0 && iseven(m) && push!(out,work[1])
    end
    out
end

"""Precision-local Stirling coefficients, with recurrence into the right half-plane."""
struct LogGammaPlan{T<:AbstractFloat}
    shift::T
    coefficients::Vector{T}
end

function LogGammaPlan(::Type{T}) where T<:AbstractFloat
    terms=cld(precision(T),4)+8
    bernoulli=bernoulli_even(terms)
    LogGammaPlan{T}(T(max(20,cld(precision(T),4))),[T(bernoulli[k]/(2k*(2k-1))) for k in 1:terms])
end

"""Complex log-Gamma from recurrence and Stirling, without an external special-function engine.

The logarithm follows the recurrence branch; exponentials and shift ratios are
independent of its integer multiples of 2pi*i. Poles on the real axis are rejected.
"""
function loggamma_complex(plan::LogGammaPlan{T},argument) where T
    z=Complex{T}(argument)
    iszero(imag(z)) && real(z)<=0 && isinteger(real(z)) && throw(DomainError(z,"Gamma pole"))
    n=max(0,ceil(Int,plan.shift-real(z)))
    w=z+n
    inverse=inv(w)
    squared=inverse^2
    power=inverse
    tail=zero(w)
    for coefficient in plan.coefficients
        tail+=coefficient*power
        power*=squared
    end
    result=(w-T(1)/2)*log(w)-w+log(2T(pi))/2+tail
    for j in 0:n-1
        result-=log(z+j)
    end
    result
end

"""Stable log(sinh(z)/z), including a convergent small-argument series."""
function log_sinhc(z::Complex{T}) where T
    if abs(z)<T(1)/2
        correction=zero(z)
        term=one(z)
        for k in 1:1000
            term*=z*z/((2k)*(2k+1))
            correction+=term
            abs(term)<=eps(T)*max(one(T),abs(correction)) && return log1p(correction)
        end
        throw(ErrorException("sinhc series did not converge"))
    end
    w=real(z)<0 ? -z : z
    real(w)<20 ? log(sinh(w)/w) : w-log(T(2))+log1p(-exp(-2w))-log(w)
end

"""Cancellation-controlled integrand of the defining Upsilon integral."""
function upsilon_integrand(b::T,x::Complex{T},t::T) where T
    d=(b+inv(b))/2-x
    iszero(t) && return -d*d
    iszero(d) && return zero(x)
    logarithm=2log_sinhc(d*t/2)-log_sinhc(Complex{T}(b*t/2))-log_sinhc(Complex{T}(t/(2b)))
    if t<1
        return d*d*(expm1(-t)-expm1(logarithm))/t
    end
    d*d*(exp(-t)-exp(logarithm))/t
end

"""Appendix C constants with explicit quadrature settings and momentum-local caches."""
mutable struct LiouvillePlan{T<:AbstractFloat}
    b::T
    Q::T
    nodes::Vector{T}
    weights::Vector{T}
    tmax::T
    gamma::LogGammaPlan{T}
    upsilon_cache::Dict{Complex{T},Complex{T}}
    sector_cache::Dict{Tuple{Bool,Complex{T}},Complex{T}}
    legs::Dict{Tuple{Bool,T},T}
end

function LiouvillePlan(::Type{T},b;order::Int=512,tmax=120) where T<:AbstractFloat
    coupling=T(b)
    coupling>0 || throw(ArgumentError("real positive Liouville coupling required"))
    nodes,weights=gauss_legendre(T,order)
    endpoint=T(tmax)
    endpoint>0 || throw(ArgumentError("positive integration endpoint required"))
    LiouvillePlan{T}(coupling,coupling+inv(coupling),(nodes .+ 1).*(endpoint/2),
        weights.*(endpoint/2),endpoint,LogGammaPlan(T),Dict(),Dict(),Dict())
end

"""Log Upsilon_b from its defining integral, with shifts into the convergence strip."""
function log_upsilon(plan::LiouvillePlan{T},argument) where T
    x=Complex{T}(argument)
    get!(plan.upsilon_cache,x) do
        b,Q=plan.b,plan.Q
        if real(x)<=Q/20
            step=min(b,inv(b))
            shift=(1-2step*x)*log(step)+loggamma_complex(plan.gamma,step*x)-loggamma_complex(plan.gamma,1-step*x)
            return log_upsilon(plan,x+step)-shift
        elseif real(x)>=19Q/20
            return log_upsilon(plan,Q-x)
        end
        value=zero(x)
        for j in eachindex(plan.nodes)
            value+=plan.weights[j]*upsilon_integrand(b,x,plan.nodes[j])
        end
        value
    end
end

"""NS/R combinations of Upsilon_b in the paper's normalization."""
function log_sector_upsilon(plan::LiouvillePlan{T},argument,ramond::Bool) where T
    x=Complex{T}(argument)
    get!(plan.sector_cache,(ramond,x)) do
        ramond ? log_upsilon(plan,(x+plan.b)/2)+log_upsilon(plan,(x+inv(plan.b))/2) :
            log_upsilon(plan,x/2)+log_upsilon(plan,(x+plan.Q)/2)
    end
end

"""Log of the positive square root of the conjugate pair of external leg factors."""
function liouville_leg(plan::LiouvillePlan{T},momentum,ramond::Bool) where T
    p=T(momentum)
    p>0 || throw(ArgumentError("continuum momenta must be positive"))
    get!(plan.legs,(ramond,p)) do
        real(log_sector_upsilon(plan,Complex{T}(0,2p),ramond))
    end
end

"""The two all-NS structure constants; the odd coefficient includes exactly 2i."""
function liouville_ns(plan::LiouvillePlan{T},momenta::NTuple{3}) where T
    p=map(T,momenta)
    total=sum(p)
    arguments=(Complex{T}(plan.Q/2,total),ntuple(j->Complex{T}(plan.Q/2,total-2p[j]),3)...)
    numerator=log_upsilon(plan,plan.b)-log(T(2))-3log(plan.b)+sum(liouville_leg(plan,x,false) for x in p)
    ntuple(j->(j==1 ? one(Complex{T}) : Complex{T}(0,2))*
        exp(numerator-sum(log_sector_upsilon(plan,x,j==2) for x in arguments)),2)
end

"""NSRR constants indexed [f+1, eta=+1 then -1], with the fixed odd-label exchange."""
function liouville_nsrr(plan::LiouvillePlan{T},momenta::NTuple{3}) where T
    p1,p2,p3=map(T,momenta)
    numerator=log_upsilon(plan,plan.b)-log(T(4))-3log(plan.b)+
        liouville_leg(plan,p1,false)+liouville_leg(plan,p2,true)+liouville_leg(plan,p3,true)
    a=(Complex{T}(plan.Q/2,p2+p3+p1),Complex{T}(plan.Q/2,p2+p3-p1))
    d=(Complex{T}(plan.Q/2,p1+p2-p3),Complex{T}(plan.Q/2,p1-p2+p3))
    out=zeros(Complex{T},2,2)
    for eta in 1:2
        denominator=sum(log_sector_upsilon(plan,x,eta==1) for x in a)+
            sum(log_sector_upsilon(plan,x,eta==2) for x in d)
        out[1,eta]=exp(numerator-denominator)
    end
    out[2,1],out[2,2]=im*out[1,2],im*out[1,1]
    out
end

"""Direct defining integrand for the independent quadrature, with its endpoint Taylor limit."""
function upsilon_reference_integrand(b::T,x::Complex{T},t::T) where T
    d=(b+inv(b))/2-x
    threshold=eps(T)^(T(1)/3)/max(one(T),abs(d),b,inv(b))
    if t<threshold
        correction=d*d/12-(b*b+inv(b*b))/24
        return d*d*(-1+(T(1)/2-correction)*t-t*t/6)
    end
    (d*d*exp(-t)-sinh(d*t/2)^2/(sinh(b*t/2)*sinh(t/(2b))))/t
end

"""Independent tanh-sinh evaluation in the convergence strip, with an explicit tail cutoff."""
function reference_log_upsilon(b::T,x::Complex{T};bits::Int=precision(T)-24) where T<:AbstractFloat
    Q=b+inv(b)
    rate=min(real(x),Q-real(x),one(T))
    rate>0 || throw(ArgumentError("reference integral requires 0<Re(x)<Q"))
    endpoint=max(T(64),(bits*log(T(2))+32)/rate)
    cuts=T[0,1,4,16,64]
    while last(cuts)<endpoint
        push!(cuts,min(4last(cuts),endpoint))
    end
    tolerance=T(2)^(-bits)
    value=zero(x)
    error=zero(T)
    evaluations=0
    for j in 1:length(cuts)-1
        result=tanh_sinh_integral(t->upsilon_reference_integrand(b,x,t),cuts[j],cuts[j+1];
            rtol=tolerance,atol=tolerance/length(cuts))
        result.converged || throw(ErrorException("reference Upsilon quadrature did not converge"))
        value+=result.value
        error+=result.error
        evaluations+=result.evaluations
    end
    (value=value,estimated_mesh_error=error,evaluations=evaluations,endpoint=endpoint)
end
