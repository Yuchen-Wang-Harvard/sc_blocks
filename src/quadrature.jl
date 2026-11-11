"""Gauss-Legendre nodes and weights, refined to the requested arithmetic precision."""
function gauss_legendre(::Type{T},n::Int) where T<:AbstractFloat
    n>0 || throw(ArgumentError("quadrature order must be positive"))
    seed=eigen(SymTridiagonal(zeros(n),[j/sqrt(4j*j-1) for j in 1:n-1]))
    T===Float64 && return seed.values,2 .* seed.vectors[1,:].^2
    nodes=T.(seed.values)
    weights=zeros(T,n)
    for j in 1:n
        x=nodes[j]
        for iteration in 1:20
            previous,value=one(T),x
            for k in 2:n
                previous,value=value,((2k-1)*x*value-(k-1)*previous)/k
            end
            derivative=n*(x*value-previous)/(x*x-1)
            change=value/derivative
            x-=change
            abs(change)<=8eps(T) && break
            iteration==20 && throw(ErrorException("Legendre root refinement failed"))
        end
        previous,value=one(T),x
        for k in 2:n
            previous,value=value,((2k-1)*x*value-(k-1)*previous)/k
        end
        derivative=n*(x*value-previous)/(x*x-1)
        nodes[j],weights[j]=x,2/((1-x*x)*derivative^2)
    end
    nodes,weights
end

"""Generalized Laguerre polynomial and its derivative at a real argument."""
function laguerre_value(n::Int,alpha::T,x::T) where T<:AbstractFloat
    n>0 || throw(ArgumentError("positive Laguerre order required"))
    previous,value=one(T),1+alpha-x
    for k in 2:n
        previous,value=value,((2k-1+alpha-x)*value-(k-1+alpha)*previous)/k
    end
    value,(n*value-(n+alpha)*previous)/x
end

"""Gauss-Laguerre quadrature for the NS/R Gaussian measures, alpha=±1/2."""
function gauss_laguerre(::Type{T},n::Int,alpha) where T<:AbstractFloat
    n>0 || throw(ArgumentError("quadrature order must be positive"))
    a=T(alpha)
    a in (T(1)/2,-T(1)/2) || throw(ArgumentError("implemented Laguerre weights require alpha=±1/2"))
    seed=eigen(SymTridiagonal([2j-1+Float64(a) for j in 1:n],
        [sqrt(j*(j+Float64(a))) for j in 1:n-1]))
    mass=a>0 ? sqrt(T(pi))/2 : sqrt(T(pi))
    T===Float64 && return seed.values,mass .* seed.vectors[1,:].^2
    nodes=T.(seed.values)
    weights=zeros(T,n)
    normalization=mass
    for j in 1:n
        normalization*=(j+a)/j
    end
    for j in 1:n
        x=nodes[j]
        for iteration in 1:20
            value,derivative=laguerre_value(n,a,x)
            change=value/derivative
            x-=change
            abs(change)<=8eps(T)*max(one(T),abs(x)) && break
            iteration==20 && throw(ErrorException("Laguerre root refinement failed"))
        end
        _,derivative=laguerre_value(n,a,x)
        nodes[j],weights[j]=x,normalization/(x*derivative^2)
    end
    nodes,weights
end

"""Independent double-exponential integration on a finite interval.

The reported error is the change on halving the mesh, not a rigorous bound.
Two consecutive changes must pass. Endpoint limits must be supplied by f.
"""
function tanh_sinh_integral(f,a::T,b::T;rtol::T=sqrt(eps(T)),atol::T=rtol,
                            max_refinements::Int=12) where T<:AbstractFloat
    a<b || throw(ArgumentError("quadrature interval must be ordered"))
    midpoint,half=(a+b)/2,(b-a)/2
    C=typeof(f(midpoint))
    cache=Dict{T,C}()
    extent=asinh((T(precision(T))*log(T(2))+32)/T(pi))
    last=zero(C)
    change=T(Inf)
    successes=0
    for refinement in 0:max_refinements
        step=inv(T(2)^refinement)
        value=zero(C)
        for k in -ceil(Int,extent/step):ceil(Int,extent/step)
            t=k*step
            term=get!(cache,t) do
                u=T(pi)/2*sinh(t)
                weight=half*T(pi)/2*cosh(t)/cosh(u)^2
                iszero(weight) ? zero(C) : weight*f(midpoint+half*tanh(u))
            end
            value+=term
        end
        value*=step
        change=abs(value-last)
        successes=refinement>=2 && change<=atol+rtol*abs(value) ? successes+1 : 0
        successes>=2 && return (value=value,error=change,refinements=refinement,evaluations=length(cache),converged=true)
        last=value
    end
    (value=last,error=change,refinements=max_refinements,evaluations=length(cache),converged=false)
end
