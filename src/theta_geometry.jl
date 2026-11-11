"""Charged free-current sewing data for the ordered (infinity,1,0) theta chart."""
struct ThetaGeometryPlan{T<:AbstractFloat}
    cutoff::Int
    covariance::Matrix{T}
    charges::Matrix{T}
end

function ThetaGeometryPlan(::Type{T},cutoff::Int=20) where T<:AbstractFloat
    cutoff>0 || throw(ArgumentError("positive oscillator cutoff required"))
    covariance=zeros(T,3cutoff,3cutoff)
    charges=zeros(T,3cutoff,2)
    for m in 1:cutoff
        i,j,k=m,cutoff+m,2cutoff+m
        charges[i,2]=charges[k,2]=-inv(sqrt(T(m)))
        charges[j,1]=parity_sign(m-1)/sqrt(T(m))
        covariance[i,k]=covariance[k,i]=-one(T)
        for n in 1:cutoff
            jn,kn=cutoff+n,2cutoff+n
            if m>=n
                covariance[i,jn]=covariance[jn,i]=-m*T(binomial(big(m-1),n-1))/sqrt(T(m*n))
            end
            value=parity_sign(m-1)*(m+n-1)*T(binomial(big(m+n-2),m-1))/sqrt(T(m*n))
            covariance[j,kn]=covariance[kn,j]=value
        end
    end
    ThetaGeometryPlan{T}(cutoff,covariance,charges)
end

"""Smallest eigenvalue of a real symmetric two-by-two matrix."""
function minimum_eigenvalue2(matrix::AbstractMatrix{T}) where T<:Real
    size(matrix)==(2,2) || throw(DimensionMismatch("two-by-two matrix required"))
    (matrix[1,1]+matrix[2,2]-sqrt((matrix[1,1]-matrix[2,2])^2+(matrix[1,2]+matrix[2,1])^2))/2
end

"""Fresh period matrix and bosonic determinant from Gaussian current sewing.

The covariance and charge pairing are bilinear: transpose, not adjoint. The
oscillator inverse metric contributes the i in each balanced sewing root.
"""
function theta_geometry(plan::ThetaGeometryPlan{T},parameters::NTuple{3}) where T
    q=map(Complex{T},parameters)
    all(x->0<abs(x)<1,q) || throw(ArgumentError("nonzero plumbing parameters inside the unit disk required"))
    roots=Complex{T}[im*exp(T(m)/2*log(qe)) for qe in q for m in 1:plan.cutoff]
    h=(roots.*plan.covariance).*transpose(roots)
    source=roots.*plan.charges
    plus=Matrix{Complex{T}}(I,3plan.cutoff,3plan.cutoff)+h
    minus=Matrix{Complex{T}}(I,3plan.cutoff,3plan.cutoff)-h
    charged=-transpose(source)*(plus\source)
    primary=Complex{T}[log(q[1])+log(q[3]) log(q[1]);
        log(q[1]) log(q[1])+log(q[2])]/2
    omega=(primary+charged)/(im*T(pi))
    log_abs_P=-(first(logabsdet(plus))+first(logabsdet(minus)))/2
    symmetry=maximum(abs,omega-transpose(omega))
    smallest=minimum_eigenvalue2(imag.(omega))
    smallest>0 || throw(DomainError(omega,"period matrix outside Siegel upper half-space"))
    (omega=omega,log_abs_P=log_abs_P,symmetry_residual=symmetry,
        minimum_imaginary_eigenvalue=smallest,mode_cutoff=plan.cutoff)
end

"""Appendix C theta characteristic for a fixed ordered tube-sign triple."""
function theta_characteristic(ramond::Bool,spins::NTuple{3,Int})
    all(s->s in (-1,1),spins) || throw(ArgumentError("tube signs must be ±1"))
    s1,s2,s3=spins
    ramond ? ((1,1),(div(1-s1*s3,2),div(1+s1*s2,2))) :
        ((0,0),(div(1+s1*s3,2),div(1+s1*s2,2)))
end

"""Theta constant with a Gaussian upper bound for the omitted lattice tail."""
function theta_constant(omega::AbstractMatrix{Complex{T}},alpha::NTuple{2,Int},beta::NTuple{2,Int};cutoff::Int=8) where T
    cutoff>=1 || throw(ArgumentError("positive theta cutoff required"))
    all(x->x in (0,1),(alpha...,beta...)) || throw(ArgumentError("binary characteristic required"))
    smallest=minimum_eigenvalue2(imag.(omega))
    smallest>0 || throw(DomainError(omega,"theta series does not converge"))
    value=zero(Complex{T})
    for n1 in -cutoff:cutoff,n2 in -cutoff:cutoff
        x,y=T(n1)+T(alpha[1])/2,T(n2)+T(alpha[2])/2
        exponent=omega[1,1]*x*x+(omega[1,2]+omega[2,1])*x*y+omega[2,2]*y*y+x*beta[1]+y*beta[2]
        value+=exp(im*T(pi)*exponent)
    end
    decay=T(pi)*smallest
    tail1=2exp(-decay*(cutoff+T(1)/2)^2)/(-expm1(-decay*(2cutoff+2)))
    bound=2*(1+sqrt(T(pi)/decay))*tail1
    (value=value,absolute_tail_bound=bound,cutoff=cutoff)
end

"""Free boson plus Majorana partition function with the manuscript's spin lift."""
function free_partition(geometry,ramond::Bool,spins::NTuple{3,Int};theta_cutoff::Int=8)
    alpha,beta=theta_characteristic(ramond,spins)
    theta=theta_constant(geometry.omega,alpha,beta;cutoff=theta_cutoff)
    imaginary=2imag.(geometry.omega)
    determinant=imaginary[1,1]*imaginary[2,2]-imaginary[1,2]*imaginary[2,1]
    determinant>0 || throw(DomainError(determinant,"nonpositive period determinant"))
    exp(3geometry.log_abs_P)*abs(theta.value)/sqrt(determinant)
end

"""Validate an integer genus-two symplectic matrix without numerical tolerance."""
function is_symplectic(matrix::AbstractMatrix{<:Integer})
    size(matrix)==(4,4) || return false
    identity2=Matrix{Int}(I,2,2)
    j=[zeros(Int,2,2) identity2;-identity2 zeros(Int,2,2)]
    matrix*j*transpose(matrix)==j
end

"""Apply the exact homology marking to a period matrix."""
function transform_period(matrix::AbstractMatrix{<:Integer},omega::AbstractMatrix)
    is_symplectic(matrix) || throw(ArgumentError("matrix is not symplectic"))
    (matrix[1:2,1:2]*omega+matrix[1:2,3:4])/(matrix[3:4,1:2]*omega+matrix[3:4,3:4])
end

"""Transport binary characteristics under the same symplectic marking."""
function transport_characteristic(matrix::AbstractMatrix{<:Integer},alpha::NTuple{2,Int},beta::NTuple{2,Int})
    is_symplectic(matrix) || throw(ArgumentError("matrix is not symplectic"))
    a,b,c,d=matrix[1:2,1:2],matrix[1:2,3:4],matrix[3:4,1:2],matrix[3:4,3:4]
    aa=mod.(d*collect(alpha)-c*collect(beta)+diag(c*transpose(d)),2)
    bb=mod.(-b*collect(alpha)+a*collect(beta)+diag(a*transpose(b)),2)
    Tuple(aa),Tuple(bb)
end

"""Recover a nearby small integer marking from freshly computed period matrices.

Enumerate isotropic lower rows (C,D), reconstruct A,B from imaginary and real
parts, and require exact integrality/symplecticity before comparing periods.
"""
function find_symplectic_map(source::AbstractMatrix{Complex{T}},target::AbstractMatrix{Complex{T}};
                             bound::Int=2,integer_tolerance=T(1)/2000) where T
    rows=[row for row in Iterators.product(ntuple(_->-bound:bound,4)...) if any(!iszero,row)]
    inverse_imag=inv(imag.(source))
    real_source=real.(source)
    w=[Complex{T}[row[1]*source[1,j]+row[2]*source[2,j]+row[2+j] for j in 1:2] for row in rows]
    best=nothing
    best_error=T(Inf)
    for (i,row0) in enumerate(rows),(j,row1) in enumerate(rows)
        row0[1]*row1[3]+row0[2]*row1[4]-row0[3]*row1[1]-row0[4]*row1[2]==0 || continue
        y=Complex{T}[target[1,1]*w[i][1]+target[1,2]*w[j][1] target[1,1]*w[i][2]+target[1,2]*w[j][2];
            target[2,1]*w[i][1]+target[2,2]*w[j][1] target[2,1]*w[i][2]+target[2,2]*w[j][2]]
        a=imag.(y)*inverse_imag
        b=real.(y)-a*real_source
        max(maximum(abs,a-round.(a)),maximum(abs,b-round.(b)))<integer_tolerance || continue
        c=Int[row0[1] row0[2];row1[1] row1[2]]
        d=Int[row0[3] row0[4];row1[3] row1[4]]
        matrix=[round.(Int,a) round.(Int,b);c d]
        is_symplectic(matrix) || continue
        error=maximum(abs,transform_period(matrix,source)-target)
        if error<best_error
            best,best_error=matrix,error
        end
    end
    isnothing(best) && throw(ErrorException("no symplectic marking found within stated search bound"))
    (matrix=best,residual=best_error,bound=bound)
end

"""Refine a printed target triple by holomorphic Newton iteration at fixed source period."""
function refine_theta_surface(plan::ThetaGeometryPlan{T},source_q::NTuple{3},target_seed::NTuple{3};
                              tolerance::T=256eps(T),max_iterations::Int=12) where T
    source=theta_geometry(plan,source_q)
    initial=theta_geometry(plan,target_seed)
    search=find_symplectic_map(ComplexF64.(source.omega),ComplexF64.(initial.omega))
    target=transform_period(search.matrix,source.omega)
    q=map(Complex{T},target_seed)
    h=eps(T)^(T(1)/3)
    components(matrix)=[matrix[1,1],matrix[1,2],matrix[2,2]]
    last_error=T(Inf)
    for iteration in 0:max_iterations
        current=theta_geometry(plan,q)
        difference=components(current.omega-target)
        error=maximum(abs,difference)
        error<=tolerance && return (source=source,target=current,source_q=map(Complex{T},source_q),target_q=q,
            matrix=search.matrix,initial_period_residual=T(search.residual),final_period_residual=error,iterations=iteration)
        iteration==max_iterations && break
        jacobian=zeros(Complex{T},3,3)
        for edge in 1:3
            plus=Base.setindex(q,q[edge]*exp(h),edge)
            minus=Base.setindex(q,q[edge]*exp(-h),edge)
            jacobian[:,edge]=components(theta_geometry(plan,plus).omega-theta_geometry(plan,minus).omega)/(2h)
        end
        step=-(jacobian\difference)
        accepted=false
        for reduction in 0:10
            candidate=ntuple(e->q[e]*exp(step[e]/2^reduction),3)
            candidate_error=maximum(abs,components(theta_geometry(plan,candidate).omega-target))
            if candidate_error<error
                q=candidate
                last_error=candidate_error
                accepted=true
                break
            end
        end
        accepted || throw(ErrorException("period refinement did not decrease the residual $error"))
    end
    throw(ErrorException("period refinement exceeded its iteration limit; residual $last_error"))
end

"""Independent numerical primitive-cycle product for the absolute bosonic Schottky factor."""
function numeric_schottky_logP(parameters::NTuple{3,Complex{T}},word_cutoff::Int) where T
    graph=theta_graph()
    domain=total_domain(3,word_cutoff)
    classes=primitive_cycles(graph,domain;first_oscillator=1)
    zero_,one_=zero(Complex{T}),one(Complex{T})
    coordinate=([zero_ one_;one_ zero_],[one_ -one_;zero_ one_],[one_ zero_;zero_ one_])
    back=(coordinate[1],[one_ one_;zero_ one_],coordinate[3])
    maps=Matrix{Complex{T}}[]
    determinants=Complex{T}[]
    for (edge,ends) in enumerate(graph.ends),direction in 0:1
        from,to=mod(ends[direction+1]-1,3)+1,mod(ends[2-direction]-1,3)+1
        push!(maps,back[to]*[zero_ parameters[edge];one_ zero_]*coordinate[from])
        push!(determinants,parity_sign(1+Int(from==1)+Int(to==1))*parameters[edge])
    end
    value=zero(T)
    for word in classes
        matrix=Matrix{Complex{T}}(I,2,2)
        determinant=one(Complex{T})
        for arc in word
            matrix=maps[arc+1]*matrix
            determinant*=determinants[arc+1]
        end
        trace=tr(matrix)
        root=sqrt(trace^2-4determinant)
        first,second=(trace+root)/2,(trace-root)/2
        large=abs(first)>abs(second) ? first : second
        multiplier=determinant/(large*large)
        abs(multiplier)<1 || throw(DomainError(multiplier,"invalid Schottky multiplier"))
        power=multiplier
        for oscillator in 1:10000
            value-=real(log1p(-power))
            power*=multiplier
            abs(power)<eps(T)/max(1,length(classes)) && break
            oscillator==10000 && throw(ErrorException("Schottky oscillator product did not converge"))
        end
    end
    (log_abs_P=value,primitive_classes=length(classes),word_cutoff=word_cutoff)
end
