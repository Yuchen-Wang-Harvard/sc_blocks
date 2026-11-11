@testset "Gaussian quadrature polynomial moments" begin
    for T in (Float64,BigFloat)
        setprecision(BigFloat,136) do
            tolerance=T===Float64 ? T(1e-12) : T(1)/T(10)^35
            x,w=gauss_legendre(T,12)
            for degree in 0:23
                expected=isodd(degree) ? zero(T) : T(2)/(degree+1)
                @test abs(sum(w.*x.^degree)-expected)<tolerance
            end
            for alpha in (-1//2,1//2)
                x,w=gauss_laguerre(T,10,alpha)
                expected=alpha<0 ? sqrt(T(pi)) : sqrt(T(pi))/2
                for degree in 0:19
                    degree>0 && (expected*=degree+T(alpha))
                    @test abs(sum(w.*x.^degree)-expected)/max(one(T),expected)<tolerance
                end
            end
        end
    end
end

setprecision(BigFloat,192) do
    @testset "Complex Gamma and independent Upsilon quadratures" begin
        gamma=SCB.LogGammaPlan(BigFloat)
        @test abs(exp(SCB.loggamma_complex(gamma,1))-1)<big"1e-50"
        @test abs(exp(SCB.loggamma_complex(gamma,big(1)/2))-sqrt(big(pi)))<big"1e-50"
        for z in (Complex{BigFloat}(2//7,3//11),Complex{BigFloat}(-5//7,2//3),Complex{BigFloat}(0,7//5))
            lhs=exp(SCB.loggamma_complex(gamma,z)+SCB.loggamma_complex(gamma,1-z))
            @test abs(lhs-big(pi)/sin(big(pi)*z))/max(big(1),abs(lhs))<big"1e-48"
        end
        b=BigFloat(7)/5
        Q=b+inv(b)
        fast=LiouvillePlan(Float64,7/5;order=512,tmax=120)
        for x in (Complex{BigFloat}(Q/4,3//7),Complex{BigFloat}(Q/2,11//5),Complex{BigFloat}(b,0))
            direct=reference_log_upsilon(b,x;bits=120)
            value=log_upsilon(fast,ComplexF64(x))
            @test abs(value-direct.value)<big"1e-10"
            println("Upsilon fast/reference error=",abs(value-direct.value),"; evaluations=",direct.evaluations)
        end
        x=Complex{BigFloat}(1//5,2//7)
        first=reference_log_upsilon(b,x;bits=120).value
        second=reference_log_upsilon(b,x+b;bits=120).value
        expected=(1-2b*x)*log(b)+SCB.loggamma_complex(gamma,b*x)-SCB.loggamma_complex(gamma,1-b*x)
        @test abs(exp(second-first-expected)-1)<big"1e-32"
        ns=liouville_ns(fast,(0.3,0.5,0.7))
        rr=liouville_nsrr(fast,(0.3,0.5,0.7))
        @test abs(imag(ns[1]))<1e-11*abs(ns[1])
        @test abs(real(ns[2]))<1e-11*abs(ns[2])
        @test rr[2,1]==im*rr[1,2] && rr[2,2]==im*rr[1,1]
        @test all(x->real(x)>0,rr[1,:])
    end
end

@testset "Fresh theta geometry and independent Schottky determinant" begin
    q=(-0.027941+0.071218im,-0.057774+0.066640im,-0.026376-0.017260im)
    seed=(-0.071250-0.042549im,-0.019696+0.005144im,0.017998-0.080979im)
    first=theta_geometry(ThetaGeometryPlan(Float64,12),q)
    second=theta_geometry(ThetaGeometryPlan(Float64,20),q)
    @test maximum(abs,first.omega-second.omega)<1e-11
    @test abs(first.log_abs_P-second.log_abs_P)<1e-11
    product=numeric_schottky_logP(q,12)
    @test abs(product.log_abs_P-second.log_abs_P)<1e-11
    refined=refine_theta_surface(ThetaGeometryPlan(Float64,20),q,seed)
    @test refined.final_period_residual<1e-12
    for (rr,ns) in (((1,-1,1),(-1,-1,-1)),((-1,-1,1),(-1,-1,1)))
        alpha,beta=theta_characteristic(true,rr)
        @test transport_characteristic(refined.matrix,alpha,beta)==theta_characteristic(false,ns)
        @test free_partition(refined.source,true,rr)>0
        @test free_partition(refined.target,false,ns)>0
    end
    println("Fresh period residual=",refined.final_period_residual,
        "; Schottky/Gaussian log-P difference=",abs(product.log_abs_P-second.log_abs_P))
end
