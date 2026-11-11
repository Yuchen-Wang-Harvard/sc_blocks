# This driver has its own executable include, so run this file separately.
using Test
isdefined(@__MODULE__,:run_cross_channel) || include(joinpath(@__DIR__,"..","checks","run_cross_channel.jl"))

setprecision(BigFloat,136) do
    @testset "Cross-channel summation order, sectors and continuum measure" begin
        q=(0.03+0.02im,-0.02+0.01im,0.01-0.04im)
        spins=((1,-1,1),(-1,-1,1))
        surface=(id="synthetic",source_q=q,target_q=q,source_spins=spins,target_spins=spins,free_ratio=(1.0,1.0))
        domain=total_domain(3,2)
        levels=[0,1,2]
        plan=cross_evaluation_plan([surface],false,domain,levels)
        ns=Dict((0,0,0)=>Complex{BigFloat}(2),(1,0,0)=>Complex{BigFloat}(0,3),
            (1,1,0)=>Complex{BigFloat}(5))
        even=evaluate_cross_series(ns,plan;form=0)
        odd=evaluate_cross_series(ns,plan;form=1)
        for spin in 1:2
            @test even[1,1,spin]==4
            @test odd[1,1,spin]==0
            @test isapprox(even[1,3,spin],abs2(2+5spins[spin][1]*spins[spin][2]*sqrt(q[1])*sqrt(q[2]));rtol=1e-14)
            @test isapprox(odd[1,3,spin],9abs(q[1]);rtol=1e-14)
        end
        rplan=cross_evaluation_plan([surface],true,domain,levels)
        coefficients=zeros(Complex{BigFloat},8)
        coefficients[1]=2
        coefficients[2]=3im
        coefficients[3]=-1+im
        rr=Dict((0,0,0)=>coefficients,(1,0,0)=>2coefficients)
        actual=evaluate_cross_series(rr,rplan)
        for spin in 1:2
            ground=2+3im*spins[spin][1]+(-1+im)*spins[spin][2]
            @test isapprox(actual[1,1,spin],abs2(ground);rtol=1e-14)
            @test isapprox(actual[1,3,spin],abs2(ground*(1+2sqrt(q[1])));rtol=1e-14)
        end
        for ramond in (false,true),mode in ("minimum","geometric")
            grid=cross_momentum_grid([surface],ramond,5,mode)
            actual=1.0
            expected=1.0
            for e in 1:3
                alpha=ramond && e>1 ? -1/2 : 1/2
                actual*=sum(grid.weights[e].*grid.nodes[e].^(2alpha+1))
                mass=alpha<0 ? sqrt(pi) : sqrt(pi)/2
                expected*=mass/(2pi*grid.reference[e]^(alpha+1))
            end
            @test isapprox(actual,expected;rtol=1e-13)
        end
        @test_throws ArgumentError cross_momentum_grid([surface],true,5,"invalid")
        reference=[0.07,0.04,0.09]
        for ramond in (false,true)
            grid=cross_momentum_grid([surface],ramond,32,"fixed";
                reference_q=reference,alpha_scheme="uniform_minus_half")
            @test grid.reference == -log.(reference)
            @test grid.alphas == fill(-1//2,3)
            # Integrate p^2 times the physical Gaussian, including the exact
            # residual relative to a different fixed reference on each edge.
            for e in 1:3
                a=grid.rates[1,e]
                actual=sum(grid.weights[e] .* grid.nodes[e].^2 .*
                    exp.(-(a-grid.reference[e]).*grid.nodes[e].^2))
                @test isapprox(actual,1/(4sqrt(pi)*a^(3/2));rtol=1e-13)
            end
        end
        @test_throws ArgumentError cross_momentum_grid([surface],true,5,"fixed")
        @test_throws ArgumentError cross_momentum_grid([surface],true,5,"fixed";reference_q=[0.1,0.2,1.0])
        @test_throws ArgumentError cross_momentum_grid([surface],true,5,"minimum";alpha_scheme="invalid")
        @test cross_order_schedule(Dict("orders"=>[1,2])) == [(1,1),(2,2)]
        @test cross_order_schedule(Dict("source_orders"=>[5,6,7],"target_orders"=>[8,9,10])) == [(5,8),(6,9),(7,10)]
        @test_throws ArgumentError cross_order_schedule(Dict("source_orders"=>[5,7]))
        @test_throws ArgumentError cross_order_schedule(Dict("source_orders"=>[5,7],"target_orders"=>[10]))
        @test_throws ArgumentError cross_order_schedule(Dict("orders"=>[0,2]))
        @test_throws ArgumentError cross_order_schedule(Dict("orders"=>[2,2]))
        @test_throws ArgumentError cross_order_schedule(Dict("orders"=>[2,1]))
        @test_throws ArgumentError cross_order_schedule(Dict("orders"=>[1],"source_orders"=>[1],"target_orders"=>[2]))
        @test cross_stability_steps(Inf,1e-4,0)==0
        @test cross_stability_steps(1e-5,1e-4,0)==1
        @test cross_stability_steps(1e-5,1e-4,1)==2
        @test cross_stability_steps(1e-3,1e-4,2)==0
        @test_throws ArgumentError cross_stability_steps(1e-5,0.0,0)
    end
end
