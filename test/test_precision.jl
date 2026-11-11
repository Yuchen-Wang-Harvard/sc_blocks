setprecision(BigFloat,192) do
    @testset "Verified span fallback and rejected out-of-span residuals" begin
        C=Complex{BigFloat}
        tiny=big"1e-35"
        columns=[Dict(1=>C(1),2=>C(1),3=>C(0)),Dict(1=>C(1),2=>C(1+tiny),3=>C(tiny))]
        solver=SpanSolver(columns)
        @test isnothing(solver.machine_factor)
        target=Dict(k=>C(3columns[1][k]-2columns[2][k]) for k in 1:3)
        coefficients,evidence=solve_span(solver,target)
        @test evidence.precise_fallback
        @test maximum(abs,coefficients-C[3,-2])<big"1e-20"
        @test evidence.residual<big"1e-48"
        # An unused row must still reject a target outside the proposed span.
        reference=SpanSolver([Dict(1=>C(1),2=>C(2),3=>C(3))])
        @test_throws ErrorException solve_span(reference,Dict(1=>C(1),2=>C(2),3=>C(4)))
        @test_throws ArgumentError solve_span(reference,Dict(4=>C(1)))
    end
end

@testset "Production coefficient precision escalation" begin
    results=Dict{Tuple{Int,Bool},Any}()
    for bits in (96,136,192),inserted in (false,true)
        results[(bits,inserted)]=setprecision(BigFloat,bits) do
            plan=RamondBlockPlan(BigFloat,theta_graph(),7//5,
                ((23//100)*im,(37//100)*im,(41//100)*im),(false,true,true),[1,1],inserted ? [1,-1] : [1,1])
            ramond_block(plan,total_domain(3,4))
        end
    end
    setprecision(BigFloat,192) do
        for inserted in (false,true),bits in (96,136)
            reference=results[(192,inserted)]
            comparison=results[(bits,inserted)]
            maximum_error=zero(BigFloat)
            for level in total_domain(3,4).indices
                a=get(reference,level,zeros(Complex{BigFloat},8))
                b=get(comparison,level,zeros(Complex{BigFloat},8))
                maximum_error=max(maximum_error,maximum(abs,a-b)/max(big(1),maximum(abs,a)))
            end
            @test maximum_error<(bits==96 ? big"1e-22" : big"1e-34")
            println("Precision escalation bits=",bits," inserted=",inserted,"; scaled error=",maximum_error)
        end
    end
end
