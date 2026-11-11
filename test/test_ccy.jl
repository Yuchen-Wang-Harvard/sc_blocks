include(joinpath(@__DIR__,"..","checks","virasoro_pbw.jl"))

setprecision(BigFloat,192) do
    @testset "Exact series algebra and graph Schottky domain" begin
        d=total_domain(3,6)
        q=Series{3,Rational{BigInt}}((0,0,0)=>1,(1,0,0)=>-1,(0,1,0)=>1//2)
        unit=series_constant(d,big(1)//1)
        @test series_multiply(q,series_inverse(q,d),d)==unit
        @test series_multiply(series_exp(Series{3,Rational{BigInt}}((1,0,0)=>1),d),
                              series_exp(Series{3,Rational{BigInt}}((1,0,0)=>-1),d),d)==unit
        bigdomain=box_domain((4,4,4))
        smalldomain=CoefficientDomain((4,3,2),7)
        full=schottky_vacuum(theta_graph(),bigdomain)
        small=schottky_vacuum(theta_graph(),smalldomain)
        @test small==filter(p->SCB.contains(smalldomain,first(p)),full)
        @test get(small,(2,0,2),0)==1
        @test get(small,(0,2,2),0)==1
        @test get(small,(2,2,0),0)==1
    end
    @testset "CCY versus independent finite-c Virasoro sewing" begin
        maxima=BigFloat[]
        for (graph,weights,external,cutoff) in (
            (theta_graph(),[17//23,19//29,31//37],Dict{Int,Rational{Int}}(),5),
            (glasses_graph(),[17//23,19//29,31//37],Dict{Int,Rational{Int}}(),5),
            (PlumbingGraph([(0,1,1)]),[19//29],Dict(1=>17//23),6),
            (PlumbingGraph([(0,0,1),(1,0,0)]),[19//29],Dict(1=>13//17,2=>17//23,5=>31//37,6=>37//41),6),
            (tetrahedron_graph(),[17//23,19//29,31//37,37//41,41//43,43//47],Dict{Int,Rational{Int}}(),3))
            domain=total_domain(length(weights),cutoff)
            plan=CCYPlan(BigFloat,graph,73//10,weights;external=external)
            reference=VirasoroPBWReference.Reference(BigFloat,graph,73//10,weights;external=external)
            answer=virasoro_block(plan,domain)
            maximum_error=zero(BigFloat)
            for level in domain.indices
                value=VirasoroPBWReference.coefficient(reference,level)
                error=abs(answer[level]-value)/max(big(1),abs(value))
                if error>=big"1e-45"
                    println("CCY mismatch: graph=",graph.slots," level=",level,
                            " direct=",value," CCY=",answer[level]," scaled=",error)
                end
                @test error<big"1e-45"
                maximum_error=max(maximum_error,error)
            end
            push!(maxima,maximum_error)
            println("Virasoro check: ",length(domain.indices)," coefficients; scaled error=",maximum_error,
                    "; transitions=",plan.transitions,"; seed terms=",plan.seed_terms)
        end
        @test maximum(maxima)<big"1e-45"
    end
end
