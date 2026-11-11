setprecision(BigFloat,192) do
    @testset "Ordinary NSRR CCY pipeline versus independent physical PBW" begin
        maximum_error=zero(BigFloat)
        comparisons=0
        for graph in (theta_graph(),glasses_graph()),p in (0,1),f in (0,1),
            eta in (graph.slots==theta_graph().slots ? ([1,1],[-1,-1]) : ([1,1],))
            plan=RamondBlockPlan(BigFloat,graph,7//5,(23//100,37//100,41//100),(false,true,true),[f,f],eta;
                primary_parities=(p,0,0))
            domain=total_domain(3,6)
            actual=ramond_block(plan,domain)
            modules=map(a->a.positive.physical,plan.actions)
            vertices=[nsrr_vertex(ntuple(j->modules[slots[j]],3)) for slots in graph.slots]
            for level in domain.indices
                expected=pbw_coefficient(graph,vertices,level,[f,f],eta)
                computed=get(actual,level,zeros(Complex{BigFloat},8))
                error=maximum(abs,computed-expected)/max(big(1),maximum(abs,expected))
                if error>=big"1e-35"
                    println("Ramond pipeline mismatch graph=",graph.slots," p=",p," f=",f," eta=",eta,
                        " level=",level," error=",error," actual=",computed," expected=",expected)
                end
                @test error<big"1e-35"
                maximum_error=max(maximum_error,error)
                comparisons+=1
            end
            @test plan.forward_residual<big"1e-35"
        end
        @test RamondBlockPlan(BigFloat,theta_graph(),7//5,(23//100,37//100,41//100),
            (false,true,true),[0,0],[1,-1]).marked==[2]
        println("Ordinary Ramond comparisons=",comparisons,"; maximum scaled error=",maximum_error)
    end
end
