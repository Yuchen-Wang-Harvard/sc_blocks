setprecision(BigFloat,192) do
    @testset "All-NS closed-branch CCY pipeline versus physical PBW" begin
        b=big(7)/5
        q=b+inv(b)
        c=3//2+3q*q
        momenta=(23//100,37//100,41//100)
        maximum_error=zero(BigFloat)
        comparisons=0
        for graph in (theta_graph(),glasses_graph()), p in ((0,0,0),(1,0,0),(0,1,1)), a in (0,1)
            forms=[a,a]
            domain=total_domain(3,6)
            plan=NSBlockPlan(BigFloat,graph,b,momenta,forms;primary_parities=p)
            answer=ns_block(plan,domain)
            modules=ntuple(e->SCAModule(BigFloat,:NS,c;h=(q*q/4-BigFloat(momenta[e])^2)/2,parity=p[e]),3)
            vertices=[ns_vertex(ntuple(j->modules[slots[j]],3)) for slots in graph.slots]
            for level in domain.indices
                direct=sum(pbw_coefficient(graph,vertices,level,forms))
                error=abs(get(answer,level,0)-direct)/max(big(1),abs(direct))
                if error>=big"1e-42"
                    println("NS pipeline mismatch graph=",graph.slots," parity=",p," a=",a,
                            " level=",level," direct=",direct," pipeline=",get(answer,level,0)," scaled=",error)
                end
                @test error<big"1e-42"
                maximum_error=max(maximum_error,error)
                comparisons+=1
            end
        end
        println("All-NS pipeline comparisons: ",comparisons,"; maximum scaled error: ",maximum_error)
    end
end
