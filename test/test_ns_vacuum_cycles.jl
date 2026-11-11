setprecision(BigFloat,192) do
    @testset "Tetrahedron first fermionic vacuum cycles versus physical PBW" begin
        graph=tetrahedron_graph()
        weights=ntuple(e->Complex{BigFloat}((2*e+1)//(11+2*e)),6)
        central=Complex{BigFloat}(29//3)
        plan=NSRecursionPlan(BigFloat,graph,central,weights)
        modules=ntuple(e->SCAModule(BigFloat,:NS,central;h=weights[e]),6)
        vertices=[ns_vertex(ntuple(j->modules[slots[j]],3)) for slots in graph.slots]
        for triangle in ((1,2,4),(1,3,5),(2,3,6),(4,5,6)),extra in (0,first(triangle))
            target=ntuple(e->(e in triangle ? 3 : 0)+2Int(e==extra),6)
            domain=box_domain(target)
            actual=ns_c_block(plan,domain)[target]
            forms=[mod(sum(target[e] for e in slots),2) for slots in graph.slots]
            mask=sum(mod(target[e],2)<<(e-1) for e in 1:6)
            expected=pbw_coefficient(graph,vertices,target,forms)[mask+1]
            @test abs(actual-expected)/max(big(1),abs(expected))<big"1e-42"
        end
    end
end
