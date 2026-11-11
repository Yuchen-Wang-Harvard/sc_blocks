setprecision(BigFloat,192) do
    @testset "Native superglobal vertices versus physical Ward forms" begin
        weights=map(Complex{BigFloat},(2//7,3//11,5//13))
        for mask in 0:7
            parities=ntuple(j->(mask>>(j-1))&1,3)
            vertex=ns_vertex(ntuple(j->SCAModule(BigFloat,:NS,29//3;h=weights[j],parity=parities[j]),3))
            for level in total_domain(3,6).indices
                states=ntuple(j->PBWState(ntuple(_->Int16(1),div(level[j],2)),isodd(level[j]) ? (Int16(1),) : ()),3)
                expected=rho(vertex,mod(sum(level),2),states)
                actual=SCB.ns_superglobal_vertex(weights,level,parities)
                @test abs(actual-expected)/max(big(1),abs(expected))<big"1e-48"
            end
        end
    end
    @testset "Independent native NS c-recursion versus physical PBW" begin
        maximum_error=zero(BigFloat)
        comparisons=0
        for graph in (theta_graph(),tetrahedron_graph()),odd in (false,true)
            count=length(graph.ends)
            weights=ntuple(e->Complex{BigFloat}((2*e+1)//(11+2*e)),count)
            parities=ntuple(e->Int(odd && e==1),count)
            central=Complex{BigFloat}(29//3)
            plan=NSRecursionPlan(BigFloat,graph,central,weights;primary_parities=parities)
            domain=total_domain(count,count==3 ? 10 : 6)
            actual=ns_c_block(plan,domain)
            modules=ntuple(e->SCAModule(BigFloat,:NS,central;h=weights[e],parity=parities[e]),count)
            vertices=[ns_vertex(ntuple(j->modules[slots[j]],3)) for slots in graph.slots]
            for level in domain.indices
                forms=[mod(sum(level[e] for e in slots),2) for slots in graph.slots]
                values=pbw_coefficient(graph,vertices,level,forms)
                mask=sum(mod(level[e]+parities[e],2)<<(e-1) for e in 1:count)
                expected=values[mask+1]
                computed=get(actual,level,zero(Complex{BigFloat}))
                discrepancy=abs(computed-expected)/max(big(1),abs(expected))
                @test discrepancy<big"1e-40"
                discrepancy<big"1e-40" || throw(ErrorException("NS recursion mismatch graph=$(graph.slots), odd=$odd, level=$level; actual=$computed expected=$expected"))
                maximum_error=max(maximum_error,discrepancy)
                comparisons+=1
            end
        end
        println("Independent NS c-recursion comparisons=",comparisons,"; maximum scaled error=",maximum_error)
    end
end
