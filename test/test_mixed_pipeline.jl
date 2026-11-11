"""Compare every requested coefficient and all fixed spins against independent physical PBW."""
function compare_mixed_block(plan,domain,vertices)
    actual=ramond_block(plan,domain)
    maximum_error=zero(BigFloat)
    for level in domain.indices
        expected=pbw_coefficient(plan.graph,vertices,level,plan.forms,plan.eta)
        computed=get(actual,level,zeros(Complex{BigFloat},2^length(level)))
        # Compare parity components and, separately, every fixed-spin evaluation.
        for (a,b) in ((computed,expected),(walsh!(copy(computed)),walsh!(copy(expected))))
            discrepancy=maximum(abs,a-b)/max(big(1),maximum(abs,b))
            @test discrepancy<big"1e-33"
            discrepancy<big"1e-33" || throw(ErrorException("mixed block discrepancy at $level; forms=$(plan.forms), eta=$(plan.eta): $discrepancy"))
            maximum_error=max(maximum_error,discrepancy)
        end
    end
    maximum_error
end

setprecision(BigFloat,192) do
    @testset "Mixed glasses and triangular-R tetrahedron versus independent PBW" begin
        maximum_error=zero(BigFloat)
        configurations=0
        for parity1 in 0:1,parity2 in 0:1
            plan=RamondBlockPlan(BigFloat,glasses_graph(),7//5,(23//100,37//100,41//100),
                (false,false,true),[0,0],[1,1];primary_parities=(parity1,parity2,0))
            modules=map(a->a.positive.physical,plan.actions)
            vertices=[ns_vertex((modules[1],modules[2],modules[2])),nsrr_vertex((modules[1],modules[3],modules[3]))]
            for f in 0:1,eta in (-1,1)
                set_vertex_labels!(plan,[f,f],[1,eta])
                maximum_error=max(maximum_error,compare_mixed_block(plan,total_domain(3,6),vertices))
                configurations+=1
            end
        end
        plan=RamondBlockPlan(BigFloat,tetrahedron_graph(),7//5,
            (23//100,37//100,41//100,43//100,47//100,53//100),
            (false,false,false,true,true,true),[0,0,0,0],[1,1,1,1])
        modules=map(a->a.positive.physical,plan.actions)
        vertices=[plan.ramond[slots[2]] ? nsrr_vertex(ntuple(j->modules[slots[j]],3)) :
            ns_vertex(ntuple(j->modules[slots[j]],3)) for slots in plan.graph.slots]
        for formmask in 0:15,etamask in 0:7
            iseven(count_ones(formmask)) || continue
            forms=[(formmask>>(j-1))&1 for j in 1:4]
            eta=[1;[SCB.parity_sign(etamask>>(j-1)) for j in 1:3]]
            set_vertex_labels!(plan,forms,eta)
            maximum_error=max(maximum_error,compare_mixed_block(plan,total_domain(6,4),vertices))
            configurations+=1
        end
        println("Mixed graph configurations=",configurations,"; maximum parity/fixed-spin error=",maximum_error)
    end
end
