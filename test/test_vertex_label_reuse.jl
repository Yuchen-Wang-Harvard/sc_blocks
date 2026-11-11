"""Check cached-label output against independent physical Ward/Gram sewing at every grade and spin."""
function check_cached_vertex_labels(plan,vertices,domain)
    actual=ramond_block(plan,domain)
    largest=zero(BigFloat)
    for level in domain.indices
        expected=pbw_coefficient(plan.graph,vertices,level,plan.forms,plan.eta)
        computed=get(actual,level,zeros(Complex{BigFloat},2^length(level)))
        for (a,b) in ((computed,expected),(walsh!(copy(computed)),walsh!(copy(expected))))
            error=maximum(abs,a-b)/max(big(1),maximum(abs,b))
            @test error<big"1e-40"
            largest=max(largest,error)
        end
    end
    largest
end

setprecision(BigFloat,192) do
    @testset "Vertex-label product reuse, insertion topology and domain invalidation" begin
        maximum_error=zero(BigFloat)
        for (graph,ramond,parities) in ((theta_graph(),(false,true,true),(1,0,0)),
                                       (glasses_graph(),(false,true,true),(0,0,0)),
                                       (glasses_graph(),(false,false,true),(1,1,0)))
            plan=RamondBlockPlan(BigFloat,graph,7//5,
                (Complex{BigFloat}(0,11//23),Complex{BigFloat}(0,13//29),Complex{BigFloat}(0,17//31)),
                ramond,[0,0],[1,1];primary_parities=parities,cache_products=true)
            modules=map(a->a.positive.physical,plan.actions)
            vertices=[ramond[slots[2]] ? nsrr_vertex(ntuple(j->modules[slots[j]],3)) :
                ns_vertex(ntuple(j->modules[slots[j]],3)) for slots in graph.slots]
            domain=total_domain(3,6)
            for f in 0:1,eta1 in (1,-1),eta2 in (1,-1)
                !ramond[graph.slots[1][2]] && eta1==-1 && continue
                set_vertex_labels!(plan,[f,f],[eta1,eta2])
                maximum_error=max(maximum_error,check_cached_vertex_labels(plan,vertices,domain))
            end
            @test !isempty(plan.products)
            @test plan.product_domain==domain.indices
            @test plan.product_domain!==domain.indices
            old_transitions=plan.ccy_transitions
            old_hits=plan.product_cache_hits
            maximum_error=max(maximum_error,check_cached_vertex_labels(plan,vertices,domain))
            @test plan.ccy_transitions==old_transitions
            # Distinct glasses labels can have disjoint branch supports. A
            # repeated identical case must reuse its existing products.
            @test plan.product_cache_hits>old_hits
            # A box with the same maximum total degree still has a different
            # multigrade domain; no old truncated product may survive.
            changed=box_domain((2,2,2))
            SCB.prepare_product_cache!(plan,changed)
            @test isempty(plan.products)
            @test plan.product_domain==changed.indices
            maximum_error=max(maximum_error,check_cached_vertex_labels(plan,vertices,changed))
            # Re-expansion after shrinking must recompute the missing grades.
            maximum_error=max(maximum_error,check_cached_vertex_labels(plan,vertices,domain))
        end
        println("Cached vertex-label physical PBW maximum error=",maximum_error)
    end
end
