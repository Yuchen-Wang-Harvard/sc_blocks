"""Independent identity-kernel sewing on the split graph, including actual descendant Gram matrices."""
function check_identity_split(plan,levels,vertices)
    graph=plan.graph
    base=plan.original
    expanded=ntuple(e->levels[plan.owners[e]],length(graph.ends))
    modules=map(a->a.positive.physical,base.actions)
    bases=[basis(modules[plan.owners[e]],expanded[e]) for e in eachindex(graph.ends)]
    inverses=[inverse_gram(modules[plan.owners[e]],expanded[e]) for e in eachindex(graph.ends)]
    tensors=[SCB.vertex_tensor(vertices[v],ntuple(j->expanded[graph.slots[v][j]],3),base.forms[v],base.eta[v])
        for v in eachindex(base.graph.slots)]
    for edge in base.marked
        g=gram(modules[edge],levels[edge])
        push!(tensors,reshape(copy(g),size(g,1),1,size(g,2)))
    end
    expected=pbw_coefficient(base.graph,vertices,levels,base.forms,base.eta)
    wrong_order_detected=false
    for mask in 0:length(expected)-1
        lifted=SCB.split_mask(plan,mask)
        selections=[findall(s->Int(state_parity(modules[plan.owners[e]],s))==((lifted>>(e-1))&1),bases[e])
            for e in eachindex(graph.ends)]
        raw=sew_tensors(graph,tensors,inverses,selections)
        @test abs(SCB.split_sewing_sign(plan,lifted)*raw-expected[mask+1])<big"1e-40"
        abs(sewing_sign(graph,lifted)*raw-expected[mask+1])>big"1e-20" && (wrong_order_detected=true)
    end
    wrong_order_detected
end

setprecision(BigFloat,192) do
    @testset "Explicit split domains, residual budgets and exact Schottky composition" begin
        @test_throws ArgumentError CoefficientDomain([(0,0),(2,0)])
        for graph in (theta_graph(),glasses_graph()),marked in ([2],[2,3])
            expanded,splits=split_graph(graph,marked)
            domain=total_domain(3,4)
            lifted=split_domain(domain,splits)
            @test !lifted.regular
            for level in box_domain(lifted.limits).indices
                collapsed=ntuple(e->haskey(splits,e) ? max(level[e],level[splits[e][2]]) : level[e],3)
                @test contains(lifted,level)==contains(domain,collapsed)
            end
            for shift in lifted.indices
                remaining=residual_domain(lifted,shift;scale=2)
                @test all(key->contains(lifted,ntuple(e->shift[e]+2key[e],length(shift))),remaining.indices)
            end
            vacuum=schottky_vacuum(graph,domain)
            expanded_vacuum=schottky_vacuum(expanded,lifted)
            expected=Dict(ntuple(e->e<=3 ? key[e] : key[marked[e-3]],length(expanded.ends))=>value for (key,value) in vacuum)
            @test expanded_vacuum==expected
        end
    end

    @testset "Inserted NSRR CCY pipeline versus independent physical PBW" begin
        maximum_error=zero(BigFloat)
        comparisons=0
        for (graph,eta) in ((theta_graph(),[1,-1]),(theta_graph(),[-1,1]),
                            (glasses_graph(),[1,-1]),(glasses_graph(),[-1,-1])),p in (0,1),f in (0,1)
            base=RamondBlockPlan(BigFloat,graph,7//5,(23//100,37//100,41//100),
                (false,true,true),[f,f],eta;primary_parities=(p,0,0))
            plan=SplitBlockPlan(base)
            # Matrix-element insertion uses a dual leg. The identity kernel
            # must preserve each original parity component before recovery.
            for mask in 0:7
                @test SCB.split_sewing_sign(plan,SCB.split_mask(plan,mask))==sewing_sign(graph,mask)
            end
            domain=total_domain(3,6)
            actual=inserted_ramond_block(plan,domain)
            modules=map(a->a.positive.physical,base.actions)
            vertices=[nsrr_vertex(ntuple(j->modules[slots[j]],3)) for slots in graph.slots]
            if f==0
                @test check_identity_split(plan,(0,0,0),vertices)
                check_identity_split(plan,(0,2,0),vertices)
            end
            for level in domain.indices
                expected=pbw_coefficient(graph,vertices,level,[f,f],eta)
                computed=get(actual,level,zeros(Complex{BigFloat},8))
                error=maximum(abs,computed-expected)/max(big(1),maximum(abs,expected))
                if error>=big"1e-34"
                    println("Inserted pipeline mismatch graph=",graph.slots," p=",p," f=",f," eta=",eta,
                        " level=",level," error=",error," actual=",computed," expected=",expected)
                end
                @test error<big"1e-34"
                error<big"1e-34" || throw(ErrorException("stop after first inserted discrepancy"))
                maximum_error=max(maximum_error,error)
                comparisons+=1
            end
            @test base.forward_residual<big"1e-34"
            @test plan.unequal_terms>0
        end
        println("Inserted Ramond comparisons=",comparisons,"; maximum scaled error=",maximum_error)
    end
end
