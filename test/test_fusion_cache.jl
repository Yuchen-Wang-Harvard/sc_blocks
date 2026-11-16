setprecision(BigFloat,192) do
    graph=PlumbingGraph([(1,3,4),(2,3,4),(1,5,6),(2,5,6)])
    central=73//10
    origin=ntuple(_->0,6)

    @testset "CCY reciprocal-pole fusion cache versus direct PBW" begin
        weights=[-2//1,-7//8,11//23,13//29,17//31,19//37]
        domain=total_domain(6,2)
        reference=VirasoroPBWReference.Reference(BigFloat,graph,central,weights)
        expected=Dict(level=>VirasoroPBWReference.coefficient(reference,level) for level in domain.indices)
        answers=[]
        for order in ((),(1,2),(2,1))
            plan=CCYPlan(BigFloat,graph,central,weights)
            # Populate the two conflicting endpoint fusion pairs in either order,
            # or leave the cache cold, before computing the complete block.
            for edge in order
                pole=SCB.ccy_pole(plan,edge,0,2,1)
                SCB.ccy_residue(plan,pole,origin,edge,2,1)
            end
            poles=ntuple(e->SCB.ccy_pole(plan,e,0,2,1),2)
            @test poles[1].central==poles[2].central==28
            @test poles[1].id==poles[2].id
            @test poles[1].x==2 && poles[2].x==1//2
            answer=virasoro_block(plan,domain)
            push!(answers,answer)
            for level in domain.indices
                @test abs(answer[level]-expected[level])/max(big(1),abs(expected[level]))<big"1e-45"
            end
        end
        @test answers[1]==answers[2]==answers[3]
        level=(0,2,0,0,0,0)
        println("Reciprocal-pole Virasoro coefficient: ",answers[1][level],
                "; direct PBW: ",expected[level],
                "; maximum scaled error: ",maximum(abs(answers[1][n]-expected[n])/max(big(1),abs(expected[n])) for n in domain.indices))
    end

    @testset "NS reciprocal-pole fusion cache versus direct PBW" begin
        weights=(-5//2,-1//1,11//23,13//29,17//31,19//37)
        domain=total_domain(6,4) # Doubled NS grades: physical total level two.
        modules=ntuple(e->SCAModule(BigFloat,:NS,central;h=weights[e]),6)
        vertices=[ns_vertex(ntuple(j->modules[slots[j]],3)) for slots in graph.slots]
        expected=Dict{NTuple{6,Int},Complex{BigFloat}}()
        for level in domain.indices
            forms=[mod(sum(level[e] for e in slots),2) for slots in graph.slots]
            mask=sum(mod(level[e],2)<<(e-1) for e in 1:6)
            expected[level]=pbw_coefficient(graph,vertices,level,forms)[mask+1]
        end
        answers=[]
        for order in ((),(1,2),(2,1))
            plan=NSRecursionPlan(BigFloat,graph,central,weights)
            for edge in order
                pole=SCB.ns_c_pole(plan,edge,0,3,1)
                # Warm both relative form channels at both incident vertices.
                for mask in (0,1<<(edge-1))
                    SCB.ns_c_residue(plan,pole,origin,edge,3,1,mask)
                end
            end
            poles=ntuple(e->SCB.ns_c_pole(plan,e,0,3,1),2)
            @test poles[1].central==poles[2].central==15
            @test poles[1].id==poles[2].id
            @test poles[1].x==2 && poles[2].x==1//2
            answer=ns_c_block(plan,domain)
            push!(answers,answer)
            for level in domain.indices
                value=get(answer,level,zero(Complex{BigFloat}))
                @test abs(value-expected[level])/max(big(1),abs(expected[level]))<big"1e-45"
            end
        end
        @test answers[1]==answers[2]==answers[3]
        level=(0,4,0,0,0,0)
        println("Reciprocal-pole NS coefficient: ",answers[1][level],
                "; direct PBW: ",expected[level],
                "; maximum scaled error: ",maximum(abs(get(answers[1],n,zero(Complex{BigFloat}))-expected[n])/max(big(1),abs(expected[n])) for n in domain.indices))
    end
end
