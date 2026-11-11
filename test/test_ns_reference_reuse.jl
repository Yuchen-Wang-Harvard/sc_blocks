setprecision(BigFloat,192) do
    @testset "Complex superglobal vertices retain the native BPZ phase" begin
        weights=(Complex{BigFloat}(2//7,1//9),Complex{BigFloat}(3//11,-2//13),Complex{BigFloat}(5//13,1//17))
        for parities in ((0,0,0),(1,0,1))
            vertex=ns_vertex(ntuple(j->SCAModule(BigFloat,:NS,29//3;h=weights[j],parity=parities[j]),3))
            for level in total_domain(3,4).indices
                states=ntuple(j->PBWState(ntuple(_->Int16(1),div(level[j],2)),isodd(level[j]) ? (Int16(1),) : ()),3)
                expected=rho(vertex,mod(sum(level),2),states)
                actual=SCB.ns_superglobal_vertex(weights,level,parities)
                @test abs(actual-expected)/max(big(1),abs(expected))<big"1e-48"
            end
        end
    end
    @testset "Compact theta reference versus complete vertex-cache path" begin
        for complex_weights in (false,true),odd in (false,true)
            weights=ntuple(e->Complex{BigFloat}((2*e+1)//(11+2*e),complex_weights ? e//19 : 0),3)
            parities=(Int(odd),0,0)
            compact=NSRecursionPlan(BigFloat,theta_graph(),29//3,weights;primary_parities=parities)
            cached=NSRecursionPlan(BigFloat,theta_graph(),29//3,weights;primary_parities=parities,compact_theta=false)
            domain=box_domain((6,6,6))
            updates=NamedTuple[]
            first=ns_c_block(compact,domain;progress=state->push!(updates,state))
            second=ns_c_block(cached,domain)
            for key in domain.indices
                @test abs(first[key]-second[key])/max(big(1),abs(second[key]))<big"1e-46"
            end
            @test compact.compact_theta && !cached.compact_theta
            @test isempty(compact.vertices) && !isempty(cached.vertices)
            @test compact.seed_terms==cached.seed_terms
            @test compact.transitions==cached.transitions
            @test count(state->state.stage=="completed_parity",updates)==8
        end
    end
end
