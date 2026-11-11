setprecision(BigFloat,192) do
    @testset "Verified L_±1/L_0 branching actions and Ramond parity reuse" begin
        maximum_error=zero(BigFloat)
        for sector in (:NS,:R),momentum in (37//100,(37//100)*im)
            plan=BranchActions(BigFloat,sector,7//5,momentum)
            ff=plan.positive
            labels=sector==:NS ? (-6,-4,-2,0,2,4,6) : (-5,-3,-1,1,3,5)
            for label in labels,alpha in (sector==:NS ? (0,) : (0,1)),mode in (-1,0,1)
                terms=branch_action(plan,label,alpha,mode)
                primary=ff_primary(ff,label,alpha)
                actual=ff_physical_apply(ff,L(mode),primary)
                reconstructed=SCB.reconstruct_action(ff,terms,alpha)
                error=scaled_vector_error(actual,reconstructed)
                @test error<big"1e-43"
                maximum_error=max(maximum_error,error)
            end
            if sector==:R
                independent=BranchActions(BigFloat,:R,7//5,momentum;parity_reuse=false)
                for label in (-5,-3,-1,1,3,5),mode in (-1,0,1)
                    reused=branch_action(plan,label,1,mode)
                    direct=branch_action(independent,label,1,mode)
                    @test length(reused)==length(direct)
                    for (a,b) in zip(reused,direct)
                        @test (a.label,a.first,a.second)==(b.label,b.first,b.second)
                        @test abs(a.coefficient-b.coefficient)/max(big(1),abs(b.coefficient))<big"1e-43"
                    end
                end
            end
        end
        println("Branch action full-vector maximum scaled error=",maximum_error)
    end
    @testset "Single-leg Virasoro factors versus independent Ward reduction" begin
        modules=ntuple(j->SCAModule(BigFloat,:NS,73//10;h=(17+6j)//29),3)
        v=ns_vertex(modules)
        weights=map(m->m.h,modules)
        for n in 0:7,word in partitions(modules[1].plan,n),slot in 1:3
            states=ntuple(j->j==slot ? PBWState(word,()) : ground_state(),3)
            actual=single_leg_vertex(slot,word,weights)
            expected=rho(v,0,states)
            @test abs(actual-expected)/max(big(1),abs(expected))<big"1e-45"
        end
    end
end
