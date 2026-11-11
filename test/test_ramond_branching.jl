setprecision(BigFloat,192) do
    @testset "Local NSRR recurrence versus direct high-label Ward forms" begin
        maximum_error=zero(BigFloat)
        comparisons=0
        second_rows=0
        for parities in (0,1)
            plan=NSRRBranching(BigFloat,7//5,(23//100,37//100,41//100);parity=parities)
            requests=NTuple{3,Int}[]
            for a in -6:2:6,b in -5:2:5,c in -5:2:5
                a*a+b*b+c*c-2<=40 && push!(requests,(a,b,c))
            end
            prepare_branching!(plan,requests)
            for labels in requests,alpha in 0:1,gamma in 0:1
                direct=SCB.direct_branch_values(plan,labels,alpha,gamma)
                f=mod(div(labels[1],2)+alpha+gamma,2)
                for (ix,eta) in enumerate((1,-1))
                    actual=branching_raw(plan,labels,alpha,gamma,f,eta)
                    error=abs(actual-direct[ix])/max(big(1),abs(actual),abs(direct[ix]))
                    if error>=big"1e-36"
                        println("Branch recurrence mismatch: ",labels," alpha=",alpha," gamma=",gamma,
                            " eta=",eta," p=",parities," error=",error)
                    end
                    @test error<big"1e-36"
                    maximum_error=max(maximum_error,error)
                    comparisons+=1
                end
            end
            second_rows+=plan.second_ward_rows
            @test plan.maximum_local_unknowns<=4
            @test plan.recursive_values>0
        end
        symmetric=NSRRBranching(BigFloat,7//5,(23//100,37//100,37//100))
        targets=[(0,5,5),(2,5,5),(-2,-5,-5),(4,3,3)]
        prepare_branching!(symmetric,targets)
        for labels in targets,alpha in 0:1,gamma in 0:1
            direct=SCB.direct_branch_values(symmetric,labels,alpha,gamma)
            solved=symmetric.tables[(alpha,gamma)][labels]
            @test maximum(abs,SCB.form_sub(direct,solved))/max(big(1),maximum(abs,direct))<big"1e-35"
        end
        @test symmetric.second_ward_rows>0
        println("NSRR recurrence comparisons=",comparisons,"; maximum scaled error=",maximum_error,
                "; alternative Ward rows=",second_rows+symmetric.second_ward_rows)
    end
end
