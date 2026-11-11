"""Enumerate two-Virasoro descendant labels at a fixed doubled Ramond grade."""
function splitting_basis(plan::SplittingPlan,n::Int)
    out=Tuple{Int,Int,Tuple,Tuple}[]
    bound=isqrt(4n+1)
    iseven(bound) && (bound-=1)
    for label in -bound:2:bound
        remaining=div(n-SCB.branch_level2(true,label),2)
        for alpha in 0:1,(a,b) in SCB.partition_pairs(plan.actions.positive.physical.plan,remaining)
            push!(out,(label,alpha,a,b))
        end
    end
    out
end

setprecision(BigFloat,192) do
    @testset "Theta psi splitting recurrence and unequal-grade Ward factorization" begin
        maximum_error=zero(BigFloat)
        for momentum in (37//100,(37//100)*im)
            plan=SplittingPlan(BranchActions(BigFloat,:R,7//5,momentum))
            ff=plan.actions.positive
            for a in (-5,-3,-1,1,3,5),b in (-5,-3,-1,1,3,5),alpha in 0:1
                direct=SCB.direct_splitting(plan,a,b,alpha)
                recursive=splitting_raw(plan,a,b,alpha)
                error=abs(direct-recursive)/max(big(1),abs(direct),abs(recursive))
                @test error<big"1e-40"
                maximum_error=max(maximum_error,error)
            end
            @test abs(splitting_raw(plan,1,-1,0)-sqrt(big(2)))<big"1e-48"
            @test abs(splitting_raw(plan,1,-1,1)-im/sqrt(big(2)))<big"1e-48"
            for n in 0:2:6,s in basis(ff.auxiliary,n)
                value=theta_psi_action(ff.auxiliary,0,s)
                @test length(value)==1
                @test abs(value[s]-SCB.parity_sign(Int(s.ground))/sqrt(big(2)))<big"1e-48"
            end
            local_wards=Dict{NTuple{3,Int},WardVertex{BigFloat,1}}()
            for leftlevel in 0:2:4,rightlevel in 0:2:4,
                left in splitting_basis(plan,leftlevel),right in splitting_basis(plan,rightlevel)
                right[2]==left[2] || continue
                out,alpha,ra,rb=right
                input,_,la,lb=left
                bra=ff_descendant(ff,out,alpha,ra,rb)
                ket=ff_descendant(ff,input,alpha,la,lb)
                actual=SCB.ff_pairing(ff,bra,theta_psi_apply(ff,leftlevel-rightlevel,ket))
                predicted=splitting_raw(plan,out,input,alpha)
                for copy in 1:2
                    ward=get!(local_wards,(out,input,copy)) do
                        weights=(branch_weight(ff.b,ff.p,out,copy),psi_weights(ff.b)[copy],branch_weight(ff.b,ff.p,input,copy))
                        ns_vertex(ntuple(j->SCAModule(BigFloat,:NS,branch_central(ff.b,copy);h=weights[j]),3))
                    end
                    first=copy==1 ? ra : rb
                    third=copy==1 ? la : lb
                    predicted*=rho(ward,0,(PBWState(first,()),ground_state(),PBWState(third,())))
                end
                error=abs(actual-predicted)/max(big(1),abs(actual),abs(predicted))
                @test error<big"1e-39"
                maximum_error=max(maximum_error,error)
            end
        end
        println("Splitting maximum scaled error=",maximum_error)
    end
end
