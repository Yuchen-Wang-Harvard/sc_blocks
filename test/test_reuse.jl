setprecision(BigFloat,192) do
    @testset "Optimized split product versus complete expanded product" begin
        for domain in (total_domain(3,6),box_domain((2,2,2)))
            base=RamondBlockPlan(BigFloat,theta_graph(),7//5,(11//23,13//29,17//31),
                (false,true,true),[0,0],[1,-1])
            @test base.vertices[1]===base.vertices[2]
            plan=SplitBlockPlan(base)
            expanded=split_domain(domain,plan.splits)
            full=split_diagonal(plan,SCB.split_enlarged_reduced(plan,expanded))
            optimized=SCB.split_enlarged_diagonal(plan,domain,expanded)
            @test plan.transpose_reuse>0
            unshared=SCB.split_enlarged_diagonal(plan,domain,expanded;reuse_transpose=false)
            for level in domain.indices
                expected=get(full,level,zeros(Complex{BigFloat},8))
                scale=max(big(1),maximum(abs,expected))
                @test maximum(abs,get(optimized,level,zeros(Complex{BigFloat},8))-expected)/scale<big"1e-45"
                @test maximum(abs,get(unshared,level,zeros(Complex{BigFloat},8))-expected)/scale<big"1e-45"
            end
            full_aux=split_diagonal(plan,SCB.split_auxiliary_block(plan,expanded))
            restricted_aux=split_diagonal(plan,SCB.split_auxiliary_block(plan,expanded;diagonal_only=true))
            @test full_aux==restricted_aux
        end
    end
end
