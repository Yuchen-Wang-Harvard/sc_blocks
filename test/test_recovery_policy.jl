setprecision(BigFloat,192) do
    @testset "Recovery diagnostics never certify a lost physical sector" begin
        C=Complex{BigFloat}
        domain=total_domain(1,0)
        auxiliary=ParitySeries{1,BigFloat}((0,)=>C[1,1])
        invalid=ParitySeries{1,BigFloat}((0,)=>C[1,-1])
        @test_throws ErrorException SCB.recover_allowed_sector(invalid,auxiliary,domain)
        diagnostics=Dict{String,Any}()
        _,residual=SCB.recover_allowed_sector(invalid,auxiliary,domain;policy=:record,diagnostics=diagnostics)
        @test !diagnostics["passed"]
        @test residual==1
        @test parse(BigFloat,diagnostics["maximum_sector_residual"])==1
        delta=big"1e-25"
        rounded=ParitySeries{1,BigFloat}((0,)=>C[1+delta,1-delta])
        @test_throws ErrorException SCB.recover_allowed_sector(rounded,auxiliary,domain;tolerance=SCB.ward_pivot_tolerance(BigFloat))
        SCB.recover_allowed_sector(rounded,auxiliary,domain;diagnostics=diagnostics)
        @test diagnostics["passed"]
        @test parse(BigFloat,diagnostics["tolerance"])==inv(BigFloat(10)^8)
        below=ParitySeries{1,BigFloat}((0,)=>C[1+big"1e-9",1-big"1e-9"])
        above=ParitySeries{1,BigFloat}((0,)=>C[1+big"1e-7",1-big"1e-7"])
        SCB.recover_allowed_sector(below,auxiliary,domain;diagnostics=diagnostics)
        @test diagnostics["passed"]
        @test_throws ErrorException SCB.recover_allowed_sector(above,auxiliary,domain)
        _,residual=SCB.recover_allowed_sector(rounded,auxiliary,domain;tolerance=big"1e-20",diagnostics=diagnostics)
        @test diagnostics["passed"]
        @test abs(residual-delta/(1+delta))<big"1e-55"
        @test_throws ErrorException SCB.recover_allowed_sector(invalid,auxiliary,domain;tolerance=big"1e-20")
    end
end
