@testset "Owned MPFR scratch precision and aliasing" begin
    for bits in (96,136,192)
        setprecision(BigFloat,bits) do
            numbers=Complex{BigFloat}[0,1,-1,im,BigFloat(2)/7+im*BigFloat(3)/11,
                BigFloat(-101)/19+im*BigFloat(31)/23]
            work=SCB.complex_scratch(Complex{BigFloat})
            for a in numbers,b in numbers
                before_a,before_b=deepcopy(a),deepcopy(b)
                expected=a*b
                actual=SCB.scratch_product!(work.product,a,b,work)
                tolerance=32eps(BigFloat)*max(big(1),abs(expected))
                @test abs(actual-expected)<=tolerance
                @test a==before_a && b==before_b
                snapshot=deepcopy(actual)
                @test abs(SCB.scratch_product!(actual,actual,b,work)-snapshot*b)<=
                    32eps(BigFloat)*max(big(1),abs(snapshot*b))
                snapshot=deepcopy(actual)
                @test abs(SCB.scratch_product!(actual,a,actual,work)-a*snapshot)<=
                    32eps(BigFloat)*max(big(1),abs(a*snapshot))
                snapshot=deepcopy(actual)
                @test abs(SCB.scratch_product!(actual,actual,actual,work)-snapshot*snapshot)<=
                    32eps(BigFloat)*max(big(1),abs(snapshot*snapshot))
                value=Complex(BigFloat(0),BigFloat(0))
                @test abs(SCB.sum_product(value,a,b,work)-expected)<=tolerance
                saved=deepcopy(value)
                SCB.scratch_product!(work.product,b,a,work)
                @test value==saved
            end
        end
    end
    @test SCB.sum_product(1.0+im,2.0+im,3.0-im,SCB.complex_scratch(ComplexF64))==8.0+2im
end
