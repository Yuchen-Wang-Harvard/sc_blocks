include(joinpath(@__DIR__,"..","checks","ns_raw_vertex.jl"))

setprecision(BigFloat,192) do
    @testset "Raw closed NS vertices versus independent chi-vector Ward forms" begin
        maximum_error=zero(BigFloat)
        comparisons=0
        b=Complex{BigFloat}(7//5)
        for imaginary in (false,true),mask in 0:7
            momenta=map(x->Complex{BigFloat}(imaginary ? im*x : x),(23//100,37//100,41//100))
            parities=ntuple(j->(mask>>(j-1))&1,3)
            reference=NSRawReference.Reference(BigFloat,b,momenta;parities=parities)
            cache=Dict{Tuple{Complex{BigFloat},Int},Complex{BigFloat}}()
            for labels in Iterators.product(-3:3,-3:3,-3:3)
                sum(x*x for x in labels)<=10 || continue
                a=mod(sum(labels),2)
                direct=NSRawReference.coefficient(reference,labels,a)
                closed=SCB.ns_branching_raw(b,momenta,labels,parities,a;ell_values=cache)
                error=abs(direct-closed)/max(big(1),abs(direct),abs(closed))
                @test error<big"1e-40"
                error<big"1e-40" || throw(ErrorException("raw NS discrepancy at $labels, parities $parities: $error"))
                maximum_error=max(maximum_error,error)
                comparisons+=1
            end
        end
        println("Raw NS vertex comparisons=",comparisons,"; maximum scaled error=",maximum_error)
    end
end
