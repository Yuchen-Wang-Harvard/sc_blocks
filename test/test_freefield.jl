include(joinpath(@__DIR__,"..","checks","tensor_pbw.jl"))

"""Sparse-vector residual normalized without discarding small coefficients."""
function scaled_vector_error(a,b)
    vector_error(a,b)/max(big(1),maximum(abs,values(a);init=big(0)),maximum(abs,values(b);init=big(0)))
end

setprecision(BigFloat,192) do
    @testset "Independent free-field and physical PBW mode actions" begin
        maximum_error=zero(BigFloat)
        for sector in (:NS,:R),parity in (sector==:NS ? (0,1) : (0,))
            ff=FreeField(BigFloat,sector,7//5,37//100;parity=parity)
            for n in 0:6,state in basis(ff.physical,n)
                vector=SCB.ff_pbw_image(ff,state)
                modes=sector==:NS ? (L(-1),L(0),L(1),L(2),G(-1),G(1),G(3)) :
                                   (L(-1),L(0),L(1),L(2),G(-2),G(0),G(2))
                for mode in modes
                    actual=ff_physical_apply(ff,mode,vector)
                    expected=Dict{UInt32,Complex{BigFloat}}()
                    for (target,coefficient) in act(ff.physical,mode,state)
                        SCB.sparse_add!(expected,SCB.ff_pbw_image(ff,target),coefficient)
                    end
                    error=scaled_vector_error(actual,expected)
                    @test error<big"1e-48"
                    maximum_error=max(maximum_error,error)
                end
            end
            reference=TensorPBWReference.Reference(BigFloat,sector,7//5,37//100;parity=parity)
            for n in 0:4,state in TensorPBWReference.tensor_basis(reference,n),copy in 1:2,mode in (-1,0,1,2)
                vector=Dict(state=>one(Complex{BigFloat}))
                fock=from_ordered_pbw(ff,vector)
                actual=ordered_pbw(ff,ff_embedded_apply(ff,copy,mode,fock))
                expected=TensorPBWReference.embedded_apply(reference,copy,mode,vector)
                error=scaled_vector_error(actual,expected)
                @test error<big"1e-46"
                maximum_error=max(maximum_error,error)
            end
        end
        println("Free-field/PBW actions: maximum scaled error=",maximum_error)
    end
    @testset "Chi-vector highest weights, current norms and reflections" begin
        maximum_highest,maximum_weight,maximum_norm=zero(BigFloat),zero(BigFloat),zero(BigFloat)
        count=0
        for momentum in (37//100,(37//100)*im),sector in (:NS,:R),
            parity in (sector==:NS ? (0,1) : (0,))
            ff=FreeField(BigFloat,sector,7//5,momentum;parity=parity)
            labels=sector==:NS ? (-6,-4,-2,0,2,4,6) : (-5,-3,-1,1,3,5)
            for n4 in labels,alpha in (sector==:NS ? (0,) : (0,1))
                vector=ff_primary(ff,n4,alpha)
                grade=sector==:NS ? div(n4*n4,4) : div(n4*n4-1,4)
                @test all(id->level2(ff.states[id])==grade,keys(vector))
                for copy in 1:2
                    weight=branch_weight(ff.b,ff.p,n4,copy)
                    expected=Dict(id=>weight*value for (id,value) in vector)
                    error=scaled_vector_error(ff_embedded_apply(ff,copy,0,vector),expected)
                    @test error<big"1e-43"
                    maximum_weight=max(maximum_weight,error)
                    for mode in 1:div(grade,2)
                        actual=ff_embedded_apply(ff,copy,mode,vector)
                        error=maximum(abs,values(actual);init=big(0))/max(big(1),maximum(abs,values(vector)))
                        @test error<big"1e-43"
                        maximum_highest=max(maximum_highest,error)
                    end
                end
                norm=ff_norm(ff,vector)
                expected=sector==:NS ? ns_branch_norm(ff.b,ff.p,div(n4,2),parity) :
                                      ramond_branch_norm(ff.b,ff.p,n4,alpha)
                error=abs(norm-expected)/max(big(1),abs(norm),abs(expected))
                @test error<big"1e-43"
                maximum_norm=max(maximum_norm,error)
                count+=1
            end
        end
        println("Free-field primary checks=",count,"; highest=",maximum_highest,
                "; weight=",maximum_weight,"; norm=",maximum_norm)
    end
end
