setprecision(BigFloat,192) do
    tolerance=big"1e-48"
    @testset "Auxiliary Clifford algebra, BPZ and master Ward" begin
        max_residual=zero(BigFloat)
        count=0
        for sectors in ((:NS,:NS,:NS),(:NS,:R,:R))
            modules=map(s->FermionModule(BigFloat,s),sectors)
            v=FermionVertex(modules)
            for levels in total_domain(3,6).indices,
                states in Iterators.product(ntuple(j->basis(modules[j],levels[j]),3)...),
                (m,n) in ((0,0),(1,0),(0,1),(-1,0),(0,-1),(1,1),(-1,1),(1,-1))
                error=last(master_ward_residual(v,states,m,n))
                max_residual=max(error,max_residual)
                @test error<tolerance
                count+=1
            end
            for module_ in modules,n in 0:6,state in basis(module_,n)
                modes=module_.ramond ? (-4,-2,0,2,4) : (-3,-1,1,3)
                for r in modes,s in modes
                    vector=Dict{FermionState,Complex{BigFloat}}()
                    for (a,b) in ((r,s),(s,r)),(x,c) in act_psi(module_,b,state)
                        SCB.sparse_add!(vector,act_psi(module_,a,x),c)
                    end
                    target=Dict(state=>Complex{BigFloat}(r+s==0 ? 1 : 0))
                    @test vector_error(vector,target)<tolerance
                end
                for r in modes,y in basis(module_,n-r)
                    lhs=sum((a*pairing(module_,x,y) for (x,a) in act_psi(module_,r,state));init=zero(Complex{BigFloat}))
                    rhs=sum((a*pairing(module_,state,x) for (x,a) in act_psi(module_,-r,y));init=zero(Complex{BigFloat}))
                    @test abs(lhs+im*SCB.parity_sign(Int(relative_parity(state)))*rhs)<tolerance
                end
            end
        end
        println("Auxiliary master Ward checks: ",count,"; scaled residual: ",max_residual)
    end
end
