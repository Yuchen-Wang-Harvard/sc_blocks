"""Unoptimized theta index sum, independent of the generic network contraction order."""
function explicit_theta_coefficient(v,levels,f,left_eta,right_eta)
    T=typeof(real(v.modules[1].c))
    bases=ntuple(j->basis(v.modules[j],levels[j]),3)
    inverses=ntuple(j->inverse_gram(v.modules[j],levels[j]),3)
    result=zeros(Complex{T},8)
    for left in Iterators.product(ntuple(j->eachindex(bases[j]),3)...),
        right in Iterators.product(ntuple(j->eachindex(bases[j]),3)...)
        ls,rs=ntuple(j->bases[j][left[j]],3),ntuple(j->bases[j][right[j]],3)
        a,b=rho(v,f,left_eta,ls),rho(v,f,right_eta,rs)
        iszero(a*b) && continue
        parity=ntuple(j->Int(state_parity(v.modules[j],ls[j])),3)
        mask=sum(parity[j]<<(j-1) for j in 1:3)
        value=a*b*prod(inverses[j][left[j],right[j]] for j in 1:3)
        result[mask+1]+=SCB.parity_sign(parity[1]*parity[2]+parity[1]*parity[3]+parity[2]*parity[3])*value
    end
    result
end

setprecision(BigFloat,192) do
    @testset "Graded graph sewing and current odd-loop convention" begin
        for graph in (theta_graph(),glasses_graph(),tetrahedron_graph())
            masks=0:2^length(graph.ends)-1
            for a in masks,b in masks
                all(v->iseven(sum((b>>(e-1))&1 for e in v)),graph.slots) || continue
                @test mixed_convolution_sign(graph,a,b)==1
            end
        end
        c=73//10
        for p in (0,1)
            modules=(SCAModule(BigFloat,:NS,c;h=83//100,parity=p),
                     SCAModule(BigFloat,:R,c;beta=(31//100)*im),
                     SCAModule(BigFloat,:R,c;beta=(47//100)*im))
            v=nsrr_vertex(modules)
            for level in total_domain(3,4).indices,f in 0:1,eta1 in (-1,1),eta2 in (-1,1)
                any(isempty,ntuple(j->basis(modules[j],level[j]),3)) && continue
                network=pbw_coefficient(theta_graph(),[v,v],level,[f,f],[eta1,eta2])
                direct=explicit_theta_coefficient(v,level,f,eta1,eta2)
                @test maximum(abs,network-direct)/max(big(1),maximum(abs,direct))<big"1e-48"
            end
            left=nsrr_vertex((modules[1],modules[2],modules[2]))
            right=nsrr_vertex((modules[1],modules[3],modules[3]))
            # The former f-dependent Ward identity gave a nonzero erased term
            # here. The no-f identity makes this eta=-1 coefficient vanish.
            coefficient=pbw_coefficient(glasses_graph(),[left,right],(1,0,0),[1,1],[-1,-1])
            @test maximum(abs,coefficient)<big"1e-48"
            # At eta=+1 the same first descendant has recoverable spin support.
            coefficient=pbw_coefficient(glasses_graph(),[left,right],(1,0,0),[1,1],[1,1])
            @test abs(evaluate_spins(coefficient,[1,-1,-1]))>big"1e-6"
            @test abs(evaluate_spins(coefficient,[1,1,1]))<big"1e-48"
        end
        physical=SCAModule(BigFloat,:R,c;beta=(31//100)*im)
        aux=FermionModule(BigFloat,:R)
        @test tensor_pairing(aux,physical,FermionState((),1),ground_state(1),FermionState((),1),ground_state(1))==-im
    end
end
