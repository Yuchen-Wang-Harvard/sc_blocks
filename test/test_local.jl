"""Maximum coefficient discrepancy between sparse state vectors."""
function vector_error(a,b)
    maximum((abs(get(a,k,0)-get(b,k,0)) for k in union(keys(a),keys(b)));init=big"0")
end

"""Direct check of a mode (anti)commutator on a supplied basis state."""
function relation_error(m,a,b,state,terms,scalar=0)
    C=typeof(m.c)
    lhs=act_vector(m,a,act(m,b,state))
    SCB.sparse_add!(lhs,act_vector(m,b,act(m,a,state)),a.odd && b.odd ? one(C) : -one(C))
    rhs=Dict(state=>C(scalar))
    for (coefficient,mode) in terms
        SCB.sparse_add!(rhs,act(m,mode,state),C(coefficient))
    end
    vector_error(lhs,rhs)
end

"""Check linear graded BPZ using mode action on each side independently."""
function bpz_error(m,mode,x,y)
    C=typeof(m.c)
    lhs=sum((a*pairing(m,s,y) for (s,a) in act(m,mode,x));init=zero(C))
    rhs=sum((a*pairing(m,x,s) for (s,a) in act(m,Mode(mode.odd,-mode.index),y));init=zero(C))
    rhs*=mode.odd ? C(0,1)*SCB.parity_sign(Int(state_parity(m,x))) : one(C)
    abs(lhs-rhs)/max(one(real(m.c)),abs(lhs),abs(rhs))
end

@testset "Combinatorial plans and ordered graphs" begin
    plan=Combinatorics()
    @test length(partitions(plan,10))==42
    @test partitions(plan,0)==[()]
    @test SCB.exact_binomial(plan,-1//2,4)==35//128
    @test length(total_domain(3,10).indices)==286
    @test length(box_domain((2,3,4)).indices)==60
    @test contains(CoefficientDomain((4,4,4),5),(2,1,2))
    @test !SCB.contains(CoefficientDomain((4,4,4),5),(2,2,2))
    theta=theta_graph()
    for mask in 0:7
        bits=ntuple(j->(mask>>(j-1))&1,3)
        @test sewing_sign(theta,mask)==SCB.parity_sign(bits[1]*bits[2]+bits[1]*bits[3]+bits[2]*bits[3])
        @test sewing_sign(glasses_graph(),mask)==1
    end
    @test insertion_edges(theta,[false,true,true],[1,1])==[]
    @test insertion_edges(theta,[false,true,true],[1,-1])==[2]
    @test insertion_edges(glasses_graph(),[false,true,true],[-1,-1])==[2,3]
    graph,splits=split_graph(theta,[2])
    @test graph.slots==[(1,2,3),(1,4,3),(4,0,2)]
    @test splits[2]==(2,4)
    @test length(graph.ends)==4
end

setprecision(BigFloat,192) do
    tolerance=big"1e-48"
    plan=Combinatorics()
    c=73//10
    @testset "SCA and graded BPZ, including odd NS primaries" begin
        for sector in (:NS,:R), p in (sector==:NS ? (0,1) : (0,))
            m=sector==:NS ? SCAModule(BigFloat,:NS,c;h=83//100,parity=p,plan=plan) :
                           SCAModule(BigFloat,:R,c;beta=(31//100)*im,plan=plan)
            @test pairing(m,ground_state(),ground_state())==1
            if sector==:NS
                @test gram(m,1)[1,1]==2im*m.h*SCB.parity_sign(p)
            else
                @test gram(m,0)==Complex{BigFloat}[1 0;0 -1]
                @test relation_error(m,G(0),G(0),ground_state(),[(2,L(0))],-m.c/12)<tolerance
            end
            for n in 0:6, state in basis(m,n)
                @test relation_error(m,L(2),L(-2),state,[(4,L(0))],m.c/2)<tolerance
                @test relation_error(m,L(1),L(-1),state,[(2,L(0))])<tolerance
                if sector==:NS
                    @test relation_error(m,G(3),G(-3),state,[(2,L(0))],2m.c/3)<tolerance
                    @test relation_error(m,G(1),G(-1),state,[(2,L(0))])<tolerance
                    @test relation_error(m,L(1),G(-3),state,[(2,G(-1))])<tolerance
                else
                    @test relation_error(m,G(2),G(-2),state,[(2,L(0))],m.c/4)<tolerance
                    @test relation_error(m,L(1),G(-2),state,[(3//2,G(0))])<tolerance
                end
                for mode in (L(-1),L(1),L(-2),G(sector==:NS ? -1 : 0),G(sector==:NS ? 3 : -2))
                    target=n-(mode.odd ? mode.index : 2mode.index)
                    for y in basis(m,target)
                        @test bpz_error(m,mode,state,y)<tolerance
                    end
                end
            end
            for n in 0:6
                isempty(basis(m,n)) && continue
                matrix=gram(m,n)
                @test maximum(abs,matrix-transpose(matrix))<tolerance
                @test maximum(abs,matrix*inverse_gram(m,n)-I)<tolerance
            end
        end
    end
    @testset "Independent master Ward identities" begin
        max_g,max_l=zero(BigFloat),zero(BigFloat)
        checks=0
        for sectors in (:NSNSNS,:NSRR), p in (0,1)
            modules=sectors==:NSNSNS ?
                ntuple(j->SCAModule(BigFloat,:NS,c;h=(17+6j)//29,parity=j==2 ? p : 0,plan=plan),3) :
                (SCAModule(BigFloat,:NS,c;h=83//100,parity=p,plan=plan),
                 SCAModule(BigFloat,:R,c;beta=(31//100)*im,plan=plan),
                 SCAModule(BigFloat,:R,c;beta=(47//100)*im,plan=plan))
            v=sectors==:NSNSNS ? ns_vertex(modules) : nsrr_vertex(modules)
            for levels in total_domain(3,4).indices,
                states in Iterators.product(ntuple(j->basis(modules[j],levels[j]),3)...),
                (m,n) in ((0,0),(1,0),(0,1),(-1,0),(0,-1),(1,1),(-1,1),(1,-1),(2,0),(0,-2))
                g=last(master_ward_residual(v,states,m,n))
                l=last(virasoro_ward_residual(v,states,m,n))
                max_g,max_l=max(max_g,g),max(max_l,l)
                @test g<tolerance
                @test l<tolerance
                checks+=1
            end
            if sectors==:NSRR
                @test rho(v,1,-1,(ground_state(),ground_state(1),ground_state()))==-im
                @test rho(v,0,1,(ground_state(),ground_state(),ground_state()))==1
            else
                @test rho(v,1,(ground_state(),PBWState((),(1,)),ground_state()))==1
            end
        end
        println("Master Ward checks: ",checks,"; scaled G residual: ",max_g,
                "; scaled L residual: ",max_l)
    end
end
