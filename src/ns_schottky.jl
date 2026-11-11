"""Fundamental directed cycles, each containing one unique non-tree edge."""
function fundamental_cycles(graph::PlumbingGraph)
    adjacency=[Tuple{Int,Int,Int}[] for _ in graph.slots]
    for (edge,ends) in enumerate(graph.ends),direction in 0:1
        source=div(ends[direction+1]-1,3)+1
        target=div(ends[2-direction]-1,3)+1
        push!(adjacency[source],(edge,target,2(edge-1)+direction))
    end
    visited=falses(length(graph.slots))
    tree=falses(length(graph.ends))
    for root in eachindex(graph.slots)
        visited[root] && continue
        visited[root]=true
        pending=[root]
        for source in pending
            for (edge,target,_) in adjacency[source]
                visited[target] && continue
                visited[target]=true
                tree[edge]=true
                push!(pending,target)
            end
        end
    end
    cycles=Tuple{Int,Vector{Int}}[]
    for edge in eachindex(graph.ends)
        tree[edge] && continue
        source,target=map(h->div(h-1,3)+1,graph.ends[edge])
        paths=Dict(target=>Int[])
        pending=[target]
        for vertex in pending
            haskey(paths,source) && break
            for (branch,next,arc) in adjacency[vertex]
                tree[branch] && !haskey(paths,next) || continue
                paths[next]=[paths[vertex];arc]
                push!(pending,next)
            end
        end
        haskey(paths,source) || throw(ArgumentError("missing spanning-tree path"))
        push!(cycles,(edge,[2(edge-1);paths[source]]))
    end
    cycles
end

"""Exact fourth-root phases of the NS spin lift, fixed by native supercurrent two-point sewing."""
function ns_spin_phases(graph::PlumbingGraph)
    N=length(graph.ends)
    zero_domain=CoefficientDomain(ntuple(_->0,N))
    Q=Rational{BigInt}
    one_=series_constant(zero_domain,one(Q))
    zero_=Series{N,Q}()
    maps=plumbing_maps(graph,zero_domain)
    determinant=[parity_sign(1+Int(mod(ends[1]-1,3)==0)+Int(mod(ends[2]-1,3)==0)) for ends in graph.ends]
    phases=zeros(Int,N)
    for (chord,word) in fundamental_cycles(graph)
        matrix=(one_,zero_,zero_,one_)
        mask,raw,sign=0,0,1
        for arc in word
            edge=div(arc,2)+1
            matrix=series_matrix_multiply(maps[arc+1],matrix,zero_domain)
            mask=xor(mask,1<<(edge-1))
            determinant[edge]<0 && (raw+=1+2mod(arc,2))
            sign*=determinant[edge]
        end
        trace=get(matrix[1],ntuple(_->0,N),zero(Q))+get(matrix[4],ntuple(_->0,N),zero(Q))
        abs(trace)==1 || throw(ErrorException("unexpected fundamental cusp trace $trace"))
        raw+=trace<0 ? 2 : 0
        desired=(sewing_sign(graph,mask)<0 ? 2 : 0)+ns_mode_phase(graph,ntuple(_->0,N),mask)
        needed=mod(desired+2-raw-(sign<0 ? 2 : 0),4)
        iseven(needed) || throw(ErrorException("NS spin lift cannot be fixed by an edge sign"))
        phases[chord]=needed
    end
    phases,determinant
end

"""Native-BPZ bosonic/NS-fermionic Schottky seed in doubled plumbing powers.

The first fermionic oscillator is 3/2. With the current full graded pairing,
the mixed convolution sign cancels on every closed even-incidence vacuum
monomial, so the seed multiplies the superglobal block by ordinary convolution.
Spin lifts are calibrated locally, before the primitive product is formed.
"""
function ns_schottky(graph::PlumbingGraph,domain::CoefficientDomain{N}) where N
    N==length(graph.ends) || throw(DimensionMismatch("graph/domain"))
    all(>(0),Iterators.flatten(graph.slots)) || throw(ArgumentError("NS seed requires a closed graph"))
    Q=Rational{BigInt}
    C=Complex{Q}
    ii=C(0,1)
    origin=ntuple(_->0,N)
    integer_domain=residual_domain(domain,origin;scale=2)
    maps=plumbing_maps(graph,integer_domain)
    phases,determinant=ns_spin_phases(graph)
    fermion=series_constant(domain,one(C))
    for word in primitive_cycles(graph,domain;first_oscillator=3)
        visits=ntuple(e->count(arc->div(arc,2)+1==e,word),N)
        residual=residual_domain(domain,map(x->3x,visits);scale=2)
        one_=series_constant(residual,one(Q))
        zero_=Series{N,Q}()
        matrix=(one_,zero_,zero_,one_)
        phase,sign=0,1
        for arc in word
            edge=div(arc,2)+1
            matrix=series_matrix_multiply(maps[arc+1],matrix,residual)
            determinant[edge]<0 && (phase+=1+2mod(arc,2))
            phase+=phases[edge]
            sign*=determinant[edge]
        end
        trace=series_add(matrix[1],matrix[4])
        constant=trace[origin]
        abs(constant)==1 || throw(ErrorException("unexpected primitive cusp trace"))
        phase+=constant<0 ? 2 : 0
        inverse=series_inverse(trace,residual)
        inverse_squared=series_multiply(inverse,inverse,residual)
        ratio=series_shift(series_scale(inverse_squared,sign),visits,residual)
        catalan,power=copy(one_),copy(one_)
        for j in 1:div(residual.total,series_order(ratio,residual))
            power=series_multiply(power,ratio,residual)
            catalan=series_add(catalan,series_scale(power,binomial(big(2j),big(j))//(j+1)))
        end
        unit=series_scale(series_multiply(inverse,catalan,residual),constant)
        multiplier=series_multiply(inverse_squared,series_multiply(catalan,catalan,residual),residual)
        u=Series{N,C}()
        k=Series{N,C}()
        for (degree,value) in unit
            target=ntuple(e->visits[e]+2degree[e],N)
            contains(domain,target) && (u[target]=ii^phase*value)
        end
        for (degree,value) in multiplier
            target=ntuple(e->2visits[e]+2degree[e],N)
            contains(domain,target) && (k[target]=C(sign*value))
        end
        kp=series_constant(domain,one(C))
        for mode in 2:div(domain.total,length(word))+1
            contains(domain,map(x->(2mode-1)*x,visits)) || break
            kp=series_multiply(kp,k,domain)
            factor=series_add(series_constant(domain,one(C)),series_scale(series_multiply(u,kp,domain),-1))
            fermion=series_multiply(fermion,factor,domain)
        end
    end
    boson=Series{N,C}(map(n->2n,key)=>C(value) for (key,value) in schottky_vacuum(graph,integer_domain))
    series_multiply(boson,fermion,domain)
end

"""Independent full all-NS SCA c-recursion, with no final phase conversion."""
function ns_c_block(plan::NSRecursionPlan{T,N},domain::CoefficientDomain{N};vacuum=nothing,progress=(state)->nothing) where {T,N}
    seed=isnothing(vacuum) ? ns_schottky(plan.graph,domain) : vacuum
    numerical=Series{N,Complex{T}}(key=>Complex{T}(value) for (key,value) in seed)
    series_multiply(ns_c_reduced(plan,domain;progress=progress),numerical,domain)
end
