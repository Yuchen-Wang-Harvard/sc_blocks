"""Lexicographic representative modulo cyclic rotation and reversal of a reduced cycle."""
function canonical_cycle(word::Vector{Int})
    reversed=[xor(a,1) for a in reverse(word)]
    best=Tuple(word)
    for sequence in (word,reversed),j in 1:length(word)
        candidate=(sequence[j:end]...,sequence[1:j-1]...)
        candidate<best && (best=candidate)
    end
    best
end

"""Primitive, cyclically reduced graph cycles whose first oscillator survives the cutoff.

Arc 2(e-1)+d starts at endpoint d+1. A cycle and its inverse share a class;
the product therefore has exponent -1, equivalent to counting both with -1/2.
"""
function primitive_cycles(graph::PlumbingGraph,d::CoefficientDomain{N};first_oscillator::Int=2) where N
    N==length(graph.ends) || throw(DimensionMismatch("graph/domain edge count"))
    first_oscillator>0 || throw(ArgumentError("nonpositive first oscillator"))
    classes=Set{Tuple{Vararg{Int}}}()
    adjacency=[Int[] for _ in graph.slots]
    for (e,ends) in enumerate(graph.ends),direction in 0:1
        push!(adjacency[div(ends[direction+1]-1,3)+1],2(e-1)+direction)
    end
    word=Int[]
    visits=zeros(Int,N)
    function visit(start,vertex)
        n=length(word)
        if n>0 && vertex==start && last(word)!=xor(first(word),1)
            primitive=all(p->n%p!=0 || any(j->word[j]!=word[mod(j-1,p)+1],1:n),1:n-1)
            primitive && push!(classes,canonical_cycle(word))
        end
        n>=div(d.total,first_oscillator) && return
        for arc in adjacency[vertex]
            n>0 && arc==xor(last(word),1) && continue
            e=div(arc,2)+1
            visits[e]+=first_oscillator
            if contains(d,Tuple(visits))
                push!(word,arc)
                endpoint=graph.ends[e][2-mod(arc,2)]
                visit(start,div(endpoint-1,3)+1)
                pop!(word)
            end
            visits[e]-=first_oscillator
        end
    end
    for v in eachindex(graph.slots)
        visit(v,v)
    end
    sort!(collect(classes))
end

"""Mobius transport matrices for the actual ordered plumbing coordinates."""
function plumbing_maps(graph::PlumbingGraph,d::CoefficientDomain{N}) where N
    Q=Rational{BigInt}
    one_=series_constant(d,one(Q))
    zero_=Series{N,Q}()
    minus=series_constant(d,-one(Q))
    coordinate=((zero_,one_,one_,zero_),(one_,minus,zero_,one_),(one_,zero_,zero_,one_))
    back=(coordinate[1],(one_,one_,zero_,one_),coordinate[3])
    maps=NTuple{4,Series{N,Q}}[]
    for (e,ends) in enumerate(graph.ends),direction in 0:1
        q=Series{N,Q}(ntuple(j->Int(j==e),N)=>one(Q))
        inversion=(zero_,q,one_,zero_)
        from,to=mod(ends[direction+1]-1,3)+1,mod(ends[2-direction]-1,3)+1
        push!(maps,series_matrix_multiply(back[to],series_matrix_multiply(inversion,coordinate[from],d),d))
    end
    maps
end

"""Exact rational bosonic Schottky vacuum on a general ordered graph.

Each primitive cycle is expanded only through the domain left after two copies
of its leading monomial. This preserves nonuniform edge budgets and avoids
expanding trace coefficients that cannot enter any requested output.
"""
function schottky_vacuum(graph::PlumbingGraph,d::CoefficientDomain{N}) where N
    Q=Rational{BigInt}
    one_=series_constant(d,one(Q))
    zero_=Series{N,Q}()
    maps=plumbing_maps(graph,d)
    logarithm=Series{N,Q}()
    for word in primitive_cycles(graph,d)
        visits=ntuple(e->count(a->div(a,2)+1==e,word),N)
        leading=map(x->2x,visits)
        residual=residual_domain(d,leading)
        matrix=(one_,zero_,zero_,one_)
        determinant_sign=1
        for arc in word
            matrix=series_matrix_multiply(maps[arc+1],matrix,residual)
            ends=graph.ends[div(arc,2)+1]
            from,to=mod(ends[mod(arc,2)+1]-1,3)+1,mod(ends[2-mod(arc,2)]-1,3)+1
            determinant_sign*=parity_sign(1+Int(from==1)+Int(to==1))
        end
        inverse=series_inverse(series_add(matrix[1],matrix[4]),residual)
        inverse_squared=series_multiply(inverse,inverse,residual)
        ratio=series_shift(series_scale(inverse_squared,determinant_sign),visits,residual)
        catalan,power=copy(one_),copy(one_)
        for j in 1:div(residual.total,series_order(ratio,residual))
            power=series_multiply(power,ratio,residual)
            coefficient=binomial(big(2j),big(j))//(j+1)
            catalan=series_add(catalan,series_scale(power,coefficient))
        end
        unit=series_multiply(inverse_squared,series_multiply(catalan,catalan,residual),residual)
        power=copy(unit)
        for s in 2:div(d.total,length(word))
            shift=map(x->s*x,visits)
            contains(d,shift) || break
            power=series_multiply(power,unit,residual)
            divisor_sum=sum(a for a in 2:s if s%a==0)
            coefficient=(big(divisor_sum)*determinant_sign^s)//s
            sparse_add!(logarithm,series_shift(power,shift,d),coefficient)
        end
    end
    prune_zeros!(logarithm)
    series_exp(logarithm,d)
end
