"""L mode has an integer index; G mode has a DOUBLED index."""
struct Mode
    odd::Bool
    index::Int16
end
L(n::Integer) = Mode(false, Int16(n))
G(r2::Integer) = Mode(true, Int16(r2))

"""Canonical L-descending/G-strictly-descending PBW state, with Ramond ground bit."""
struct PBWState
    l::Tuple{Vararg{Int16}}
    g::Tuple{Vararg{Int16}}
    ground::Bool
end
PBWState(l::Tuple, g::Tuple, ground::Integer=0) =
    PBWState(map(Int16,l), map(Int16,g), Bool(ground))
ground_state(bit::Integer=0) = PBWState((), (), bit)
level2(s::PBWState) = 2sum(Int, s.l; init=0) + sum(Int, s.g; init=0)
relative_parity(s::PBWState) = xor(isodd(length(s.g)), s.ground)
pbw_word(s::PBWState) = (map(n->L(-n),s.l)..., map(r->G(-r),s.g)...)

"""An SCA Verma module with momentum-local caches and explicit NS ground parity."""
mutable struct SCAModule{T<:AbstractFloat}
    ramond::Bool
    c::Complex{T}
    h::Complex{T}
    beta::Complex{T}
    primary_parity::Bool
    plan::Combinatorics
    bases::Dict{Int,Vector{PBWState}}
    indices::Dict{Int,Dict{PBWState,Int}}
    normal_cache::Dict{Tuple{Tuple{Vararg{Mode}},Bool},Dict{PBWState,Complex{T}}}
    actions::Dict{Tuple{Mode,PBWState},Dict{PBWState,Complex{T}}}
    pairings::Dict{Tuple{PBWState,PBWState},Complex{T}}
    grams::Dict{Int,Matrix{Complex{T}}}
    factors::Dict{Int,LU{Complex{T},Matrix{Complex{T}},Vector{Int}}}
    inverses::Dict{Int,Matrix{Complex{T}}}
end

"""Construct a module without converting inputs through machine precision."""
function SCAModule(::Type{T}, sector::Symbol, c; h=nothing, beta=nothing,
                   parity::Integer=0, plan=Combinatorics()) where T<:AbstractFloat
    sector in (:NS,:R) || throw(ArgumentError("sector must be :NS or :R"))
    parity in (0,1) || throw(ArgumentError("primary parity must be 0 or 1"))
    C = Complex{T}
    cc = C(c)
    r = sector == :R
    !r && isnothing(h) && throw(ArgumentError("NS module requires h"))
    r && parity != 0 && throw(ArgumentError("R ground parity is carried by its ground bit"))
    isnothing(h) && isnothing(beta) && throw(ArgumentError("weight or R beta required"))
    bb = r ? (isnothing(beta) ? sqrt(cc/24-C(h)) : C(beta)) : zero(C)
    hh = isnothing(h) ? cc/24-bb^2 : C(h)
    SCAModule{T}(r, cc, hh, bb, Bool(parity), plan,
                 Dict(), Dict(), Dict(), Dict(), Dict(), Dict(), Dict(), Dict())
end
state_parity(m::SCAModule, s::PBWState) = xor(relative_parity(s), m.primary_parity)

"""Enumerate the module's canonical basis at a doubled descendant level."""
function basis(m::SCAModule, n::Integer)
    get!(m.bases, Int(n)) do
        out = PBWState[]
        (n < 0 || (m.ramond && isodd(n))) && return out
        for total_l in 0:div(n,2), lp in partitions(m.plan,total_l),
            gp in fermion_partitions(m.plan,n-2total_l,m.ramond),
            bit in 0:Int(m.ramond)
            push!(out, PBWState(lp,gp,bit))
        end
        out
    end
end

"""Cached inverse lookup for a module basis."""
function basis_index(m::SCAModule, n::Integer)
    get!(m.indices, Int(n)) do
        Dict(s=>i for (i,s) in enumerate(basis(m,n)))
    end
end

"""Accumulate a sparse vector; no tolerance-based loss of small coefficients."""
function sparse_add!(out::Dict{K,V}, source::Dict{K,V}, scale=one(V)) where {K,V}
    iszero(scale) && return out
    for (s,value) in source
        out[s] = get(out,s,zero(V)) + scale*value
    end
    out
end
prune_zeros!(out::Dict) = filter!(p->!iszero(last(p)), out)

"""SCA supercommutator, with `nothing` denoting a central scalar."""
function commutator(m::SCAModule{T}, a::Mode, b::Mode) where T
    C = Complex{T}
    x,y = Int(a.index),Int(b.index)
    out = Tuple{C,Union{Mode,Nothing}}[]
    if !a.odd && !b.odd
        push!(out,(C(x-y),L(x+y)))
        x+y == 0 && push!(out,(m.c*(x^3-x)/12,nothing))
    elseif !a.odd
        push!(out,(C(x-y)/2,G(2x+y)))
    elseif !b.odd
        push!(out,(C(x-y)/2,G(x+2y)))
    else
        push!(out,(C(2),L(div(x+y,2))))
        x+y == 0 && push!(out,(m.c*(x*x-1)/12,nothing))
    end
    out
end

"""Normal order an operator word on a primary, sharing all intermediate words."""
function normal_order(m::SCAModule{T}, word::Tuple{Vararg{Mode}}, bit::Bool) where T
    C = Complex{T}
    key = (word,bit)
    haskey(m.normal_cache,key) && return m.normal_cache[key]
    out = Dict{PBWState,C}()
    if isempty(word)
        out[ground_state(Int(bit))] = one(C)
    elseif last(word).index >= 0
        mode = last(word)
        if mode.index == 0
            scale = m.h
            target = bit
            if mode.odd
                m.ramond || throw(ArgumentError("G0 in an NS module"))
                scale = C(0,1)*m.beta*C(1,bit ? 1 : -1)/sqrt(T(2))
                target = !bit
            end
            sparse_add!(out,normal_order(m,Base.front(word),target),scale)
        end
    else
        ordered = true
        for k in 1:length(word)-1
            a,b = word[k],word[k+1]
            equalodd = a == b && a.odd
            inversion = (a.index >= 0 && b.index < 0) ||
                (a.index < 0 && b.index < 0 &&
                 (Int(a.odd),a.index) > (Int(b.odd),b.index))
            (equalodd || inversion) || continue
            prefix,suffix = word[1:k-1],word[k+2:end]
            if !equalodd
                sparse_add!(out,normal_order(m,(prefix...,b,a,suffix...),bit),
                            a.odd && b.odd ? -one(C) : one(C))
            end
            for (coef,mode) in commutator(m,a,b)
                replacement = isnothing(mode) ? () : (mode,)
                sparse_add!(out,normal_order(m,(prefix...,replacement...,suffix...),bit),
                            equalodd ? coef/2 : coef)
            end
            ordered = false
            break
        end
        if ordered
            s = PBWState(Tuple(-x.index for x in word if !x.odd),
                         Tuple(-x.index for x in word if x.odd),Int(bit))
            out[s] = one(C)
        end
    end
    prune_zeros!(out)
    m.normal_cache[key] = out
    out
end

"""Cached action of one SCA mode on a PBW basis state."""
function act(m::SCAModule{T}, mode::Mode, state::PBWState) where T
    mode.odd && (iseven(mode.index) != m.ramond) &&
        throw(ArgumentError("G mode has wrong sector"))
    get!(m.actions, (mode,state)) do
        shift = mode.odd ? Int(mode.index) : 2Int(mode.index)
        shift > level2(state) && return Dict{PBWState,Complex{T}}()
        normal_order(m,(mode,pbw_word(state)...),state.ground)
    end
end

"""Linear action on a sparse state, reusing basis actions."""
function act_vector(m::SCAModule{T}, mode::Mode,
                    vector::Dict{PBWState,Complex{T}}) where T
    out = Dict{PBWState,Complex{T}}()
    for (s,value) in vector
        sparse_add!(out,act(m,mode,s),value)
    end
    prune_zeros!(out)
end

"""Dense action matrix between the appropriate graded PBW spaces."""
function action_matrix(m::SCAModule{T}, mode::Mode, n::Integer) where T
    target = n - (mode.odd ? Int(mode.index) : 2Int(mode.index))
    rows,cols = basis(m,target),basis(m,n)
    ix = basis_index(m,target)
    result = zeros(Complex{T},length(rows),length(cols))
    for (j,s) in enumerate(cols), (dest,value) in act(m,mode,s)
        result[ix[dest],j] += value
    end
    result
end

"""Current linear graded BPZ pairing, including odd NS primary parity."""
function pairing(m::SCAModule{T}, left::PBWState, right::PBWState) where T
    C = Complex{T}
    (level2(left) == level2(right) && state_parity(m,left) == state_parity(m,right)) ||
        return zero(C)
    get!(m.pairings,(left,right)) do
        k = length(left.g)
        factor = C(0,1)^k * parity_sign(div(k*(k-1),2)+(Int(left.ground)+Int(m.primary_parity))*k)
        vector = Dict(right=>one(C))
        for mode in pbw_word(left)
            vector = act_vector(m,Mode(mode.odd,-mode.index),vector)
            isempty(vector) && return zero(C)
        end
        factor * (m.ramond ? parity_sign(Int(left.ground)) : 1) *
            get(vector,ground_state(Int(left.ground)),zero(C))
    end
end

"""Cached full BPZ Gram matrix. This is bilinear, not Hermitian."""
function gram(m::SCAModule{T}, n::Integer) where T
    get!(m.grams,Int(n)) do
        states = basis(m,n)
        C = Complex{T}
        C[pairing(m,a,b) for a in states, b in states]
    end
end
"""Cache one LU factorization per module level for all subsequent solves."""
gram_factor(m::SCAModule,n::Integer) = get!(() -> lu(gram(m,n)),m.factors,Int(n))
"""Inverse BPZ Gram for direct sewing checks, computed from the cached LU."""
function inverse_gram(m::SCAModule{T},n::Integer) where T
    get!(m.inverses,Int(n)) do
        gram_factor(m,n) \ Matrix{Complex{T}}(I,length(basis(m,n)),length(basis(m,n)))
    end
end
