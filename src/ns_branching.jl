"""Appendix-B ell function, including its negative-index reflection rule."""
function ell(x::Complex{T},m::Integer,b::Complex{T}) where T<:AbstractFloat
    q=b+inv(b)
    if m<0
        value=ell(q-x,-m,b)
        return isodd(m) ? value : parity_sign(div(-m,2))*value
    end
    value=isodd(m) ? Complex{T}(sqrt(sqrt(sqrt(T(2))))) : one(Complex{T})
    for r in 0:m-1,s in 0:m-r-1
        iseven(r+s-m) && (value*=x+r*b+s/b)
    end
    value
end

"""The two Virasoro central charges in the enlarged module."""
function branch_central(b::Complex{T},copy::Integer) where T
    copy in (1,2) || throw(ArgumentError("Virasoro copy must be 1 or 2"))
    t=copy==1 ? b : inv(b)
    1+3(b+inv(b))^2/(1-t*t)
end

"""Embedded weight h_n^(copy), with n4=4n an exact integer label."""
function branch_weight(b::Complex{T},p::Complex{T},n4::Integer,copy::Integer) where T
    copy in (1,2) || throw(ArgumentError("Virasoro copy must be 1 or 2"))
    t=copy==1 ? b : inv(b)
    ((b+inv(b))^2/4-(p+n4*t/2)^2)/(2(1-t*t))
end

"""Current NS chi-vector BPZ norm for n2=2n, including intrinsic primary parity."""
function ns_branch_norm(b::Complex{T},p::Complex{T},n2::Integer,parity::Integer=0) where T
    parity in (0,1) || throw(ArgumentError("NS primary parity must be zero or one"))
    n2<0 && return ns_branch_norm(b,-p,-n2,parity)
    n2==0 && return one(Complex{T})
    Complex{T}(0,1)^mod(n2,2)*parity_sign(n2*(1+parity))*Complex{T}(2)^n2*
        ell(2p,2n2,b)/ell(b+inv(b)+2p,2n2,b)
end

"""Current Ramond chi-vector BPZ norm for n4=4n and total branch parity alpha."""
function ramond_branch_norm(b::Complex{T},p::Complex{T},n4::Integer,alpha::Integer) where T
    isodd(n4) && alpha in (0,1) || throw(ArgumentError("invalid Ramond branch label"))
    n4<0 && return ramond_branch_norm(b,-p,-n4,alpha)
    m=div(n4-1,2)
    prefactor=alpha==0 ? Complex{T}(2)^(2fld(m,2)+1) : Complex{T}(0,-1)*Complex{T}(2)^(2cld(m,2))
    prefactor*ell(2p,n4,b)/ell(b+inv(b)+2p,n4,b)
end

"""Closed normalized all-NS branching coefficient in the manuscript's signed-epsilon convention."""
function ns_branching(b::Complex{T},momenta::NTuple{3,Complex{T}},n2::NTuple{3,Int},
                      parities::NTuple{3,Int},a::Integer;ell_values=nothing) where T
    a in (0,1) && all(p->p in (0,1),parities) || throw(ArgumentError("invalid NS parity label"))
    mod(sum(n2),2)==a || return zero(Complex{T})
    q=b+inv(b)
    evaluate_ell(x,m)=isnothing(ell_values) ? ell(x,m,b) : get!(() -> ell(x,m,b),ell_values,(x,m))
    # Signed integer differences, not their mod-two representatives.
    epsilon=ntuple(j->mod(parities[j]+n2[j],2)-parities[j],3)
    eighth_root=Complex{T}(1,1)/sqrt(T(2))
    value=parity_sign(n2[3]*(1+parities[2]))*Complex{T}(0,-1)^a*
          eighth_root^(epsilon[1]-epsilon[2]-epsilon[3])
    for sign1 in (-1,1),sign3 in (-1,1)
        value*=evaluate_ell(q/2+sign1*momenta[1]+momenta[2]+sign3*momenta[3],
                           sign1*n2[1]+n2[2]+sign3*n2[3])
    end
    for j in 1:3
        value/=sqrt(evaluate_ell(2momenta[j],2n2[j])*evaluate_ell(q+2momenta[j],2n2[j]))
    end
    value
end

"""Closed all-NS vertex in the same unnormalized chi basis as the Ramond recurrence.

For m=abs(2n), the raw-vector root relative to the displayed closed formula is
(-1)^(m(m+1)/2) exp(-i*pi*(m mod 2)/4) (-i)^(p*(m mod 2))
times 2^(m/2)*sqrt(ell(2P,4n)*ell(Q+2P,4n))/ell(Q+2sign(n)P,2m).
The sign comes from the ordered chi word and the norm's BPZ phase. Multiplying
the three roots into the normalized vertex cancels all square roots exactly.
This formula performs that cancellation analytically, with reflected labels
and intrinsic NS parity retained. No numerical branch-root choice is needed.
"""
function ns_branching_raw(b::Complex{T},momenta::NTuple{3,Complex{T}},n2::NTuple{3,Int},
                          parities::NTuple{3,Int},a::Integer;ell_values=nothing) where T
    a in (0,1) && all(p->p in (0,1),parities) || throw(ArgumentError("invalid NS parity label"))
    mod(sum(n2),2)==a || return zero(Complex{T})
    q=b+inv(b)
    evaluate_ell(x,m)=isnothing(ell_values) ? ell(x,m,b) : get!(() -> ell(x,m,b),ell_values,(x,m))
    magnitudes=map(abs,n2)
    bits=map(n->mod(n,2),n2)
    exponent=n2[3]*(1+parities[2])+parities[1]*bits[1]+sum(div(m*(m+1),2) for m in magnitudes)
    value=parity_sign(exponent)*Complex{T}(0,-1)^(a+bits[2]+bits[3])*sqrt(T(2))^sum(magnitudes)
    for sign1 in (-1,1),sign3 in (-1,1)
        value*=evaluate_ell(q/2+sign1*momenta[1]+momenta[2]+sign3*momenta[3],
            sign1*n2[1]+n2[2]+sign3*n2[3])
    end
    for j in 1:3
        value/=evaluate_ell(q+2sign(n2[j])*momenta[j],2magnitudes[j])
    end
    value
end
