"""Private MPFR work buffers; never alias cached coefficients or returned series values."""
struct BigComplexScratch
    product::Complex{BigFloat}
    parts::NTuple{4,BigFloat}
end

BigComplexScratch()=BigComplexScratch(Complex(BigFloat(0),BigFloat(0)),ntuple(_->BigFloat(0),4))

"""Use reusable MPFR buffers for BigFloat and ordinary arithmetic for other types."""
complex_scratch(::Type{Complex{BigFloat}})=BigComplexScratch()
complex_scratch(::Type{Complex{T}}) where T=nothing

"""Multiply into an owned complex buffer; either input may alias the destination.

All four real products are computed before either destination component is
written. MPFR rounding follows Julia's current BigFloat rounding mode. The
library is the MPFR already shipped with Julia, not an external backend.
"""
function scratch_product!(destination::Complex{BigFloat},a::Complex{BigFloat},b::Complex{BigFloat},work::BigComplexScratch)
    ac,bd,ad,bc=work.parts
    mode=Base.MPFR.rounding_raw(BigFloat)
    # The owned destination has distinct real/imaginary buffers. Write its
    # imaginary part first so full-input aliasing remains safe in either path.
    if iszero(imag(b))
        ccall((:mpfr_mul,Base.MPFR.libmpfr),Cint,(Ref{BigFloat},Ref{BigFloat},Ref{BigFloat},Base.MPFR.MPFRRoundingMode),imag(destination),imag(a),real(b),mode)
        ccall((:mpfr_mul,Base.MPFR.libmpfr),Cint,(Ref{BigFloat},Ref{BigFloat},Ref{BigFloat},Base.MPFR.MPFRRoundingMode),real(destination),real(a),real(b),mode)
        return destination
    elseif iszero(imag(a))
        ccall((:mpfr_mul,Base.MPFR.libmpfr),Cint,(Ref{BigFloat},Ref{BigFloat},Ref{BigFloat},Base.MPFR.MPFRRoundingMode),imag(destination),real(a),imag(b),mode)
        ccall((:mpfr_mul,Base.MPFR.libmpfr),Cint,(Ref{BigFloat},Ref{BigFloat},Ref{BigFloat},Base.MPFR.MPFRRoundingMode),real(destination),real(a),real(b),mode)
        return destination
    end
    ccall((:mpfr_mul,Base.MPFR.libmpfr),Cint,(Ref{BigFloat},Ref{BigFloat},Ref{BigFloat},Base.MPFR.MPFRRoundingMode),ac,real(a),real(b),mode)
    ccall((:mpfr_mul,Base.MPFR.libmpfr),Cint,(Ref{BigFloat},Ref{BigFloat},Ref{BigFloat},Base.MPFR.MPFRRoundingMode),bd,imag(a),imag(b),mode)
    ccall((:mpfr_mul,Base.MPFR.libmpfr),Cint,(Ref{BigFloat},Ref{BigFloat},Ref{BigFloat},Base.MPFR.MPFRRoundingMode),ad,real(a),imag(b),mode)
    ccall((:mpfr_mul,Base.MPFR.libmpfr),Cint,(Ref{BigFloat},Ref{BigFloat},Ref{BigFloat},Base.MPFR.MPFRRoundingMode),bc,imag(a),real(b),mode)
    ccall((:mpfr_sub,Base.MPFR.libmpfr),Cint,(Ref{BigFloat},Ref{BigFloat},Ref{BigFloat},Base.MPFR.MPFRRoundingMode),real(destination),ac,bd,mode)
    ccall((:mpfr_add,Base.MPFR.libmpfr),Cint,(Ref{BigFloat},Ref{BigFloat},Ref{BigFloat},Base.MPFR.MPFRRoundingMode),imag(destination),ad,bc,mode)
    destination
end

"""Add into an owned complex accumulator without mutating its input."""
function scratch_add!(destination::Complex{BigFloat},value::Complex{BigFloat})
    mode=Base.MPFR.rounding_raw(BigFloat)
    ccall((:mpfr_add,Base.MPFR.libmpfr),Cint,(Ref{BigFloat},Ref{BigFloat},Ref{BigFloat},Base.MPFR.MPFRRoundingMode),real(destination),real(destination),real(value),mode)
    ccall((:mpfr_add,Base.MPFR.libmpfr),Cint,(Ref{BigFloat},Ref{BigFloat},Ref{BigFloat},Base.MPFR.MPFRRoundingMode),imag(destination),imag(destination),imag(value),mode)
    destination
end

"""Accumulate a product; the BigFloat method mutates only its owned accumulator/workspace."""
sum_product(value,a,b,::Nothing)=value+a*b
function sum_product(value::Complex{BigFloat},a::Complex{BigFloat},b::Complex{BigFloat},work::BigComplexScratch)
    scratch_product!(work.product,a,b,work)
    scratch_add!(value,work.product)
end
