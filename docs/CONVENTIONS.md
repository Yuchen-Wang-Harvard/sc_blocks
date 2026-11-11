# Mathematical conventions

These are the mathematical conventions implemented by the package. Primary
normalizations are fixed; the returned blocks require no convention conversion.

We store twice every descendant level and twice each G-mode index as integers.
Branching labels can be stored as 4n, so NS and Ramond quarter shifts are exact.
Intrinsic NS parity is explicit. All signs in Ward identities and sewing use
total state parity; form selection uses relative descendant/R-ground parity.

BPZ is bilinear, not Hermitian: numerical contractions use transpose, never
complex adjoint. L_n^dagger=L_-n, G_r^dagger=iG_-r, psi_r^dagger=-i psi_-r;
(AB)^dagger=(-1)^(|A||B|)B^dagger A^dagger. The NS ground norm is one for
either parity. R grounds satisfy G_0 w^±=i beta exp(∓i pi/4) w^∓ and have
Gram diag(1,-1). Auxiliary R grounds have Gram diag(1,-i).

The NSRR master identity has origin sign `-(-1)^(epsilon_1+epsilon_2+n)`,
with no form-parity f. Its ground forms are
rho_0(++ )=1, rho_0(--)=eta, rho_1(+-)=1, rho_1(-+)=i eta. Auxiliary-first
tensor states use the full graded BPZ pairing

```
<a tensor x, b tensor y> = (-1)^(|a||x|) <a,b>_F <x,y>_SCA.
```

The ordered enlarged vertex is

```
rho_hat(a1 tensor x1, a2 tensor x2, a3 tensor x3)
  = (-1)^(|a1||x1| + |x2||a3|) rho_F(a1,a2,a3) rho(x1,x2,x3).
```

Sewing uses the inverse full tensor Gram without an additional edge-parity
factor. A physical odd NS ground is never rescaled by i.

All Ramond edges form cycles in a closed trivalent graph with NSRR vertices.
A cycle requires an insertion precisely when the product of vertex eta's is
-1. Construct both inserted blocks on the split graph and extract equal final
left/right powers before recovering physical coefficients. Divide only within
the allowed parity sector; a zero auxiliary component is not invertible.

The insertion is a bra/ket operator kernel. Its dual-input ordering differs
from an ordinary ordered three-point tensor at (right,external,left) by one
interchange of equal-parity legs. `split_sewing_sign` includes this sign in
both auxiliary and enlarged assembly. This convention is fixed by requiring
an identity kernel between two inverse Grams to reproduce one original
propagator. It leaves the local splitting matrix element unchanged and is
never applied as a conversion of a recovered block.

Antiholomorphic labels obey tilde a=a+sum_i(epsilon(phi_i)+epsilon(tilde phi_i))
and tilde f=f+epsilon(phi_NS)+epsilon(tilde phi_NS), modulo two. The residual
nonchiral sewing sign remains sum_(v<w) p_v p_w. There is no extra eta flip.

The current Appendix B signed epsilon_i is an ordinary integer difference.
Choose odd-primary norm roots to be (-i)^(2n mod2) times even-primary roots.
In Ramond reflection, reflect beta along with P, then negate components
ending in w-minus. The Liouville specialization is P=i p, beta=i p/sqrt(2),
with C^(1,eta)=i C^(0,-eta). Schottky products over both inverse conjugacy
classes carry the compensating factor 1/2 in their exponent.

## Numerical recovery settings

The default normalized auxiliary-sector residual limit is **1e-8**. It is independent of the
working precision and does not change the local Ward/pivot tolerances. Actual
sector and forward-convolution residuals are always recorded. The optional
`recovery_tolerance` input explicitly overrides this numerical limit.

`recovery_policy="record"` is a diagnostic mode: it retains coefficients and
residuals even after a failed sector check, marks the output as unvalidated,
and is rejected by the benchmark harness. Ordinary production uses
`recovery_policy="error"`. A failed numerical call writes a `.failed.toml`
record with its input, source hashes, elapsed time and available stage data.
