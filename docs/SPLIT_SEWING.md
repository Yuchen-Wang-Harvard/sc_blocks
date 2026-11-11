# Split-edge matrix elements and parity ordering

The insertion is the matrix element `<right|Theta psi|left>`. It is evaluated
with the full current tensor BPZ metric. The local splitting recurrence returns
exactly that raw matrix element; it contains no normalization-root adjustment.

An operator kernel is not an ordinary ordered trilinear form. In a tensor
contraction its input is a dual leg. The graph stores the insertion axes as
`(right, external, left)`, so identifying the kernel with those axes reverses
the two equal-parity internal legs. That contributes `(-1)^parity` at each
insertion. `split_sewing_sign` applies this local ordering rule in both the
auxiliary and enlarged graph assembly, before any division.

There is an independent way to fix the sign. Replace the inserted operator by
the identity. Its matrix element is the BPZ Gram matrix, so two split inverse
Grams contract to the original inverse Gram. The original parity component
must therefore be unchanged. For the endpoint order in `split_graph`, the
generic ordinary-vertex permutation has

```
K_split(lift(mask)) = K_original(mask) + sum_marked parity(edge)  (mod 2).
```

Treating the matrix kernel as another ordinary trilinear vertex would thus
insert a parity operator even when the inserted operator is the identity.
The dual-leg interchange cancels that spurious operator. This is why the code
uses `K_split + sum_marked parity(edge)` for a matrix-element insertion. It
does not alter BPZ norms, the local Ward identities, or recovered block
coefficients by a final conversion.

The generic ordinary-vertex `sewing_sign` is used for ordinary vertices;
`split_sewing_sign` handles the insertion-kernel ordering described above.

# Unequal intermediate powers

For a requested original coefficient domain D, the split calculation uses
the explicit downward-closed domain whose membership criterion is

```
(unsplit powers, max(left_1,right_1), ...) in D.
```

For each branching primary, subtract its independent left/right leading
powers before dividing the remaining budgets by two for the integer Virasoro
levels. Both CCY factors are evaluated on this full domain. Their product is
formed before retaining equal final split powers and equal split parity bits.

The insertion sphere's zero/infinity maps compose to the original sewing map
with `q = q_left*q_right`. Its Schottky vacuum therefore contains equal powers
only. Exact rational comparisons on one- and two-split graphs check this map
composition in the implemented series arithmetic. The common vacuum square
can consequently be restored once after diagonal extraction and auxiliary
division, while unequal intermediate CCY levels must still be retained.
