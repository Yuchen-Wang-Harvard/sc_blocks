# Code file and function reference

This reference describes the numerical library, command-line drivers,
independent checks and their input configurations.

## Package and project files

- `Project.toml`: package identity and explicit Julia/standard-library dependencies.
- `README.md`: package layout, standalone invocation and check commands.
- `docs/CONVENTIONS.md`: mathematical signs, parities and normalization.
- `docs/SPLIT_SEWING.md`: insertion-kernel ordering, identity check and unequal-grade domains.
- `.gitignore`: excludes generated outputs and local build/cache files.

## Source files

### `src/SuperconformalBlocks.jl`

Package entry point, standard-library imports, convention identifier and public
exports. It includes the numerical files in dependency order. It defines no
numerical functions itself.

### `src/combinatorics.jl`

Momentum-independent combinatorics and reusable downward-closed level domains.

| Function/type | Purpose |
| --- | --- |
| `Combinatorics` | Own partition, fermionic-partition and exact-binomial caches. |
| `partitions` | Enumerate descending unrestricted integer partitions. |
| `fermion_partitions` | Enumerate strict partitions into odd NS or even R doubled modes. |
| `exact_binomial` | Evaluate generalized integer/half-integer binomials as exact big rationals. |
| `CoefficientDomain` | Enumerate a box intersected with a total-degree bound, or validate an explicit downward-closed domain; cache index lookup. |
| `CoefficientDomain`'s `visit` | Recursively enumerate the allowed exponent tuples. |
| `total_domain`, `box_domain` | Construct pure total or independent-edge cutoffs. |
| `contains` | Test edge/total limits for regular domains, or exact membership for explicit domains. |
| `parity_sign` | Evaluate `(-1)^n` safely for signed integers. |
| `rising` | Rising factorial without machine-precision constants. |

### `src/graphs.jl`

Ordered graph incidence, Koszul signs, split-edge construction and loop selection.

| Function/type | Purpose |
| --- | --- |
| `PlumbingGraph` | Validate two endpoints per internal edge and retain their order. |
| Constructor-local `edgeof` | Map a half-edge ID to its internal edge during inversion counting. |
| `theta_graph`, `glasses_graph`, `tetrahedron_graph` | Construct the standard ordered closed graphs. |
| `sewing_sign` | Evaluate the edge-pair-to-vertex Koszul permutation. |
| `split_graph` | Add insertion vertices while retaining separate left/right edges and endpoint orders. |
| `ramond_cycles` | Find Ramond cycle components and validate NSRR ordering. |
| `insertion_edges` | Select one edge per cycle with product of vertex eta equal to minus one. |
| `mixed_convolution_sign` | Include full tensor-Gram and ordered-vertex signs for ordinary closed-graph convolution. Explicitly rejects external insertions. |

### `src/sca.jl`

Independent physical PBW algebra and current linear graded BPZ metric. This is
the direct check engine and supplies local boundary data; high-level production
blocks must use the branching/CCY algorithm.

| Function/type | Purpose |
| --- | --- |
| `Mode`, `L`, `G` | Represent integer L indices and doubled G indices. |
| `PBWState`, `ground_state` | Store canonical L/G words and the Ramond ground bit. |
| `level2` | Return doubled descendant grade. |
| `relative_parity` | Count lowering-operator parity plus the R ground bit. |
| `pbw_word` | Expand a canonical state into its ordered lowering modes. |
| `SCAModule` | Own momentum, intrinsic NS parity, basis/action/Gram caches. Inputs are converted directly to the requested precision. |
| `state_parity` | Add intrinsic primary parity to relative parity. |
| `basis`, `basis_index` | Cache states and inverse state lookup at each grade. |
| `sparse_add!`, `prune_zeros!` | Accumulate sparse vectors and remove exact zeros only. |
| `commutator` | SCA supercommutators including central terms. |
| `normal_order` | Reduce mode words on the primary, caching intermediate words. |
| `act`, `act_vector`, `action_matrix` | Apply modes to basis states, sparse vectors or complete grades. |
| `pairing` | Linear graded BPZ pairing, including odd NS primary signs. |
| `gram`, `gram_factor`, `inverse_gram` | Share one Gram matrix, LU factorization and inverse per grade. |

### `src/ward.jl`

Independent ordered physical forms. NSRR eta signs are evaluated together;
relative form labels and total state parities remain distinct.

| Function/type | Purpose |
| --- | --- |
| `WardVertex`, `ns_vertex`, `nsrr_vertex` | Construct/validate local physical forms and cache vertices. |
| `form_zero`, `form_add`, `form_scale`, `form_sub` | Small fixed-tuple arithmetic for batched eta values. |
| `mode_value` | Insert a mode; use translation directly for middle/origin L_-1. |
| `ward_sum` | Truncate each binomial contour sum using annihilation above the state grade. |
| `reduce_l` | Remove one L lowering mode by the Virasoro Ward identity. |
| `reduce_g` | Remove one G lowering mode in the NSNSNS or current no-f NSRR identity. |
| `ward_values` | Cache normalized form values in the states' relative-parity channel. |
| `rho` | Select relative form label a/f and, for NSRR, eta. |
| `vertex_tensor` | Cache complete local three-point tensors for repeated contractions. |
| `master_ward_residual` | Test arbitrary master G contour choices independently of recursive pull choices. |
| `virasoro_ward_residual` | Test arbitrary master L contour choices. |
| `residual_pair` | Absolute and max(1,magnitudes)-scaled discrepancy without precision loss. |

### `src/fermion.jl`

Independent auxiliary Clifford algebra and even three-point Ward forms.

| Function/type | Purpose |
| --- | --- |
| `FermionState` | Distinct descending doubled modes and auxiliary ground bit. |
| `level2`, `relative_parity`, `state_parity` | Auxiliary grade and parity. |
| `FermionModule`, `basis` | Own sector, shared combinatorics, basis and action caches. |
| `act_psi` | Creation, annihilation and normalized Ramond zero-mode action. |
| `fermion_norm`, `pairing`, `gram`, `inverse_gram` | Auxiliary BPZ norm and diagonal Gram data; these are distinct from ordered ground vertices. |
| `FermionVertex` | Validate NSNSNS/NSRR slots and own the auxiliary vertex cache. |
| `psi_value`, `psi_sum` | Mode insertions and finite binomial sums. |
| `rho` | Compute the even auxiliary form from its own Ward identity. |
| `master_ward_residual` | Check arbitrary auxiliary master contours. |
| `tensor_pairing` | Full auxiliary-first graded tensor BPZ pairing. |
| `enlarged_rho` | Ordered enlarged vertex with both required local signs. |

### `src/series.jl`

Sparse multivariate arithmetic, used with exact rational Schottky coefficients
and numerical recursion coefficients.

| Function/type | Purpose |
| --- | --- |
| `Series` | Dictionary from exponent tuples to generic scalar coefficients. |
| `series_constant`, `series_add`, `series_scale` | Construct and combine sparse coefficient series. |
| `series_multiply` | Truncated multivariate convolution. |
| `series_order` | Lowest nonzero total degree. |
| `series_shift` | Multiply by a monomial, respecting the original domain. |
| `residual_domain` | Domain remaining after a leading monomial has been removed, optionally dividing exponents by an integer scale. Preserves explicit split-domain shapes. |
| `series_inverse` | Formal inverse with a nonzero constant coefficient. |
| `series_exp` | Formal exponential of a zero-constant logarithm. |
| `series_matrix_multiply` | 2-by-2 Mobius matrix product over truncated series. |

### `src/schottky.jl`

Exact bosonic vacuum factor for the actual ordered graph. Inverse cycles are
identified and assigned exponent -1, equivalent to counting both with -1/2.

| Function/type | Purpose |
| --- | --- |
| `canonical_cycle` | Minimize a cycle over rotations and reversal. |
| `primitive_cycles` | Enumerate cyclically reduced primitive classes with cutoff pruning. |
| Enumeration-local `visit` | Traverse directed arcs, forbidding immediate reversal and pruning by first oscillator degree. |
| `plumbing_maps` | Build exact local-coordinate and edge inversion matrices. |
| `schottky_vacuum` | Expand the exact multiplier product; prune each trace by its remaining domain before exponentiation. |

### `src/numeric_scratch.jl`

Owned work buffers for Julia's bundled MPFR library. These buffers never alias
cached vertices, input coefficients or another returned coefficient.

| Function/type | Purpose |
| --- | --- |
| `BigComplexScratch` | Allocate a private complex product and four real MPFR work buffers at the current precision. |
| `complex_scratch` | Select the MPFR workspace for BigFloat or the ordinary-arithmetic fallback for other types. |
| `scratch_product!` | Multiply into an owned buffer, specializing exact-real factors; safe when either or both inputs are that buffer. |
| `scratch_add!` | Add to an owned complex accumulator using the current BigFloat rounding mode. |
| `sum_product` | Accumulate a product, with MPFR buffer reuse or ordinary arithmetic as appropriate. |

### `src/ccy.jl`

Production Virasoro central-charge recursion. It never uses finite-c Gram blocks
as its seed. Cache ownership is local to a numerical plan/worker.

| Function/type | Purpose |
| --- | --- |
| `CCYPole` | Store central pole, squared coupling, universal residue and center ID. |
| `CCYPlan` | Validate edge/external weights and own numerical caches and counters. |
| `ccy_pole` | Cached Kac pole geometry, with removable x±1 factors canceled before numerical evaluation. |
| `ccy_fusion` | Cached fusion polynomial, pairing factors to avoid extra square-root choices. |
| `ccy_denominator` | Cache inverse pole differences and reject unresolved collisions. |
| `global_norm` | Cache n!(2h)_n factors. |
| `global_vertex_core` | Evaluate the part independent of the middle descendant level, cached across all such levels; use exact integer combinatorial factors and real arithmetic when all weights are exactly real. |
| `global_vertex` | Closed L_-1 three-point factor using reused rising/falling products. |
| `ccy_vertex` | Cache complete local factors independently of spectator edges; identical ordered vertices share one representative. |
| `ccy_normalized_vertex` | Cache a vertex together with its assigned inverse edge norms; each edge is assigned exactly once, including self-loops. |
| `ccy_global` | Weight-dependent graph global seed. |
| `ccy_seed_coefficient` | Attach all residue amplitudes at a target grade, reusing MPFR work buffers for BigFloat and using ordinary arithmetic for other types. |
| `ccy_residue` | Ordered endpoint fusion factors; self-loops currently require slots 2 and 3. |
| `ccy_reduced` | Forward batched recursion, merging accumulated shifts and incoming poles, then attaching global seeds. |
| `virasoro_block` | Restore the exact Schottky factor once after the reduced sum. |

### `src/contractions.jl`

Bilinear tensor contractions shared by direct reference calculations. No
branching, pole, or production coefficient data are used here.

| Function/type | Purpose |
| --- | --- |
| `NetworkTensor` | Validate a dense tensor and its distinct half-edge axis labels. |
| `contract_pair` | Contract all shared labels with ordinary matrix multiplication. |
| `contraction_cost` | Estimate contraction work in log units without overflow. |
| `contract_network` | Greedily eliminate shared indices and close disconnected components. |
| `sew_tensors` | Assemble local forms and inverse Grams with optional parity selections. |
| `pbw_coefficient` | Physical closed-graph PBW coefficient in every edge-parity sector. Validates module agreement across each edge. |
| `evaluate_spins` | Evaluate parity coefficients at specified edge spin signs. |

### `src/ns_branching.jl`

Current Appendix-B closed formulas and double-Virasoro branch parameters.

| Function | Purpose |
| --- | --- |
| `ell` | Exact-index product and reflection rule, evaluated at the chosen numeric precision. |
| `branch_central`, `branch_weight` | Central charges and weights of the two commuting Virasoro factors. |
| `ns_branch_norm` | Current NS chi-vector BPZ norm, with explicit intrinsic primary parity. |
| `ramond_branch_norm` | Current R chi-vector BPZ norm, independently checked against free-field/PBW pairings. |
| `ns_branching` | Normalized all-NS branching coefficient with signed integer epsilon differences. |
| `ns_branching_raw` | Closed unnormalized chi-vector vertex; cancels all normalization square roots analytically and retains reflected labels and intrinsic NS parity. |
| Local `evaluate_ell` | Optionally reuse ell products within one fixed-coupling numerical plan. |

### `src/ns_pipeline.jl`

Complete ordinary all-NS path. All plumbing exponents passed into this file
are doubled levels.

| Function/type | Purpose |
| --- | --- |
| `NSBlockPlan` | Own closed-graph data, intrinsic parities, local branching/auxiliary caches and counters. |
| `ns_local_branch` | Cache a local branching coefficient across spectator edge choices. |
| `ns_branch_weight` | Cache embedded weights by edge, branch label and Virasoro copy. |
| `ns_auxiliary_tensor` | Cache the direct free-fermion tensor at each incident grade triple. |
| `ns_auxiliary_block` | Direct auxiliary Fock sewing with its own inverse Grams and Koszul sign. |
| `recover_unit` | Triangular division when the auxiliary constant term is nonzero; rejects zero divisors. |
| `ns_block` | Sum allowed branch labels with remaining-domain CCY computations, divide by the auxiliary block, then restore the two universal Schottky factors once. |
| `ns_spin_coefficient` | Apply edge spin signs to one coefficient, including intrinsic NS primary parities. |

### `src/freefield.jl`

Native boson/fermion realization in the current chi-vector convention. Sparse
oscillator states use interned integer IDs; actions strip spectator factors so
the same result is reused across auxiliary/physical states. Ordered PBW maps
are internal local basis changes, never final-block phase conversions.

| Function/type | Purpose |
| --- | --- |
| `OscillatorState`, `physical_level2`, `level2` | Store boson, physical-fermion and auxiliary words; compute physical or total doubled grades. |
| `FockMapping` | Cache the physical PBW-to-oscillator matrix and its LU factorization at one grade. |
| `FreeField` | Own momentum/sector data, interned states, actions, reflected chart, primary and descendant caches. |
| `intern!`, `ff_add!` | Intern states and accumulate coefficients by stable integer ID. |
| `boson_action`, `eta_action` | Apply free boson/physical Clifford modes, including the fixed Ramond zero-mode phase. |
| `ff_physical_action`, `ff_physical_apply` | Compute/cache physical L/G oscillator actions and apply them to sparse vectors. |
| `ff_auxiliary_action` | Compute/cache auxiliary Virasoro action without physical spectators. |
| `ff_embedded_action`, `ff_embedded_apply` | Apply either commuting Virasoro algebra, including graded psi-G signs. |
| `chi_apply` | Apply psi-i eta, including the reflected zero mode used in Ramond primaries. |
| `ff_pbw_image` | Embed physical PBW words, sharing suffix computations. |
| `fock_mapping`, `ordered_pbw`, `from_ordered_pbw` | Share one physical map per grade and convert local vectors in either direction. |
| `ff_primary` | Construct raw chi highest vectors; negative branches reflect P and beta and identify the w-minus component. |
| `ff_descendant` | Build both Virasoro descendant words with a shared suffix cache. |
| `ff_pairing`, `ff_norm` | Evaluate full graded tensor BPZ pairings/norms through the independent PBW metric. |

### `src/span_solve.jl`

Sparse-span coefficient extraction with precision verified on every oscillator
row. Machine arithmetic chooses/preconditions a square subsystem; all residuals,
refinement and final coefficients use the requested working precision.

| Function/type | Purpose |
| --- | --- |
| `span_tolerance` | Precision-scaled residual threshold. |
| `SpanSolver` | Normalize columns, select rows by pivoted QR, and cache a machine or high-precision LU. |
| `solve_span` | Iteratively refine at working precision; fall back to precise LU and verify the full rectangular system. Return coefficients and conditioning/residual evidence. |

### `src/branch_actions.jl`

| Function/type | Purpose |
| --- | --- |
| `BranchActionTerm` | One target branch and pair of Virasoro words with its raw coefficient. |
| `BranchActions` | Share positive/reflected free-field charts and physical L_-1, L_0, L_1 decompositions. |
| `partition_pairs` | Enumerate two Virasoro words with fixed combined integer level. |
| `action_descriptors` | Enumerate the exact neighboring-branch support for a physical L action. |
| `direct_branch_action` | Fit that support with the verified sparse-span solver. |
| `branch_action` | Cache reflected actions and reuse Ramond alpha=0 fits for alpha=1 with the independently checked raw-vector factor. |
| `reconstruct_action` | Reconstruct a full oscillator vector for independent action validation. |
| `single_leg_vertex` | Closed Virasoro Ward reduction with descendants on one slot. |

### `src/ramond_branching.jl`

Local physical-L Ward recurrence for raw NSRR branching vertices. Both eta
values share each solve. Direct physical-G Ward evaluation is confined to a
fixed low-label boundary box in production; higher-label direct values appear
only in independent checks.

| Function/type | Purpose |
| --- | --- |
| `BranchWardRow` | Coefficients plus absolute product bounds used to detect unresolved cancellation. |
| `NSRRBranching` | Own or share three edge action plans, local Ward forms, boundary values and recurrence caches. |
| `branch_height`, `branch_boundary` | Dependency ordering and fixed direct-anchor box. |
| `branch_group`, `branch_group_members` | Group the coupled low-label ambiguities into systems of at most four unknowns. |
| `branching_ordered` | Cache raw primaries in ordered tensor PBW coordinates. |
| `direct_branch_values` | Evaluate both eta values by physical and auxiliary Ward forms. |
| `known_branch_value` | Retrieve solved values; allow new direct evaluations only inside the fixed boundary box. |
| `branching_weight` | Cache embedded branch weights. |
| `single_leg_bound`, `branch_term_factor` | Descendant factors and absolute bounds for Ward cancellation tests. |
| `branch_ward_row` | First physical L_1 identity or alternative L_-1/L_0/L_1 identity. |
| Row-local `add_terms` | Accumulate changed-label coefficients and absolute product bounds for one acted-on leg. |
| `ward_pivot_tolerance` | Precision-scaled cancellation threshold. |
| `solve_branch_group` | Solve a small row-scaled system for both eta values by full pivoting; report unresolved rank. |
| `branch_row_residual` | Check every retained Ward equation after solving. |
| `branch_support` | Build the complete dependency closure and coupled boundary groups. |
| `prepare_branching!` | Solve requested labels in dependency order, adding the alternative Ward rows when necessary. |
| `branching_raw` | Select f and eta from the cached raw coefficients without normalization-root choices. |

### `src/parity_series.jl`

| Function/type | Purpose |
| --- | --- |
| `ParitySeries` | Sparse multigrade coefficients with a vector indexed by edge-parity mask. |
| `walsh!` | In-place parity/fixed-spin Walsh transform. |
| `parity_product` | Plumbing convolution with XOR of parity masks. |
| `parity_scalar_product` | Multiply all parity components by a universal scalar series. |
| `recover_allowed_sector` | Divide triangularly in the allowed sector, using the original default 1e-8 residual guard; record actual sector/forward residuals. Explicit diagnostic mode never certifies its output. Caller inserts on odd cycles first. |

### `src/ramond_pipeline.jl`

Closed graphs allow NSNSNS and ordered NSRR vertices. Raw local vertices and
one inverse norm per edge avoid arbitrary square-root gauges.

| Function/type | Purpose |
| --- | --- |
| `RamondBlockPlan` | Validate graph data, detect odd cycles and share edge actions, local vertices, norms, weights and auxiliary tensors. Optional `cache_products=true` retains momentum-local Virasoro products across form/eta choices; it is off for single-block benchmarks. |
| `record_ccy_stages!` | Aggregate pole and seed timings as substages already included in the overall CCY time. |
| `set_vertex_labels!` | Change form/eta choices while retaining all momentum-local caches; recompute the odd-cycle insertion choices. |
| `prepare_product_cache!` | Clear optional products when the exact requested multigrade domain changes; retain a separate copy of the domain indices. Products are keyed by insertion topology and branch labels. |
| `record_stage!` | Accumulate wall time for branching, CCY, products, auxiliary recovery and vacuum restoration without altering arithmetic. |
| `branch_level2` | Exact NS/R primary grade from the integer label 4n. |
| `ramond_branch_labels` | Enumerate labels whose leading powers fit the requested domain. |
| `ramond_plan_norm`, `ramond_plan_weight` | Cache raw BPZ norms and embedded Virasoro weights. |
| `ramond_auxiliary_tensor`, `ramond_auxiliary_block` | Cache local auxiliary vertices and sew the ordinary auxiliary parity series. |
| `ramond_branch_coefficients` | Sew raw branching vertices and norms in all parity sectors. |
| `graph_branch_value` | Dispatch to a cached raw closed NS vertex or the local NSRR recurrence; caches remain reusable across form/eta choices. |
| `ramond_block` | Dispatch odd cycles to split assembly; otherwise sum two CCY factors per label, reuse them across Ramond parities, recover physical coefficients and restore the vacuum once. |

### `src/splitting.jl`

| Function/type | Purpose |
| --- | --- |
| `theta_psi_action`, `theta_psi_apply` | Apply the even Theta-psi operator on auxiliary states or full oscillator vectors. |
| `psi_weights` | External weights in the two embedded Virasoro factors. |
| `SplittingPlan` | Share the edge's physical-L actions and own raw splitting/factor caches. |
| `direct_splitting` | Independent raw highest-vector matrix element, used only for ground anchors in production. |
| `splitting_factor` | Cache one-leg Virasoro factors at the insertion. |
| `splitting_zero_pivot` | Alternative physical L_0 identity when the L_1 pivot is unresolved. |
| `splitting_raw` | Physical-L recurrence with exact neighboring-branch support, current BPZ norms and ground anchors. |

### `src/split_pipeline.jl`

| Function/type | Purpose |
| --- | --- |
| `SplitBlockPlan` | Expanded graph, original-edge ownership, shared splitting plans and insertion/domain caches. |
| `split_domain` and local `visit` | Enumerate the downward closure of diagonal targets by requiring each pair's maximum to fit the original domain. |
| `split_mask` | Lift one original parity to both halves of every marked edge. |
| `split_sewing_sign` | Include dual-input ordering for inserted bra/ket kernels; the identity kernel must preserve the original propagator. Applied before recovery in both auxiliary and enlarged sewing. |
| `split_diagonal` | Extract equal final split powers and parities after expanded assembly. |
| `split_branch_labels` | Enumerate independent left/right branch labels with the exact splitting support. |
| `split_branch_coefficients` | Multiply raw ordinary vertices, splitting matrix elements, full inverse norms and the operator-kernel sewing sign. |
| `auxiliary_insertion_tensor` | Direct auxiliary bra/ket matrix element at independent left/right grades. |
| `split_auxiliary_block` | Sew the auxiliary block on the expanded graph; production requests only equal final split powers, while the diagnostic route evaluates the full domain. |
| `split_diagonal_product` | Convolve all unequal intermediate powers for each required shifted diagonal output coefficient. |
| `split_transpose_supported` | Identify the ordered theta graph symmetry that permits exchanging the inserted edge halves. |
| `split_enlarged_diagonal` | Production diagonal assembly with shared local branching requests and reuse of the two-factor product under the theta transposition. Raw branching coefficients are still evaluated separately for both orientations. |
| `split_enlarged_reduced` | Diagnostic complete expanded product, retaining unequal output coefficients for comparison with the optimized diagonal route. |
| `inserted_ramond_block` | Extract the diagonal, divide in the allowed sector, verify the forward convolution and restore the Schottky square. |

### `src/ns_recursion.jl`

Independent native-BPZ all-NS SCA central-charge recursion. Its superglobal
vertices and null transport include the current phases locally; it does not
convert a completed block. The current reference supports closed graphs
without self-loops, including the required theta and tetrahedron checks.

| Function/type | Purpose |
| --- | --- |
| `NSPole` | NS central pole, regularized residue, squared coupling and center ID. |
| `NSRecursionPlan` | Own physical weights, intrinsic parities, current sewing phases and numerical caches. |
| `ns_mode_phase` | Exact fourth-root exponent supplied by current local vertices and inverse BPZ norms. |
| `ns_c_pole` | Stable NS Kac geometry and residue with removable x±1 factors canceled algebraically. |
| `ns_c_fusion` | Paired NS fusion polynomial in either relative form channel. |
| `ns_c_denominator` | Cache inverse central-charge differences, rejecting unresolved collisions. |
| `ns_superglobal_vertex` | Independently reduce osp(1\|2) descendants with current ordered-primary normalization. |
| `ns_c_vertex`, `ns_c_global` | Cache local superglobal vertices/norms and sew the graph seed. |
| `ns_c_residue` | Include native local BPZ, vertex and Koszul phases in each null-state transition. |
| `ns_c_reduced` and local `maskof` | Batch by output parity; merge null shifts and incoming poles before attaching superglobal seeds. |

### `src/ns_schottky.jl`

| Function | Purpose |
| --- | --- |
| `fundamental_cycles` | Build a spanning forest and one directed cycle per non-tree edge. |
| `ns_spin_phases` | Fix edge spin lifts from native supercurrent two-point sewing at the cusp. |
| `ns_schottky` | Exact complex-rational bosonic/fermionic primitive product; the first fermionic oscillator is 3/2. |
| `ns_c_block` | Restore this seed by ordinary convolution in the current BPZ convention. Accepts a precomputed momentum-independent seed. |

## Command-line driver

### `bin/blocks.jl`

Standalone coefficient CLI for all-NS and mixed NSNSNS/NSRR graphs, including automatic
odd-cycle insertions. No dependencies on git, parent folders or external
numerical engines. The `NSRR` driver route takes an explicit `ramond` edge list.

| Function | Purpose |
| --- | --- |
| `parse_real` | Parse integer, quoted decimal or quoted rational parameters at working precision. |
| `parse_complex` | Parse a real parameter or a two-component complex pair. |
| `production_hashes` | Hash production sources and drivers for portable provenance. |
| `run_blocks` | Read configuration, compute/timestamp coefficients and atomically save TOML; persist failed-call diagnostics and distinguish diagnostic recovery from validated output. |

`examples/ns_theta.toml` demonstrates a total-level-3 all-NS theta computation
at 136 working bits. It is not a matched paper-scale runtime configuration.
`examples/nsrr_theta.toml` and `examples/nsrr_theta_inserted.toml` demonstrate
ordinary and inserted NSRR theta computations at the same cutoff. Ramond output
stores all parity components: array index `1+mask`, edge e at bit e-1. This
includes intrinsic NS parity; a fixed-spin value is the corresponding Walsh sum.
`examples/benchmark_theta_ordinary.toml` and
`examples/benchmark_theta_inserted.toml` specify independent
edge-level-10, 136-working-bit benchmarks at b=7/5 and
P=(11/23,13/29,17/31).
The corresponding `*_level5.toml` files specify the smaller matched
performance-comparison inputs at the same precision and momentum point.

## Checks and tests

### `checks/virasoro_pbw.jl`

Independent finite-c Virasoro reference module. It shares elementary L-mode
algebra and generic contractions, but uses no CCY, fusion or Schottky functions.

| Function/type | Purpose |
| --- | --- |
| `Reference` | Construct graph modules and local reference caches. |
| `vir_basis` | Pure L-descendant partitions at integer levels. |
| `vir_inverse` | Direct restricted Virasoro Gram inverse. |
| `vir_tensor` | Descendant tensor from L Ward reductions. |
| `coefficient` | Direct finite-c graph sewing at one multidegree. |

### `checks/run_foundations.jl`

Reproducible driver for the implemented subset. It does not claim full release
validation. `source_hashes` hashes every Julia source/CLI/check/test file;
`run_foundations` runs the suite, records hardware/runtime/precision and actual
timing/allocation data, and saves success or failure in TOML. A supplied output
path permits retaining each run separately; existing evidence is never overwritten.

### `checks/tensor_pbw.jl`

Independent tensor-PBW realization of both Virasoro embeddings. It does not
use the free-field oscillator action or its state-selection rules.

| Function/type | Purpose |
| --- | --- |
| `Reference` | Own physical SCA/auxiliary modules and independent action/basis caches. |
| `tensor_basis` | Enumerate full auxiliary-first tensor PBW bases. |
| `auxiliary_virasoro` | Evaluate auxiliary L modes using the complete bounded Clifford sum. |
| `embedded_action`, `embedded_apply` | Apply either embedding to states or sparse vectors using physical PBW actions and explicit graded signs. |

### `checks/ns_raw_vertex.jl`

Independent raw all-NS local reference. `Reference` owns three free-field
modules, physical-G/auxiliary Ward forms and ordered-primary caches.
`coefficient` contracts raw chi vectors directly, without calling the closed
branching formula or choosing normalization roots.

### `checks/run_pbw_matrix.jl`

Checkpointing full-cutoff comparison driver. Families are `ns_theta`,
`ns_glasses`, `nsrr_theta`, `mixed_glasses`, `double_r_glasses` and
`tetrahedron`. Cutoffs are per edge except the tetrahedron's total level.
Optional arguments specify the intrinsic NS primary-parity mask and working
bits. Each source/reference timing includes its measured compilation time.

| Function | Purpose |
| --- | --- |
| `measurement_fields` | Retain timing/allocation/compilation metadata without numerical coefficient data. |
| `save_matrix_report` | Atomically checkpoint complete, partial or failed matrices. |
| `pbw_matrix_setup` | Construct production and separate physical-PBW modules from fixed input parameters; enumerate every allowed form/eta choice. |
| `compare_pbw_matrix` | Independently sew every PBW coefficient and compare both parity and all fixed-spin evaluations. |
| `run_pbw_matrix` | Run the requested family/cutoff, reuse local caches across labels, time both algorithms separately and checkpoint each case. |

### `checks/benchmark_blocks.jl`

`benchmark_blocks` runs each input repeatedly with a fresh numerical plan.
The first run includes cold compilation; later runs reuse compiled Julia
methods only. Each report records its source hashes and both numerical and
complete-driver timings. Existing benchmark outputs are never overwritten.
Diagnostic recovery results are rejected as unvalidated benchmarks.

### `checks/compare_saved_blocks.jl`

`compare_saved_blocks` checks identical inputs/conventions/domains, reads every
saved complex parity coefficient at 256 bits, and records the largest scaled
discrepancy with file hashes. It rejects nonfinite results and discrepancies
above its stated tolerance. This regression check supplements independent PBW
verification; it does not establish correctness independently. An explicit
keyword permits a precision difference while still requiring all other inputs
to match.

`examples/benchmark_theta_inserted_level5_192.toml` repeats the 136-bit timing
point at 192 bits to diagnose cancellation and summation-order sensitivity.

### `test/runtests.jl`

Standalone test entry point and includes; no function definitions.

### `test/test_local.jl`

Tests exact combinatorics, graph signs, physical SCA, BPZ and master Ward
identities. `vector_error` compares sparse states; `relation_error` checks
commutators by independent mode applications; `bpz_error` applies the mode on
each side of the pairing independently.

### `test/test_fermion.jl`

Checks Clifford anticommutators, auxiliary BPZ, and independent master psi
identities. It defines no reusable functions.

### `test/test_ccy.jl`

Checks exact series algebra, Schottky cutoff consistency, and CCY against the
independent finite-c reference on five graph types. It defines no functions.

### `test/test_sewing.jl`

Checks ordinary convolution signs, physical parity sectors, tensor ground norm,
and the corrected glasses example. `explicit_theta_coefficient` evaluates an
unoptimized six-index sum independent of the generic contraction order.

### `test/test_ns_pipeline.jl`

Compares the closed-branch all-NS production path against independent physical
PBW sewing for theta/glasses through total level 3, both form labels and three
intrinsic-parity assignments. It defines no functions.

### Additional numerical tests

- `test/test_numeric_scratch.jl`: precision and input/output aliasing checks at
  96, 136 and 192 bits, plus the ordinary-arithmetic fallback. No functions.
- `test/test_freefield.jl`: oscillator/PBW intertwiners, both embedded algebras,
  reflected highest vectors, weights and norms. `scaled_vector_error` compares
  sparse vectors on the union of their supports.
- `test/test_branch_actions.jl`: full-vector L actions, independent fits of both
  Ramond parities, and closed one-leg factors versus Ward evaluation. No functions.
- `test/test_ramond_branching.jl`: higher-label recurrence versus direct physical-G
  Ward forms, both NS parities and eta signs, including alternative Ward pivots.
  No functions.
- `test/test_ramond_pipeline.jl`: ordinary theta/glasses versus physical PBW through
  total level 3 and forward-convolution checks. No functions.
- `test/test_splitting.jl`: splitting recurrence versus direct matrix elements,
  ground anchors, and unequal-grade Virasoro factorization. `splitting_basis`
  enumerates both Virasoro descendant words at a fixed doubled Ramond grade.
- `test/test_split_pipeline.jl`: explicit-domain completeness, exact split Schottky
  composition, identity-kernel ordering and inserted theta/glasses PBW comparisons.
  `check_identity_split` contracts actual descendant Grams with an identity
  kernel and includes a negative control with the incorrect ordinary-vertex order.
- `test/test_reuse.jl`: compare diagonal-target and transposed-reuse assembly
  against the complete expanded product on box and total domains; compare
  auxiliary diagonal restriction with full expanded auxiliary sewing. No functions.
- `test/test_ns_raw.jl`: raw closed NS formula versus independent local Ward
  values for all intrinsic parities, reflected labels and real/imaginary momenta.
  No functions.
- `test/test_mixed_pipeline.jl`: mixed glasses and all 64 allowed tetrahedron
  form/eta choices at preliminary cutoffs. `compare_mixed_block` compares all
  parity coefficients and all Walsh-transformed fixed-spin values against PBW.
- `test/test_precision.jl`: deficient machine-span fallback, rejected missing or
  inconsistent rows, and independent 96/136/192-bit ordinary/inserted blocks.
  No reusable functions.
- `test/test_ns_recursion.jl`: native superglobal vertices versus physical Ward
  evaluation and independent SCA c-recursion versus full physical PBW. Current
  graph cutoffs are theta total level 5 and tetrahedron total level 3, with even
  and odd intrinsic primaries. No reusable functions.

### `src/quadrature.jl`

| Function | Purpose |
| --- | --- |
| `gauss_legendre` | Compute Gaussian nodes/weights and refine machine seeds at the requested precision. |
| `laguerre_value` | Evaluate a generalized Laguerre polynomial and its derivative by recurrence. |
| `gauss_laguerre` | Gaussian quadrature for alpha=±1/2, with precision-refined roots and analytic normalization. |
| `tanh_sinh_integral` | Independent double-exponential integration; report mesh convergence, evaluations and the observed change. |

### `src/liouville.jl`

Current Appendix C constants. These functions use the displayed normalization
directly, including the odd NS factor 2i and the NSRR odd-label exchange.

| Function/type | Purpose |
| --- | --- |
| `bernoulli_even` | Generate exact even Bernoulli numbers. |
| `LogGammaPlan` | Prepare precision-local Stirling coefficients and the recurrence shift. |
| `loggamma_complex` | Evaluate complex log-Gamma by recurrence and Stirling. |
| `log_sinhc` | Stable logarithm of sinh(z)/z, with a small-argument series. |
| `upsilon_integrand` | Evaluate the defining integral with controlled small-t cancellation. |
| `LiouvillePlan` | Own quadrature nodes, coupling and special-function/leg caches. |
| `log_upsilon` | Defining Gaussian quadrature plus analytic shifts/reflection into its convergence strip. |
| `log_sector_upsilon` | Cache the NS/R combinations in Appendix C. |
| `liouville_leg` | Positive conjugate-pair external-leg normalization in logarithmic form. |
| `liouville_ns` | Return C^(0), C^(1), including the exact factor 2i. |
| `liouville_nsrr` | Return the four constants with the relative half and C^(1,eta)=i C^(0,-eta). |
| `upsilon_reference_integrand` | Direct defining integrand, with a separate endpoint Taylor evaluation. |
| `reference_log_upsilon` | Independent tanh-sinh integration on explicit panels with a stated tail cutoff. |

### `src/theta_geometry.jl`

Fresh charged-current geometry and exact homology bookkeeping; no archived
periods, free-frame factors or fitted marking matrices are read.

| Function/type | Purpose |
| --- | --- |
| `ThetaGeometryPlan` | Prepare ordered free-current covariances and charge sources. |
| `minimum_eigenvalue2` | Smallest eigenvalue of a real symmetric two-by-two matrix. |
| `theta_geometry` | Gaussian sewing determinant, period matrix, symmetry residual and positivity diagnostic. |
| `theta_characteristic` | Convert the paper's fixed tube signs to binary characteristics. |
| `theta_constant` | Lattice sum with an explicit Gaussian bound on omitted terms. |
| `free_partition` | Free boson plus Majorana partition function from the current period and determinant. |
| `is_symplectic` | Verify an integer genus-two marking exactly. |
| `transform_period` | Apply the fractional-linear period transformation. |
| `transport_characteristic` | Transport a spin characteristic using the same integer marking. |
| `find_symplectic_map` | Search bounded lower rows, reconstruct upper rows, and test the resulting exact symplectic map. |
| `refine_theta_surface` and local `components` | Refine printed target plumbing seeds by holomorphic Newton iteration at fixed source period. |
| `numeric_schottky_logP` | Independent numerical primitive-cycle determinant product, with determinants tracked separately to avoid subtraction loss. |

### Additional analytic and diagnostic tests

- `test/test_recovery_policy.jl`: genuine lost-sector rejection, explicit
  diagnostic non-certification, precision-specific overrides, and the original
  default 1e-8 threshold. No functions.
- `test/test_analytic_components.jl`: Gaussian polynomial moments, independent
  Gamma/Upsilon checks, determinant/primitive-product comparison, fresh period
  matching and exact spin transport. No functions; included in the full suite.
- `test/test_ns_vacuum_cycles.jl`: independent physical PBW checks of each
  tetrahedron triangle's first fermionic vacuum term and its first even-level
  correction. No functions; included in the full suite.
- `examples/diagnostic_theta_ordinary_level10.toml`: explicitly unvalidated
  recovery-residual diagnostic; it is not an accepted benchmark input.
- `examples/cross_channel.toml`: the ten rounded Appendix C seed pairs and
  spins, with explicit integration and cutoff settings; no saved periods,
  block coefficients or comparison answers.

### `checks/run_ns_recursion_matrix.jl`

| Function | Purpose |
| --- | --- |
| `save_ns_matrix` | Atomically checkpoint a partial, failed or completed all-NS comparison. |
| `run_ns_recursion_matrix` | Compute an independent SCA reference once and compare every vertex-form sector from double Virasoro at the explicit theta per-edge or tetrahedron total cutoff, with real or physical imaginary momenta. Reports coefficient comparisons and the single-parity argument certifying all spin lifts separately. |

### `checks/run_cross_channel.jl`

Fresh Appendix C integration driver. Its source channel uses the production
algorithm; its target may use independent NS recursion or double Virasoro.
The input specifies separate source/target grid orders, working precision,
Laguerre weights and reference Gaussian scales. `required_successive_steps`
controls quadrature convergence; `ratio_tolerance` controls the channel-ratio
comparison. The driver records these two criteria separately.

| Function | Purpose |
| --- | --- |
| `cross_encode` | Encode complex scalars and collections as portable real/imaginary pairs. |
| `save_cross_report` | Atomically preserve partial, failed and completed integration reports. |
| `prepare_cross_surfaces` | Recompute/refine geometry, compare oscillator cutoffs, transport both spins and compute free-frame ratios. |
| `cross_evaluation_plan` | Cache plumbing monomials and spin characters at block precision. |
| `evaluate_cross_series` | Separate NS forms, sum Ramond parity components before squaring, and evaluate multiple total cutoffs. |
| `cross_momentum_grid` | Build fixed or surface-derived generalized Laguerre grids, including original uniform alpha=-1/2, with the actual dp/pi measure and residual-Gaussian compensation. |
| `cross_order_schedule` | Validate either legacy common orders or separate strictly increasing source/target grids. |
| `cross_stability_steps` | Count successive stable changes and reset after a failed or unavailable comparison. |
| `cross_integral_records` | Serialize integral values with explicit surface, cutoff and spin labels. |
| `integrate_cross_channel` | Compute fresh momentum-local blocks/constants, reuse all surface/spin evaluations, and record separate component times. |
| `read_cross_piece` | Reconstruct and validate every labelled integral entry from a freshly completed worker. |
| `dispatch_cross_channel` | Partition the grid among isolated Julia processes using explicit local job/result files, preserve worker checkpoints and reduce in deterministic chunk order. |
| `run_cross_channel` and local `progress` | Run all requested orders, checkpoint progress, test quadrature/ratio tolerances and retain an explicit incomplete release-validation flag. |

`test/test_cross_driver.jl` separately checks form separation, coherent spin
summation, total cutoffs and continuum-measure normalization. It defines no
additional reusable functions and includes the executable integration driver.

### `checks/compare_cross_channel.jl`

| Function | Purpose |
| --- | --- |
| `cross_comparison_values` | Read integral/ratio values keyed by surface, cutoff and spin. |
| `compare_cross_channel` | Check compatible completed runs under different precision, quadrature, reference Gaussian or target-algorithm settings; hash both inputs and preserve failed convergence. |

`checks/run_cross_smoke.jl` defines `run_cross_smoke`: a separately timed driver
validation of the small-grid serial/parallel reductions on all geometries and
a fresh independent target-method comparison. Its smaller cutoffs are explicit;
it does not certify the paper's integration check.

`checks/cross_worker.jl` is the local-process worker executable. Its
`worker_cross_surfaces` decodes fresh geometry/spin inputs; `run_cross_worker`
and its local `worker_progress` compute one disjoint momentum chunk and retain
partial or failed evidence. Workers use atomic TOML files, without sockets or
shared mutable numerical caches.

`checks/run_liouville_precision.jl` defines `run_liouville_precision`: compare
the Float64/order-512 constants with 192-bit/order-1024 constants at four
momentum triples, then independently integrate oscillatory Upsilon arguments
by tanh-sinh. It preserves timings, actual errors, settings and failures.

### NS reference memory and arithmetic checks

`ns_superglobal_core` in `src/ns_recursion.jl` evaluates the independent
osp(1|2) descendant polynomial before its native local phase.
`ns_superglobal_vertex` uses real arithmetic only for exactly real weights and
retains the general complex path. `ns_c_inverse_norm` caches each inverse BPZ
norm, avoiding repeated division in seed sums. `NSRecursionPlan` enables the
original compact-theta strategy: evaluate one of the two identical ordered
vertices, square it bilinearly, and do not retain vertex values that have no
later reuse. Setting `compact_theta=false` retains the complete cache path.
`ns_c_reduced` and `ns_c_block` accept a progress callback; the matrix driver's
local `reference_progress` checkpoints each parity stage without changing the
arithmetic.

`test/test_ns_reference_reuse.jl` compares compact and full-cache paths on real
and complex weights with both intrinsic NS parities, checks equal path counts,
and independently checks complex local vertices against physical Ward forms.
It defines no additional reusable functions.

### Vertex-label product reuse

`test/test_vertex_label_reuse.jl` defines `check_cached_vertex_labels`, which
compares cached-label physical outputs with independent Ward/Gram PBW sewing,
including all fixed spins. The test covers theta, mixed glasses and double-R
glasses, odd intrinsic NS parities, physical imaginary momenta, changing
insertion topology, and shrinking then expanding the requested domain.

The optional `cache_products` argument to `run_pbw_matrix` enables this path in
full-cutoff independent comparisons. The cross-channel driver enables it by
default; `reuse_vertex_label_products=false` requests the uncached comparison.
Products remain local to one momentum plan and are discarded with that plan.
The source integration report records actual cache hits and CCY transitions;
multi-process reports sum these counts over the disjoint momentum workers.
`run_cross_smoke` also compares freshly computed source integrals with and
without product retention, recording discrepancies, transition counts and
timings separately from the independent target-method comparison.
