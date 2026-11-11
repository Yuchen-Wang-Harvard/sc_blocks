# SuperconformalBlocks.jl

Julia implementation of arbitrary-plumbing N=1 superconformal blocks. The
production algorithm uses branching recurrences, Virasoro CCY recursion,
Schottky/global seeds and auxiliary-fermion recovery. Odd Ramond loops use
split-edge insertions. Independent PBW sewing and SCA c-recursion provide
checks of the resulting coefficients.

## Package contents

- `src/`: numerical library.
- `bin/blocks.jl`: command-line block computation.
- `checks/`: independent comparison and timing drivers.
- `test/`: algebra, Ward-identity and component tests.
- `examples/`: runnable block and check input configurations.
- `docs/`: code/function reference and mathematical implementation details.
- `Project.toml`: Julia package dependencies and compatibility.

The package uses Julia standard libraries only. Install Julia separately;
Julia 1.13.3 was used for the checks. Run the commands below from this directory.
Choose unused output paths: the drivers preserve existing output files.

## Compute blocks

```sh
julia --project=. bin/blocks.jl examples/ns_theta.toml results/ns_theta.toml
julia --project=. bin/blocks.jl examples/nsrr_theta.toml results/nsrr_theta.toml
julia --project=. bin/blocks.jl examples/nsrr_theta_inserted.toml results/nsrr_theta_inserted.toml
```

The example blocks have total descendant level at most 3. Benchmark inputs
use explicit per-edge cutoffs. Couplings and momenta are quoted decimal or
rational strings; complex parameters are pairs of such strings. Intrinsic NS
primary parities are explicit input data.

Output `level2` indices are twice the descendant levels. NSRR coefficients
retain every edge-parity component: array index `1+mask`, with edge `e` at
bit `e-1`. The driver determines insertion edges from the vertex eta signs.
Working precision and achieved numerical accuracy are distinct; output records
include the precision, recovery diagnostics, timings and source hashes.

## Run checks

```sh
julia --project=. test/runtests.jl
julia --project=. checks/run_foundations.jl results/foundations.toml
julia --project=. checks/run_pbw_matrix.jl nsrr_theta 5 results/pbw_theta.toml 0 192
julia --project=. checks/run_ns_recursion_matrix.jl theta 10 results/ns_recursion_theta.toml 192 0 imaginary
julia --project=. checks/run_cross_channel.jl examples/cross_channel.toml results/cross_channel.toml
julia --project=. checks/benchmark_blocks.jl examples/benchmark_theta_ordinary.toml results/ordinary 1
```

The PBW driver supports all-NS theta/glasses, NSRR theta, mixed/double-R
glasses and the tetrahedron with a Ramond face. Its genus-two cutoffs are per
edge; the tetrahedron cutoff is total. The NS c-recursion driver uses per-edge
cutoffs for theta and total cutoffs for tetrahedron. See the file/function
reference for all check drivers and their roles.

The cross-channel input specifies source/target total cutoffs, separate momentum
grids, precision and convergence controls. Its quadrature and channel-ratio
checks are reported separately. Generated outputs are not part of the package.

## Documentation

- [Code files and functions](docs/CODE_MAP.md)
- [Mathematical conventions and numerical settings](docs/CONVENTIONS.md)
- [Split-edge ordering and unequal intermediate grades](docs/SPLIT_SEWING.md)
