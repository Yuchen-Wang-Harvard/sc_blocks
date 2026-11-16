using Test
using LinearAlgebra
using SuperconformalBlocks
const SCB = SuperconformalBlocks

include("test_local.jl")
include("test_fermion.jl")
include("test_numeric_scratch.jl")
include("test_ccy.jl")
include("test_sewing.jl")
include("test_ns_pipeline.jl")
include("test_freefield.jl")
include("test_branch_actions.jl")
include("test_ramond_branching.jl")
include("test_ramond_pipeline.jl")
include("test_splitting.jl")
include("test_split_pipeline.jl")
include("test_reuse.jl")
include("test_ns_raw.jl")
include("test_mixed_pipeline.jl")
include("test_vertex_label_reuse.jl")
include("test_precision.jl")
include("test_recovery_policy.jl")
include("test_ns_recursion.jl")
include("test_fusion_cache.jl")
include("test_ns_reference_reuse.jl")
include("test_ns_vacuum_cycles.jl")
include("test_analytic_components.jl")
