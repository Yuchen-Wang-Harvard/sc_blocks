module SuperconformalBlocks

using LinearAlgebra
using Dates
using SHA
using TOML
using Printf
import Base: contains

const CONVENTION = "bpz_native_november_2026_no_f_ward_unit_ns_ground"

include("combinatorics.jl")
include("graphs.jl")
include("sca.jl")
include("ward.jl")
include("fermion.jl")
include("series.jl")
include("numeric_scratch.jl")
include("schottky.jl")
include("ccy.jl")
include("contractions.jl")
include("ns_branching.jl")
include("ns_pipeline.jl")
include("freefield.jl")
include("span_solve.jl")
include("branch_actions.jl")
include("ramond_branching.jl")
include("parity_series.jl")
include("ramond_pipeline.jl")
include("splitting.jl")
include("split_pipeline.jl")
include("ns_recursion.jl")
include("ns_schottky.jl")
include("quadrature.jl")
include("liouville.jl")
include("theta_geometry.jl")

export CONVENTION, Combinatorics, partitions, fermion_partitions,
       CoefficientDomain, total_domain, box_domain, contains,
       PlumbingGraph, theta_graph, glasses_graph, tetrahedron_graph,
       sewing_sign, split_graph, ramond_cycles, insertion_edges, mixed_convolution_sign,
       Mode, L, G, PBWState, ground_state, level2, relative_parity,
       SCAModule, state_parity, basis, act, act_vector, action_matrix,
       pairing, gram, gram_factor, inverse_gram,
       WardVertex, ns_vertex, nsrr_vertex, rho, ward_values,
       mode_value, master_ward_residual, virasoro_ward_residual,
       FermionState, FermionModule, FermionVertex, act_psi, fermion_norm,
       psi_value, tensor_pairing, enlarged_rho,
       Series, series_constant, series_add, series_scale, series_multiply,
       series_inverse, series_exp, residual_domain,
       primitive_cycles, schottky_vacuum,
       CCYPlan, ccy_reduced, virasoro_block, global_vertex,
       NetworkTensor, contract_network, sew_tensors, pbw_coefficient, evaluate_spins,
       ell, branch_central, branch_weight, ns_branch_norm, ramond_branch_norm,
       ns_branching, ns_branching_raw, NSBlockPlan, ns_block, ns_spin_coefficient,
       OscillatorState, FreeField, intern!, ff_physical_apply, ff_embedded_apply,
       ff_primary, ff_descendant, ff_norm, ordered_pbw, from_ordered_pbw,
       SpanSolver, solve_span, BranchActions, BranchActionTerm, branch_action,
       single_leg_vertex, NSRRBranching, prepare_branching!, branching_raw,
       ParitySeries, walsh!, parity_product, RamondBlockPlan, ramond_block, set_vertex_labels!,
       theta_psi_action, theta_psi_apply, psi_weights, SplittingPlan, splitting_raw,
       SplitBlockPlan, split_domain, split_diagonal, inserted_ramond_block,
       NSRecursionPlan, ns_c_reduced, ns_schottky, ns_c_block,
       gauss_legendre, gauss_laguerre, tanh_sinh_integral,
       LiouvillePlan, log_upsilon, log_sector_upsilon, liouville_ns, liouville_nsrr,
       reference_log_upsilon, ThetaGeometryPlan, theta_geometry, theta_characteristic,
       theta_constant, free_partition, find_symplectic_map, refine_theta_surface,
       transform_period, transport_characteristic, numeric_schottky_logP

end
