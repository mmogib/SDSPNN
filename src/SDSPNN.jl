module SDSPNN

using LinearAlgebra
using OrdinaryDiffEq
using Dates
using SHA
using JSON3
using SQLite
using DBInterface
using CSV
using DataFrames

const JCODE_ROOT = normpath(joinpath(@__DIR__, ".."))

include("types.jl")
include("metrics.jl")
include("problems.jl")
include("projection.jl")
include("active_projection.jl")
include("solver.jl")
include("residuals.jl")
include("classifier.jl")
include("io_utils.jl")

export
    DISK_RADIUS, ALPHA, LAMBDA, MU, K_OPERATOR, F_MAX, KM1_BOUND,
    ANNULUS_PERIOD, PROJECTION_ROOT_TOL, PROJECTION_GATE_TOL,
    ProjectionResult, FlowResult, ClassifierConfig,
    chi_radius, metric_tau, metric_bounds_tau, metric_hessian_control,
    metric_diagonal_control, metric_manufactured_control, control_constants,
    operator_matrix, operator_A, analytic_constants, rho_tau, cg_tau,
    cg_root, residual_norm_tau,
    weighted_disk_projection, weighted_disk_projection_spectral,
    active_projection_constants, active_projection_point, active_boundary_geometry,
    active_clamp_signature, active_sample_angles,
    active_projection_joint_bounds,
    flow_map, flow_rhs!, solve_flow, exact_annular_state, weighted_radius,
    euclidean_radius_for_weighted,
    central_residual, richardson_residual, target_inverse_metric,
    analytic_target_residual,
    classifier_config, classify_trajectory,
    JCODE_ROOT, TeeIO, setup_logging, teardown_logging, figure_output_path,
    canonical_json, config_hash, selection_hash,
    percentage_from_counts, validate_percentage, validate_accounting_rows,
    open_artifact_db, close_artifact_db, init_schema!, begin_artifact!,
    register_expected_run!, record_run!, insert_artifact_row!,
    record_percentage!, register_plot_series!, record_gate!,
    validate_artifact!, finalize_artifact!, artifact_status, export_artifact!

end
