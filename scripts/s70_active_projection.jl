# Static active weighted-disk projection artifact. No ODE is integrated.
# Usage: julia --project=. scripts/s70_active_projection.jl [--force]

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using SDSPNN
using DBInterface
using LinearAlgebra
using Printf

const ARTIFACT_ID = "AP1_R2_20260802"
const SAMPLE_COUNT = 8
const TAU_SPECS = [
    (tau=1.0, run_id="tau_1_primary", role="primary"),
    (tau=0.1, run_id="tau_0p1_secondary", role="secondary"),
]

function point_payload(point, anchor, constants, geometry, signature,
                       tau_role, point_role, sample_index, phi)
    projection = point.projection
    payload = Dict{String,Any}(
        "tau" => constants.tau,
        "tau_role" => tau_role,
        "point_role" => point_role,
        "sample_index" => sample_index,
        "phi" => phi,
        "x1" => point.x[1],
        "x2" => point.x[2],
        "y1" => point.y[1],
        "y2" => point.y[2],
        "P1" => projection.z[1],
        "P2" => projection.z[2],
        "nu" => projection.nu,
        "active" => projection.active,
        "active_margin" => point.active_margin,
        "projection_displacement" => point.projection_displacement,
        "map_displacement" => point.map_displacement,
        "displacement_bound" => constants.displacement_bound,
        "displacement_ratio" => point.map_displacement / constants.displacement_bound,
        "clamp1" => point.clamp[1],
        "clamp2" => point.clamp[2],
        "weighted_vs_clamp_gap" => point.clamp_gap,
        "wrong_metric_map_gap" => point.wrong_metric_map_gap,
        "spectral_nu" => point.spectral.nu,
        "spectral_nu_gap" => abs(point.spectral.nu - projection.nu),
        "spectral_projection_gap" => norm(point.spectral.z - projection.z),
        "kkt_primal" => projection.primal_residual,
        "kkt_dual" => projection.dual_residual,
        "kkt_stationarity" => projection.stationarity_residual,
        "kkt_complementarity" => projection.complementarity_residual,
        "L" => constants.L,
        "KM_bound" => constants.KM_bound,
        "rho" => constants.rho,
        "cg" => constants.cg,
        "arc_lo" => geometry.arc_lo,
        "arc_hi" => geometry.arc_hi,
        "radial_threshold" => geometry.radial_threshold,
        "clamp_equilibrium_phi" => signature.phi,
        "clamp_equilibrium_distance" => signature.distance_from_solution,
    )
    if point_role == "equilibrium_anchor"
        merge!(payload, Dict(
            "pair_rule" => "not_applicable_anchor_row",
            "joint_lhs" => nothing,
            "joint_rhs_exact" => nothing,
            "joint_ratio_exact" => nothing,
            "joint_rhs_L" => nothing,
            "joint_ratio_L" => nothing,
            "joint_rhs_L2" => nothing,
            "joint_ratio_L2" => nothing,
            "joint_argument_gap" => nothing,
            "joint_metric_gap" => nothing,
            "joint_anchor" => nothing,
        ))
    else
        bounds = active_projection_joint_bounds(point, anchor, constants.L)
        merge!(payload, Dict(
            "pair_rule" => "sample is (Q1,z1); equilibrium is (Q2,z2)",
            "joint_lhs" => bounds.lhs,
            "joint_rhs_exact" => bounds.rhs_exact,
            "joint_ratio_exact" => bounds.ratio_exact,
            "joint_rhs_L" => bounds.rhs_L,
            "joint_ratio_L" => bounds.ratio_L,
            "joint_rhs_L2" => bounds.rhs_L2,
            "joint_ratio_L2" => bounds.ratio_L2,
            "joint_argument_gap" => bounds.argument_gap,
            "joint_metric_gap" => bounds.metric_gap,
            "joint_anchor" => bounds.anchor,
        ))
    end
    return payload
end

function generate_run(spec)
    tau = spec.tau
    constants = active_projection_constants(tau)
    geometry = active_boundary_geometry(tau)
    signature = active_clamp_signature(tau)
    xbar = [DISK_RADIUS, 0.0]
    anchor = active_projection_point(xbar, tau)
    rows = Dict{String,Any}[
        point_payload(anchor, anchor, constants, geometry, signature,
                      spec.role, "equilibrium_anchor", 0, 0.0),
    ]
    for (index, phi) in enumerate(active_sample_angles(tau, SAMPLE_COUNT))
        x = DISK_RADIUS .* [cos(phi), sin(phi)]
        point = active_projection_point(x, tau)
        push!(rows, point_payload(point, anchor, constants, geometry, signature,
                                  spec.role, "active_arc_sample", index, phi))
    end
    return (constants=constants, geometry=geometry, signature=signature,
            anchor=anchor, rows=rows)
end

function validate_stage1!(generated)
    formulas = Dict(
        "L" => "1/(1-((1+sqrt(5))/4)*tau)",
        "KM_bound" => "3.571*tau",
        "rho" => "sqrt(1-2*alpha*mu/L+alpha^2*K^2)",
        "cg" => "1-rho-(alpha/2)*L^2*KM_bound*Fmax",
        "mu" => "1",
        "K" => "1",
        "Fmax" => "3.1",
    )
    g0_values = Dict{String,Any}()
    for spec in TAU_SPECS
        run = generated[spec.run_id]
        constants = run.constants
        L_formula = inv(1 - ((1 + sqrt(5.0)) / 4) * spec.tau)
        KM_formula = 3.571 * spec.tau
        rho_formula = sqrt(1 - 2ALPHA / L_formula + ALPHA^2)
        cg_formula = 1 - rho_formula - (ALPHA / 2) * L_formula^2 * KM_formula * 3.1
        constants.L == L_formula || error("G0 L formula failed at tau=$(spec.tau)")
        constants.KM_bound == KM_formula || error("G0 KM formula failed at tau=$(spec.tau)")
        constants.rho == rho_formula || error("G0 rho formula failed at tau=$(spec.tau)")
        constants.cg == cg_formula || error("G0 cg formula failed at tau=$(spec.tau)")
        g0_values[spec.run_id] = Dict(
            "tau" => spec.tau, "L" => constants.L,
            "KM_bound" => constants.KM_bound, "rho" => constants.rho,
            "cg" => constants.cg, "displacement_bound" => constants.displacement_bound,
        )
    end

    all_rows = reduce(vcat, (generated[spec.run_id].rows for spec in TAU_SPECS))
    anchors = [generated[spec.run_id].anchor for spec in TAU_SPECS]
    all(abs(anchor.projection.nu - 9 / 220) <= 5e-13 for anchor in anchors) ||
        error("G1 equilibrium multiplier failed")
    all(anchor.map_displacement <= 5e-13 for anchor in anchors) ||
        error("G2 equilibrium fixed-point gate failed")
    max_nu_gap = maximum(row["spectral_nu_gap"] for row in all_rows)
    max_projection_gap = maximum(row["spectral_projection_gap"] for row in all_rows)
    max_nu_gap <= 2e-13 || error("G3 independent multiplier comparison failed")
    max_projection_gap <= 2e-13 || error("G3 independent projection comparison failed")
    all(row["active"] && row["active_margin"] > 0 for row in all_rows) ||
        error("G4 inactive point entered the stage-1 sample")
    all(row["map_displacement"] <= row["displacement_bound"] + 2e-12 for row in all_rows) ||
        error("G4 displacement bound failed")
    pair_rows = [row for row in all_rows if row["point_role"] == "active_arc_sample"]
    all(row["joint_lhs"] <= row["joint_rhs_exact"] + 2e-12 &&
        row["joint_lhs"] <= row["joint_rhs_L"] + 2e-12 &&
        row["joint_lhs"] <= row["joint_rhs_L2"] + 2e-12 for row in pair_rows) ||
        error("G4 joint-projection bound failed")
    length(all_rows) == 2 * (SAMPLE_COUNT + 1) || error("G5 row count failed")
    all(count(row -> row["point_role"] == "active_arc_sample",
              generated[spec.run_id].rows) == SAMPLE_COUNT for spec in TAU_SPECS) ||
        error("G5 per-parameter sample count failed")
    return (formulas=formulas, g0_values=g0_values, all_rows=all_rows,
            max_nu_gap=max_nu_gap, max_projection_gap=max_projection_gap)
end

function print_table(io, generated)
    println(io, "tau role idx phi nu P1 P2 ||P-y|| joint_lhs rhs_exact ratio_exact rhs_L ratio_L ratio_L2 ||T-x|| clamp_gap")
    for spec in TAU_SPECS
        for row in generated[spec.run_id].rows
            lhs = row["joint_lhs"] === nothing ? NaN : row["joint_lhs"]
            rhs_exact = row["joint_rhs_exact"] === nothing ? NaN : row["joint_rhs_exact"]
            ratio_exact = row["joint_ratio_exact"] === nothing ? NaN : row["joint_ratio_exact"]
            rhs_L = row["joint_rhs_L"] === nothing ? NaN : row["joint_rhs_L"]
            ratio_L = row["joint_ratio_L"] === nothing ? NaN : row["joint_ratio_L"]
            ratio_L2 = row["joint_ratio_L2"] === nothing ? NaN : row["joint_ratio_L2"]
            @printf(io, "%.1f %-9s %2d % .6f %.12f % .8f % .8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n",
                    row["tau"], row["tau_role"], row["sample_index"], row["phi"],
                    row["nu"], row["P1"], row["P2"], row["projection_displacement"],
                    lhs, rhs_exact, ratio_exact, rhs_L, ratio_L, ratio_L2,
                    row["map_displacement"],
                    row["weighted_vs_clamp_gap"])
        end
    end
end

function main()
    force = "--force" in ARGS
    logpath, tee = setup_logging("s70_active_projection")
    db = open_artifact_db()
    try
        generated = Dict(spec.run_id => generate_run(spec) for spec in TAU_SPECS)
        validation = validate_stage1!(generated)
        expected_rows = 2 * (SAMPLE_COUNT + 1)
        config = Dict(
            "stage" => 1,
            "tau_order" => [spec.tau for spec in TAU_SPECS],
            "tau_roles" => Dict(string(spec.tau) => spec.role for spec in TAU_SPECS),
            "sample_count_per_tau" => SAMPLE_COUNT,
            "sample_rule" => "eight equal interior fractions of each exact active boundary arc",
            "pair_rule" => "each active sample is (Q1,z1); its same-tau equilibrium is (Q2,z2)",
            "formulas" => validation.formulas,
        )
        protocol = Dict(
            "kind" => "static active weighted-disk projection",
            "ODE_integration" => "none",
            "RNG" => "none",
            "active_test" => "norm(x-alpha*M_tau(x)*(x-c))-1.1 > 0",
            "projection_primary" => "matrix-solve secular equation with certified bracket",
            "projection_independent" => "eigendecomposition secular equation with separate bisection",
            "joint_rhs_exact" => "sqrt(lambda_max(Q1)/lambda_min(Q1))*norm(z1-z2)+norm(Q1-Q2)/lambda_min(Q1)*anchor",
            "joint_rhs_L" => "L*norm(z1-z2)+L*norm(Q1-Q2)*anchor",
            "joint_rhs_L2" => "L^2*norm(z1-z2)+L*norm(Q1-Q2)*anchor",
            "aggregation_rules" => Dict(
                "G3" => "maximum absolute discrepancy across all 18 per-point rows",
                "G4" => "universal check across all 18 rows and all 16 non-anchor pairs",
                "G5" => "exact row count, nine rows for each of two tau values",
            ),
        )
        created = begin_artifact!(db, ARTIFACT_ID, "AP1",
                                  "scripts/s70_active_projection.jl", config;
                                  expected_runs=length(TAU_SPECS), expected_rows=expected_rows,
                                  protocol=protocol, force=force,
                                  notes="Stage 1 only. Static computation. No ODE and no figure.")
        if !created
            println(tee, "$ARTIFACT_ID already exists; use --force to replace")
            return
        end
        for spec in TAU_SPECS
            run_config = Dict(
                "tau" => spec.tau, "role" => spec.role,
                "sample_count" => SAMPLE_COUNT,
                "sample_rule" => config["sample_rule"],
            )
            register_expected_run!(db, ARTIFACT_ID, spec.run_id, run_config;
                                   expected_samples=SAMPLE_COUNT + 1)
        end

        DBInterface.execute(db, "BEGIN IMMEDIATE TRANSACTION")
        try
            for spec in TAU_SPECS
                run = generated[spec.run_id]
                run_config = Dict(
                    "tau" => spec.tau, "role" => spec.role,
                    "sample_count" => SAMPLE_COUNT,
                    "sample_rule" => config["sample_rule"],
                )
                max_kkt = maximum(maximum((row["kkt_primal"], row["kkt_dual"],
                                           row["kkt_stationarity"], row["kkt_complementarity"]))
                                   for row in run.rows)
                record_run!(db, ARTIFACT_ID, spec.run_id, run_config;
                            run_status="completed", terminal_success=true,
                            expected_samples=length(run.rows), actual_samples=length(run.rows),
                            details=Dict(
                                "arc" => [run.geometry.arc_lo, run.geometry.arc_hi],
                                "radial_threshold" => run.geometry.radial_threshold,
                                "clamp_signature" => Dict(
                                    "phi" => run.signature.phi,
                                    "distance_from_solution" => run.signature.distance_from_solution,
                                ),
                                "max_kkt_residual" => max_kkt,
                                "max_kkt_aggregation_rule" => "maximum of four diagnostics over nine per-point rows",
                            ))
                for (index, row) in enumerate(run.rows)
                    insert_artifact_row!(db, ARTIFACT_ID, spec.run_id, index, row)
                end
            end
            DBInterface.execute(db, "COMMIT")
        catch
            DBInterface.execute(db, "ROLLBACK")
            rethrow()
        end
        query = "SELECT run_id,row_index,payload_json FROM artifact_rows WHERE artifact_id='$ARTIFACT_ID' ORDER BY CASE run_id WHEN 'tau_1_primary' THEN 1 ELSE 2 END,row_index"
        register_plot_series!(db, ARTIFACT_ID, "active_projection_stage1_table", query)
        record_gate!(db, "AP_G0", "formula_recomputation", true,
                     Dict("formulas" => validation.formulas, "values" => validation.g0_values))
        record_gate!(db, "AP_G1", "equilibrium_multiplier", true,
                     Dict("target" => "9/220", "tolerance" => 5e-13))
        record_gate!(db, "AP_G2", "equilibrium_fixed_point", true,
                     Dict("tolerance" => 5e-13,
                          "wrong_metric_gaps" => Dict(spec.run_id =>
                              generated[spec.run_id].anchor.wrong_metric_map_gap for spec in TAU_SPECS)))
        record_gate!(db, "AP_G3", "independent_spectral_projection", true,
                     Dict("max_nu_gap" => validation.max_nu_gap,
                          "max_projection_gap" => validation.max_projection_gap,
                          "aggregation_rule" => "maximum across all 18 per-point rows"))
        record_gate!(db, "AP_G4", "all_points_active_and_bounds_hold", true,
                     Dict("rows" => expected_rows, "pairs" => 2SAMPLE_COUNT,
                          "aggregation_rule" => "universal check across every stored row and pair"))
        record_gate!(db, "AP_G5", "per_point_rows", true,
                     Dict("rows_per_tau" => SAMPLE_COUNT + 1,
                          "sample_rows_per_tau" => SAMPLE_COUNT,
                          "aggregation_rule" => "exact counts, not a summary statistic"))
        print_table(tee, generated)
        println(tee, "Stage 1 is static; no ODE was integrated and no figure was generated.")
        @printf(tee, "The exact active arcs are (%.6f, %.6f) at tau=1 and (%.6f, %.6f) at tau=0.1.\n",
                generated["tau_1_primary"].geometry.arc_lo,
                generated["tau_1_primary"].geometry.arc_hi,
                generated["tau_0p1_secondary"].geometry.arc_lo,
                generated["tau_0p1_secondary"].geometry.arc_hi)
        @printf(tee, "At the equilibrium, the weighted-versus-radial-clamp gap is %.8f at tau=1 and %.8f at tau=0.1.\n",
                generated["tau_1_primary"].anchor.clamp_gap,
                generated["tau_0p1_secondary"].anchor.clamp_gap)
        println(tee, "Pending artifact $ARTIFACT_ID has $expected_rows per-point rows; run s90 to validate, finalize, and export it.")
    finally
        close_artifact_db(db)
        teardown_logging(tee, logpath)
    end
end

main()
