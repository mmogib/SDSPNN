# F5 analytic certified band versus true convergence band.
# Usage: julia --project=. scripts/s60_f5_variation_band.jl [--force]

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using SDSPNN
using DBInterface
using Printf

const ARTIFACT_ID = "F5_R2_20260801"

function band_payload(tau, grid_kind)
    c = analytic_constants(tau)
    return Dict(
        "grid_kind" => grid_kind, "tau" => tau, "L" => c.L,
        "KM_bound" => c.KM_bound, "rho" => c.rho, "cg" => c.cg,
        "certified" => c.cg > 0, "analytic_decay_rate" => ALPHA * (1 - tau),
        "residual_norm_r_0p95" => residual_norm_tau(tau, 0.95),
        "analytic_flow_class" => tau < 1 ? "convergent" : "periodic_endpoint",
    )
end

function main()
    force = "--force" in ARGS
    logpath, tee = setup_logging("s60_f5_variation_band")
    db = open_artifact_db()
    try
        plot_grid = collect(range(0.0, 1.0; length=1001))
        refinement_grid = collect(range(0.139, 0.140; length=101))
        root = cg_root(0.139, 0.140; tol=1e-14)
        cg_tau(0.139) > 0 > cg_tau(0.140) || error("F5 analytic root bracket failed")
        specs = [
            (run_id="plot_grid", samples=length(plot_grid)),
            (run_id="root_refinement", samples=length(refinement_grid)),
            (run_id="root_summary", samples=1),
        ]
        expected_rows = sum(s.samples for s in specs)
        config = Dict(
            "tau_domain" => [0.0, 1.0], "plot_grid_points" => length(plot_grid),
            "analytic_root_bracket" => [0.139, 0.140],
            "refinement_points" => length(refinement_grid), "root_tolerance" => 1e-14,
            "residual_radius" => 0.95,
        )
        protocol = Dict(
            "cg" => "analytic condition (19) with KM upper bound 3.571*tau",
            "decay_rate" => "analytic alpha*(1-tau); no fixed-horizon classifier used",
            "residual" => "analytic tau*sqrt(5)/(2*det(D_tau)*0.95)",
            "root_method" => "analytic sign bracket followed by bisection",
            "RNG" => "none", "simulation" => "none; analytic rate is the source of the claim",
        )
        created = begin_artifact!(db, ARTIFACT_ID, "F5", "scripts/s60_f5_variation_band.jl",
                                  config; expected_runs=length(specs), expected_rows=expected_rows,
                                  protocol=protocol, force=force)
        if !created
            println(tee, "$ARTIFACT_ID already exists; use --force to replace")
            return
        end
        for spec in specs
            run_config = Dict("run_id" => spec.run_id, "samples" => spec.samples,
                              "root_bracket" => [0.139, 0.140])
            register_expected_run!(db, ARTIFACT_ID, spec.run_id, run_config;
                                   expected_samples=spec.samples)
        end

        DBInterface.execute(db, "BEGIN IMMEDIATE TRANSACTION")
        try
            for (run_id, grid, kind) in (("plot_grid", plot_grid, "plot"),
                                         ("root_refinement", refinement_grid, "refinement"))
                run_config = Dict("run_id" => run_id, "samples" => length(grid),
                                  "root_bracket" => [0.139, 0.140])
                record_run!(db, ARTIFACT_ID, run_id, run_config;
                            run_status="completed", terminal_success=true,
                            expected_samples=length(grid), actual_samples=length(grid))
                for (j, tau) in enumerate(grid)
                    insert_artifact_row!(db, ARTIFACT_ID, run_id, j, band_payload(tau, kind))
                end
            end
            factor = 1 / root
            summary_config = Dict("run_id" => "root_summary", "samples" => 1,
                                  "root_bracket" => [0.139, 0.140])
            record_run!(db, ARTIFACT_ID, "root_summary", summary_config;
                        run_status="completed", terminal_success=true,
                        expected_samples=1, actual_samples=1,
                        details=Dict("root" => root, "true_band_endpoint" => 1.0,
                                     "conservatism_factor" => factor))
            insert_artifact_row!(db, ARTIFACT_ID, "root_summary", 1, Dict(
                "grid_kind" => "root_summary", "tau" => root, "cg" => cg_tau(root),
                "bracket_lo" => 0.139, "bracket_hi" => 0.140,
                "c_g_lo" => cg_tau(0.139), "c_g_hi" => cg_tau(0.140),
                "true_band_endpoint" => 1.0, "conservatism_factor" => factor,
                "residual_norm_r_0p95" => residual_norm_tau(root, 0.95),
                "analytic_decay_rate" => ALPHA * (1 - root),
            ))
            DBInterface.execute(db, "COMMIT")
        catch
            DBInterface.execute(db, "ROLLBACK")
            rethrow()
        end
        record_percentage!(db, ARTIFACT_ID, "completed", length(specs), length(specs))
        query = "SELECT run_id,row_index,payload_json FROM artifact_rows WHERE artifact_id='$ARTIFACT_ID' ORDER BY run_id,row_index"
        register_plot_series!(db, ARTIFACT_ID, "variation_band", query)
        @printf(tee, "Pending F5: root=%.12f c_g(0.1)=%.10f residual(0.1)=%.10f factor=%.6f\n",
                root, cg_tau(0.1), residual_norm_tau(0.1, 0.95), 1/root)
    finally
        close_artifact_db(db)
        teardown_logging(tee, logpath)
    end
end

main()
