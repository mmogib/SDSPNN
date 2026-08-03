# F2 certified weighted radius, theorem envelope, and analytic basin marker.
# Usage: julia --project=. scripts/s40_f2_certified_radius.jl [--force]

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using SDSPNN
using DBInterface
using LinearAlgebra
using Printf

const ARTIFACT_ID = "F2_R2_20260801"

function main()
    force = "--force" in ARGS
    logpath, tee = setup_logging("s40_f2_certified_radius")
    db = open_artifact_db()
    try
        c = analytic_constants(1.0)
        abs(c.CL - 42.71491868) <= 2e-8 || error("fixed C_L mismatch")
        Rclean = c.R_clean
        eta = c.eta_R
        target = Rclean / 2
        weighted_basin_radius = weighted_radius([0.8, 0.0]; tau=1.0)
        abs(weighted_basin_radius - 0.8sqrt(2.0)) <= 2e-15 ||
            error("weighted basin boundary mismatch")
        r_small = euclidean_radius_for_weighted(target, 0.0; tau=1.0)
        abs(weighted_radius([r_small, 0.0]; tau=1.0) - target) <= 2e-18 ||
            error("weighted initial radius construction failed")
        save_grid = collect(range(0.0, 2000.0; length=2001))
        specs = [
            (run_id="certified_baseline", tol=1e-10, radius=r_small, samples=2001),
            (run_id="certified_strict", tol=1e-12, radius=r_small, samples=2001),
            (run_id="outer_r0_0p5", tol=1e-10, radius=0.5, samples=2001),
            (run_id="radius_markers", tol=0.0, radius=0.0, samples=2),
        ]
        config = Dict("CL_fixed" => 42.71491868, "R_clean" => Rclean,
                      "eta_R" => eta, "weighted_initial_radius" => target,
                      "weighted_basin_boundary_radius" => weighted_basin_radius,
                      "outer_initial_radius" => 0.5, "horizon" => 2000.0,
                      "save_points" => 2001, "envelope_slack" => 5e-11)
        protocol = Dict(
            "algorithm" => "Rodas5P(autodiff=AutoFiniteDiff())",
            "baseline_tolerance" => 1e-10, "strict_tolerance" => 1e-12,
            "dtmax" => 0.25, "maxiters" => 10_000_000,
            "save_grid" => "uniform 2001 points on [0,2000]",
            "weighted_radius" => "sqrt(x' * inv(M(x)) * x)",
            "event_handling" => "no state reset; weighted projection active set evaluated at every RHS call",
            "failure_policy" => "envelope excess above 5e-11 or non-success retcode fails",
            "RNG" => "none",
        )
        expected_rows = sum(s.samples for s in specs)
        created = begin_artifact!(db, ARTIFACT_ID, "F2", "scripts/s40_f2_certified_radius.jl",
                                  config; expected_runs=length(specs), expected_rows=expected_rows,
                                  protocol=protocol, force=force)
        if !created
            println(tee, "$ARTIFACT_ID already exists; use --force to replace")
            return
        end
        for spec in specs
            run_config = Dict("run_id" => spec.run_id, "tol" => spec.tol,
                              "initial_radius" => spec.radius, "samples" => spec.samples)
            register_expected_run!(db, ARTIFACT_ID, spec.run_id, run_config;
                                   expected_samples=spec.samples)
        end

        solutions = Dict{String,FlowResult}()
        for spec in specs[1:3]
            solutions[spec.run_id] = solve_flow([spec.radius, 0.0], (0.0, 2000.0);
                tau=1.0, abstol=spec.tol, reltol=spec.tol, save_grid=save_grid,
                dtmax=0.25, maxiters=10_000_000)
        end
        baseline_weighted = [weighted_radius(view(solutions["certified_baseline"].xs, :, j); tau=1.0)
                             for j in axes(solutions["certified_baseline"].xs, 2)]
        strict_weighted = [weighted_radius(view(solutions["certified_strict"].xs, :, j); tau=1.0)
                           for j in axes(solutions["certified_strict"].xs, 2)]
        convergence_difference = maximum(abs.(baseline_weighted .- strict_weighted))

        DBInterface.execute(db, "BEGIN IMMEDIATE TRANSACTION")
        try
            for spec in specs[1:3]
                sol = solutions[spec.run_id]
                is_certified = startswith(spec.run_id, "certified")
                wradii = [weighted_radius(view(sol.xs, :, j); tau=1.0) for j in axes(sol.xs, 2)]
                envelope = is_certified ? [target * exp(-eta*t) for t in sol.ts] : fill(nothing, length(sol.ts))
                max_excess = is_certified ? maximum(max(wradii[j] - envelope[j], 0.0) for j in eachindex(wradii)) : 0.0
                max_excess <= 5e-11 || error("F2 envelope gate failed for $(spec.run_id)")
                run_config = Dict("run_id" => spec.run_id, "tol" => spec.tol,
                                  "initial_radius" => spec.radius, "samples" => spec.samples)
                record_run!(db, ARTIFACT_ID, spec.run_id, run_config;
                            run_status=sol.run_status, terminal_success=sol.run_status == "completed",
                            expected_samples=spec.samples, actual_samples=length(sol.ts),
                            details=Dict("max_envelope_excess" => max_excess,
                                         "convergence_difference_baseline_vs_strict" => convergence_difference,
                                         "weighted_initial_radius" => wradii[1],
                                         "feasibility_drift" => sol.feasibility_drift,
                                         "retcode" => sol.retcode))
                for (j, t) in enumerate(sol.ts)
                    x = view(sol.xs, :, j)
                    insert_artifact_row!(db, ARTIFACT_ID, spec.run_id, j, Dict(
                        "series" => spec.run_id, "t" => t, "x1" => x[1], "x2" => x[2],
                        "euclidean_radius" => norm(x), "weighted_radius" => wradii[j],
                        "theorem_envelope" => envelope[j],
                        "projection_active" => sol.projection_active[j],
                    ))
                end
            end
            ratio = weighted_basin_radius / Rclean
            decades = log10(ratio)
            marker_config = Dict("run_id" => "radius_markers", "tol" => 0.0,
                                 "initial_radius" => 0.0, "samples" => 2)
            record_run!(db, ARTIFACT_ID, "radius_markers", marker_config;
                        run_status="completed", terminal_success=true,
                        expected_samples=2, actual_samples=2,
                        details=Dict("weighted_basin_to_certificate_ratio" => ratio,
                                     "decades" => decades,
                                     "radius_kind" => "weighted"))
            insert_artifact_row!(db, ARTIFACT_ID, "radius_markers", 1, Dict(
                "series" => "radius_marker", "label" => "R_half", "radius" => Rclean,
                "radius_kind" => "weighted", "ratio" => ratio, "decades" => decades))
            insert_artifact_row!(db, ARTIFACT_ID, "radius_markers", 2, Dict(
                "series" => "radius_marker", "label" => "weighted_analytic_basin_boundary",
                "radius" => weighted_basin_radius, "euclidean_boundary_radius" => 0.8,
                "radius_kind" => "weighted", "ratio" => ratio, "decades" => decades))
            DBInterface.execute(db, "COMMIT")
        catch
            DBInterface.execute(db, "ROLLBACK")
            rethrow()
        end
        record_percentage!(db, ARTIFACT_ID, "completed", length(specs), length(specs))
        query = "SELECT run_id,row_index,payload_json FROM artifact_rows WHERE artifact_id='$ARTIFACT_ID' ORDER BY run_id,row_index"
        register_plot_series!(db, ARTIFACT_ID, "certified_radius", query)
        @printf(tee, "Pending F2: R_half=%.12e weighted_boundary=%.12e ratio=%.6e decades=%.6f strict_diff=%.3e\n",
                Rclean, weighted_basin_radius, weighted_basin_radius/Rclean,
                log10(weighted_basin_radius/Rclean), convergence_difference)
    finally
        close_artifact_db(db)
        teardown_logging(tee, logpath)
    end
end

main()
