# F1 annulus trajectories and boundary curves.
# Usage: julia --project=. scripts/s30_f1_annulus.jl [--force]

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using SDSPNN
using DBInterface
using LinearAlgebra
using Printf

const ARTIFACT_ID = "F1_R2_20260801"

function main()
    force = "--force" in ARGS
    logpath, tee = setup_logging("s30_f1_annulus")
    db = open_artifact_db()
    try
        radii = [0.2, 0.5, 0.75, 0.85, 0.95, 1.05]
        save_grid = collect(range(0.0, 520.0; length=2001))
        boundary_angles = collect(range(0.0, 2pi; length=361))
        run_specs = [(run_id="r_$(replace(string(r), "." => "p"))", kind="trajectory",
                      radius=r, samples=length(save_grid)) for r in radii]
        append!(run_specs, [
            (run_id="disk_boundary", kind="boundary", radius=DISK_RADIUS, samples=length(boundary_angles)),
            (run_id="basin_boundary", kind="boundary", radius=0.8, samples=length(boundary_angles)),
        ])
        expected_rows = sum(spec.samples for spec in run_specs)
        config = Dict("radii" => radii, "initial_angle" => 0.0, "horizon" => 520.0,
                      "save_points" => 2001, "boundary_points" => 361,
                      "period" => ANNULUS_PERIOD)
        protocol = Dict(
            "algorithm" => "Rodas5P(autodiff=AutoFiniteDiff())",
            "abstol" => 1e-10, "reltol" => 1e-10, "dtmax" => 0.25,
            "maxiters" => 10_000_000, "save_grid" => "uniform 2001 points on [0,520]",
            "event_handling" => "no state reset; weighted projection active set evaluated at every RHS call",
            "failure_policy" => "non-success retcode, nonfinite state, KKT failure, or feasibility drift fails",
            "RNG" => "none",
        )
        created = begin_artifact!(db, ARTIFACT_ID, "F1", "scripts/s30_f1_annulus.jl",
                                  config; expected_runs=length(run_specs), expected_rows=expected_rows,
                                  protocol=protocol, force=force)
        if !created
            println(tee, "$ARTIFACT_ID already exists; use --force to replace")
            return
        end

        for spec in run_specs
            run_config = Dict("run_id" => spec.run_id, "kind" => spec.kind,
                              "radius" => spec.radius, "samples" => spec.samples)
            register_expected_run!(db, ARTIFACT_ID, spec.run_id, run_config;
                                   expected_samples=spec.samples)
        end

        DBInterface.execute(db, "BEGIN IMMEDIATE TRANSACTION")
        try
            for radius in radii
                run_id = "r_$(replace(string(radius), "." => "p"))"
                run_config = Dict("run_id" => run_id, "kind" => "trajectory",
                                  "radius" => radius, "samples" => length(save_grid))
                sol = solve_flow([radius, 0.0], (0.0, 520.0); tau=1.0,
                                 abstol=1e-10, reltol=1e-10, save_grid=save_grid,
                                 dtmax=0.25, maxiters=10_000_000)
                analytic_class = radius < 0.8 ? "convergent" : "periodic"
                radial_error = radius >= 0.8 ?
                    maximum(abs(norm(view(sol.xs, :, j)) - radius) for j in axes(sol.xs, 2)) : nothing
                phase_error = radius >= 0.8 ? maximum(
                    abs(atan(det(hcat(exact_annular_state([radius, 0.0], t), view(sol.xs, :, j))),
                             dot(exact_annular_state([radius, 0.0], t), view(sol.xs, :, j))))
                    for (j, t) in enumerate(sol.ts)) : nothing
                record_run!(db, ARTIFACT_ID, run_id, run_config;
                            run_status=sol.run_status, terminal_success=sol.run_status == "completed",
                            expected_samples=length(save_grid), actual_samples=length(sol.ts),
                            details=Dict("analytic_class" => analytic_class,
                                         "radial_error_if_annular" => radial_error,
                                         "phase_error_if_annular" => phase_error,
                                         "feasibility_drift" => sol.feasibility_drift,
                                         "retcode" => sol.retcode))
                for (j, t) in enumerate(sol.ts)
                    x = view(sol.xs, :, j)
                    insert_artifact_row!(db, ARTIFACT_ID, run_id, j, Dict(
                        "series" => "trajectory", "start_radius" => radius,
                        "analytic_class" => analytic_class, "t" => t,
                        "x1" => x[1], "x2" => x[2], "euclidean_radius" => norm(x),
                        "weighted_radius" => weighted_radius(x; tau=1.0),
                        "projection_active" => sol.projection_active[j],
                    ))
                end
            end
            for (run_id, radius) in (("disk_boundary", DISK_RADIUS), ("basin_boundary", 0.8))
                run_config = Dict("run_id" => run_id, "kind" => "boundary",
                                  "radius" => radius, "samples" => length(boundary_angles))
                record_run!(db, ARTIFACT_ID, run_id, run_config;
                            run_status="completed", terminal_success=true,
                            expected_samples=length(boundary_angles), actual_samples=length(boundary_angles))
                for (j, theta) in enumerate(boundary_angles)
                    insert_artifact_row!(db, ARTIFACT_ID, run_id, j, Dict(
                        "series" => "boundary", "boundary_name" => run_id,
                        "theta" => theta, "x1" => radius*cos(theta),
                        "x2" => radius*sin(theta), "euclidean_radius" => radius,
                    ))
                end
            end
            DBInterface.execute(db, "COMMIT")
        catch
            DBInterface.execute(db, "ROLLBACK")
            rethrow()
        end
        record_percentage!(db, ARTIFACT_ID, "completed", length(run_specs), length(run_specs))
        query = "SELECT run_id,row_index,payload_json FROM artifact_rows WHERE artifact_id='$ARTIFACT_ID' ORDER BY run_id,row_index"
        register_plot_series!(db, ARTIFACT_ID, "phase_portrait", query)
        @printf(tee, "Pending F1 artifact: %s, runs=%d, rows=%d, horizon/period=%.6f\n",
                ARTIFACT_ID, length(run_specs), expected_rows, 520 / ANNULUS_PERIOD)
    finally
        close_artifact_db(db)
        teardown_logging(tee, logpath)
    end
end

main()
