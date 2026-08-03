# F4 Richardson integrability residual with zero/nonzero controls.
# Usage: julia --project=. scripts/s50_f4_integrability.jl [--force]

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using SDSPNN
using DBInterface
using LinearAlgebra
using Statistics
using Printf

const ARTIFACT_ID = "F4_R2_20260801"

function payload_for(case_name, Hfun, x, analytic)
    result = richardson_residual(Hfun, x; h=1e-4)
    error_vector = result.extrapolated - analytic
    return Dict(
        "case" => case_name, "x1" => x[1], "x2" => x[2],
        "radius" => norm(x), "theta" => norm(x) == 0 ? 0.0 : atan(x[2], x[1]),
        "raw_h" => result.raw_h, "raw_h2" => result.raw_h2, "raw_h4" => result.raw_h4,
        "residual1" => result.extrapolated[1], "residual2" => result.extrapolated[2],
        "residual_norm" => norm(result.extrapolated),
        "analytic1" => analytic[1], "analytic2" => analytic[2],
        "analytic_norm" => norm(analytic), "vector_error" => norm(error_vector),
        "observed_order" => result.observed_order,
        "truncation_estimate" => result.truncation_estimate, "h" => 1e-4,
    )
end

function main()
    force = "--force" in ARGS
    logpath, tee = setup_logging("s50_f4_integrability")
    db = open_artifact_db()
    try
        unit_axis = collect(range(0.0, 1.0; length=21))
        square_points = [[x, y] for x in unit_axis for y in unit_axis]
        disk_axis = collect(range(-1.1, 1.1; length=21))
        disk_points = [[x, y] for x in disk_axis for y in disk_axis if hypot(x, y) <= 1.1 + 1e-14]
        annulus_radii = [0.85, 0.90, 0.95, 1.00, 1.05]
        annulus_angles = collect(range(0.0, 2pi; length=33))[1:end-1]
        target_points = [r .* [cos(theta), sin(theta)] for r in annulus_radii for theta in annulus_angles]
        specs = [
            (run_id="C1_separable", points=square_points),
            (run_id="C2_hessian", points=disk_points),
            (run_id="C3_manufactured", points=square_points),
            (run_id="T_rotating_annulus", points=target_points),
        ]
        expected_rows = sum(length(s.points) for s in specs)
        config = Dict(
            "h" => 1e-4, "levels" => [1.0, 0.5, 0.25],
            "square_grid" => [21, 21], "disk_cartesian_grid" => [21, 21],
            "disk_points_after_filter" => length(disk_points),
            "target_radii" => annulus_radii, "target_angles" => length(annulus_angles),
        )
        protocol = Dict(
            "central_differences" => true, "richardson" => "fourth-order extrapolation",
            "observed_order" => "log2(norm(R_h-R_h2)/norm(R_h2-R_h4))",
            "truncation_estimate" => "norm(Rich_h2h4-Rich_hh2)/15",
            "join_policy" => "target differentiated only for radii 0.85 through 1.05; never across r=0.8",
            "RNG" => "none",
        )
        created = begin_artifact!(db, ARTIFACT_ID, "F4", "scripts/s50_f4_integrability.jl",
                                  config; expected_runs=length(specs), expected_rows=expected_rows,
                                  protocol=protocol, force=force)
        if !created
            println(tee, "$ARTIFACT_ID already exists; use --force to replace")
            return
        end
        for spec in specs
            run_config = Dict("run_id" => spec.run_id, "samples" => length(spec.points),
                              "h" => 1e-4)
            register_expected_run!(db, ARTIFACT_ID, spec.run_id, run_config;
                                   expected_samples=length(spec.points))
        end

        generated = Dict{String,Vector{Dict{String,Any}}}()
        generated["C1_separable"] = [payload_for("C1", metric_diagonal_control, x, [0.0, 0.0])
                                     for x in square_points]
        generated["C2_hessian"] = [payload_for("C2", metric_hessian_control, x, [0.0, 0.0])
                                   for x in disk_points]
        generated["C3_manufactured"] = [payload_for("C3", metric_manufactured_control, x, [0.0, x[2]/2])
                                        for x in square_points]
        generated["T_rotating_annulus"] = [payload_for("T", target_inverse_metric, x,
                                                        analytic_target_residual(x)) for x in target_points]

        max_c1 = maximum(row["residual_norm"] for row in generated["C1_separable"])
        max_c2 = maximum(row["residual_norm"] for row in generated["C2_hessian"])
        max_c3_error = maximum(row["vector_error"] for row in generated["C3_manufactured"])
        max_target_error = maximum(row["vector_error"] for row in generated["T_rotating_annulus"])
        target_orders = [row["observed_order"] for row in generated["T_rotating_annulus"]
                         if row["observed_order"] isa Real && isfinite(row["observed_order"])]
        median_order = median(target_orders)
        max_truncation = maximum(row["truncation_estimate"] for row in generated["T_rotating_annulus"])
        min_target_norm = minimum(row["residual_norm"] for row in generated["T_rotating_annulus"])
        max_c1 <= 1e-8 || error("F4 C1 floor gate failed")
        max_c2 <= 1e-8 || error("F4 C2 floor gate failed")
        max_c3_error <= 2e-8 || error("F4 C3 nonzero gate failed")
        max_target_error <= 2e-8 || error("F4 target closed-form gate failed")
        1.7 <= median_order <= 2.3 || error("F4 observed-order gate failed")
        min_target_norm >= 4.2 || error("F4 target norm has an artificial node")

        DBInterface.execute(db, "BEGIN IMMEDIATE TRANSACTION")
        try
            for spec in specs
                rows = generated[spec.run_id]
                run_config = Dict("run_id" => spec.run_id, "samples" => length(spec.points), "h" => 1e-4)
                record_run!(db, ARTIFACT_ID, spec.run_id, run_config;
                            run_status="completed", terminal_success=true,
                            expected_samples=length(rows), actual_samples=length(rows),
                            details=Dict("max_residual_norm" => maximum(row["residual_norm"] for row in rows),
                                         "max_vector_error" => maximum(row["vector_error"] for row in rows)))
                for (j, row) in enumerate(rows)
                    insert_artifact_row!(db, ARTIFACT_ID, spec.run_id, j, row)
                end
            end
            DBInterface.execute(db, "COMMIT")
        catch
            DBInterface.execute(db, "ROLLBACK")
            rethrow()
        end
        record_percentage!(db, ARTIFACT_ID, "completed", length(specs), length(specs))
        query = "SELECT run_id,row_index,payload_json FROM artifact_rows WHERE artifact_id='$ARTIFACT_ID' ORDER BY run_id,row_index"
        register_plot_series!(db, ARTIFACT_ID, "integrability_residual", query)
        record_gate!(db, "G3", "F4_full_grid", true, Dict(
            "C1_max_norm" => max_c1, "C2_max_norm" => max_c2,
            "C3_max_vector_error" => max_c3_error,
            "target_max_vector_error" => max_target_error,
            "target_median_order" => median_order,
            "target_max_truncation_estimate" => max_truncation,
            "target_min_norm" => min_target_norm,
        ))
        @printf(tee, "Pending F4: rows=%d C1=%.3e C2=%.3e C3err=%.3e Terr=%.3e order=%.6f\n",
                expected_rows, max_c1, max_c2, max_c3_error, max_target_error, median_order)
    finally
        close_artifact_db(db)
        teardown_logging(tee, logpath)
    end
end

main()
