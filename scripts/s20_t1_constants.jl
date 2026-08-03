# T1 constants and hypothesis table artifact.
# Usage: julia --project=. scripts/s20_t1_constants.jl [--force]

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using SDSPNN
using DBInterface
using Printf

const ARTIFACT_ID = "T1_R2_20260801"

function tau_row(instance, tau; provenance_label="analytic here")
    c = analytic_constants(tau)
    return Dict(
        "instance" => instance, "tau" => tau,
        "mu" => c.mu, "K" => c.K, "L" => c.L,
        "KM" => c.KM_bound, "Fmax" => c.Fmax, "alpha" => c.alpha,
        "window" => c.window, "inside_window" => c.inside_window,
        "rho" => c.rho, "CL" => c.CL, "R" => c.R_clean,
        "cg" => c.cg,
        "integrability_status" => tau == 0 ? "integrable" : "nonintegrable",
        "provenance" => Dict(
            "mu" => "manuscript/exact", "K" => "manuscript/exact",
            "L" => "$provenance_label/exact", "KM" => "$provenance_label/certified_bound",
            "Fmax" => "manuscript/exact", "alpha" => "manuscript/exact",
            "window" => "analytic here/exact from declared bounds",
            "rho" => "analytic here/exact from declared bounds",
            "CL" => "analytic here/exact from declared bounds",
            "R" => "analytic here/exact from declared bounds",
            "cg" => "analytic here/certified lower bound",
        ),
        "exact_or_bound" => Dict(
            "mu" => "exact", "K" => "exact", "L" => "exact",
            "KM" => "upper_bound", "Fmax" => "exact", "alpha" => "exact",
            "window" => "bound_derived", "rho" => "bound_derived",
            "CL" => "bound_derived", "R" => "bound_derived", "cg" => "lower_bound",
        ),
    )
end

function control_row(instance, L, KM, spectral, integrability, provenance)
    na = "not_applicable_F4_control_only"
    return Dict(
        "instance" => instance, "tau" => na, "mu" => na, "K" => na,
        "L" => L, "KM" => KM, "Fmax" => na, "alpha" => na,
        "window" => na, "inside_window" => na, "rho" => na,
        "CL" => na, "R" => na, "cg" => na,
        "spectral_bounds" => spectral,
        "integrability_status" => integrability,
        "provenance" => provenance,
        "exact_or_bound" => Dict("L" => "certified", "KM" => "certified_bound"),
    )
end

function main()
    force = "--force" in ARGS
    logpath, tee = setup_logging("s20_t1_constants")
    db = open_artifact_db()
    try
        rows = [
            tau_row("A_counterexample", 1.0; provenance_label="manuscript"),
            tau_row("B_tau_0.1", 0.1),
            tau_row("B_tau_0.1394", 0.1394),
            tau_row("B_tau_1", 1.0),
            control_row("C_bounded_hessian", control_constants(:C).L, control_constants(:C).KM_bound,
                        [control_constants(:C).spectral_lower, control_constants(:C).spectral_upper],
                        control_constants(:C).integrability,
                        Dict("metric" => "revision 2", "L" => "analytic/exact", "KM" => "analytic/certified_bound")),
            control_row("D_separable_diagonal", control_constants(:D).L, control_constants(:D).KM_bound,
                        [control_constants(:D).spectral_lower, control_constants(:D).spectral_upper],
                        control_constants(:D).integrability,
                        Dict("metric" => "revision 2", "L" => "analytic/exact", "KM" => "analytic/exact_on_domain")),
            control_row("E_manufactured_nonzero", control_constants(:E).L, control_constants(:E).KM_bound,
                        [control_constants(:E).spectral_lower, control_constants(:E).spectral_upper],
                        control_constants(:E).integrability,
                        Dict("metric" => "revision 2", "L" => "analytic/exact", "KM" => "analytic/exact_on_domain")),
        ]
        config = Dict("instances" => [row["instance"] for row in rows], "revision" => 2)
        protocol = Dict("kind" => "analytic constants", "RNG" => "none",
                        "flow_fields_for_controls" => "not_applicable")
        created = begin_artifact!(db, ARTIFACT_ID, "T1", "scripts/s20_t1_constants.jl",
                                  config; expected_runs=length(rows), expected_rows=length(rows),
                                  protocol=protocol, force=force)
        if !created
            println(tee, "$ARTIFACT_ID already exists; use --force to replace")
            return
        end
        DBInterface.execute(db, "BEGIN IMMEDIATE TRANSACTION")
        try
            for (index, row) in enumerate(rows)
                run_id = String(row["instance"])
                run_config = Dict("instance" => run_id, "revision" => 2)
                register_expected_run!(db, ARTIFACT_ID, run_id, run_config; expected_samples=1)
                record_run!(db, ARTIFACT_ID, run_id, run_config;
                            run_status="completed", terminal_success=true,
                            expected_samples=1, actual_samples=1)
                insert_artifact_row!(db, ARTIFACT_ID, run_id, 1, row)
            end
            DBInterface.execute(db, "COMMIT")
        catch
            DBInterface.execute(db, "ROLLBACK")
            rethrow()
        end
        record_percentage!(db, ARTIFACT_ID, "completed", length(rows), length(rows))
        query = "SELECT run_id,row_index,payload_json FROM artifact_rows WHERE artifact_id='$ARTIFACT_ID' ORDER BY run_id,row_index"
        register_plot_series!(db, ARTIFACT_ID, "constants_table", query)
        println(tee, "Pending T1 artifact: $ARTIFACT_ID, rows=$(length(rows))")
    finally
        close_artifact_db(db)
        teardown_logging(tee, logpath)
    end
end

main()
