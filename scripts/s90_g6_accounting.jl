# G6 accounting, atomic final-manifest promotion, then export.
# Usage: julia --project=. scripts/s90_g6_accounting.jl

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using SDSPNN
using Printf

const ARTIFACT_IDS = [
    "T1_R2_20260801", "F1_R2_20260801", "F2_R2_20260801",
    "F4_R2_20260801", "F5_R2_20260801", "AP1_R2_20260802",
]

function main()
    logpath, tee = setup_logging("s90_g6_accounting")
    db = open_artifact_db()
    try
        reports = Dict{String,Any}()
        # Validate every artifact before finalizing or exporting any artifact.
        # Accept `final` to resume safely if logging/export was interrupted after G6.
        for artifact_id in ARTIFACT_IDS
            artifact_status(db, artifact_id) in ("pending", "final") ||
                error("G6 requires pending or final artifact $artifact_id")
            report = validate_artifact!(db, artifact_id)
            reports[artifact_id] = report
            record_gate!(db, "G6", "$(artifact_id)_accounting", true, Dict(
                "N" => report.N, "expected_rows" => report.expected_rows,
                "actual_rows" => report.actual_rows, "status_counts" => report.status_counts,
            ))
            @printf(tee, "G6 validation PASS %-18s N=%d rows=%d\n",
                    artifact_id, report.N, report.actual_rows)
        end

        hashes = Dict{String,String}()
        for artifact_id in ARTIFACT_IDS
            hashes[artifact_id] = finalize_artifact!(db, artifact_id)
            println(tee, "Finalized $artifact_id hash=$(hashes[artifact_id])")
        end

        exports = Dict{String,Any}()
        for artifact_id in ARTIFACT_IDS
            exports[artifact_id] = export_artifact!(db, artifact_id)
            println(tee, "Exported $artifact_id only after G6: ", exports[artifact_id])
        end
        record_gate!(db, "G6", "all_artifacts_final", true,
                     Dict("artifact_hashes" => hashes, "exports" => exports))
    finally
        close_artifact_db(db)
        teardown_logging(tee, logpath)
    end
end

main()
