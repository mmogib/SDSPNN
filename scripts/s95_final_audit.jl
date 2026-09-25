# Read-only final audit. It validates final artifacts and prints referee-facing evidence.

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using SDSPNN
using DBInterface
using DataFrames
using JSON3

const IDS = [
    "T1_R2_20260801", "F1_R2_20260801", "F2_R2_20260801",
    "F4_R2_20260801", "F5_R2_20260801", "AP1_R2_20260802",
]

function json_rows(db, query, params=())
    frame = DataFrame(DBInterface.execute(db, query, params))
    return [JSON3.read(String(value)) for value in frame.payload_json]
end

function tau_star_consistency(db; atol=1e-14)
    t1 = only(json_rows(db, """SELECT payload_json FROM artifact_rows
        WHERE artifact_id='T1_R2_20260801' AND run_id='B_tau_star'"""))
    f5 = only(json_rows(db, """SELECT payload_json FROM artifact_rows
        WHERE artifact_id='F5_R2_20260801' AND run_id='root_summary'"""))
    t1_tau = Float64(t1.tau)
    f5_tau = Float64(f5.tau)
    difference = abs(t1_tau - f5_tau)
    isapprox(t1_tau, f5_tau; atol=atol, rtol=0.0) ||
        error("T1/F5 tau-star mismatch: T1=$t1_tau, F5=$f5_tau, atol=$atol")
    return Dict(
        "T1_tau_star" => t1_tau,
        "F5_tau_star" => f5_tau,
        "absolute_difference" => difference,
        "tolerance" => atol,
    )
end

function main()
    db = open_artifact_db()
    try
        println("=== FINAL MANIFESTS ===")
        for artifact_id in IDS
            report = validate_artifact!(db, artifact_id)
            artifact_status(db, artifact_id) == "final" || error("$artifact_id is not final")
            manifest = DataFrame(DBInterface.execute(db, """SELECT actual_rows,
                final_manifest_hash FROM artifact_manifests WHERE artifact_id=?""", (artifact_id,)))
            series = DataFrame(DBInterface.execute(db, """SELECT selection_query,
                selection_query_hash, final_manifest_hash FROM plot_series WHERE artifact_id=?""", (artifact_id,)))
            for row in eachrow(series)
                selection_hash(String(row.selection_query)) == String(row.selection_query_hash) ||
                    error("selection hash mismatch for $artifact_id")
                String(row.final_manifest_hash) == String(manifest.final_manifest_hash[1]) ||
                    error("series manifest hash mismatch for $artifact_id")
            end
            println(canonical_json(Dict(
                "artifact_id" => artifact_id, "N" => report.N,
                "rows" => report.actual_rows,
                "manifest_hash" => String(manifest.final_manifest_hash[1]),
                "selection_hashes" => String.(series.selection_query_hash),
            )))
        end

        println("=== CROSS-ARTIFACT CONSISTENCY ===")
        println(canonical_json(tau_star_consistency(db)))

        println("=== GATE EVIDENCE ===")
        gates = DataFrame(DBInterface.execute(db, """SELECT gate_name, check_name,
            passed, value_json FROM gate_results ORDER BY gate_name,check_name"""))
        for row in eachrow(gates)
            println(canonical_json(Dict("gate" => row.gate_name, "check" => row.check_name,
                                        "passed" => row.passed, "value" => JSON3.read(String(row.value_json)))))
        end

        println("=== F1 RUN DETAILS ===")
        f1 = DataFrame(DBInterface.execute(db, """SELECT run_id, details_json FROM runs
            WHERE artifact_id='F1_R2_20260801' ORDER BY run_id"""))
        for row in eachrow(f1)
            println(canonical_json(Dict("run_id" => row.run_id,
                                        "details" => JSON3.read(String(row.details_json)))))
        end

        println("=== F2 RUN DETAILS AND MARKERS ===")
        f2 = DataFrame(DBInterface.execute(db, """SELECT run_id, details_json FROM runs
            WHERE artifact_id='F2_R2_20260801' ORDER BY run_id"""))
        for row in eachrow(f2)
            println(canonical_json(Dict("run_id" => row.run_id,
                                        "details" => JSON3.read(String(row.details_json)))))
        end
        for payload in json_rows(db, """SELECT payload_json FROM artifact_rows
            WHERE artifact_id='F2_R2_20260801' AND run_id='radius_markers' ORDER BY row_index""")
            println(canonical_json(payload))
        end

        println("=== F5 ROOT SUMMARY ===")
        for payload in json_rows(db, """SELECT payload_json FROM artifact_rows
            WHERE artifact_id='F5_R2_20260801' AND run_id='root_summary'""")
            println(canonical_json(payload))
        end

        println("=== T1 ROWS ===")
        for payload in json_rows(db, """SELECT payload_json FROM artifact_rows
            WHERE artifact_id='T1_R2_20260801' ORDER BY run_id""")
            println(canonical_json(payload))
        end

        println("=== AP1 STAGE-1 ROWS ===")
        for payload in json_rows(db, """SELECT payload_json FROM artifact_rows
            WHERE artifact_id='AP1_R2_20260802' ORDER BY
            CASE run_id WHEN 'tau_1_primary' THEN 1 ELSE 2 END,row_index""")
            println(canonical_json(payload))
        end
    finally
        close_artifact_db(db)
    end
end

main()
