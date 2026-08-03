# Shared read-only loader for pNN scripts. It never computes numerical quantities.

using SDSPNN
using SQLite
using DBInterface
using DataFrames
using JSON3

function load_final_artifact(artifact_id::AbstractString, series_name::AbstractString)
    db = open_artifact_db()
    try
        manifest_df = DataFrame(DBInterface.execute(db,
            "SELECT status, final_manifest_hash FROM artifact_manifests WHERE artifact_id = ?",
            (artifact_id,)))
        nrow(manifest_df) == 1 || error("missing manifest for $artifact_id")
        manifest_df.status[1] == "final" || error("artifact $artifact_id is not final")
        final_hash = String(manifest_df.final_manifest_hash[1])

        series_df = DataFrame(DBInterface.execute(db, """SELECT selection_query,
            selection_query_hash, final_manifest_hash FROM plot_series
            WHERE artifact_id = ? AND series_name = ?""", (artifact_id, series_name)))
        nrow(series_df) == 1 || error("missing stored plot series $series_name")
        query = String(series_df.selection_query[1])
        query_hash = String(series_df.selection_query_hash[1])
        selection_hash(query) == query_hash || error("selection query hash mismatch")
        String(series_df.final_manifest_hash[1]) == final_hash ||
            error("plot series does not name the final manifest hash")

        rows = DataFrame(DBInterface.execute(db, query))
        payloads = [JSON3.read(String(payload)) for payload in rows.payload_json]
        return (
            run_ids=String.(rows.run_id), row_indices=Int.(rows.row_index),
            payloads=payloads, artifact_id=String(artifact_id),
            final_manifest_hash=final_hash, selection_query=query,
            selection_query_hash=query_hash,
        )
    finally
        close_artifact_db(db)
    end
end

function write_figure_provenance(pdf_path::AbstractString, artifact)
    sidecar = replace(pdf_path, r"\.pdf$" => "_provenance.txt")
    open(sidecar, "w") do io
        println(io, "artifact_id: ", artifact.artifact_id)
        println(io, "final_manifest_hash: ", artifact.final_manifest_hash)
        println(io, "selection_query_hash: ", artifact.selection_query_hash)
        println(io, "selection_query: ", artifact.selection_query)
    end
    return sidecar
end
