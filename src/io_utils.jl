struct TeeIO <: IO
    console::IO
    logfile::IO
end

function Base.unsafe_write(tee::TeeIO, p::Ptr{UInt8}, n::UInt)
    Base.unsafe_write(tee.console, p, n)
    Base.unsafe_write(tee.logfile, p, n)
    return n
end
Base.flush(tee::TeeIO) = (flush(tee.console); flush(tee.logfile))
Base.isopen(tee::TeeIO) = isopen(tee.console) && isopen(tee.logfile)
function Base.write(tee::TeeIO, byte::UInt8)
    write(tee.console, byte)
    write(tee.logfile, byte)
    return 1
end
function Base.write(tee::TeeIO, bytes::Vector{UInt8})
    write(tee.console, bytes)
    write(tee.logfile, bytes)
    return length(bytes)
end
function Base.write(tee::TeeIO, bytes::SubArray{UInt8,1})
    write(tee.console, bytes)
    write(tee.logfile, bytes)
    return length(bytes)
end

function utc_timestamp()
    return Dates.format(Dates.now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sssZ")
end

function setup_logging(script_name::AbstractString)
    logdir = joinpath(JCODE_ROOT, "results", "logs")
    mkpath(logdir)
    stamp = Dates.format(Dates.now(Dates.UTC), dateformat"yyyymmdd_HHMMSS")
    logpath = joinpath(logdir, "$(script_name)_$(stamp)_UTC.log")
    logfile = open(logpath, "w")
    tee = TeeIO(stdout, logfile)
    println(tee, "Log: $logpath")
    println(tee, "Started UTC: $(utc_timestamp())")
    return logpath, tee
end

"""Return the contained output path for a top-level PDF figure filename."""
function figure_output_path(filename::AbstractString)
    name = String(filename)
    (isempty(name) || occursin('/', name) || occursin('\\', name) || basename(name) != name) &&
        throw(ArgumentError("figure filename must be a basename"))
    endswith(lowercase(name), ".pdf") ||
        throw(ArgumentError("figure filename must end in .pdf"))
    return joinpath(JCODE_ROOT, "results", "figures", name)
end

function teardown_logging(tee::TeeIO, logpath::AbstractString)
    println(tee, "Finished UTC: $(utc_timestamp())")
    flush(tee)
    close(tee.logfile)
    println("Log saved to: $logpath")
    return nothing
end

function canonical_json(value)
    if value isa NamedTuple
        return canonical_json(Dict(string(k) => v for (k, v) in pairs(value)))
    elseif value isa AbstractDict
        entries = String[]
        for key in sort!(collect(keys(value)); by=string)
            push!(entries, String(JSON3.write(string(key))) * ":" * canonical_json(value[key]))
        end
        return "{" * join(entries, ",") * "}"
    elseif value isa Tuple || value isa AbstractVector
        return "[" * join((canonical_json(v) for v in value), ",") * "]"
    elseif value isa Symbol
        return String(JSON3.write(string(value)))
    elseif value isa AbstractFloat && !isfinite(value)
        return String(JSON3.write(string(value)))
    else
        return String(JSON3.write(value))
    end
end

_sha256_hex(value::AbstractString) = bytes2hex(SHA.sha256(Vector{UInt8}(codeunits(value))))
config_hash(config) = _sha256_hex(canonical_json(config))
selection_hash(query::AbstractString) = _sha256_hex(String(query))

function percentage_from_counts(n::Integer, N::Integer)
    0 <= n <= N || throw(ArgumentError("percentage counts require 0 <= n <= N"))
    N > 0 || throw(ArgumentError("percentage denominator must be positive"))
    return 100.0 * n / N
end

validate_percentage(n::Integer, N::Integer, displayed::Real; atol::Real=1e-12) =
    isapprox(Float64(displayed), percentage_from_counts(n, N); atol=atol, rtol=0)

const _RUN_STATUSES = Set([
    "completed", "solver_failure", "timeout", "evaluation_cap",
    "nonfinite_state", "projection_failure", "feasibility_failure",
])

"""Pure G6 validator used both by tests and the SQLite gate."""
function validate_accounting_rows(expected_runs, actual_runs,
                                  actual_row_counts::AbstractDict)
    expected_ids = [String(r.run_id) for r in expected_runs]
    actual_ids = [String(r.run_id) for r in actual_runs]
    length(unique(expected_ids)) == length(expected_ids) ||
        error("G6: duplicate expected run_id")
    length(unique(actual_ids)) == length(actual_ids) ||
        error("G6: duplicate actual run_id")
    Set(expected_ids) == Set(actual_ids) ||
        error("G6: expected and actual run sets differ")

    actual_keys = [(String(r.config_hash), String(r.replicate_key)) for r in actual_runs]
    length(unique(actual_keys)) == length(actual_keys) ||
        error("G6: duplicate config hash without a distinct replicate key")

    expected_by_id = Dict(String(r.run_id) => r for r in expected_runs)
    actual_by_id = Dict(String(r.run_id) => r for r in actual_runs)
    expected_rows = 0
    actual_rows = 0
    status_counts = Dict{String,Int}()
    for run_id in expected_ids
        exp = expected_by_id[run_id]
        act = actual_by_id[run_id]
        String(exp.config_hash) == String(act.config_hash) ||
            error("G6: config hash mismatch for $run_id")
        String(exp.replicate_key) == String(act.replicate_key) ||
            error("G6: replicate key mismatch for $run_id")
        status = String(act.run_status)
        status in _RUN_STATUSES || error("G6: invalid exclusive run_status $status")
        status_counts[status] = get(status_counts, status, 0) + 1
        expected_samples = Int(exp.expected_samples)
        actual_samples = get(actual_row_counts, run_id, -1)
        actual_samples == expected_samples ||
            error("G6: row count mismatch for $run_id: $actual_samples != $expected_samples")
        expected_rows += expected_samples
        actual_rows += actual_samples
    end
    sum(values(status_counts)) == length(expected_ids) ||
        error("G6: exclusive status counts do not sum to N")
    return (N=length(expected_ids), expected_rows=expected_rows,
            actual_rows=actual_rows, status_counts=status_counts)
end

function init_schema!(db::SQLite.DB)
    statements = [
        "PRAGMA foreign_keys = ON",
        """CREATE TABLE IF NOT EXISTS artifact_manifests (
            artifact_id TEXT PRIMARY KEY,
            deliverable TEXT NOT NULL,
            compute_script TEXT NOT NULL,
            status TEXT NOT NULL CHECK(status IN ('pending','final','failed')),
            created_utc TEXT NOT NULL,
            finalized_utc TEXT,
            config_json TEXT NOT NULL,
            config_hash TEXT NOT NULL,
            protocol_json TEXT NOT NULL,
            expected_runs INTEGER NOT NULL,
            expected_rows INTEGER NOT NULL,
            actual_rows INTEGER,
            final_manifest_hash TEXT,
            notes TEXT NOT NULL DEFAULT ''
        )""",
        """CREATE TABLE IF NOT EXISTS expected_runs (
            artifact_id TEXT NOT NULL,
            run_id TEXT NOT NULL,
            config_json TEXT NOT NULL,
            config_hash TEXT NOT NULL,
            replicate_key TEXT NOT NULL DEFAULT '',
            expected_samples INTEGER NOT NULL,
            PRIMARY KEY (artifact_id, run_id),
            FOREIGN KEY (artifact_id) REFERENCES artifact_manifests(artifact_id) ON DELETE CASCADE
        )""",
        """CREATE TABLE IF NOT EXISTS runs (
            artifact_id TEXT NOT NULL,
            run_id TEXT NOT NULL,
            config_json TEXT NOT NULL,
            config_hash TEXT NOT NULL,
            replicate_key TEXT NOT NULL DEFAULT '',
            run_status TEXT NOT NULL,
            terminal_success INTEGER NOT NULL,
            confirmed_hit INTEGER NOT NULL,
            hit_right_censored INTEGER NOT NULL,
            t_hit REAL,
            expected_samples INTEGER NOT NULL,
            actual_samples INTEGER NOT NULL,
            details_json TEXT NOT NULL,
            PRIMARY KEY (artifact_id, run_id),
            UNIQUE (artifact_id, config_hash, replicate_key),
            FOREIGN KEY (artifact_id) REFERENCES artifact_manifests(artifact_id) ON DELETE CASCADE
        )""",
        """CREATE TABLE IF NOT EXISTS artifact_rows (
            artifact_id TEXT NOT NULL,
            run_id TEXT NOT NULL,
            row_index INTEGER NOT NULL,
            payload_json TEXT NOT NULL,
            row_hash TEXT NOT NULL,
            PRIMARY KEY (artifact_id, run_id, row_index),
            FOREIGN KEY (artifact_id, run_id) REFERENCES runs(artifact_id, run_id) ON DELETE CASCADE
        )""",
        """CREATE TABLE IF NOT EXISTS displayed_percentages (
            artifact_id TEXT NOT NULL,
            label TEXT NOT NULL,
            numerator INTEGER NOT NULL,
            denominator INTEGER NOT NULL,
            percentage REAL NOT NULL,
            PRIMARY KEY (artifact_id, label),
            FOREIGN KEY (artifact_id) REFERENCES artifact_manifests(artifact_id) ON DELETE CASCADE
        )""",
        """CREATE TABLE IF NOT EXISTS plot_series (
            artifact_id TEXT NOT NULL,
            series_name TEXT NOT NULL,
            selection_query TEXT NOT NULL,
            selection_query_hash TEXT NOT NULL,
            final_manifest_hash TEXT NOT NULL DEFAULT '',
            PRIMARY KEY (artifact_id, series_name),
            FOREIGN KEY (artifact_id) REFERENCES artifact_manifests(artifact_id) ON DELETE CASCADE
        )""",
        """CREATE TABLE IF NOT EXISTS gate_results (
            gate_name TEXT NOT NULL,
            check_name TEXT NOT NULL,
            passed INTEGER NOT NULL,
            value_json TEXT NOT NULL,
            checked_utc TEXT NOT NULL,
            PRIMARY KEY (gate_name, check_name)
        )""",
    ]
    for statement in statements
        DBInterface.execute(db, statement)
    end
    return db
end

function open_artifact_db(path::AbstractString=joinpath(JCODE_ROOT, "results", "experiments.db"))
    path != ":memory:" && mkpath(dirname(path))
    db = SQLite.DB(path)
    init_schema!(db)
    return db
end

close_artifact_db(db::SQLite.DB) = SQLite.close(db)

function _materialize_rows(iterator)
    rows = NamedTuple[]
    for row in iterator
        names = Tuple(propertynames(row))
        values = Tuple(getproperty(row, name) for name in names)
        push!(rows, NamedTuple{names}(values))
    end
    return rows
end

function _clear_artifact!(db, artifact_id)
    for table in ("artifact_rows", "runs", "expected_runs", "displayed_percentages",
                  "plot_series", "artifact_manifests")
        DBInterface.execute(db, "DELETE FROM $table WHERE artifact_id = ?", (artifact_id,))
    end
end

function begin_artifact!(db::SQLite.DB, artifact_id::AbstractString,
                         deliverable::AbstractString, compute_script::AbstractString,
                         config; expected_runs::Integer, expected_rows::Integer,
                         protocol, force::Bool=false, notes::AbstractString="")
    existing = _materialize_rows(DBInterface.execute(db,
        "SELECT status FROM artifact_manifests WHERE artifact_id = ?", (artifact_id,)))
    if !isempty(existing)
        force || return false
        _clear_artifact!(db, artifact_id)
    end
    cfg_json = canonical_json(config)
    DBInterface.execute(db, """INSERT INTO artifact_manifests
        (artifact_id, deliverable, compute_script, status, created_utc,
         config_json, config_hash, protocol_json, expected_runs, expected_rows, notes)
        VALUES (?, ?, ?, 'pending', ?, ?, ?, ?, ?, ?, ?)""",
        (artifact_id, deliverable, compute_script, utc_timestamp(), cfg_json,
         config_hash(config), canonical_json(protocol), Int(expected_runs),
         Int(expected_rows), String(notes)))
    return true
end

function register_expected_run!(db::SQLite.DB, artifact_id::AbstractString,
                                run_id::AbstractString, config;
                                replicate_key::AbstractString="",
                                expected_samples::Integer)
    cfg_json = canonical_json(config)
    DBInterface.execute(db, """INSERT INTO expected_runs
        (artifact_id, run_id, config_json, config_hash, replicate_key, expected_samples)
        VALUES (?, ?, ?, ?, ?, ?)""",
        (artifact_id, run_id, cfg_json, config_hash(config), replicate_key, Int(expected_samples)))
    return nothing
end

function record_run!(db::SQLite.DB, artifact_id::AbstractString,
                     run_id::AbstractString, config;
                     replicate_key::AbstractString="", run_status::AbstractString,
                     terminal_success::Bool=false, confirmed_hit::Bool=false,
                     hit_right_censored::Bool=false, t_hit=nothing,
                     expected_samples::Integer, actual_samples::Integer,
                     details=Dict{String,Any}())
    run_status in _RUN_STATUSES || throw(ArgumentError("invalid run_status $run_status"))
    cfg_json = canonical_json(config)
    DBInterface.execute(db, """INSERT INTO runs
        (artifact_id, run_id, config_json, config_hash, replicate_key, run_status,
         terminal_success, confirmed_hit, hit_right_censored, t_hit,
         expected_samples, actual_samples, details_json)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        (artifact_id, run_id, cfg_json, config_hash(config), replicate_key, run_status,
         Int(terminal_success), Int(confirmed_hit), Int(hit_right_censored), t_hit,
         Int(expected_samples), Int(actual_samples), canonical_json(details)))
    return nothing
end

function insert_artifact_row!(db::SQLite.DB, artifact_id::AbstractString,
                              run_id::AbstractString, row_index::Integer, payload)
    payload_json = canonical_json(payload)
    DBInterface.execute(db, """INSERT INTO artifact_rows
        (artifact_id, run_id, row_index, payload_json, row_hash)
        VALUES (?, ?, ?, ?, ?)""",
        (artifact_id, run_id, Int(row_index), payload_json, _sha256_hex(payload_json)))
    return nothing
end

function record_percentage!(db::SQLite.DB, artifact_id::AbstractString,
                            label::AbstractString, n::Integer, N::Integer)
    percentage = percentage_from_counts(n, N)
    DBInterface.execute(db, """INSERT INTO displayed_percentages
        (artifact_id, label, numerator, denominator, percentage) VALUES (?, ?, ?, ?, ?)""",
        (artifact_id, label, Int(n), Int(N), percentage))
    return percentage
end

function register_plot_series!(db::SQLite.DB, artifact_id::AbstractString,
                               series_name::AbstractString, query::AbstractString)
    DBInterface.execute(db, """INSERT INTO plot_series
        (artifact_id, series_name, selection_query, selection_query_hash)
        VALUES (?, ?, ?, ?)""",
        (artifact_id, series_name, String(query), selection_hash(query)))
    return nothing
end

function record_gate!(db::SQLite.DB, gate_name::AbstractString,
                      check_name::AbstractString, passed::Bool, value)
    DBInterface.execute(db, """INSERT OR REPLACE INTO gate_results
        (gate_name, check_name, passed, value_json, checked_utc)
        VALUES (?, ?, ?, ?, ?)""",
        (gate_name, check_name, Int(passed), canonical_json(value), utc_timestamp()))
    return nothing
end

function artifact_status(db::SQLite.DB, artifact_id::AbstractString)
    rows = _materialize_rows(DBInterface.execute(db,
        "SELECT status FROM artifact_manifests WHERE artifact_id = ?", (artifact_id,)))
    isempty(rows) && return "missing"
    return String(rows[1].status)
end

function validate_artifact!(db::SQLite.DB, artifact_id::AbstractString)
    manifest_rows = _materialize_rows(DBInterface.execute(db,
        "SELECT expected_runs, expected_rows FROM artifact_manifests WHERE artifact_id = ?",
        (artifact_id,)))
    length(manifest_rows) == 1 || error("G6: missing manifest for $artifact_id")
    manifest = manifest_rows[1]

    expected = [(run_id=String(r.run_id), config_hash=String(r.config_hash),
                 replicate_key=String(r.replicate_key), expected_samples=Int(r.expected_samples))
                for r in DBInterface.execute(db,
                    "SELECT * FROM expected_runs WHERE artifact_id = ? ORDER BY run_id", (artifact_id,))]
    actual = [(run_id=String(r.run_id), config_hash=String(r.config_hash),
               replicate_key=String(r.replicate_key), run_status=String(r.run_status))
              for r in DBInterface.execute(db,
                  "SELECT * FROM runs WHERE artifact_id = ? ORDER BY run_id", (artifact_id,))]
    row_counts = Dict{String,Int}()
    for r in DBInterface.execute(db, """SELECT run_id, COUNT(*) AS n
        FROM artifact_rows WHERE artifact_id = ? GROUP BY run_id""", (artifact_id,))
        row_counts[String(r.run_id)] = Int(r.n)
    end
    report = validate_accounting_rows(expected, actual, row_counts)
    report.N == Int(manifest.expected_runs) || error("G6: manifest expected_runs mismatch")
    report.expected_rows == Int(manifest.expected_rows) || error("G6: manifest expected_rows mismatch")

    for p in DBInterface.execute(db,
        "SELECT numerator, denominator, percentage FROM displayed_percentages WHERE artifact_id = ?", (artifact_id,))
        validate_percentage(Int(p.numerator), Int(p.denominator), Float64(p.percentage)) ||
            error("G6: displayed percentage is not 100n/N")
    end
    series = _materialize_rows(DBInterface.execute(db,
        "SELECT selection_query, selection_query_hash FROM plot_series WHERE artifact_id = ?", (artifact_id,)))
    isempty(series) && error("G6: no stored selection query for $artifact_id")
    for s in series
        selection_hash(String(s.selection_query)) == String(s.selection_query_hash) ||
            error("G6: selection query hash mismatch")
    end
    return report
end

function finalize_artifact!(db::SQLite.DB, artifact_id::AbstractString)
    report = validate_artifact!(db, artifact_id)
    manifest = only(_materialize_rows(DBInterface.execute(db,
        "SELECT * FROM artifact_manifests WHERE artifact_id = ?", (artifact_id,))))
    run_hashes = [String(r.config_hash) for r in DBInterface.execute(db,
        "SELECT config_hash FROM runs WHERE artifact_id = ? ORDER BY run_id", (artifact_id,))]
    row_hashes = [String(r.row_hash) for r in DBInterface.execute(db,
        "SELECT row_hash FROM artifact_rows WHERE artifact_id = ? ORDER BY run_id,row_index", (artifact_id,))]
    query_hashes = [String(r.selection_query_hash) for r in DBInterface.execute(db,
        "SELECT selection_query_hash FROM plot_series WHERE artifact_id = ? ORDER BY series_name", (artifact_id,))]
    core = Dict(
        "artifact_id" => artifact_id,
        "deliverable" => String(manifest.deliverable),
        "compute_script" => String(manifest.compute_script),
        "created_utc" => String(manifest.created_utc),
        "config_hash" => String(manifest.config_hash),
        "expected_runs" => Int(manifest.expected_runs),
        "expected_rows" => Int(manifest.expected_rows),
        "actual_rows" => report.actual_rows,
        "run_config_hashes" => run_hashes,
        "row_hashes" => row_hashes,
        "selection_query_hashes" => query_hashes,
    )
    final_hash = config_hash(core)
    finalized = utc_timestamp()
    DBInterface.execute(db, """UPDATE artifact_manifests
        SET status='final', finalized_utc=?, actual_rows=?, final_manifest_hash=?
        WHERE artifact_id=?""", (finalized, report.actual_rows, final_hash, artifact_id))
    DBInterface.execute(db, """UPDATE plot_series SET final_manifest_hash=?
        WHERE artifact_id=?""", (final_hash, artifact_id))
    return final_hash
end

function export_artifact!(db::SQLite.DB, artifact_id::AbstractString;
                          export_dir::AbstractString=joinpath(JCODE_ROOT, "results", "exports"))
    artifact_status(db, artifact_id) == "final" ||
        error("G6: export is forbidden before finalization")
    mkpath(export_dir)
    rows = DataFrame(DBInterface.execute(db, """SELECT r.artifact_id, r.run_id,
        r.row_index, r.payload_json, r.row_hash, u.run_status, u.config_hash
        FROM artifact_rows r JOIN runs u USING (artifact_id, run_id)
        WHERE r.artifact_id = ? ORDER BY r.run_id, r.row_index""", (artifact_id,)))
    csv_path = joinpath(export_dir, "$(lowercase(artifact_id))_rows.csv")
    CSV.write(csv_path, rows)
    manifest = only(_materialize_rows(DBInterface.execute(db,
        "SELECT * FROM artifact_manifests WHERE artifact_id = ?", (artifact_id,))))
    manifest_dict = Dict(string(name) => getproperty(manifest, name) for name in propertynames(manifest))
    manifest_path = joinpath(export_dir, "$(lowercase(artifact_id))_manifest.json")
    open(manifest_path, "w") do io
        write(io, canonical_json(manifest_dict))
        write(io, "\n")
    end
    return (csv_path=csv_path, manifest_path=manifest_path)
end
