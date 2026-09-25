# Plot F1 from final artifact F1_R2_20260801. No ODE solve or projection occurs here.
# Mohammed runs: julia --project=. scripts/p30_f1_annulus.jl

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using SDSPNN
using Plots
using LaTeXStrings
include(joinpath(@__DIR__, "plot_artifact_utils.jl"))

function main()
    artifact = load_final_artifact("F1_R2_20260801", "phase_portrait")
    plt = plot(; aspect_ratio=:equal, xlabel=L"x_1", ylabel=L"x_2",
               title="Periodic annulus and analytic basin boundary", legend=:outerright)

    for run_id in ("disk_boundary", "basin_boundary")
        indices = findall(==(run_id), artifact.run_ids)
        payloads = artifact.payloads[indices]
        plot!(plt, [p.x1 for p in payloads], [p.x2 for p in payloads];
              color=:black, linestyle=:dash, linewidth=1.2,
              label=run_id == "disk_boundary" ? L"\|x\|=1.1" : L"\|x\|=0.8")
    end

    palette = [:blue, :cyan, :green, :orange, :red, :purple]
    for (color, run_id) in zip(palette, ("r_0p2", "r_0p5", "r_0p75", "r_0p85", "r_0p95", "r_1p05"))
        indices = findall(==(run_id), artifact.run_ids)
        payloads = artifact.payloads[indices]
        plot!(plt, [p.x1 for p in payloads], [p.x2 for p in payloads];
              color=color, linewidth=1.4, label=replace(run_id, "r_" => "s(0)=", "p" => "."))
    end

    output = figure_output_path("f1_annulus_r2.pdf")
    mkpath(dirname(output))
    savefig(plt, output)
    sidecar = write_figure_provenance(output, artifact)
    println("Wrote $output")
    println("Wrote $sidecar")
end

main()
