# Plot F2 from final artifact F2_R2_20260801. No ODE solve or derived value is computed.
# Mohammed runs: julia --project=. scripts/p40_f2_certified_radius.jl

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using SDSPNN
using Plots
using LaTeXStrings
using Printf
include(joinpath(@__DIR__, "plot_artifact_utils.jl"))

function main()
    artifact = load_final_artifact("F2_R2_20260801", "certified_radius")

    marker_indices = findall(==("radius_markers"), artifact.run_ids)
    markers = artifact.payloads[marker_indices]
    certified_marker = only(filter(p -> String(p.label) == "R_half", markers))
    analytic_marker = only(filter(
        p -> String(p.label) == "weighted_analytic_basin_boundary", markers))
    comparison_x = sqrt(certified_marker.radius * analytic_marker.radius)
    panel_a = plot([certified_marker.radius, analytic_marker.radius], [1.0, 1.0];
        xscale=:log10, ylims=(0.86, 1.14), yticks=false, linewidth=1.1,
        color=:gray, label=false, xlabel="weighted radius",
        title="Local certificate and analytic basin boundary", legend=:top)
    scatter!(panel_a, [certified_marker.radius], [1.0];
             markersize=7, marker=:circle, color=:blue, label=L"R_{1/2}")
    scatter!(panel_a, [analytic_marker.radius], [1.0];
             markersize=7, marker=:diamond, color=:red,
             label=L"r(\|x\|=0.8)=0.8\sqrt{2}")
    annotate!(panel_a, comparison_x, 0.91,
              text(@sprintf("%.4f decades", analytic_marker.decades), 9, :center))

    panel_b = plot(; yscale=:log10, xlabel="t", ylabel=L"r(x(t))",
                   title="Weighted-radius envelope", legend=:topright)
    styles = Dict(
        "certified_baseline" => (:blue, :solid, "certified start, 1e-10"),
        "certified_strict" => (:cyan, :dash, "certified start, 1e-12"),
        "outer_r0_0p5" => (:red, :solid, L"\mathrm{start}\ \|x(0)\|=0.5"),
    )
    for run_id in ("certified_baseline", "certified_strict", "outer_r0_0p5")
        indices = findall(==(run_id), artifact.run_ids)
        payloads = artifact.payloads[indices]
        color, linestyle, label = styles[run_id]
        plot!(panel_b, [p.t for p in payloads], [p.weighted_radius for p in payloads];
              color=color, linestyle=linestyle, linewidth=1.4, label=label)
    end
    baseline_indices = findall(==("certified_baseline"), artifact.run_ids)
    baseline = artifact.payloads[baseline_indices]
    plot!(panel_b, [p.t for p in baseline], [p.theorem_envelope for p in baseline];
          color=:black, linestyle=:dot, linewidth=1.5, label="theorem envelope")

    layout = @layout [comparison{0.27h}; envelope]
    figure = plot(panel_a, panel_b; layout=layout, size=(720, 780))
    output = figure_output_path("f2_certified_radius_r2.pdf")
    mkpath(dirname(output))
    savefig(figure, output)
    sidecar = write_figure_provenance(output, artifact)
    println("Wrote $output")
    println("Wrote $sidecar")
end

main()
