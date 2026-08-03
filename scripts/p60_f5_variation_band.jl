# Plot F5 from final artifact F5_R2_20260801. All curves are stored analytic values.
# Mohammed runs: julia --project=. scripts/p60_f5_variation_band.jl

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using SDSPNN
using Plots
using LaTeXStrings
include(joinpath(@__DIR__, "plot_artifact_utils.jl"))

function main()
    artifact = load_final_artifact("F5_R2_20260801", "variation_band")
    grid = artifact.payloads[findall(==("plot_grid"), artifact.run_ids)]
    summary = only(artifact.payloads[findall(==("root_summary"), artifact.run_ids)])

    tau = [p.tau for p in grid]
    panel_top = plot(tau, [p.cg for p in grid]; color=:blue, linewidth=1.6,
                     xlabel=L"\tau", ylabel=L"c_g(\tau)", label=L"c_g",
                     title="Certified variation band")
    hline!(panel_top, [0.0]; color=:black, linestyle=:dash, label="zero")
    vline!(panel_top, [summary.tau]; color=:red, linestyle=:dot,
           label="certified endpoint")

    panel_bottom = plot(tau, [p.analytic_decay_rate for p in grid];
                        color=:green, linewidth=1.6, xlabel=L"\tau",
                        ylabel="analytic radial decay rate", label=L"\alpha(1-\tau)",
                        title="True band and non-integrability")
    right_axis = twinx(panel_bottom)
    plot!(right_axis, tau, [p.residual_norm_r_0p95 for p in grid];
          color=:purple, linestyle=:dash, linewidth=1.4,
          ylabel=L"\|\mathcal{R}_\tau\|\ (r=0.95)", label="residual norm")

    figure = plot(panel_top, panel_bottom; layout=(2, 1), size=(760, 800))
    output = figure_output_path("f5_variation_band_r2.pdf")
    mkpath(dirname(output))
    savefig(figure, output)
    sidecar = write_figure_provenance(output, artifact)
    println("Wrote $output")
    println("Wrote $sidecar")
end

main()
