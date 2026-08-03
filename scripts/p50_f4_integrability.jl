# Plot F4 from final artifact F4_R2_20260801. No differentiation occurs here.
# Mohammed runs: julia --project=. scripts/p50_f4_integrability.jl

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using SDSPNN
using Plots
using LaTeXStrings
include(joinpath(@__DIR__, "plot_artifact_utils.jl"))

function case_panel(artifact, run_id, title_text; color_by_residual=true)
    indices = findall(==(run_id), artifact.run_ids)
    payloads = artifact.payloads[indices]
    common = (; markersize=4, markerstrokewidth=0,
              xlabel=L"x_1", ylabel=L"x_2", title=title_text,
              aspect_ratio=:equal, legend=false)
    if color_by_residual
        return scatter([p.x1 for p in payloads], [p.x2 for p in payloads];
                       marker_z=[p.residual_norm for p in payloads], color=:viridis,
                       colorbar=true, common...)
    end
    return scatter([p.x1 for p in payloads], [p.x2 for p in payloads];
                   color=:gray, colorbar=false, common...)
end

function main()
    artifact = load_final_artifact("F4_R2_20260801", "integrability_residual")
    panels = [
        case_panel(artifact, "C1_separable", "separable diagonal\nexact zero";
                   color_by_residual=false),
        case_panel(artifact, "C2_hessian", "bounded Hessian\nnumerical floor";
                   color_by_residual=false),
        case_panel(artifact, "C3_manufactured", "manufactured nonzero"),
        case_panel(artifact, "T_rotating_annulus", "rotating inverse metric"),
    ]
    figure = plot(panels...; layout=(2, 2), size=(900, 760),
                  plot_title=L"\|\mathcal{R}(x)\|")
    output = figure_output_path("f4_integrability_residual_r2.pdf")
    mkpath(dirname(output))
    savefig(figure, output)
    sidecar = write_figure_provenance(output, artifact)
    println("Wrote $output")
    println("Wrote $sidecar")
end

main()
