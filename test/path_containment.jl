using Test

const TEST_JCODE_ROOT = normpath(joinpath(@__DIR__, ".."))
push!(LOAD_PATH, joinpath(TEST_JCODE_ROOT, "src"))
using SDSPNN

const PLOT_SCRIPTS = [
    joinpath(TEST_JCODE_ROOT, "scripts", "p30_f1_annulus.jl"),
    joinpath(TEST_JCODE_ROOT, "scripts", "p40_f2_certified_radius.jl"),
    joinpath(TEST_JCODE_ROOT, "scripts", "p50_f4_integrability.jl"),
    joinpath(TEST_JCODE_ROOT, "scripts", "p60_f5_variation_band.jl"),
]
const PRODUCTION_JULIA_FILES = vcat(
    filter(path -> endswith(path, ".jl"), readdir(joinpath(TEST_JCODE_ROOT, "src"); join=true)),
    filter(path -> endswith(path, ".jl"), readdir(joinpath(TEST_JCODE_ROOT, "scripts"); join=true)),
)

@testset "jcode self-containment" begin
    @test figure_output_path("probe.pdf") ==
          joinpath(JCODE_ROOT, "results", "figures", "probe.pdf")
    @test_throws ArgumentError figure_output_path("../escape.pdf")
    @test_throws ArgumentError figure_output_path("nested/escape.pdf")
    @test_throws ArgumentError figure_output_path("probe.png")

    for script in PLOT_SCRIPTS
        source = read(script, String)
        @test !occursin(r"paper[\\/]imgs", source)
        @test !occursin("joinpath(JCODE_ROOT, \"..\"", source)
        @test occursin("figure_output_path", source)
    end

    f2_source = read(joinpath(TEST_JCODE_ROOT, "scripts", "p40_f2_certified_radius.jl"), String)
    @test occursin("scatter!(panel_a", f2_source)
    @test occursin("marker=:diamond", f2_source)
    @test occursin("weighted_analytic_basin_boundary", f2_source)
    @test occursin("R_{1/2}", f2_source)
    @test occursin("start", f2_source) && occursin("x(0)", f2_source)
    @test !occursin("R_{\\mathrm{clean}}", f2_source)
    @test !occursin("Euclidean r(0)=0.5", f2_source)

    f2_compute_source = read(joinpath(TEST_JCODE_ROOT, "scripts", "s40_f2_certified_radius.jl"), String)
    @test occursin("weighted_analytic_basin_boundary", f2_compute_source)
    @test occursin("weighted_basin_radius", f2_compute_source)

    f4_source = read(joinpath(TEST_JCODE_ROOT, "scripts", "p50_f4_integrability.jl"), String)
    @test occursin("color_by_residual=false", f4_source)
    @test occursin("numerical floor", f4_source)
    @test occursin("\"separable diagonal\\nexact zero\"", f4_source)
    @test occursin("\"bounded Hessian\\nnumerical floor\"", f4_source)
    @test occursin("\"manufactured nonzero\"", f4_source)
    @test occursin("\"rotating inverse metric\"", f4_source)
    @test !occursin("\"C1:", f4_source)
    @test !occursin("\"C2:", f4_source)
    @test !occursin("\"C3:", f4_source)
    @test !occursin("\"T:", f4_source)

    for source_path in PRODUCTION_JULIA_FILES
        source = read(source_path, String)
        @test !occursin(r"paper[\\/]", source)
        @test !occursin(r"channels[\\/]", source)
        @test !occursin(r"joinpath\(JCODE_ROOT,\s*\"\.\.\"", source)
    end
end
