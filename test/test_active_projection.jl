using Test
using LinearAlgebra

push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..", "src")))
using SDSPNN

@testset "active projection G0 formulas" begin
    expected_cg = Dict(
        0.1 => 0.012989725292288541,
        1.0 => -7.579225540623604,
    )
    for tau in (0.1, 1.0)
        constants = active_projection_constants(tau)
        L_formula = inv(1 - ((1 + sqrt(5.0)) / 4) * tau)
        KM_formula = 3.571 * tau
        rho_formula = sqrt(1 - 2ALPHA / L_formula + ALPHA^2)
        cg_formula = 1 - rho_formula - (ALPHA / 2) * L_formula^2 * KM_formula * 3.1

        @test constants.mu == 1.0
        @test constants.K == 1.0
        @test constants.Fmax == 3.1
        @test constants.L == L_formula
        @test constants.KM_bound == KM_formula
        @test constants.rho == rho_formula
        @test constants.cg == cg_formula
        @test constants.cg ≈ expected_cg[tau] atol=5e-15
        @test constants.displacement_bound == ALPHA * L_formula * 3.1
    end
end

@testset "stage 1 compute script boundary" begin
    script = normpath(joinpath(@__DIR__, "..", "scripts", "s70_active_projection.jl"))
    @test isfile(script)
    source = read(script, String)
    @test occursin("function main()", source)
    @test occursin("main()", source)
    @test !occursin("solve_flow", source)
    @test !occursin("ODEProblem", source)
    @test !occursin("using Plots", source)
    @test !occursin("using Random", source)
    @test !occursin("rand(", source)
    @test !occursin("paper/", source)
    @test !occursin("channels/", source)
    @test occursin("rhs_exact ratio_exact rhs_L ratio_L ratio_L2", source)
end

@testset "joint projection bounds on active pairs" begin
    xbar = [DISK_RADIUS, 0.0]
    for tau in (0.1, 1.0)
        constants = active_projection_constants(tau)
        anchor = active_projection_point(xbar, tau)
        for phi in active_sample_angles(tau, 8)
            x = DISK_RADIUS .* [cos(phi), sin(phi)]
            point = active_projection_point(x, tau)
            bounds = active_projection_joint_bounds(point, anchor, constants.L)

            @test bounds.argument_gap > 0
            @test bounds.metric_gap > 0
            @test bounds.anchor == anchor.projection_displacement
            @test bounds.lhs <= bounds.rhs_exact + 2e-12
            @test bounds.lhs <= bounds.rhs_L + 2e-12
            @test bounds.lhs <= bounds.rhs_L2 + 2e-12
            @test 0 <= bounds.ratio_exact <= 1 + 2e-12
            @test 0 <= bounds.ratio_L <= 1 + 2e-12
            @test 0 <= bounds.ratio_L2 <= bounds.ratio_L + 2e-12
            @test bounds.rhs_L2 >= bounds.rhs_L
        end
    end
end

@testset "deterministic per-parameter arc samples" begin
    for tau in (0.1, 1.0)
        geometry = active_boundary_geometry(tau)
        angles = active_sample_angles(tau, 8)
        @test length(angles) == 8
        @test issorted(angles)
        @test any(<(0), angles)
        @test any(>(0), angles)
        @test all(geometry.arc_lo < phi < geometry.arc_hi for phi in angles)
        @test all(active_projection_point(
            DISK_RADIUS .* [cos(phi), sin(phi)], tau).projection.active for phi in angles)
        expected = [geometry.arc_lo + j * (geometry.arc_hi - geometry.arc_lo) / 9
                    for j in 1:8]
        @test angles ≈ expected atol=2e-15
    end
end

@testset "radial-clamp failure signatures" begin
    expected = Dict(
        0.1 => (-0.0236716758769695, 0.026038235518088218),
        1.0 => (-0.3859513874437237, 0.42191643637545007),
    )
    xbar = [DISK_RADIUS, 0.0]
    for tau in (0.1, 1.0)
        signature = active_clamp_signature(tau)
        phi, distance = expected[tau]
        @test signature.phi ≈ phi atol=5e-13
        @test signature.distance_from_solution ≈ distance atol=5e-13
        x = DISK_RADIUS .* [cos(signature.phi), sin(signature.phi)]
        point = active_projection_point(x, tau)
        @test point.projection.active
        @test point.clamp ≈ x atol=5e-13
        @test norm(x - xbar) ≈ signature.distance_from_solution atol=5e-13
    end
end

@testset "exact active boundary geometry" begin
    expected = Dict(
        0.1 => (-0.9723992517440125, 1.076470791733393, 1.0526288949968579),
        1.0 => (-0.5937570949301803, 1.5230104348896472, 1.052363002380429),
    )
    for tau in (0.1, 1.0)
        geometry = active_boundary_geometry(tau)
        lo, hi, radial = expected[tau]
        @test geometry.arc_lo ≈ lo atol=5e-12
        @test geometry.arc_hi ≈ hi atol=5e-12
        @test geometry.radial_threshold ≈ radial atol=5e-12
        @test geometry.arc_lo < 0 < geometry.arc_hi

        for phi in (geometry.arc_lo, geometry.arc_hi)
            x = DISK_RADIUS .* [cos(phi), sin(phi)]
            @test abs(active_projection_point(x, tau).active_margin) <= 2e-12
        end
        @test active_projection_point([DISK_RADIUS, 0.0], tau).active_margin > 0
    end
end

@testset "active projection equilibrium anchors" begin
    xbar = [DISK_RADIUS, 0.0]
    expected_clamp_gap = Dict(
        0.1 => 0.0021615689223323667,
        1.0 => 0.021612591227205737,
    )
    for tau in (0.1, 1.0)
        point = active_projection_point(xbar, tau)
        @test point.projection.active
        @test point.projection.nu ≈ 9 / 220 atol=5e-13
        @test point.projection.z ≈ xbar atol=5e-13
        @test point.map_displacement <= 5e-13
        @test point.projection_displacement > 0
        @test point.active_margin > 0
        @test point.spectral.active
        @test point.spectral.nu ≈ point.projection.nu atol=2e-13
        @test point.spectral.z ≈ point.projection.z atol=2e-13
        @test point.clamp_gap ≈ expected_clamp_gap[tau] atol=5e-14
        @test point.wrong_metric_map_gap > 1e-4
        @test maximum((point.projection.primal_residual,
                       point.projection.dual_residual,
                       point.projection.stationarity_residual,
                       point.projection.complementarity_residual)) <= PROJECTION_GATE_TOL
    end
end

@testset "independent spectral projection check" begin
    y = [1.8, 0.9]
    Q = [2.0 0.4; 0.4 0.8]
    reference = weighted_disk_projection(y, Q, DISK_RADIUS)
    spectral = weighted_disk_projection_spectral(y, Q, DISK_RADIUS)

    @test spectral.active
    @test spectral.nu ≈ reference.nu atol=2e-13
    @test spectral.z ≈ reference.z atol=2e-13
    @test norm(spectral.z) ≈ DISK_RADIUS atol=2e-13

    inside = weighted_disk_projection_spectral([0.2, -0.1], Q, DISK_RADIUS)
    @test !inside.active
    @test inside.nu == 0.0
    @test inside.z == [0.2, -0.1]
end
