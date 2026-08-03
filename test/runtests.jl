using Test
using LinearAlgebra

push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..", "src")))
using SDSPNN

include("path_containment.jl")
include("test_active_projection.jl")

const G1_POINTS = [
    [0.0, 0.0], [0.2, 0.0], [0.4, 0.3], [0.7999, 0.0],
    [0.8, 0.0], [0.8001, 0.0], [0.95, 0.0], [0.6, 0.7],
]

@testset "G1 static unit tests" begin
    @testset "operator and declared constants" begin
        A = operator_matrix()
        @test A == [1.0 -2.0; 2.0 1.0]
        @test minimum(eigvals(Symmetric((A + A') / 2))) == 1.0
        @test opnorm(A) ≈ sqrt(5.0) atol=2e-15

        c1 = analytic_constants(1.0)
        @test c1.mu == 1.0
        @test c1.K ≈ sqrt(5.0) atol=2e-15
        @test c1.Fmax ≈ 1.1sqrt(5.0) atol=2e-15
        @test c1.L ≈ 3.0 + sqrt(5.0) atol=2e-14
        @test c1.KM_bound == 3.571
        @test c1.window ≈ 2 / ((3 + sqrt(5.0)) * 5) atol=2e-15
        @test c1.alpha < c1.window
        @test c1.rho ≈ 0.996695389493 atol=5e-13
        @test c1.CL ≈ 42.71491868 atol=2e-8
        @test c1.R_clean ≈ 3.8682158467e-5 rtol=2e-9
        @test c1.eta_R ≈ 1.65230525328e-3 rtol=2e-9
    end

    @testset "tau metric regularity and bounds" begin
        for tau in (0.0, 0.05, 0.1, 0.139, 0.14, 0.2, 1.0), x in G1_POINTS
            M = metric_tau(x, tau)
            @test M ≈ M' atol=2e-15
            lo, hi, L = metric_bounds_tau(tau)
            ev = eigvals(Symmetric(M))
            @test minimum(ev) >= lo - 2e-12
            @test maximum(ev) <= hi + 2e-12
            @test L >= max(maximum(ev), inv(minimum(ev))) - 2e-12
        end

        @test metric_tau([0.0, 0.0], 1.0) == Matrix{Float64}(I, 2, 2)
        @test chi_radius(0.0) == 0.0
        @test chi_radius(0.8) == 1.0
        @test chi_radius(1.0) == 1.0

        h = 1e-5
        left1 = (chi_radius(0.8) - chi_radius(0.8 - h)) / h
        right1 = (chi_radius(0.8 + h) - chi_radius(0.8)) / h
        left2 = (chi_radius(0.8) - 2chi_radius(0.8 - h) + chi_radius(0.8 - 2h)) / h^2
        right2 = (chi_radius(0.8 + 2h) - 2chi_radius(0.8 + h) + chi_radius(0.8)) / h^2
        @test abs(left1 - right1) <= 2e-7
        # One-sided second differences converge only linearly at this patched join.
        @test abs(left2 - right2) <= 2e-3

        eps_fd = 1e-7
        directions = ([1.0, 0.0], [0.0, 1.0], normalize([1.0, 1.0]))
        for tau in (0.1, 1.0), x in G1_POINTS[2:end], d in directions
            if norm(x + eps_fd*d) <= DISK_RADIUS
                quotient = opnorm(metric_tau(x + eps_fd*d, tau) - metric_tau(x, tau)) / eps_fd
                @test quotient <= 3.571tau + 3e-4
            end
        end
    end

    @testset "analytic tau-family values" begin
        expected = Dict(
            0.05 => (1.04216, 0.17855, 0.957364, 0.03071),
            0.10 => (1.08802, 0.35710, 0.959474, 0.01453),
            0.139 => (1.12670, 0.49637, 0.961117, 0.00014),
            0.140 => (1.12773, 0.49994, 0.961159, -0.00026),
            0.20 => (1.19304, 0.71420, 0.963681, -0.02619),
        )
        for (tau, values) in expected
            c = analytic_constants(tau)
            @test c.L ≈ values[1] atol=5e-6
            # The comparison values in revision 2 are displayed to five decimals.
            @test c.KM_bound ≈ values[2] atol=1.1e-6
            @test c.rho ≈ values[3] atol=6e-7
            @test c.cg ≈ values[4] atol=6e-6
        end
        @test cg_tau(0.139) > 0
        @test cg_tau(0.140) < 0
        root = cg_root(0.139, 0.140; tol=1e-14)
        @test 0.1393 < root < 0.1394
        @test residual_norm_tau(0.1, 0.95) ≈ 0.1242087475 atol=2e-10
    end

    @testset "control metrics" begin
        controls = Dict(
            :C => (0.5, 1.35, 2.0, 0.600),
            :D => (1.0, 5.0, 5.0, 16.0),
            :E => (1.0, 1.25, 1.25, 0.5),
        )
        for (name, expected) in controls
            constants = control_constants(name)
            @test constants.spectral_lower == expected[1]
            @test constants.spectral_upper == expected[2]
            @test constants.L == expected[3]
            @test constants.KM_bound == expected[4]
            @test constants.spectral_lower >= inv(constants.L)
            @test constants.spectral_upper <= constants.L
        end
        for x in ([0.0, 0.0], [0.5, 0.25], [1.0, 1.0])
            C = metric_hessian_control(x)
            @test C ≈ C' atol=2e-15
            @test minimum(eigvals(Symmetric(C))) >= 0.5 - 2e-15
            @test maximum(eigvals(Symmetric(C))) <= 1.35 + 2e-15

            D = metric_diagonal_control(x)
            @test D ≈ D' atol=2e-15
            @test minimum(eigvals(Symmetric(D))) >= 1.0 - 2e-15
            @test maximum(eigvals(Symmetric(D))) <= 5.0 + 2e-15

            E = metric_manufactured_control(x)
            @test E ≈ E' atol=2e-15
            @test minimum(eigvals(Symmetric(E))) >= 1.0 - 2e-15
            @test maximum(eigvals(Symmetric(E))) <= 1.25 + 2e-15
        end
    end

    @testset "weighted disk projection KKT" begin
        Q = [2.0 0.4; 0.4 0.8]
        inside = weighted_disk_projection([0.2, -0.1], Q, DISK_RADIUS)
        @test !inside.active
        @test inside.z ≈ [0.2, -0.1] atol=2e-15
        @test inside.nu == 0.0

        outside = weighted_disk_projection([1.8, 0.9], Q, DISK_RADIUS)
        @test outside.active
        @test outside.bracket_initial_lo == 0.0
        @test outside.bracket_initial_hi ≈ norm(Q * [1.8, 0.9]) / DISK_RADIUS atol=2e-15
        @test outside.primal_residual <= PROJECTION_GATE_TOL
        @test outside.dual_residual <= PROJECTION_GATE_TOL
        @test outside.stationarity_residual <= PROJECTION_GATE_TOL
        @test outside.complementarity_residual <= PROJECTION_GATE_TOL
        @test norm(outside.z) ≈ DISK_RADIUS atol=2e-12
        @test abs(det(hcat(outside.z, [1.8, 0.9]))) > 1e-3

        radial = weighted_disk_projection([2.0, -1.0], Matrix{Float64}(I, 2, 2), DISK_RADIUS)
        @test radial.z ≈ DISK_RADIUS * [2.0, -1.0] / sqrt(5.0) atol=2e-12
    end
end

@testset "G2 integrator verification" begin
    x0 = [0.95, 0.0]
    horizon = 2 * ANNULUS_PERIOD
    save_grid = range(0.0, horizon; length=1001)
    errors = Dict{Float64,Tuple{Float64,Float64}}()
    for tol in (1e-10, 1e-12)
        sol = solve_flow(x0, (0.0, horizon); tau=1.0, abstol=tol, reltol=tol,
                         save_grid=save_grid, dtmax=0.25, maxiters=10_000_000)
        @test sol.run_status == "completed"
        radial_error = maximum(abs.(vec(sqrt.(sum(abs2, sol.xs; dims=1))) .- 0.95))
        phase_error = 0.0
        for (j, t) in enumerate(sol.ts)
            exact = exact_annular_state(x0, t)
            phase_error = max(phase_error,
                abs(atan(det(hcat(exact, sol.xs[:, j])), dot(exact, sol.xs[:, j]))))
        end
        errors[tol] = (radial_error, phase_error)
    end
    @test errors[1e-10][1] <= 2e-8
    @test errors[1e-10][2] <= 2e-8
    @test errors[1e-12][1] <= 2e-10
    @test errors[1e-12][2] <= 2e-10
    @test errors[1e-12][1] <= errors[1e-10][1] + 2e-13
    @test errors[1e-12][2] <= errors[1e-10][2] + 2e-13
end

@testset "G3 residual validation" begin
    points_square = ([0.2, 0.2], [0.5, 0.7], [0.9, 1.0])
    for x in points_square
        c1 = richardson_residual(metric_diagonal_control, x; h=1e-4)
        @test norm(c1.extrapolated) <= 1e-8

        c3 = richardson_residual(metric_manufactured_control, x; h=1e-4)
        @test c3.extrapolated ≈ [0.0, x[2] / 2] atol=2e-10
    end

    for x in ([0.2, 0.1], [0.6, -0.3], [1.0, 0.2])
        c2 = richardson_residual(metric_hessian_control, x; h=1e-4)
        @test norm(c2.extrapolated) <= 1e-8
    end

    angles = (0.0, 0.4, 0.9, 1.3, 2.0)
    refs = (2.353, 4.000, 5.148, 5.163, 3.301)
    for (theta, ref1) in zip(angles, refs)
        x = 0.85 .* [cos(theta), sin(theta)]
        target = richardson_residual(target_inverse_metric, x; h=1e-4)
        exact = analytic_target_residual(x)
        @test target.extrapolated ≈ exact atol=2e-8
        # The spec's three-decimal spot values differ from its closed form by up to 1.13e-3.
        @test exact[1] ≈ ref1 atol=1.2e-3
        @test norm(target.extrapolated) ≈ 2sqrt(5.0) / 0.85 atol=2e-8
        @test 1.7 <= target.observed_order <= 2.3
        @test target.truncation_estimate <= 2e-8
    end
end

@testset "G5 classifier validation" begin
    cfg = classifier_config()
    @test cfg.convergence_ratio == 0.1
    @test cfg.convergence_slope == -1e-4
    @test cfg.periodic_relative_span == 1e-6
    @test cfg.periodic_min_turns == 0.9

    ts_c = collect(range(0.0, 100.0; length=1001))
    xs_c = hcat(([exp(-0.05t) * cos(0.1t), exp(-0.05t) * sin(0.1t)] for t in ts_c)...)
    @test classify_trajectory(ts_c, xs_c; config=cfg).label == "convergent"

    ts_p = collect(range(0.0, ANNULUS_PERIOD; length=1001))
    xs_p = hcat((exact_annular_state([0.95, 0.0], t) for t in ts_p)...)
    @test classify_trajectory(ts_p, xs_p; config=cfg).label == "periodic"

    ts_r = collect(range(0.0, 400.0; length=1001))
    rate = ALPHA * (1 - 0.99)
    xs_r = hcat(([0.95exp(-rate*t) * cos(0.025t),
                   0.95exp(-rate*t) * sin(0.025t)] for t in ts_r)...)
    @test classify_trajectory(ts_r, xs_r; config=cfg).label == "right_censored"
end

@testset "G6 accounting primitives" begin
    @test config_hash(Dict("b" => 2, "a" => [1, 3])) ==
          config_hash(Dict("a" => [1, 3], "b" => 2))
    @test selection_hash("SELECT * FROM rows ORDER BY row_index") ==
          selection_hash("SELECT * FROM rows ORDER BY row_index")

    expected = [
        (run_id="r1", config_hash="h1", replicate_key="", expected_samples=3),
        (run_id="r2", config_hash="h2", replicate_key="", expected_samples=2),
    ]
    actual = [
        (run_id="r1", config_hash="h1", replicate_key="", run_status="completed"),
        (run_id="r2", config_hash="h2", replicate_key="", run_status="completed"),
    ]
    report = validate_accounting_rows(expected, actual, Dict("r1" => 3, "r2" => 2))
    @test report.N == 2
    @test report.expected_rows == 5
    @test report.actual_rows == 5
    @test percentage_from_counts(2, 5) == 40.0
    @test validate_percentage(2, 5, 40.0)
    @test !validate_percentage(2, 5, 39.9)

    @test_throws ErrorException validate_accounting_rows(
        expected, actual[1:1], Dict("r1" => 3))
    duplicate = [actual[1], merge(actual[2], (config_hash="h1",))]
    @test_throws ErrorException validate_accounting_rows(
        expected, duplicate, Dict("r1" => 3, "r2" => 2))
    @test_throws ErrorException validate_accounting_rows(
        expected, actual, Dict("r1" => 3, "r2" => 1))

    db = open_artifact_db(":memory:")
    @test begin_artifact!(db, "TEST", "T0", "test/runtests.jl",
                          Dict("grid" => [1, 2]); expected_runs=1,
                          expected_rows=2, protocol=Dict("kind" => "fixture"), force=true)
    register_expected_run!(db, "TEST", "r1", Dict("case" => 1); expected_samples=2)
    record_run!(db, "TEST", "r1", Dict("case" => 1);
                run_status="completed", expected_samples=2, actual_samples=2)
    insert_artifact_row!(db, "TEST", "r1", 1, Dict("value" => 10))
    insert_artifact_row!(db, "TEST", "r1", 2, Dict("value" => 20))
    record_percentage!(db, "TEST", "completed", 1, 1)
    query = "SELECT payload_json FROM artifact_rows WHERE artifact_id='TEST' ORDER BY run_id,row_index"
    register_plot_series!(db, "TEST", "fixture", query)
    db_report = validate_artifact!(db, "TEST")
    @test db_report.N == 1
    @test db_report.actual_rows == 2
    manifest_hash = finalize_artifact!(db, "TEST")
    @test length(manifest_hash) == 64
    @test artifact_status(db, "TEST") == "final"
    close_artifact_db(db)
end
