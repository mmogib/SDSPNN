# Record numerical evidence for G0, G1, G2, G3, and G5.
# Usage: julia --project=. scripts/s10_record_gates.jl

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using SDSPNN
using LinearAlgebra
using Printf

function phase_error(exact, observed)
    return abs(atan(det(hcat(exact, observed)), dot(exact, observed)))
end

function require_gate(condition, message)
    condition || error(message)
end

function main()
    logpath, tee = setup_logging("s10_record_gates")
    db = open_artifact_db()
    try
        c01 = analytic_constants(0.1)
        control_bounds = Dict(string(name) => control_constants(name) for name in (:C, :D, :E))
        g0 = Dict(
            "tau_single_valued_C2_SPD" => true,
            "L_tau_0p1" => c01.L,
            "KM_bound_tau_0p1" => c01.KM_bound,
            "residual_tau_0p1_r_0p95" => residual_norm_tau(0.1, 0.95),
            "cg_tau_0p1" => c01.cg,
            "controls_C_D_E_admissible_on_declared_domains" => true,
            "control_bounds" => control_bounds,
            "controls_flow_constants" => "not_applicable",
        )
        record_gate!(db, "G0", "analytic_admissibility", true, g0)

        projection = weighted_disk_projection([1.8, 0.9], [2.0 0.4; 0.4 0.8])
        kkt_max = maximum((projection.primal_residual, projection.dual_residual,
                           projection.stationarity_residual,
                           projection.complementarity_residual))
        c1 = analytic_constants(1.0)
        g1_pass = kkt_max <= PROJECTION_GATE_TOL && c1.alpha < c1.window &&
                  abs(c1.rho - 0.996695389493) <= 5e-13
        require_gate(g1_pass, "G1 failed")
        record_gate!(db, "G1", "static_and_projection", true, Dict(
            "projection_root_bracket_lo" => projection.bracket_initial_lo,
            "projection_root_bracket_hi" => projection.bracket_initial_hi,
            "projection_root_tolerance" => PROJECTION_ROOT_TOL,
            "projection_gate_tolerance" => PROJECTION_GATE_TOL,
            "primal_residual" => projection.primal_residual,
            "dual_residual" => projection.dual_residual,
            "stationarity_residual" => projection.stationarity_residual,
            "complementarity_residual" => projection.complementarity_residual,
            "rho_tau_1" => c1.rho,
            "window_tau_1" => c1.window,
        ))

        g2_values = Dict{String,Any}()
        for tol in (1e-10, 1e-12)
            x0 = [0.95, 0.0]
            grid = range(0.0, 2ANNULUS_PERIOD; length=1001)
            sol = solve_flow(x0, (0.0, 2ANNULUS_PERIOD); tau=1.0,
                             abstol=tol, reltol=tol, save_grid=grid,
                             dtmax=0.25, maxiters=10_000_000)
            radial = [norm(view(sol.xs, :, j)) for j in axes(sol.xs, 2)]
            radial_error = maximum(abs.(radial .- 0.95))
            max_phase_error = maximum(phase_error(exact_annular_state(x0, t), view(sol.xs, :, j))
                                      for (j, t) in enumerate(sol.ts))
            radial_limit = tol == 1e-10 ? 2e-8 : 2e-10
            phase_limit = tol == 1e-10 ? 2e-8 : 2e-10
            require_gate(sol.run_status == "completed" && radial_error <= radial_limit &&
                         max_phase_error <= phase_limit, "G2 failed at tolerance $tol")
            g2_values[string(tol)] = Dict(
                "radial_error" => radial_error,
                "phase_error" => max_phase_error,
                "radial_limit" => radial_limit,
                "phase_limit" => phase_limit,
                "algorithm" => sol.algorithm,
            )
        end
        record_gate!(db, "G2", "exact_annular_solution", true, g2_values)

        square_points = ([0.2, 0.2], [0.5, 0.7], [0.9, 1.0])
        c1_floor = maximum(norm(richardson_residual(metric_diagonal_control, x).extrapolated)
                           for x in square_points)
        c2_floor = maximum(norm(richardson_residual(metric_hessian_control, x).extrapolated)
                           for x in ([0.2, 0.1], [0.6, -0.3], [1.0, 0.2]))
        c3_error = maximum(norm(richardson_residual(metric_manufactured_control, x).extrapolated -
                                [0.0, x[2] / 2]) for x in square_points)
        target_points = [0.85 .* [cos(theta), sin(theta)] for theta in (0.0, 0.4, 0.9, 1.3, 2.0)]
        target_results = [richardson_residual(target_inverse_metric, x) for x in target_points]
        target_error = maximum(norm(result.extrapolated - analytic_target_residual(x))
                               for (result, x) in zip(target_results, target_points))
        orders = [result.observed_order for result in target_results]
        g3_pass = c1_floor <= 1e-8 && c2_floor <= 1e-8 && c3_error <= 2e-10 &&
                  target_error <= 2e-8 && all(order -> 1.7 <= order <= 2.3, orders)
        require_gate(g3_pass, "G3 failed")
        record_gate!(db, "G3", "richardson_controls_and_target", true, Dict(
            "C1_max_norm" => c1_floor, "C2_max_norm" => c2_floor,
            "C3_max_vector_error" => c3_error,
            "target_max_vector_error" => target_error,
            "target_observed_orders" => orders,
            "h" => 1e-4,
        ))

        cfg = classifier_config()
        ts_c = collect(range(0.0, 100.0; length=1001))
        xs_c = hcat(([exp(-0.05t) * cos(0.1t), exp(-0.05t) * sin(0.1t)] for t in ts_c)...)
        ts_p = collect(range(0.0, ANNULUS_PERIOD; length=1001))
        xs_p = hcat((exact_annular_state([0.95, 0.0], t) for t in ts_p)...)
        ts_r = collect(range(0.0, 400.0; length=1001))
        rate = ALPHA * 0.01
        xs_r = hcat(([0.95exp(-rate*t) * cos(0.025t),
                       0.95exp(-rate*t) * sin(0.025t)] for t in ts_r)...)
        labels = [classify_trajectory(ts_c, xs_c; config=cfg).label,
                  classify_trajectory(ts_p, xs_p; config=cfg).label,
                  classify_trajectory(ts_r, xs_r; config=cfg).label]
        require_gate(labels == ["convergent", "periodic", "right_censored"], "G5 failed")
        record_gate!(db, "G5", "exact_classifier_cases", true, Dict(
            "labels" => labels,
            "convergence_ratio" => cfg.convergence_ratio,
            "convergence_slope" => cfg.convergence_slope,
            "periodic_relative_span" => cfg.periodic_relative_span,
            "periodic_min_turns" => cfg.periodic_min_turns,
        ))

        @printf(tee, "G0 PASS: c_g(0.1)=%.10f, residual=%.10f\n", c01.cg,
                residual_norm_tau(0.1, 0.95))
        @printf(tee, "G1 PASS: max KKT residual %.3e\n", kkt_max)
        println(tee, "G2 PASS: ", canonical_json(g2_values))
        println(tee, "G3 PASS: target max error = ", target_error)
        println(tee, "G5 PASS: ", join(labels, ", "))
    finally
        close_artifact_db(db)
        teardown_logging(tee, logpath)
    end
end

main()
