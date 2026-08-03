function flow_map(x::AbstractVector; tau::Real=1.0)
    M = metric_tau(x, tau)
    forward = Float64.(x) .- ALPHA .* (M * operator_A(x))
    Q = inv(Symmetric(M))
    return weighted_disk_projection(forward, Q, DISK_RADIUS)
end

function flow_rhs!(du, u, tau, t)
    projection = flow_map(u; tau=tau)
    @. du = LAMBDA * (projection.z - u)
    return nothing
end

"""Integrate the SD-SPNN flow with the fully declared revision-2 protocol."""
function solve_flow(x0::AbstractVector, tspan::Tuple{<:Real,<:Real};
                    tau::Real=1.0, abstol::Real=1e-10, reltol::Real=1e-10,
                    save_grid=range(tspan[1], tspan[2]; length=2001),
                    dtmax::Real=0.25, maxiters::Int=10_000_000)
    prob = ODEProblem(flow_rhs!, Float64.(x0), (Float64(tspan[1]), Float64(tspan[2])), Float64(tau))
    algorithm = Rodas5P(autodiff=AutoFiniteDiff())
    sol = solve(prob, algorithm; abstol=Float64(abstol), reltol=Float64(reltol),
                saveat=collect(Float64, save_grid), dtmax=Float64(dtmax),
                maxiters=maxiters, dense=false)
    xs = Matrix{Float64}(Array(sol))
    ts = Float64.(sol.t)
    active = [flow_map(view(xs, :, j); tau=tau).active for j in axes(xs, 2)]
    drift = maximum(max(norm(view(xs, :, j)) - DISK_RADIUS, 0.0) for j in axes(xs, 2))
    status = successful_retcode(sol) ? "completed" : "solver_failure"
    return FlowResult(
        ts, xs, active, status, string(sol.retcode), drift,
        "Rodas5P(autodiff=AutoFiniteDiff())", Float64(abstol), Float64(reltol), Float64(dtmax), maxiters,
        "No state reset; the continuous weighted-projection active set is evaluated at every RHS call.",
        "Any non-success retcode, nonfinite state, projection KKT failure, or feasibility drift above tolerance fails the run.",
    )
end

function exact_annular_state(x0::AbstractVector, t::Real)
    angle = -ALPHA * Float64(t) / 2
    R = [cos(angle) -sin(angle); sin(angle) cos(angle)]
    return R * Float64.(x0)
end

function weighted_radius(x::AbstractVector; tau::Real=1.0)
    M = metric_tau(x, tau)
    return sqrt(dot(x, M \ x))
end

function euclidean_radius_for_weighted(target::Real, angle::Real=0.0; tau::Real=1.0)
    target >= 0 || throw(ArgumentError("target weighted radius must be nonnegative"))
    target == 0 && return 0.0
    direction = [cos(Float64(angle)), sin(Float64(angle))]
    lo, hi = 0.0, DISK_RADIUS
    weighted_radius(hi .* direction; tau=tau) >= target ||
        throw(ArgumentError("target weighted radius lies outside the disk"))
    for _ in 1:100
        mid = (lo + hi) / 2
        if weighted_radius(mid .* direction; tau=tau) < target
            lo = mid
        else
            hi = mid
        end
    end
    return (lo + hi) / 2
end
