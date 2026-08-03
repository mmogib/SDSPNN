const ACTIVE_CENTER = [2.0, 0.0]
const ACTIVE_MU = 1.0
const ACTIVE_K = 1.0
const ACTIVE_FMAX = 3.1

_active_operator(x::AbstractVector) = Float64.(x) - ACTIVE_CENTER

function _active_forward_point(x::AbstractVector, tau::Real)
    xx = Float64.(x)
    return xx - ALPHA * (metric_tau(xx, tau) * _active_operator(xx))
end

function _sign_bisection(f::Function, lo::Real, hi::Real; tol::Real=1e-14)
    a, b = Float64(lo), Float64(hi)
    fa, fb = f(a), f(b)
    fa == 0 && return a
    fb == 0 && return b
    fa * fb < 0 || throw(ArgumentError("bisection bracket must change sign"))
    while b - a > tol
        mid = (a + b) / 2
        fm = f(mid)
        if fa * fm <= 0
            b = mid
        else
            a = mid
            fa = fm
        end
    end
    return (a + b) / 2
end

function weighted_disk_projection_spectral(y::AbstractVector, Q::AbstractMatrix,
                                           radius::Real=DISK_RADIUS;
                                           root_tol::Real=1e-14)
    length(y) == 2 || throw(DimensionMismatch("spectral projection requires a 2-vector"))
    size(Q) == (2, 2) || throw(DimensionMismatch("Q must be 2 by 2"))
    radius > 0 || throw(ArgumentError("radius must be positive"))
    root_tol > 0 || throw(ArgumentError("root tolerance must be positive"))

    yy = Float64.(y)
    Qs = Matrix(Symmetric(Matrix{Float64}(Q)))
    decomposition = eigen(Symmetric(Qs))
    minimum(decomposition.values) > 0 || throw(ArgumentError("Q must be SPD"))
    if norm(yy) <= radius
        return (z=yy, nu=0.0, active=false, iterations=0)
    end

    coordinates = decomposition.vectors' * yy
    z_at(nu) = decomposition.vectors *
        ((decomposition.values ./ (decomposition.values .+ nu)) .* coordinates)
    lo = 0.0
    hi = norm(Qs * yy) / radius
    initial_hi = hi
    norm(z_at(hi)) <= radius + 10eps(Float64) ||
        throw(ErrorException("spectral multiplier bracket failed"))
    iterations = 0
    while hi - lo > root_tol * max(1.0, initial_hi)
        iterations += 1
        iterations <= 256 ||
            throw(ErrorException("spectral multiplier bisection exceeded 256 iterations"))
        mid = (lo + hi) / 2
        if norm(z_at(mid)) > radius
            lo = mid
        else
            hi = mid
        end
    end
    nu = (lo + hi) / 2
    return (z=Vector{Float64}(z_at(nu)), nu=nu, active=true, iterations=iterations)
end

function active_projection_constants(tau::Real)
    0 <= tau <= 1 || throw(ArgumentError("tau must lie in [0,1]"))
    L = inv(1 - ((1 + sqrt(5.0)) / 4) * tau)
    KM = KM1_BOUND * tau
    rho = sqrt(1 - 2ALPHA * ACTIVE_MU / L + ALPHA^2 * ACTIVE_K^2)
    cg = 1 - rho - (ALPHA / 2) * L^2 * KM * ACTIVE_FMAX
    return (
        tau=Float64(tau), mu=ACTIVE_MU, K=ACTIVE_K, Fmax=ACTIVE_FMAX,
        L=L, KM_bound=KM, rho=rho, cg=cg,
        window=2ACTIVE_MU / (L * ACTIVE_K^2),
        displacement_bound=ALPHA * L * ACTIVE_FMAX,
    )
end

function active_projection_point(x::AbstractVector, tau::Real)
    length(x) == 2 || throw(DimensionMismatch("active projection point requires a 2-vector"))
    xx = Float64.(x)
    F = _active_operator(xx)
    M = metric_tau(xx, tau)
    Q = Matrix(inv(Symmetric(M)))
    y = _active_forward_point(xx, tau)
    projection = weighted_disk_projection(y, Q, DISK_RADIUS)
    spectral = weighted_disk_projection_spectral(y, Q, DISK_RADIUS)
    clamp = norm(y) <= DISK_RADIUS ? copy(y) : DISK_RADIUS .* y ./ norm(y)
    wrong_metric = weighted_disk_projection(y, M, DISK_RADIUS)
    return (
        x=xx, F=F, M=M, Q=Q, y=y,
        projection=projection, spectral=spectral, clamp=clamp,
        active_margin=norm(y) - DISK_RADIUS,
        projection_displacement=norm(projection.z - y),
        map_displacement=norm(projection.z - xx),
        clamp_gap=norm(projection.z - clamp),
        wrong_metric_map_gap=norm(wrong_metric.z - xx),
    )
end

function active_boundary_geometry(tau::Real)
    0 <= tau <= 1 || throw(ArgumentError("tau must lie in [0,1]"))
    margin(phi) = norm(_active_forward_point(
        DISK_RADIUS .* [cos(phi), sin(phi)], tau)) - DISK_RADIUS
    grid = collect(range(-pi, pi; length=4097))
    roots = Float64[]
    previous_phi = grid[1]
    previous_value = margin(previous_phi)
    for phi in grid[2:end]
        value = margin(phi)
        if previous_value * value < 0
            push!(roots, _sign_bisection(margin, previous_phi, phi))
        end
        previous_phi = phi
        previous_value = value
    end
    length(roots) == 2 ||
        throw(ErrorException("expected exactly two active-arc endpoints"))
    roots[1] < 0 < roots[2] ||
        throw(ErrorException("active boundary arc must contain the equilibrium angle"))
    margin(0.0) > 0 || throw(ErrorException("equilibrium must lie in the active arc"))

    radial_margin(r) = norm(_active_forward_point([r, 0.0], tau)) - DISK_RADIUS
    radial = _sign_bisection(radial_margin, 0.8, DISK_RADIUS)
    return (arc_lo=roots[1], arc_hi=roots[2], radial_threshold=radial)
end

function active_clamp_signature(tau::Real)
    0 <= tau <= 1 || throw(ArgumentError("tau must lie in [0,1]"))
    amplitude = hypot(tau, 2 - tau)
    phase = atan(tau, 2 - tau)
    rhs = tau * DISK_RADIUS / (2amplitude)
    abs(rhs) <= 1 || throw(ErrorException("radial-clamp signature has no real root"))
    phi = asin(rhs) - phase
    distance = 2DISK_RADIUS * abs(sin(phi / 2))
    return (phi=phi, distance_from_solution=distance,
            amplitude=amplitude, phase=phase, normalized_rhs=rhs)
end

function active_sample_angles(tau::Real, count::Integer=8)
    count >= 1 || throw(ArgumentError("sample count must be positive"))
    geometry = active_boundary_geometry(tau)
    width = geometry.arc_hi - geometry.arc_lo
    return [geometry.arc_lo + j * width / (count + 1) for j in 1:count]
end

function active_projection_joint_bounds(point, anchor, L::Real)
    L >= 1 || throw(ArgumentError("common spectral bound must be at least one"))
    Q1 = point.Q
    Q2 = anchor.Q
    spectrum = eigvals(Symmetric(Q1))
    lambda_min = minimum(spectrum)
    lambda_max = maximum(spectrum)
    argument_gap = norm(point.y - anchor.y)
    metric_gap = opnorm(Q1 - Q2)
    anchor_displacement = anchor.projection_displacement
    lhs = norm(point.projection.z - anchor.projection.z)
    rhs_exact = sqrt(lambda_max / lambda_min) * argument_gap +
                metric_gap * anchor_displacement / lambda_min
    rhs_L = L * argument_gap + L * metric_gap * anchor_displacement
    rhs_L2 = L^2 * argument_gap + L * metric_gap * anchor_displacement
    rhs_exact > 0 || throw(ErrorException("joint-bound denominator must be positive"))
    rhs_L > 0 || throw(ErrorException("L-bound denominator must be positive"))
    rhs_L2 > 0 || throw(ErrorException("L-squared-bound denominator must be positive"))
    return (
        lhs=lhs, rhs_exact=rhs_exact, ratio_exact=lhs / rhs_exact,
        rhs_L=rhs_L, ratio_L=lhs / rhs_L,
        rhs_L2=rhs_L2, ratio_L2=lhs / rhs_L2,
        argument_gap=argument_gap, metric_gap=metric_gap,
        anchor=anchor_displacement,
        lambda_min_Q1=lambda_min, lambda_max_Q1=lambda_max,
    )
end
