const PROJECTION_ROOT_TOL = 1e-13
const PROJECTION_GATE_TOL = 5e-11

function _projection_diagnostics(z, y, Q, radius, nu)
    primal = max(norm(z) - radius, 0.0)
    dual = max(-nu, 0.0)
    stationarity = norm(Q * (z - y) + nu * z)
    complementarity = abs(nu * (dot(z, z) - radius^2))
    return primal, dual, stationarity, complementarity
end

"""
    weighted_disk_projection(y, Q, radius; root_tol=1e-13, gate_tol=5e-11)

Project `y` onto the Euclidean disk in the SPD metric `Q`. Outside the disk,
solve `norm((Q+nu*I)\\(Q*y))=radius` by bisection. The certified initial
bracket is `[0, norm(Q*y)/radius]`.
"""
function weighted_disk_projection(y::AbstractVector, Q::AbstractMatrix,
                                  radius::Real=DISK_RADIUS;
                                  root_tol::Real=PROJECTION_ROOT_TOL,
                                  gate_tol::Real=PROJECTION_GATE_TOL)
    length(y) == 2 || throw(DimensionMismatch("weighted disk projection requires a 2-vector"))
    size(Q) == (2, 2) || throw(DimensionMismatch("Q must be 2 by 2"))
    radius > 0 || throw(ArgumentError("radius must be positive"))
    Qs = Matrix(Symmetric(Matrix{Float64}(Q)))
    minimum(eigvals(Symmetric(Qs))) > 0 || throw(ArgumentError("Q must be SPD"))
    yy = Float64.(y)

    if norm(yy) <= radius
        primal, dual, stationarity, comp = _projection_diagnostics(yy, yy, Qs, radius, 0.0)
        return ProjectionResult(yy, 0.0, false, 0, 0.0, 0.0, 0.0, 0.0,
                                primal, dual, stationarity, comp)
    end

    qy = Qs * yy
    lo = 0.0
    hi = norm(qy) / radius
    initial_hi = hi
    z_at(nu) = (Qs + nu * I) \ qy
    norm(z_at(hi)) <= radius + 10eps(Float64) ||
        throw(ErrorException("certified multiplier bracket failed"))

    iterations = 0
    while hi - lo > root_tol * max(1.0, initial_hi)
        iterations += 1
        iterations <= 256 || throw(ErrorException("multiplier bisection exceeded 256 iterations"))
        mid = (lo + hi) / 2
        if norm(z_at(mid)) > radius
            lo = mid
        else
            hi = mid
        end
    end
    nu = (lo + hi) / 2
    z = Vector{Float64}(z_at(nu))
    primal, dual, stationarity, comp = _projection_diagnostics(z, yy, Qs, radius, nu)
    maximum((primal, dual, stationarity, comp)) <= gate_tol ||
        throw(ErrorException("weighted-disk KKT gate failed"))
    return ProjectionResult(z, nu, true, iterations, 0.0, initial_hi, lo, hi,
                            primal, dual, stationarity, comp)
end
