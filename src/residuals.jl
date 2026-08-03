function central_residual(Hfun::Function, x::AbstractVector, h::Real)
    h > 0 || throw(ArgumentError("finite-difference step must be positive"))
    xx = Float64.(x)
    e1 = [1.0, 0.0]
    e2 = [0.0, 1.0]
    derivative(i, j, e) = (Hfun(xx .+ h .* e)[i, j] - Hfun(xx .- h .* e)[i, j]) / (2h)
    return [derivative(2, 2, e1) - derivative(1, 2, e2),
            derivative(1, 1, e2) - derivative(1, 2, e1)]
end

function richardson_residual(Hfun::Function, x::AbstractVector; h::Real=1e-4)
    raw_h = central_residual(Hfun, x, h)
    raw_h2 = central_residual(Hfun, x, h / 2)
    raw_h4 = central_residual(Hfun, x, h / 4)
    rich_coarse = (4raw_h2 - raw_h) / 3
    rich_fine = (4raw_h4 - raw_h2) / 3
    d1 = norm(raw_h - raw_h2)
    d2 = norm(raw_h2 - raw_h4)
    order = d1 > 0 && d2 > 0 ? log2(d1 / d2) : NaN
    truncation = norm(rich_fine - rich_coarse) / 15
    return (
        raw_h=raw_h, raw_h2=raw_h2, raw_h4=raw_h4,
        richardson_coarse=rich_coarse, extrapolated=rich_fine,
        observed_order=order, truncation_estimate=truncation,
    )
end

target_inverse_metric(x::AbstractVector) = inv(Symmetric(metric_tau(x, 1.0)))

function analytic_target_residual(x::AbstractVector)
    r = norm(x)
    r > 0 || throw(ArgumentError("target residual is singular at the origin"))
    c, s = x[1] / r, x[2] / r
    return [2(c + 2s) / r, 2(s - 2c) / r]
end
