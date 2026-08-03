const DISK_RADIUS = 1.1
const ALPHA = 0.05
const LAMBDA = 1.0
const MU = 1.0
const K_OPERATOR = sqrt(5.0)
const F_MAX = 1.1sqrt(5.0)
const KM1_BOUND = 3.571
const ANNULUS_PERIOD = 4pi / ALPHA

const _B = [1.0 -0.5; -0.5 0.5]
const _A_HESSIAN = [1.0 0.0; 0.0 0.5]

"""C2 smootherstep cutoff used by the rotating metric."""
function chi_radius(r::Real)
    r <= 0 && return 0.0
    r >= 0.8 && return 1.0
    u = Float64(r) / 0.8
    return 10u^3 - 15u^4 + 6u^5
end

"""The revision-2 interpolation M_tau = I + tau(M_1-I)."""
function metric_tau(x::AbstractVector, tau::Real)
    length(x) == 2 || throw(DimensionMismatch("metric_tau requires a 2-vector"))
    0 <= tau <= 1 || throw(ArgumentError("tau must lie in [0,1]"))
    r = norm(x)
    c = chi_radius(r)
    if r == 0 || tau == 0 || c == 0
        return Matrix{Float64}(I, 2, 2)
    end
    er = Float64.(x) / r
    et = [-er[2], er[1]]
    R = hcat(er, et)
    return (1 - tau*c) .* Matrix{Float64}(I, 2, 2) .+ tau*c .* (R * _B * R')
end

function metric_bounds_tau(tau::Real)
    0 <= tau <= 1 || throw(ArgumentError("tau must lie in [0,1]"))
    lower = 1 - ((1 + sqrt(5.0)) / 4) * tau
    upper = 1 + ((sqrt(5.0) - 1) / 4) * tau
    return lower, upper, inv(lower)
end

"""Globally bounded Hessian control from revision 2."""
function metric_hessian_control(x::AbstractVector)
    length(x) == 2 || throw(DimensionMismatch("control metric requires a 2-vector"))
    xx = Float64.(x)
    s2 = 1 + dot(xx, xx)
    return _A_HESSIAN .+ 0.35 .* (inv(sqrt(s2)) .* Matrix{Float64}(I, 2, 2) .-
                                  inv(s2^(3/2)) .* (xx * xx'))
end

"""Separable diagonal control on [0,1]^2."""
function metric_diagonal_control(x::AbstractVector)
    length(x) == 2 || throw(DimensionMismatch("control metric requires a 2-vector"))
    m(t) = 1 + 16t * (1 - t)
    return [m(Float64(x[1])) 0.0; 0.0 m(Float64(x[2]))]
end

"""Manufactured non-integrable control with residual (0,x2/2)."""
function metric_manufactured_control(x::AbstractVector)
    length(x) == 2 || throw(DimensionMismatch("control metric requires a 2-vector"))
    return [1 + Float64(x[2])^2 / 4 0.0; 0.0 1.0]
end

function control_constants(name::Symbol)
    if name === :C
        return (spectral_lower=0.5, spectral_upper=1.35, L=2.0,
                KM_bound=0.600, integrability="integrable")
    elseif name === :D
        return (spectral_lower=1.0, spectral_upper=5.0, L=5.0,
                KM_bound=16.0, integrability="integrable_on_[0,1]^2")
    elseif name === :E
        return (spectral_lower=1.0, spectral_upper=1.25, L=1.25,
                KM_bound=0.5, integrability="nonintegrable_except_x2=0")
    end
    throw(ArgumentError("unknown control metric $name"))
end

function residual_norm_tau(tau::Real, r::Real)
    r > 0 || throw(ArgumentError("r must be positive"))
    detD = 1 - tau / 2 - tau^2 / 4
    return tau * sqrt(5.0) / (2detD * r)
end
