operator_matrix() = [1.0 -2.0; 2.0 1.0]

function operator_A(x::AbstractVector)
    length(x) == 2 || throw(DimensionMismatch("operator_A requires a 2-vector"))
    return operator_matrix() * x
end

function rho_tau(tau::Real)
    _, _, L = metric_bounds_tau(tau)
    return sqrt(1 - 2ALPHA * MU / L + ALPHA^2 * K_OPERATOR^2)
end

function cg_tau(tau::Real)
    _, _, L = metric_bounds_tau(tau)
    rho = rho_tau(tau)
    return 1 - rho - (ALPHA / 2) * L^2 * (KM1_BOUND * tau) * F_MAX
end

function analytic_constants(tau::Real)
    lower, upper, L = metric_bounds_tau(tau)
    KM = KM1_BOUND * tau
    rho = rho_tau(tau)
    window = 2MU / (L * K_OPERATOR^2)
    CL = 0.5 * L^(3/2) * KM * (1 + rho)
    R_clean = KM == 0 ? Inf : (1 - rho) / (L^(3/2) * KM * (1 + rho))
    return (
        tau=Float64(tau), mu=MU, K=K_OPERATOR, Fmax=F_MAX,
        alpha=ALPHA, lambda=LAMBDA, spectral_lower=lower,
        spectral_upper=upper, L=L, KM_bound=KM, window=window,
        inside_window=ALPHA < window, rho=rho, CL=CL,
        R_clean=R_clean, eta_R=(1 - rho) / 2, cg=cg_tau(tau),
    )
end

function cg_root(lo::Real=0.139, hi::Real=0.140; tol::Real=1e-14)
    flo = cg_tau(lo)
    fhi = cg_tau(hi)
    flo > 0 && fhi < 0 || throw(ArgumentError("root bracket must have c_g(lo)>0>c_g(hi)"))
    a, b = Float64(lo), Float64(hi)
    while b - a > tol
        mid = (a + b) / 2
        if cg_tau(mid) > 0
            a = mid
        else
            b = mid
        end
    end
    return (a + b) / 2
end
