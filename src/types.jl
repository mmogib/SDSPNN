struct ProjectionResult
    z::Vector{Float64}
    nu::Float64
    active::Bool
    iterations::Int
    bracket_initial_lo::Float64
    bracket_initial_hi::Float64
    bracket_final_lo::Float64
    bracket_final_hi::Float64
    primal_residual::Float64
    dual_residual::Float64
    stationarity_residual::Float64
    complementarity_residual::Float64
end

struct FlowResult
    ts::Vector{Float64}
    xs::Matrix{Float64}
    projection_active::Vector{Bool}
    run_status::String
    retcode::String
    feasibility_drift::Float64
    algorithm::String
    abstol::Float64
    reltol::Float64
    dtmax::Float64
    maxiters::Int
    event_handling::String
    failure_policy::String
end

Base.@kwdef struct ClassifierConfig
    convergence_ratio::Float64 = 0.1
    convergence_slope::Float64 = -1e-4
    periodic_relative_span::Float64 = 1e-6
    periodic_min_turns::Float64 = 0.9
end
