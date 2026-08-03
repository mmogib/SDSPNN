classifier_config() = ClassifierConfig()

function _angular_travel(xs::AbstractMatrix)
    angles = atan.(view(xs, 2, :), view(xs, 1, :))
    travel = 0.0
    for j in 2:length(angles)
        delta = angles[j] - angles[j-1]
        travel += atan(sin(delta), cos(delta))
    end
    return travel
end

function classify_trajectory(ts::AbstractVector, xs::AbstractMatrix;
                             config::ClassifierConfig=classifier_config())
    length(ts) == size(xs, 2) || throw(DimensionMismatch("one state is required per time"))
    length(ts) >= 2 || throw(ArgumentError("at least two samples are required"))
    radii = vec(sqrt.(sum(abs2, xs; dims=1)))
    ratio = radii[end] / radii[1]
    slope = (log(max(radii[end], floatmin(Float64))) -
             log(max(radii[1], floatmin(Float64)))) / (ts[end] - ts[1])
    relative_span = (maximum(radii) - minimum(radii)) / max(sum(radii) / length(radii), eps())
    turns = abs(_angular_travel(xs)) / (2pi)

    label = if ratio <= config.convergence_ratio && slope <= config.convergence_slope
        "convergent"
    elseif relative_span <= config.periodic_relative_span && turns >= config.periodic_min_turns
        "periodic"
    else
        "right_censored"
    end
    return (label=label, radius_ratio=ratio, log_slope=slope,
            periodic_relative_span=relative_span, turns=turns)
end
