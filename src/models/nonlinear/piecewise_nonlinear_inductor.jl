export PIECEWISE_NONLINEAR_INDUCTOR_TYPE,
       PSEUDO_NONLINEAR_INDUCTOR_TYPE,
       piecewise_nonlinear_inductor_network_step

const PIECEWISE_NONLINEAR_INDUCTOR_TYPE = 93
const PSEUDO_NONLINEAR_INDUCTOR_TYPE = 98

_is_pseudo_nonlinear_inductor_type(type_code::Integer) =
    Int(type_code) in (PSEUDO_NONLINEAR_INDUCTOR_TYPE, SATURATED_TRANSFORMER_NONLINEAR_TYPE)

function _piecewise_nonlinear_inductor_segment(
    current::Float64,
    table_start::Int,
    table_end::Int,
    currents_a::AbstractVector{<:Real},
)
    current <= Float64(currents_a[table_start]) && return table_start
    current >= Float64(currents_a[table_end]) && return table_end - 1
    return clamp(
        searchsortedlast(
            @view(currents_a[table_start:table_end]),
            current,
            by = Float64,
        ) + table_start - 1,
        table_start,
        table_end - 1,
    )
end

"""Solve one or more network-coupled true nonlinear inductors at one accepted step."""
function piecewise_nonlinear_inductor_network_step(
    predictor_flux_wb::AbstractVector{<:Real},
    base_branch_voltage_v::AbstractVector{<:Real},
    prior_segments::AbstractVector{Int},
    table_starts::AbstractVector{Int},
    table_ends::AbstractVector{Int},
    currents_a::AbstractVector{<:Real},
    fluxes_wb::AbstractVector{<:Real},
    branch_voltage_response_ohm::AbstractMatrix{<:Real};
    half_timestep_s::Real,
    max_iterations::Int=32,
)
    owner_count = length(predictor_flux_wb)
    length(base_branch_voltage_v) == owner_count == length(prior_segments) ==
        length(table_starts) == length(table_ends) ||
        throw(ArgumentError("nonlinear-inductor owner vectors must have equal length"))
    size(branch_voltage_response_ohm) == (owner_count, owner_count) ||
        throw(ArgumentError("nonlinear-inductor network response must be square"))
    owner_count > 0 || throw(ArgumentError("at least one nonlinear inductor is required"))
    max_iterations > 0 || throw(ArgumentError("max_iterations must be positive"))
    delta2 = Float64(half_timestep_s)
    isfinite(delta2) && delta2 > 0.0 ||
        throw(ArgumentError("half_timestep_s must be finite and positive"))

    predictor = Float64.(predictor_flux_wb)
    base_voltage = Float64.(base_branch_voltage_v)
    response = Matrix{Float64}(branch_voltage_response_ohm)
    all(isfinite, predictor) && all(isfinite, base_voltage) && all(isfinite, response) ||
        throw(ArgumentError("nonlinear-inductor network state must be finite"))
    for owner in 1:owner_count
        table_start = table_starts[owner]
        table_end = table_ends[owner]
        1 <= table_start < table_end <= length(currents_a) ||
            throw(ArgumentError("nonlinear-inductor table range is invalid"))
        table_end <= length(fluxes_wb) ||
            throw(ArgumentError("nonlinear-inductor flux table range is invalid"))
        for index in table_start:table_end
            current_value = Float64(currents_a[index])
            flux_value = Float64(fluxes_wb[index])
            isfinite(current_value) && isfinite(flux_value) ||
                throw(ArgumentError("nonlinear-inductor characteristic must be finite"))
            index == table_start && continue
            current_value > Float64(currents_a[index - 1]) &&
                flux_value > Float64(fluxes_wb[index - 1]) ||
                throw(ArgumentError("nonlinear-inductor current and flux must increase"))
        end
    end
    target_flux = predictor + delta2 .* base_voltage
    segments = copy(prior_segments)
    current = zeros(Float64, owner_count)
    extrapolated = falses(owner_count)
    converged = false
    iteration_count = 0
    for iteration in 1:max_iterations
        iteration_count = iteration
        jacobian = -delta2 .* response
        right_hand_side = copy(target_flux)
        for owner in 1:owner_count
            table_start = table_starts[owner]
            table_end = table_ends[owner]
            1 <= table_start < table_end <= length(currents_a) ||
                throw(ArgumentError("nonlinear-inductor table range is invalid"))
            table_end <= length(fluxes_wb) ||
                throw(ArgumentError("nonlinear-inductor flux table range is invalid"))
            segment = clamp(segments[owner], table_start, table_end - 1)
            current_delta = Float64(currents_a[segment + 1]) -
                            Float64(currents_a[segment])
            flux_delta = Float64(fluxes_wb[segment + 1]) -
                         Float64(fluxes_wb[segment])
            current_delta > 0.0 && flux_delta > 0.0 ||
                throw(ArgumentError("nonlinear-inductor current and flux must increase"))
            slope_h = flux_delta / current_delta
            intercept_wb = Float64(fluxes_wb[segment]) -
                           slope_h * Float64(currents_a[segment])
            jacobian[owner, owner] += slope_h
            right_hand_side[owner] -= intercept_wb
        end
        current = jacobian \ right_hand_side
        all(isfinite, current) ||
            throw(ArgumentError("nonlinear-inductor coupled solve produced non-finite current"))
        new_segments = similar(segments)
        for owner in 1:owner_count
            new_segments[owner] = _piecewise_nonlinear_inductor_segment(
                current[owner],
                table_starts[owner],
                table_ends[owner],
                currents_a,
            )
            extrapolated[owner] =
                current[owner] < Float64(currents_a[table_starts[owner]]) ||
                current[owner] > Float64(currents_a[table_ends[owner]])
        end
        if new_segments == segments
            converged = true
            break
        end
        segments = new_segments
    end
    converged || throw(ArgumentError("nonlinear-inductor active-set solve did not converge"))
    any(extrapolated) &&
        throw(ArgumentError("nonlinear-inductor current lies outside its characteristic"))

    accepted_voltage = base_voltage + response * current
    accepted_flux = predictor + delta2 .* accepted_voltage
    next_predictor = accepted_flux + delta2 .* accepted_voltage
    return (
        accepted_current_a = current,
        accepted_flux_wb = accepted_flux,
        predictor_flux_wb = next_predictor,
        accepted_branch_voltage_v = accepted_voltage,
        active_segments = segments,
        segment_change_count = count(identity, segments .!= prior_segments),
        extrapolated = extrapolated,
        iteration_count = iteration_count,
        converged = converged,
        base_branch_voltage_v = base_voltage,
    )
end
