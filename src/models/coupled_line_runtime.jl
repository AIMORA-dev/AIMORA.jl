module CoupledLineRuntime

using ..CoupledLineFitting
using LinearAlgebra
using SHA
using TOML

export CoupledLineRuntimeSettings,
       CoupledLineRuntimePreparation,
       CoupledLineRuntimeState,
       CoupledLineRuntimeSnapshot,
       prepare_coupled_line_runtime,
       coupled_line_runtime_state,
       initialize_coupled_line_runtime_deenergized!,
       initialize_coupled_line_runtime_sinusoidal!,
       accept_coupled_line_runtime_step!,
       coupled_line_runtime_discrete_response,
       coupled_line_runtime_terminal_admittance,
       coupled_line_runtime_snapshot,
       restore_coupled_line_runtime_snapshot!,
       write_coupled_line_runtime_snapshot,
       read_coupled_line_runtime_snapshot,
       coupled_line_runtime_diagnostics,
       coupled_line_runtime_report_text

const COUPLED_LINE_RUNTIME_SCHEMA_VERSION = 1
const COUPLED_LINE_RUNTIME_MAXIMUM_PHASE_COUNT = 12

function _runtime_sha256(write_payload)
    io = IOBuffer()
    write_payload(io)
    return bytes2hex(sha256(take!(io)))
end

function _write_runtime_value(io, value::Float64)
    println(io, bitstring(value))
end

function _write_runtime_value(io, value::Integer)
    println(io, value)
end

function _write_runtime_value(io, value::Bool)
    println(io, value)
end

function _write_runtime_value(io, value::Symbol)
    println(io, String(value))
end

function _write_runtime_value(io, value::AbstractString)
    println(io, value)
end

function _write_runtime_values(io, values)
    for value in values
        _write_runtime_value(io, value)
    end
end

function _checked_runtime_signature(value::AbstractString, label::AbstractString)
    signature = lowercase(String(value))
    occursin(r"^[0-9a-f]{64}$", signature) ||
        throw(ArgumentError("$label must be lowercase SHA-256"))
    return signature
end

function _checked_positive_runtime_value(value::Real, label::AbstractString)
    checked = Float64(value)
    isfinite(checked) && checked > 0.0 ||
        throw(ArgumentError("$label must be finite and positive"))
    return checked
end

struct CoupledLineRuntimeSettings
    timestep_s::Float64
    prewarp_frequency_hz::Union{Nothing,Float64}
    maximum_condition_number::Float64
    maximum_relative_frequency_warp::Float64
    passivity_tolerance::Float64
    reciprocity_relative_tolerance::Float64
    kcl_absolute_tolerance_a::Float64
    energy_absolute_tolerance_j::Float64
    deterministic_signature_sha256::String
end

function CoupledLineRuntimeSettings(;
    timestep_s::Real,
    prewarp_frequency_hz::Union{Nothing,Real}=nothing,
    maximum_condition_number::Real=1.0e12,
    maximum_relative_frequency_warp::Real=0.10,
    passivity_tolerance::Real=1.0e-8,
    reciprocity_relative_tolerance::Real=1.0e-8,
    kcl_absolute_tolerance_a::Real=1.0e-9,
    energy_absolute_tolerance_j::Real=1.0e-9,
)
    timestep = _checked_positive_runtime_value(timestep_s, "coupled line runtime timestep_s")
    prewarp = prewarp_frequency_hz === nothing ? nothing :
        _checked_positive_runtime_value(
            prewarp_frequency_hz,
            "coupled line runtime prewarp_frequency_hz",
        )
    prewarp === nothing ||
        prewarp * timestep < 0.5 - 64.0 * eps(Float64) ||
        throw(ArgumentError(
            "coupled line runtime prewarp frequency must be below Nyquist",
        ))
    condition_limit = _checked_positive_runtime_value(
        maximum_condition_number,
        "coupled line runtime maximum_condition_number",
    )
    condition_limit > 1.0 ||
        throw(ArgumentError("coupled line runtime maximum_condition_number must exceed one"))
    warp_limit = _checked_positive_runtime_value(
        maximum_relative_frequency_warp,
        "coupled line runtime maximum_relative_frequency_warp",
    )
    passivity_limit = _checked_positive_runtime_value(
        passivity_tolerance,
        "coupled line runtime passivity_tolerance",
    )
    reciprocity_limit = _checked_positive_runtime_value(
        reciprocity_relative_tolerance,
        "coupled line runtime reciprocity_relative_tolerance",
    )
    kcl_limit = _checked_positive_runtime_value(
        kcl_absolute_tolerance_a,
        "coupled line runtime kcl_absolute_tolerance_a",
    )
    energy_limit = _checked_positive_runtime_value(
        energy_absolute_tolerance_j,
        "coupled line runtime energy_absolute_tolerance_j",
    )
    signature = _runtime_sha256() do io
        println(io, COUPLED_LINE_RUNTIME_SCHEMA_VERSION)
        _write_runtime_value(io, timestep)
        println(io, prewarp === nothing ? "no_prewarp" : bitstring(prewarp))
        _write_runtime_values(
            io,
            (
                condition_limit,
                warp_limit,
                passivity_limit,
                reciprocity_limit,
                kcl_limit,
                energy_limit,
            ),
        )
    end
    return CoupledLineRuntimeSettings(
        timestep,
        prewarp,
        condition_limit,
        warp_limit,
        passivity_limit,
        reciprocity_limit,
        kcl_limit,
        energy_limit,
        signature,
    )
end

struct CoupledLineRuntimePreparation
    schema_version::Int
    source_signature_sha256::String
    response_signature_sha256::String
    fit_signature_sha256::String
    uncertainty_alternative_fit_signatures_sha256::Vector{String}
    uncertainty_set_complete::Bool
    model_signature_sha256::String
    phase_order::Vector{Symbol}
    port_order::Vector{Symbol}
    reference_impedance_ohm::Vector{Float64}
    source_frequency_band_hz::NTuple{2,Float64}
    settings::CoupledLineRuntimeSettings
    bilinear_alpha_per_s::Float64
    state_transition::Matrix{Float64}
    endpoint_input::Matrix{Float64}
    output_matrix::Matrix{Float64}
    continuous_direct_term::Matrix{Float64}
    discrete_direct_term::Matrix{Float64}
    history_state_output::Matrix{Float64}
    history_incident_output::Matrix{Float64}
    incident_from_voltage::Matrix{Float64}
    incident_from_history::Matrix{Float64}
    history_to_terminal_current::Matrix{Float64}
    companion_admittance_s::Matrix{Float64}
    reference_impedance_sqrt_ohm_sqrt::Vector{Float64}
    reference_impedance_inverse_sqrt_per_ohm_sqrt::Vector{Float64}
    bilinear_matrix_condition_number::Float64
    wave_matrix_condition_number::Float64
    maximum_discrete_pole_magnitude::Float64
    maximum_relative_frequency_warp::Float64
    sampled_maximum_scattering_singular_value::Float64
    companion_reciprocity_relative_error::Float64
    minimum_companion_conductance_eigenvalue_s::Float64
    deterministic_signature_sha256::String
end

function _runtime_bilinear_alpha(settings::CoupledLineRuntimeSettings)
    settings.prewarp_frequency_hz === nothing &&
        return 2.0 / settings.timestep_s
    angular_frequency = 2.0 * pi * something(settings.prewarp_frequency_hz)
    half_angle = 0.5 * angular_frequency * settings.timestep_s
    half_angle < 0.5 * pi ||
        throw(ArgumentError("coupled line runtime prewarp frequency must be below Nyquist"))
    tangent = tan(half_angle)
    abs(tangent) > eps(Float64) ||
        throw(ArgumentError("coupled line runtime prewarp mapping is numerically singular"))
    return angular_frequency / tangent
end

function _runtime_frequency_warp(
    frequency_hz::Float64,
    timestep_s::Float64,
    alpha_per_s::Float64,
)
    frequency_hz == 0.0 && return 0.0
    angle = pi * frequency_hz * timestep_s
    angle < 0.5 * pi ||
        throw(ArgumentError("coupled line runtime frequency band reaches or exceeds Nyquist"))
    mapped_angular_frequency = alpha_per_s * tan(angle)
    physical_angular_frequency = 2.0 * pi * frequency_hz
    return abs(mapped_angular_frequency / physical_angular_frequency - 1.0)
end

function _runtime_discrete_response(
    state_transition::AbstractMatrix{Float64},
    endpoint_input::AbstractMatrix{Float64},
    output_matrix::AbstractMatrix{Float64},
    continuous_direct_term::AbstractMatrix{Float64},
    frequency_hz::Float64,
    timestep_s::Float64,
)
    angle = 2.0 * pi * frequency_hz * timestep_s
    angle < pi ||
        throw(ArgumentError("coupled line runtime response frequency must be below Nyquist"))
    z = cis(angle)
    state_count = size(state_transition, 1)
    identity_state = Matrix{ComplexF64}(I, state_count, state_count)
    state_response = (z .* identity_state .- state_transition) \
        (ComplexF64.(endpoint_input) .* (z + 1.0))
    return ComplexF64.(continuous_direct_term) .+
        ComplexF64.(output_matrix) * state_response
end

function coupled_line_runtime_discrete_response(
    preparation::CoupledLineRuntimePreparation,
    frequency_hz::Real,
)
    frequency = Float64(frequency_hz)
    isfinite(frequency) && frequency >= 0.0 ||
        throw(ArgumentError("coupled line runtime frequency must be finite and nonnegative"))
    frequency <= last(preparation.source_frequency_band_hz) ||
        throw(ArgumentError("coupled line runtime frequency exceeds its accepted source band"))
    return _runtime_discrete_response(
        preparation.state_transition,
        preparation.endpoint_input,
        preparation.output_matrix,
        preparation.continuous_direct_term,
        frequency,
        preparation.settings.timestep_s,
    )
end

function coupled_line_runtime_terminal_admittance(
    preparation::CoupledLineRuntimePreparation,
    frequency_hz::Real,
)
    scattering = coupled_line_runtime_discrete_response(preparation, frequency_hz)
    return coupled_line_scattering_to_admittance(
        scattering,
        preparation.reference_impedance_ohm,
    )
end

function _coupled_line_runtime_preparation(
    model::CoupledLineRationalModel,
    settings::CoupledLineRuntimeSettings;
    source_signature_sha256::AbstractString,
    response_signature_sha256::AbstractString,
    fit_signature_sha256::AbstractString,
    phase_order,
    frequencies_hz,
    continuous_passivity_passed::Bool,
    uncertainty_alternative_fit_signatures_sha256=String[],
    uncertainty_set_complete::Bool=false,
)
    continuous_passivity_passed ||
        throw(ArgumentError("coupled line runtime requires an accepted passive L205 fit"))
    source_signature = _checked_runtime_signature(
        source_signature_sha256,
        "coupled line runtime source signature",
    )
    response_signature = _checked_runtime_signature(
        response_signature_sha256,
        "coupled line runtime response signature",
    )
    fit_signature = _checked_runtime_signature(
        fit_signature_sha256,
        "coupled line runtime fit signature",
    )
    uncertainty_signatures = _checked_runtime_signature.(
        String.(uncertainty_alternative_fit_signatures_sha256),
        Ref("coupled line runtime uncertainty alternative fit signature"),
    )
    length(unique(uncertainty_signatures)) == length(uncertainty_signatures) ||
        throw(ArgumentError(
            "coupled line runtime uncertainty alternative signatures must be unique",
        ))
    uncertainty_set_complete && isempty(uncertainty_signatures) && throw(ArgumentError(
        "a complete coupled line runtime uncertainty set must contain an alternative",
    ))
    fit_signature in uncertainty_signatures && throw(ArgumentError(
        "coupled line runtime nominal fit cannot also be an uncertainty alternative",
    ))
    model_signature = _checked_runtime_signature(
        model.deterministic_signature_sha256,
        "coupled line runtime model signature",
    )
    phases = Symbol.(phase_order)
    phase_count = length(phases)
    1 <= phase_count <= COUPLED_LINE_RUNTIME_MAXIMUM_PHASE_COUNT ||
        throw(ArgumentError("coupled line runtime requires 1 through 12 active phases"))
    length(unique(phases)) == phase_count ||
        throw(ArgumentError("coupled line runtime phase identities must be unique"))
    port_count = 2 * phase_count
    size(model.direct_term) == (port_count, port_count) ||
        throw(ArgumentError("coupled line runtime model port count must equal twice its phase count"))
    length(model.port_order) == port_count ||
        throw(ArgumentError("coupled line runtime model port order is incomplete"))
    expected_port_order = Symbol[
        (Symbol(:sending_, phase) for phase in phases)...,
        (Symbol(:receiving_, phase) for phase in phases)...,
    ]
    model.port_order == expected_port_order ||
        throw(ArgumentError("coupled line runtime model port order does not match its phase order"))
    frequencies = Float64.(frequencies_hz)
    length(frequencies) >= 2 &&
        all(value -> isfinite(value) && value > 0.0, frequencies) &&
        issorted(frequencies) && all(diff(frequencies) .> 0.0) ||
        throw(ArgumentError("coupled line runtime source frequencies must be positive and strictly increasing"))
    maximum_frequency = last(frequencies)
    maximum_frequency * settings.timestep_s < 0.5 ||
        throw(ArgumentError("coupled line runtime timestep does not resolve its source frequency band"))

    state = model.state_matrix_per_s
    input = model.input_matrix
    output = model.output_matrix_per_s
    direct = model.direct_term
    state_count = size(state, 1)
    state_count > 0 && size(state) == (state_count, state_count) ||
        throw(ArgumentError("coupled line runtime state matrix must be nonempty and square"))
    size(input) == (state_count, port_count) &&
        size(output) == (port_count, state_count) ||
        throw(DimensionMismatch("coupled line runtime state input/output dimensions disagree"))
    all(matrix -> all(isfinite, matrix), (state, input, output, direct)) ||
        throw(ArgumentError("coupled line runtime matrices must be finite"))

    alpha = _runtime_bilinear_alpha(settings)
    identity_state = Matrix{Float64}(I, state_count, state_count)
    bilinear_matrix = alpha .* identity_state .- state
    bilinear_condition = cond(bilinear_matrix)
    isfinite(bilinear_condition) &&
        bilinear_condition <= settings.maximum_condition_number ||
        throw(ArgumentError("coupled line runtime bilinear state solve is ill-conditioned"))
    bilinear_factorization = lu(bilinear_matrix)
    state_transition = bilinear_factorization \
        (alpha .* identity_state .+ state)
    endpoint_input = bilinear_factorization \ input
    discrete_direct = direct + output * endpoint_input
    history_state_output = output * state_transition
    history_incident_output = output * endpoint_input

    identity_ports = Matrix{Float64}(I, port_count, port_count)
    wave_matrix = identity_ports + discrete_direct
    wave_condition = cond(wave_matrix)
    isfinite(wave_condition) && wave_condition <= settings.maximum_condition_number ||
        throw(ArgumentError("coupled line runtime wave companion solve is ill-conditioned"))
    wave_factorization = lu(wave_matrix)
    reference_sqrt = sqrt.(model.reference_impedance_ohm)
    reference_inverse_sqrt = inv.(reference_sqrt)
    reference_inverse_matrix = Diagonal(reference_inverse_sqrt)
    incident_from_voltage = wave_factorization \ Matrix(reference_inverse_matrix)
    incident_from_history = -(wave_factorization \ identity_ports)
    history_to_current = reference_inverse_matrix * (
        (identity_ports - discrete_direct) * incident_from_history -
        identity_ports
    )
    raw_companion = reference_inverse_matrix *
        (identity_ports - discrete_direct) * incident_from_voltage
    reciprocity_scale = max(opnorm(raw_companion), eps(Float64))
    reciprocity_error = opnorm(raw_companion - transpose(raw_companion)) /
        reciprocity_scale
    reciprocity_error <= settings.reciprocity_relative_tolerance ||
        throw(ArgumentError("coupled line runtime companion is not reciprocal"))
    companion = 0.5 .* (raw_companion .+ transpose(raw_companion))
    minimum_conductance = minimum(eigvals(Symmetric(companion)); init=Inf)
    conductance_floor = settings.passivity_tolerance *
        max(opnorm(companion), 1.0)
    minimum_conductance >= -conductance_floor ||
        throw(ArgumentError("coupled line runtime companion is active"))

    discrete_pole_magnitude = maximum(abs, eigvals(state_transition); init=0.0)
    discrete_pole_magnitude < 1.0 ||
        throw(ArgumentError("coupled line runtime discrete realization is unstable"))
    maximum_warp = maximum(
        frequency -> _runtime_frequency_warp(
            frequency,
            settings.timestep_s,
            alpha,
        ),
        frequencies;
        init=0.0,
    )
    maximum_warp <= settings.maximum_relative_frequency_warp ||
        throw(ArgumentError("coupled line runtime frequency warping exceeds policy"))
    sampled_maximum_singular_value = maximum(
        frequency -> opnorm(_runtime_discrete_response(
            state_transition,
            endpoint_input,
            output,
            direct,
            frequency,
            settings.timestep_s,
        )),
        vcat(0.0, frequencies);
        init=opnorm(direct),
    )
    sampled_maximum_singular_value <= 1.0 + settings.passivity_tolerance ||
        throw(ArgumentError("coupled line runtime discrete response is not passive"))

    preparation_signature = _runtime_sha256() do io
        println(io, COUPLED_LINE_RUNTIME_SCHEMA_VERSION)
        println(io, source_signature)
        println(io, response_signature)
        println(io, fit_signature)
        _write_runtime_values(io, uncertainty_signatures)
        _write_runtime_value(io, uncertainty_set_complete)
        println(io, model_signature)
        println(io, settings.deterministic_signature_sha256)
        _write_runtime_values(io, phases)
        _write_runtime_values(io, model.port_order)
        _write_runtime_values(io, model.reference_impedance_ohm)
        _write_runtime_values(io, frequencies)
        _write_runtime_value(io, alpha)
        for matrix in (
            state_transition,
            endpoint_input,
            output,
            direct,
            discrete_direct,
            companion,
            history_to_current,
        )
            _write_runtime_values(io, matrix)
        end
    end
    return CoupledLineRuntimePreparation(
        COUPLED_LINE_RUNTIME_SCHEMA_VERSION,
        source_signature,
        response_signature,
        fit_signature,
        uncertainty_signatures,
        uncertainty_set_complete,
        model_signature,
        phases,
        copy(model.port_order),
        copy(model.reference_impedance_ohm),
        (first(frequencies), last(frequencies)),
        settings,
        alpha,
        state_transition,
        endpoint_input,
        copy(output),
        copy(direct),
        discrete_direct,
        history_state_output,
        history_incident_output,
        incident_from_voltage,
        incident_from_history,
        history_to_current,
        companion,
        reference_sqrt,
        reference_inverse_sqrt,
        bilinear_condition,
        wave_condition,
        discrete_pole_magnitude,
        maximum_warp,
        sampled_maximum_singular_value,
        reciprocity_error,
        minimum_conductance,
        preparation_signature,
    )
end

function prepare_coupled_line_runtime(
    fit::CoupledLineFitResult,
    settings::CoupledLineRuntimeSettings,
)
    fit.certificate_after_enforcement.continuous_passivity_passed ||
        throw(ArgumentError("coupled line runtime requires a continuously passive L205 fit"))
    fit.source_signature_sha256 == fit.source_response.source_signature_sha256 ||
        throw(ArgumentError("coupled line runtime fit source signature is inconsistent"))
    fit.response_signature_sha256 == fit.source_response.response_signature_sha256 ||
        throw(ArgumentError("coupled line runtime fit response signature is inconsistent"))
    fit.model.port_order == fit.source_response.port_order ||
        throw(ArgumentError("coupled line runtime fit model and response port orders differ"))
    fit.model.reference_impedance_ohm ==
        fit.source_response.reference_impedance_ohm ||
        throw(ArgumentError("coupled line runtime fit model and response references differ"))
    return _coupled_line_runtime_preparation(
        fit.model,
        settings;
        source_signature_sha256=fit.source_signature_sha256,
        response_signature_sha256=fit.response_signature_sha256,
        fit_signature_sha256=fit.deterministic_signature_sha256,
        phase_order=fit.source_response.phase_order,
        frequencies_hz=fit.source_response.frequencies_hz,
        continuous_passivity_passed=true,
        uncertainty_alternative_fit_signatures_sha256=
            fit.uncertainty === nothing ? String[] :
            fit.uncertainty.alternative_fit_signatures_sha256,
        uncertainty_set_complete=
            fit.uncertainty === nothing ? false : fit.uncertainty.complete_set,
    )
end

mutable struct CoupledLineRuntimeState
    preparation::CoupledLineRuntimePreparation
    rational_state::Vector{Float64}
    previous_incident_wave::Vector{Float64}
    incident_wave::Vector{Float64}
    outgoing_wave::Vector{Float64}
    terminal_voltage_v::Vector{Float64}
    terminal_current_a::Vector{Float64}
    history_current_a::Vector{Float64}
    state_workspace::Vector{Float64}
    wave_history_workspace::Vector{Float64}
    current_workspace::Vector{Float64}
    incident_workspace::Vector{Float64}
    outgoing_workspace::Vector{Float64}
    previous_terminal_power_w::Float64
    terminal_power_w::Float64
    cumulative_supplied_energy_j::Float64
    minimum_cumulative_supplied_energy_j::Float64
    maximum_kcl_residual_a::Float64
    maximum_state_magnitude::Float64
    accepted_time_s::Float64
    accepted_step_count::Int
    initialization_kind::Symbol
    initialization_frequency_hz::Float64
end

function coupled_line_runtime_state(preparation::CoupledLineRuntimePreparation)
    state_count = size(preparation.state_transition, 1)
    port_count = length(preparation.port_order)
    state = CoupledLineRuntimeState(
        preparation,
        zeros(state_count),
        zeros(port_count),
        zeros(port_count),
        zeros(port_count),
        zeros(port_count),
        zeros(port_count),
        zeros(port_count),
        zeros(state_count),
        zeros(port_count),
        zeros(port_count),
        zeros(port_count),
        zeros(port_count),
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0,
        :deenergized,
        0.0,
    )
    return initialize_coupled_line_runtime_deenergized!(state)
end

function _runtime_history_output!(
    state::CoupledLineRuntimeState,
)
    preparation = state.preparation
    mul!(
        state.wave_history_workspace,
        preparation.history_state_output,
        state.rational_state,
    )
    mul!(
        state.wave_history_workspace,
        preparation.history_incident_output,
        state.previous_incident_wave,
        1.0,
        1.0,
    )
    mul!(
        state.history_current_a,
        preparation.history_to_terminal_current,
        state.wave_history_workspace,
    )
    return state.history_current_a
end

function initialize_coupled_line_runtime_deenergized!(
    state::CoupledLineRuntimeState,
)
    fill!(state.rational_state, 0.0)
    fill!(state.previous_incident_wave, 0.0)
    fill!(state.incident_wave, 0.0)
    fill!(state.outgoing_wave, 0.0)
    fill!(state.terminal_voltage_v, 0.0)
    fill!(state.terminal_current_a, 0.0)
    fill!(state.history_current_a, 0.0)
    fill!(state.state_workspace, 0.0)
    fill!(state.wave_history_workspace, 0.0)
    fill!(state.current_workspace, 0.0)
    fill!(state.incident_workspace, 0.0)
    fill!(state.outgoing_workspace, 0.0)
    state.previous_terminal_power_w = 0.0
    state.terminal_power_w = 0.0
    state.cumulative_supplied_energy_j = 0.0
    state.minimum_cumulative_supplied_energy_j = 0.0
    state.maximum_kcl_residual_a = 0.0
    state.maximum_state_magnitude = 0.0
    state.accepted_time_s = 0.0
    state.accepted_step_count = 0
    state.initialization_kind = :deenergized
    state.initialization_frequency_hz = 0.0
    return state
end

function initialize_coupled_line_runtime_sinusoidal!(
    state::CoupledLineRuntimeState,
    terminal_voltage_phasor_v::AbstractVector{<:Complex},
    frequency_hz::Real,
)
    preparation = state.preparation
    port_count = length(preparation.port_order)
    length(terminal_voltage_phasor_v) == port_count ||
        throw(ArgumentError("coupled line runtime initialization voltage count must match ports"))
    voltage_phasor = ComplexF64.(terminal_voltage_phasor_v)
    all(value -> isfinite(real(value)) && isfinite(imag(value)), voltage_phasor) ||
        throw(ArgumentError("coupled line runtime initialization phasors must be finite"))
    frequency = _checked_positive_runtime_value(
        frequency_hz,
        "coupled line runtime initialization frequency_hz",
    )
    frequency <= last(preparation.source_frequency_band_hz) ||
        throw(ArgumentError("coupled line runtime initialization frequency exceeds source band"))
    scattering = coupled_line_runtime_discrete_response(preparation, frequency)
    terminal_admittance = coupled_line_scattering_to_admittance(
        scattering,
        preparation.reference_impedance_ohm,
    )
    current_phasor = terminal_admittance * voltage_phasor
    incident_phasor = 0.5 .* (
        preparation.reference_impedance_inverse_sqrt_per_ohm_sqrt .*
            voltage_phasor .+
        preparation.reference_impedance_sqrt_ohm_sqrt .*
            current_phasor
    )
    angle = 2.0 * pi * frequency * preparation.settings.timestep_s
    z = cis(angle)
    state_count = size(preparation.state_transition, 1)
    identity_state = Matrix{ComplexF64}(I, state_count, state_count)
    state_phasor = (z .* identity_state .- preparation.state_transition) \
        (ComplexF64.(preparation.endpoint_input) * ((z + 1.0) .* incident_phasor))
    outgoing_phasor = ComplexF64.(preparation.output_matrix) * state_phasor .+
        ComplexF64.(preparation.continuous_direct_term) * incident_phasor

    state.rational_state .= real.(state_phasor)
    state.previous_incident_wave .= real.(incident_phasor)
    state.incident_wave .= real.(incident_phasor)
    state.outgoing_wave .= real.(outgoing_phasor)
    state.terminal_voltage_v .= real.(voltage_phasor)
    state.terminal_current_a .= real.(current_phasor)
    state.previous_terminal_power_w = dot(
        state.terminal_voltage_v,
        state.terminal_current_a,
    )
    state.terminal_power_w = state.previous_terminal_power_w
    state.cumulative_supplied_energy_j = 0.0
    state.minimum_cumulative_supplied_energy_j = 0.0
    state.maximum_kcl_residual_a = 0.0
    state.maximum_state_magnitude = maximum(abs, state.rational_state; init=0.0)
    state.accepted_time_s = 0.0
    state.accepted_step_count = 0
    state.initialization_kind = :sinusoidal_discrete_operating_point
    state.initialization_frequency_hz = frequency
    _runtime_history_output!(state)
    return state
end

function accept_coupled_line_runtime_step!(
    state::CoupledLineRuntimeState,
    terminal_voltage_v::AbstractVector{<:Real},
)
    preparation = state.preparation
    port_count = length(preparation.port_order)
    length(terminal_voltage_v) == port_count ||
        throw(ArgumentError("coupled line runtime terminal voltage count must match ports"))
    voltage = Float64.(terminal_voltage_v)
    all(isfinite, voltage) ||
        throw(ArgumentError("coupled line runtime terminal voltages must be finite"))

    _runtime_history_output!(state)
    mul!(
        state.incident_workspace,
        preparation.incident_from_voltage,
        voltage,
    )
    mul!(
        state.incident_workspace,
        preparation.incident_from_history,
        state.wave_history_workspace,
        1.0,
        1.0,
    )
    mul!(
        state.state_workspace,
        preparation.state_transition,
        state.rational_state,
    )
    mul!(
        state.state_workspace,
        preparation.endpoint_input,
        state.incident_workspace,
        1.0,
        1.0,
    )
    mul!(
        state.state_workspace,
        preparation.endpoint_input,
        state.previous_incident_wave,
        1.0,
        1.0,
    )
    mul!(
        state.outgoing_workspace,
        preparation.output_matrix,
        state.state_workspace,
    )
    mul!(
        state.outgoing_workspace,
        preparation.continuous_direct_term,
        state.incident_workspace,
        1.0,
        1.0,
    )
    @inbounds for port in 1:port_count
        state.current_workspace[port] =
            preparation.reference_impedance_inverse_sqrt_per_ohm_sqrt[port] *
            (state.incident_workspace[port] - state.outgoing_workspace[port])
    end
    mul!(
        state.wave_history_workspace,
        preparation.companion_admittance_s,
        voltage,
    )
    state.wave_history_workspace .+= state.history_current_a
    kcl_residual = maximum(
        abs,
        state.current_workspace .- state.wave_history_workspace;
        init=0.0,
    )
    all(isfinite, state.state_workspace) &&
        all(isfinite, state.incident_workspace) &&
        all(isfinite, state.outgoing_workspace) &&
        all(isfinite, state.current_workspace) ||
        throw(ArgumentError("coupled line runtime accepted state became nonfinite"))
    kcl_residual <= preparation.settings.kcl_absolute_tolerance_a +
        128.0 * eps(Float64) *
        max(maximum(abs, state.current_workspace; init=0.0), 1.0) ||
        throw(ArgumentError("coupled line runtime companion current failed its KCL reconstruction"))

    power = dot(voltage, state.current_workspace)
    energy_increment = 0.5 * preparation.settings.timestep_s *
        (state.previous_terminal_power_w + power)
    cumulative_energy = state.cumulative_supplied_energy_j + energy_increment
    energy_roundoff = 256.0 * eps(Float64) * max(
        abs(state.cumulative_supplied_energy_j),
        abs(energy_increment),
        preparation.settings.timestep_s * abs(power),
        1.0,
    )
    state.initialization_kind == :deenergized &&
        cumulative_energy <
            -preparation.settings.energy_absolute_tolerance_j - energy_roundoff &&
        throw(ArgumentError(
            "coupled line runtime passive energy balance became negative",
        ))

    state.terminal_voltage_v .= voltage
    state.incident_wave .= state.incident_workspace
    state.outgoing_wave .= state.outgoing_workspace
    state.terminal_current_a .= state.current_workspace
    state.cumulative_supplied_energy_j = cumulative_energy
    state.minimum_cumulative_supplied_energy_j = min(
        state.minimum_cumulative_supplied_energy_j,
        state.cumulative_supplied_energy_j,
    )
    state.previous_terminal_power_w = power
    state.terminal_power_w = power
    state.maximum_kcl_residual_a = max(
        state.maximum_kcl_residual_a,
        kcl_residual,
    )
    state.maximum_state_magnitude = max(
        state.maximum_state_magnitude,
        maximum(abs, state.state_workspace; init=0.0),
    )
    state.rational_state .= state.state_workspace
    state.previous_incident_wave .= state.incident_workspace
    state.accepted_step_count += 1
    state.accepted_time_s =
        state.accepted_step_count * preparation.settings.timestep_s
    _runtime_history_output!(state)
    return state
end

struct CoupledLineRuntimeSnapshot
    schema_version::Int
    preparation_signature_sha256::String
    rational_state::Vector{Float64}
    previous_incident_wave::Vector{Float64}
    incident_wave::Vector{Float64}
    outgoing_wave::Vector{Float64}
    terminal_voltage_v::Vector{Float64}
    terminal_current_a::Vector{Float64}
    history_current_a::Vector{Float64}
    previous_terminal_power_w::Float64
    terminal_power_w::Float64
    cumulative_supplied_energy_j::Float64
    minimum_cumulative_supplied_energy_j::Float64
    maximum_kcl_residual_a::Float64
    maximum_state_magnitude::Float64
    accepted_time_s::Float64
    accepted_step_count::Int
    initialization_kind::Symbol
    initialization_frequency_hz::Float64
    deterministic_signature_sha256::String
end

function _coupled_line_runtime_snapshot_signature(
    preparation_signature,
    rational_state,
    previous_incident_wave,
    incident_wave,
    outgoing_wave,
    terminal_voltage,
    terminal_current,
    history_current,
    scalars,
    accepted_step_count,
    initialization_kind,
)
    return _runtime_sha256() do io
        println(io, COUPLED_LINE_RUNTIME_SCHEMA_VERSION)
        println(io, preparation_signature)
        for values in (
            rational_state,
            previous_incident_wave,
            incident_wave,
            outgoing_wave,
            terminal_voltage,
            terminal_current,
            history_current,
        )
            _write_runtime_values(io, values)
        end
        _write_runtime_values(io, scalars)
        _write_runtime_value(io, accepted_step_count)
        _write_runtime_value(io, initialization_kind)
    end
end

function coupled_line_runtime_snapshot(state::CoupledLineRuntimeState)
    scalars = (
        state.previous_terminal_power_w,
        state.terminal_power_w,
        state.cumulative_supplied_energy_j,
        state.minimum_cumulative_supplied_energy_j,
        state.maximum_kcl_residual_a,
        state.maximum_state_magnitude,
        state.accepted_time_s,
        state.initialization_frequency_hz,
    )
    signature = _coupled_line_runtime_snapshot_signature(
        state.preparation.deterministic_signature_sha256,
        state.rational_state,
        state.previous_incident_wave,
        state.incident_wave,
        state.outgoing_wave,
        state.terminal_voltage_v,
        state.terminal_current_a,
        state.history_current_a,
        scalars,
        state.accepted_step_count,
        state.initialization_kind,
    )
    return CoupledLineRuntimeSnapshot(
        COUPLED_LINE_RUNTIME_SCHEMA_VERSION,
        state.preparation.deterministic_signature_sha256,
        copy(state.rational_state),
        copy(state.previous_incident_wave),
        copy(state.incident_wave),
        copy(state.outgoing_wave),
        copy(state.terminal_voltage_v),
        copy(state.terminal_current_a),
        copy(state.history_current_a),
        state.previous_terminal_power_w,
        state.terminal_power_w,
        state.cumulative_supplied_energy_j,
        state.minimum_cumulative_supplied_energy_j,
        state.maximum_kcl_residual_a,
        state.maximum_state_magnitude,
        state.accepted_time_s,
        state.accepted_step_count,
        state.initialization_kind,
        state.initialization_frequency_hz,
        signature,
    )
end

function _coupled_line_runtime_snapshot_dictionary(
    snapshot::CoupledLineRuntimeSnapshot,
)
    return Dict{String,Any}(
        "schema" => "aimora.coupled_line_runtime_snapshot.v1",
        "schema_version" => snapshot.schema_version,
        "preparation_signature_sha256" => snapshot.preparation_signature_sha256,
        "rational_state" => snapshot.rational_state,
        "previous_incident_wave" => snapshot.previous_incident_wave,
        "incident_wave" => snapshot.incident_wave,
        "outgoing_wave" => snapshot.outgoing_wave,
        "terminal_voltage_v" => snapshot.terminal_voltage_v,
        "terminal_current_a" => snapshot.terminal_current_a,
        "history_current_a" => snapshot.history_current_a,
        "previous_terminal_power_w" => snapshot.previous_terminal_power_w,
        "terminal_power_w" => snapshot.terminal_power_w,
        "cumulative_supplied_energy_j" => snapshot.cumulative_supplied_energy_j,
        "minimum_cumulative_supplied_energy_j" =>
            snapshot.minimum_cumulative_supplied_energy_j,
        "maximum_kcl_residual_a" => snapshot.maximum_kcl_residual_a,
        "maximum_state_magnitude" => snapshot.maximum_state_magnitude,
        "accepted_time_s" => snapshot.accepted_time_s,
        "accepted_step_count" => snapshot.accepted_step_count,
        "initialization_kind" => String(snapshot.initialization_kind),
        "initialization_frequency_hz" => snapshot.initialization_frequency_hz,
        "deterministic_signature_sha256" =>
            snapshot.deterministic_signature_sha256,
    )
end

function _coupled_line_runtime_snapshot_from_dictionary(data)
    get(data, "schema", "") == "aimora.coupled_line_runtime_snapshot.v1" ||
        throw(ArgumentError("coupled line runtime snapshot schema is unsupported"))
    try
        return CoupledLineRuntimeSnapshot(
            Int(data["schema_version"]),
            String(data["preparation_signature_sha256"]),
            Float64.(data["rational_state"]),
            Float64.(data["previous_incident_wave"]),
            Float64.(data["incident_wave"]),
            Float64.(data["outgoing_wave"]),
            Float64.(data["terminal_voltage_v"]),
            Float64.(data["terminal_current_a"]),
            Float64.(data["history_current_a"]),
            Float64(data["previous_terminal_power_w"]),
            Float64(data["terminal_power_w"]),
            Float64(data["cumulative_supplied_energy_j"]),
            Float64(data["minimum_cumulative_supplied_energy_j"]),
            Float64(data["maximum_kcl_residual_a"]),
            Float64(data["maximum_state_magnitude"]),
            Float64(data["accepted_time_s"]),
            Int(data["accepted_step_count"]),
            Symbol(data["initialization_kind"]),
            Float64(data["initialization_frequency_hz"]),
            String(data["deterministic_signature_sha256"]),
        )
    catch error
        error isa ArgumentError && rethrow()
        throw(ArgumentError("coupled line runtime snapshot is malformed"))
    end
end

"""Write one public, portable, integrity-bound coupled-line runtime snapshot."""
function write_coupled_line_runtime_snapshot(
    path::AbstractString,
    snapshot::CoupledLineRuntimeSnapshot,
)
    output_path = abspath(path)
    mkpath(dirname(output_path))
    mktemp(dirname(output_path)) do temporary_path, io
        TOML.print(io, _coupled_line_runtime_snapshot_dictionary(snapshot); sorted=true)
        close(io)
        mv(temporary_path, output_path; force=true)
    end
    return output_path
end

"""Read one public coupled-line snapshot; integrity is rechecked on restore."""
function read_coupled_line_runtime_snapshot(path::AbstractString)
    isfile(path) || throw(ArgumentError("coupled line runtime snapshot file does not exist"))
    data = try
        TOML.parsefile(path)
    catch
        throw(ArgumentError("coupled line runtime snapshot is not valid TOML"))
    end
    return _coupled_line_runtime_snapshot_from_dictionary(data)
end

function restore_coupled_line_runtime_snapshot!(
    state::CoupledLineRuntimeState,
    snapshot::CoupledLineRuntimeSnapshot,
)
    snapshot.schema_version == COUPLED_LINE_RUNTIME_SCHEMA_VERSION ||
        throw(ArgumentError("coupled line runtime snapshot schema is unsupported"))
    snapshot.preparation_signature_sha256 ==
        state.preparation.deterministic_signature_sha256 ||
        throw(ArgumentError("coupled line runtime snapshot preparation is stale"))
    expected_state_count = length(state.rational_state)
    expected_port_count = length(state.previous_incident_wave)
    length(snapshot.rational_state) == expected_state_count ||
        throw(ArgumentError("coupled line runtime snapshot state count is incompatible"))
    for values in (
        snapshot.previous_incident_wave,
        snapshot.incident_wave,
        snapshot.outgoing_wave,
        snapshot.terminal_voltage_v,
        snapshot.terminal_current_a,
        snapshot.history_current_a,
    )
        length(values) == expected_port_count ||
            throw(ArgumentError("coupled line runtime snapshot port count is incompatible"))
    end
    all(values -> all(isfinite, values), (
        snapshot.rational_state,
        snapshot.previous_incident_wave,
        snapshot.incident_wave,
        snapshot.outgoing_wave,
        snapshot.terminal_voltage_v,
        snapshot.terminal_current_a,
        snapshot.history_current_a,
    )) || throw(ArgumentError("coupled line runtime snapshot contains nonfinite state"))
    scalars = (
        snapshot.previous_terminal_power_w,
        snapshot.terminal_power_w,
        snapshot.cumulative_supplied_energy_j,
        snapshot.minimum_cumulative_supplied_energy_j,
        snapshot.maximum_kcl_residual_a,
        snapshot.maximum_state_magnitude,
        snapshot.accepted_time_s,
        snapshot.initialization_frequency_hz,
    )
    all(isfinite, scalars) ||
        throw(ArgumentError("coupled line runtime snapshot contains nonfinite diagnostics"))
    snapshot.accepted_step_count >= 0 ||
        throw(ArgumentError("coupled line runtime snapshot step count is negative"))
    expected_time = snapshot.accepted_step_count *
        state.preparation.settings.timestep_s
    abs(snapshot.accepted_time_s - expected_time) <=
        64.0 * eps(Float64) * max(abs(expected_time), 1.0) ||
        throw(ArgumentError("coupled line runtime snapshot time and step count disagree"))
    expected_signature = _coupled_line_runtime_snapshot_signature(
        snapshot.preparation_signature_sha256,
        snapshot.rational_state,
        snapshot.previous_incident_wave,
        snapshot.incident_wave,
        snapshot.outgoing_wave,
        snapshot.terminal_voltage_v,
        snapshot.terminal_current_a,
        snapshot.history_current_a,
        scalars,
        snapshot.accepted_step_count,
        snapshot.initialization_kind,
    )
    snapshot.deterministic_signature_sha256 == expected_signature ||
        throw(ArgumentError("coupled line runtime snapshot integrity signature does not match"))

    state.rational_state .= snapshot.rational_state
    state.previous_incident_wave .= snapshot.previous_incident_wave
    state.incident_wave .= snapshot.incident_wave
    state.outgoing_wave .= snapshot.outgoing_wave
    state.terminal_voltage_v .= snapshot.terminal_voltage_v
    state.terminal_current_a .= snapshot.terminal_current_a
    state.history_current_a .= snapshot.history_current_a
    state.previous_terminal_power_w = snapshot.previous_terminal_power_w
    state.terminal_power_w = snapshot.terminal_power_w
    state.cumulative_supplied_energy_j = snapshot.cumulative_supplied_energy_j
    state.minimum_cumulative_supplied_energy_j =
        snapshot.minimum_cumulative_supplied_energy_j
    state.maximum_kcl_residual_a = snapshot.maximum_kcl_residual_a
    state.maximum_state_magnitude = snapshot.maximum_state_magnitude
    state.accepted_time_s = snapshot.accepted_time_s
    state.accepted_step_count = snapshot.accepted_step_count
    state.initialization_kind = snapshot.initialization_kind
    state.initialization_frequency_hz = snapshot.initialization_frequency_hz
    fill!(state.state_workspace, 0.0)
    fill!(state.wave_history_workspace, 0.0)
    fill!(state.current_workspace, 0.0)
    fill!(state.incident_workspace, 0.0)
    fill!(state.outgoing_workspace, 0.0)
    return state
end

function coupled_line_runtime_diagnostics(state::CoupledLineRuntimeState)
    preparation = state.preparation
    return (
        source_signature_sha256=preparation.source_signature_sha256,
        response_signature_sha256=preparation.response_signature_sha256,
        fit_signature_sha256=preparation.fit_signature_sha256,
        uncertainty_alternative_fit_signatures_sha256=
            copy(preparation.uncertainty_alternative_fit_signatures_sha256),
        uncertainty_set_complete=preparation.uncertainty_set_complete,
        unknown_uncertainty_explicit=!preparation.uncertainty_set_complete,
        model_signature_sha256=preparation.model_signature_sha256,
        settings_signature_sha256=
            preparation.settings.deterministic_signature_sha256,
        runtime_signature_sha256=preparation.deterministic_signature_sha256,
        schema_version=preparation.schema_version,
        phase_order=copy(preparation.phase_order),
        port_order=copy(preparation.port_order),
        timestep_s=preparation.settings.timestep_s,
        bilinear_alpha_per_s=preparation.bilinear_alpha_per_s,
        state_count=length(state.rational_state),
        port_count=length(state.terminal_voltage_v),
        accepted_time_s=state.accepted_time_s,
        accepted_step_count=state.accepted_step_count,
        initialization_kind=state.initialization_kind,
        initialization_frequency_hz=state.initialization_frequency_hz,
        terminal_power_w=state.terminal_power_w,
        cumulative_supplied_energy_j=state.cumulative_supplied_energy_j,
        minimum_cumulative_supplied_energy_j=
            state.minimum_cumulative_supplied_energy_j,
        maximum_kcl_residual_a=state.maximum_kcl_residual_a,
        maximum_state_magnitude=state.maximum_state_magnitude,
        maximum_discrete_pole_magnitude=
            preparation.maximum_discrete_pole_magnitude,
        maximum_relative_frequency_warp=
            preparation.maximum_relative_frequency_warp,
        sampled_maximum_scattering_singular_value=
            preparation.sampled_maximum_scattering_singular_value,
        minimum_companion_conductance_eigenvalue_s=
            preparation.minimum_companion_conductance_eigenvalue_s,
        passive_energy_balance_passed=
            state.initialization_kind != :deenergized ||
            state.minimum_cumulative_supplied_energy_j >=
                -preparation.settings.energy_absolute_tolerance_j,
        coupled_phase_domain_runtime_executed=state.accepted_step_count > 0,
        ulm_compatible_architecture=true,
        ulm_file_compatibility_claimed=false,
        atp_or_pscad_equivalence_claimed=false,
    )
end

function coupled_line_runtime_report_text(state::CoupledLineRuntimeState)
    diagnostics = coupled_line_runtime_diagnostics(state)
    io = IOBuffer()
    println(io, "AIMORA coupled frequency-dependent line runtime")
    for key in keys(diagnostics)
        value = getfield(diagnostics, key)
        rendered = value isa AbstractVector ? join(value, ",") : string(value)
        println(io, key, "=", rendered)
    end
    return String(take!(io))
end

end
