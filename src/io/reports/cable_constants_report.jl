using ..CableConstantsStudy:
    CableConstantsStudyResult
using ..Lines:
    NestedCableTransientLineState

export CableConstantsReportArtifact,
       cable_constants_report_artifact,
       cable_constants_report_text,
       write_cable_constants_report_text,
       NestedCableTransientReportArtifact,
       nested_cable_transient_report_artifact,
       nested_cable_transient_report_text,
       write_nested_cable_transient_report_text

struct CableConstantsReportArtifact
    title::String
    result::CableConstantsStudyResult
end

struct NestedCableTransientReportArtifact
    title::String
    line_length_m::Float64
    dt_s::Float64
    target_frequency_hz::Float64
    fit_max_abs_error::Float64
    fit_relative_max_abs_error::Float64
    pole_term_count::Int
    history_update_count::Int
    retained_current_injection_a::Vector{Float64}
    nodal_rhs_a::Vector{Float64}
    solved_voltage_v::Vector{Float64}
    sending_phase_current_a::Vector{ComplexF64}
    receiving_phase_current_a::Vector{ComplexF64}
    modal_history_current_a::Vector{ComplexF64}
    nodal_kcl_max_abs_error_a::Float64
    real_current_projection_max_imag_abs_a::Float64
    physical_checks_passed::Bool
    deferred_effects::Tuple{Vararg{Symbol}}
end

function nested_cable_transient_report_artifact(
    line_state::NestedCableTransientLineState,
    timestep_result::NamedTuple,
    nodal_conductance_s::AbstractVector;
    title::AbstractString = "AIMORA nested cable transient report",
)
    line_result = timestep_result.frequency_dependent_line_result
    line_result === nothing &&
        throw(ArgumentError("timestep result contains no frequency-dependent line update"))
    line_result.nested_cable_frequency_state_consumed ||
        throw(ArgumentError("timestep result did not consume a nested cable frequency state"))
    length(line_result.recursive_convolution_updates) == 1 ||
        throw(ArgumentError("nested cable transient report currently requires one line state"))
    update = only(line_result.recursive_convolution_updates)
    length(line_result.phase_current_injections) == 1 ||
        throw(ArgumentError("nested cable transient report requires one retained current injection"))
    injection = only(line_result.phase_current_injections)
    conductance = Float64.(nodal_conductance_s)
    length(conductance) == length(timestep_result.voltage_pu) ||
        throw(ArgumentError("nodal conductance count must match solved voltage count"))
    all(value -> isfinite(value) && value > 0.0, conductance) ||
        throw(ArgumentError("nodal conductances must be finite and positive"))
    rhs = Float64.(line_result.rhs_after_values)
    voltage = Float64.(timestep_result.voltage_pu)
    kcl_error = maximum(abs.(conductance .* voltage .- rhs); init = 0.0)
    return NestedCableTransientReportArtifact(
        String(title),
        line_state.line_length_m,
        line_state.dt_s,
        only(line_result.target_frequency_hz_values),
        line_state.response_fit.max_abs_error,
        line_state.fit_relative_max_abs_error,
        line_state.response_fit.term_count,
        update.recursive_convolution_update_count,
        Float64.(injection.rhs_after_values .- injection.rhs_before_values),
        rhs,
        voltage,
        ComplexF64.(update.sending_phase_current),
        ComplexF64.(update.receiving_phase_current),
        ComplexF64.(update.convolution_modal_current),
        kcl_error,
        injection.real_current_projection_max_imag_abs,
        line_state.physical_checks_passed &&
            line_result.state_mutated &&
            isempty(line_result.deferred_effects),
        line_result.deferred_effects,
    )
end

function cable_constants_report_artifact(
    result::CableConstantsStudyResult;
    title::AbstractString = "AIMORA cable constants report",
)
    return CableConstantsReportArtifact(String(title), result)
end

function _write_cable_report_values(io, label::AbstractString, values)
    print(io, label)
    for value in values
        @printf(io, " %.12e", value)
    end
    println(io)
    return io
end

function _write_cable_report_complex_matrix(io, label::AbstractString, matrix)
    println(io, label)
    for row in axes(matrix, 1)
        print(io, "row ", row)
        for value in matrix[row, :]
            @printf(io, " %.12e %.12e", real(value), imag(value))
        end
        println(io)
    end
    return io
end

function _write_cable_constants_report(io, artifact::CableConstantsReportArtifact)
    result = artifact.result
    case = result.case
    state = result.pipe_sheath_state
    println(io, artifact.title)
    println(io, "source ", result.source)
    println(io, "cable_kind_code ", case.cable_kind_code)
    println(io, "surface_position_code ", case.surface_position_code)
    println(io, "phase_count ", case.phase_count)
    println(io, "selected_grounded_conductor_count ", state.selected_grounded_conductor_count)
    println(io, "geometry_conductor_count ", result.geometry.conductor_count)
    for phase in 1:case.phase_count
        println(io, "phase ", phase, " layer_count ", case.layer_counts[phase])
        _write_cable_report_values(io, "boundary_radii_m", case.boundary_radii_m[phase, :])
        _write_cable_report_values(io, "resistivity_ohm_m", case.resistivity_ohm_m[phase, :])
        _write_cable_report_values(
            io,
            "conductor_relative_permeability",
            case.conductor_relative_permeability[phase, :],
        )
        _write_cable_report_values(
            io,
            "insulation_relative_permeability",
            case.insulation_relative_permeability[phase, :],
        )
        _write_cable_report_values(
            io,
            "insulation_relative_permittivity",
            case.insulation_relative_permittivity[phase, :],
        )
        _write_cable_report_values(
            io,
            "dielectric_wave_speed_m_per_s",
            state.layer_wave_speeds_m_per_s[phase, :],
        )
        @printf(
            io,
            "position_m depth %.12e horizontal %.12e\n",
            case.depths_m[phase],
            case.horizontal_positions_m[phase],
        )
    end
    for (index, schedule) in enumerate(result.frequency_schedules)
        println(io, "frequency_schedule ", index)
        @printf(io, "earth_resistivity_ohm_m %.12e\n", schedule.final_earth_resistivity_ohm_m)
        @printf(io, "start_frequency_hz %.12e\n", schedule.start_frequency_hz)
        println(io, "decade_count ", schedule.decade_count)
        println(io, "points_per_decade ", schedule.points_per_decade)
        @printf(io, "distance_m %.12e\n", schedule.initial_distance_m)
        _write_cable_report_values(io, "frequencies_hz", schedule.frequencies_hz)
    end
    for (index, frequency_state) in enumerate(result.frequency_states)
        println(io, "frequency_state ", index)
        @printf(io, "frequency_hz %.12e\n", frequency_state.frequency_hz)
        _write_cable_report_complex_matrix(
            io,
            "series_impedance_matrix_ohm_per_m",
            frequency_state.series_impedance_matrix_ohm_per_m,
        )
        _write_cable_report_complex_matrix(
            io,
            "shunt_admittance_matrix_s_per_m",
            frequency_state.shunt_admittance_matrix_s_per_m,
        )
        for mode in eachindex(frequency_state.modal_series_impedance_ohm_per_m)
            @printf(
                io,
                "mode %d attenuation_db_per_km %.12e velocity_m_per_s %.12e modal_series_ohm_per_m %.12e %.12e modal_shunt_s_per_m %.12e %.12e characteristic_impedance_ohm %.12e %.12e characteristic_admittance_s %.12e %.12e\n",
                mode,
                frequency_state.modal_attenuation_db_per_km[mode],
                frequency_state.modal_velocity_m_per_s[mode],
                real(frequency_state.modal_series_impedance_ohm_per_m[mode]),
                imag(frequency_state.modal_series_impedance_ohm_per_m[mode]),
                real(frequency_state.modal_shunt_admittance_s_per_m[mode]),
                imag(frequency_state.modal_shunt_admittance_s_per_m[mode]),
                real(frequency_state.modal_characteristic_impedance_ohm[mode]),
                imag(frequency_state.modal_characteristic_impedance_ohm[mode]),
                real(frequency_state.modal_characteristic_admittance_s[mode]),
                imag(frequency_state.modal_characteristic_admittance_s[mode]),
            )
        end
    end
    println(io, "physical_checks_passed ", result.physical_checks_passed)
    println(io, "deferred_effects ", join(string.(result.deferred_effects), ","))
    println(io, "END CABLE CONSTANTS")
    return io
end

function cable_constants_report_text(artifact::CableConstantsReportArtifact)
    io = IOBuffer()
    _write_cable_constants_report(io, artifact)
    return String(take!(io))
end

cable_constants_report_text(result::CableConstantsStudyResult; kwargs...) =
    cable_constants_report_text(cable_constants_report_artifact(result; kwargs...))

function write_cable_constants_report_text(
    path::AbstractString,
    artifact::CableConstantsReportArtifact,
)
    open(path, "w") do io
        _write_cable_constants_report(io, artifact)
    end
    return String(path)
end

function write_cable_constants_report_text(
    path::AbstractString,
    result::CableConstantsStudyResult;
    kwargs...,
)
    return write_cable_constants_report_text(
        path,
        cable_constants_report_artifact(result; kwargs...),
    )
end

function _write_nested_cable_transient_report(
    io,
    artifact::NestedCableTransientReportArtifact,
)
    println(io, artifact.title)
    @printf(io, "line_length_m %.12e\n", artifact.line_length_m)
    @printf(io, "dt_s %.12e\n", artifact.dt_s)
    @printf(io, "target_frequency_hz %.12e\n", artifact.target_frequency_hz)
    @printf(io, "fit_max_abs_error %.12e\n", artifact.fit_max_abs_error)
    @printf(io, "fit_relative_max_abs_error %.12e\n", artifact.fit_relative_max_abs_error)
    println(io, "pole_term_count ", artifact.pole_term_count)
    println(io, "history_update_count ", artifact.history_update_count)
    _write_cable_report_values(
        io,
        "retained_current_injection_a",
        artifact.retained_current_injection_a,
    )
    _write_cable_report_values(io, "nodal_rhs_a", artifact.nodal_rhs_a)
    _write_cable_report_values(io, "solved_voltage_v", artifact.solved_voltage_v)
    _write_cable_report_values(
        io,
        "sending_phase_current_re_a",
        real.(artifact.sending_phase_current_a),
    )
    _write_cable_report_values(
        io,
        "sending_phase_current_im_a",
        imag.(artifact.sending_phase_current_a),
    )
    _write_cable_report_values(
        io,
        "receiving_phase_current_re_a",
        real.(artifact.receiving_phase_current_a),
    )
    _write_cable_report_values(
        io,
        "receiving_phase_current_im_a",
        imag.(artifact.receiving_phase_current_a),
    )
    _write_cable_report_values(
        io,
        "modal_history_current_re_a",
        real.(artifact.modal_history_current_a),
    )
    _write_cable_report_values(
        io,
        "modal_history_current_im_a",
        imag.(artifact.modal_history_current_a),
    )
    @printf(io, "nodal_kcl_max_abs_error_a %.12e\n", artifact.nodal_kcl_max_abs_error_a)
    @printf(
        io,
        "real_current_projection_max_imag_abs_a %.12e\n",
        artifact.real_current_projection_max_imag_abs_a,
    )
    println(io, "physical_checks_passed ", artifact.physical_checks_passed)
    println(io, "deferred_effects ", join(string.(artifact.deferred_effects), ","))
    println(io, "END NESTED CABLE TRANSIENT")
    return io
end

function nested_cable_transient_report_text(artifact::NestedCableTransientReportArtifact)
    io = IOBuffer()
    _write_nested_cable_transient_report(io, artifact)
    return String(take!(io))
end

function write_nested_cable_transient_report_text(
    path::AbstractString,
    artifact::NestedCableTransientReportArtifact,
)
    open(path, "w") do io
        _write_nested_cable_transient_report(io, artifact)
    end
    return String(path)
end
