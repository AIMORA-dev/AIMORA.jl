module CableConstantsStudy

using ..DeckParser:
    DeckCableConstantsCase,
    DeckParseResult,
    assert_deck_valid!,
    deck_cable_constants_cases
using ..Lines:
    CableFrequencyScanLoopSchedule,
    CableGeometryConductor,
    CableGeometryConstants,
    LineStepResponseExponentialFitResult,
    NestedCableFrequencyState,
    CablePipeSheathDerivedState,
    cable_frequency_scan_loop_schedule,
    cable_geometry_constants,
    nested_cable_semlyen_frequency_dependent_line_from_fit,
    nested_cable_frequency_state,
    cable_pipe_sheath_derived_state

export CableConstantsStudyResult,
       build_semlyen_frequency_dependent_line,
       run_cable_constants_study

struct CableConstantsStudyResult
    source::String
    case::DeckCableConstantsCase
    pipe_sheath_state::CablePipeSheathDerivedState
    geometry::CableGeometryConstants
    frequency_schedules::Vector{CableFrequencyScanLoopSchedule}
    frequency_states::Vector{NestedCableFrequencyState}
    physical_checks_passed::Bool
    deferred_effects::Tuple{Vararg{Symbol}}
end

function _cable_study_frequency_states(
    case::DeckCableConstantsCase,
    state::CablePipeSheathDerivedState,
    geometry::CableGeometryConstants,
    schedules::AbstractVector{CableFrequencyScanLoopSchedule},
)
    rows = NestedCableFrequencyState[]
    for schedule in schedules, frequency in schedule.frequencies_hz
        push!(
            rows,
            nested_cable_frequency_state(
                state,
                geometry,
                case.boundary_radii_m,
                case.resistivity_ohm_m,
                case.conductor_relative_permeability,
                case.insulation_relative_permeability,
                case.insulation_relative_permittivity,
                frequency,
                schedule.final_earth_resistivity_ohm_m;
                layer_counts = case.layer_counts,
            ),
        )
    end
    return rows
end

function _cable_study_pipe_sheath_state(case::DeckCableConstantsCase)
    conductor_count = sum(case.layer_counts)
    phase_and_sheath_count = case.phase_count + count(>=(2), case.layer_counts)
    return cable_pipe_sheath_derived_state(
        cable_kind_code = case.cable_kind_code,
        surface_position_code = case.surface_position_code,
        conductor_count = conductor_count,
        active_phase_count = case.phase_count,
        pipe_count = case.pipe_count,
        active_phase_count_without_pipe = phase_and_sheath_count,
        requested_grounded_count = case.grounding_selector,
        grounded_flags = ones(Int, conductor_count),
        layer_counts = case.layer_counts,
        boundary_radii_m = case.boundary_radii_m,
        resistivity_ohm_m = case.resistivity_ohm_m,
        relative_permeability = case.conductor_relative_permeability,
        relative_permittivity = case.insulation_relative_permittivity,
        pipe_radii_m = ones(3),
        pipe_resistivity_ohm_m = 1.0,
        pipe_relative_permeability = 1.0,
        pipe_inner_insulator_relative_permittivity = 1.0,
        pipe_outer_insulator_relative_permittivity = 1.0,
        conductor_depths_m = case.depths_m,
        conductor_distances_m = case.horizontal_positions_m,
        conductor_pipe_center_distances_m = zeros(case.phase_count),
        conductor_angles_rad = zeros(case.phase_count),
    )
end

function _cable_study_geometry(case::DeckCableConstantsCase)
    conductors = CableGeometryConductor[
        CableGeometryConductor(
            case.boundary_radii_m[phase, 7],
            case.horizontal_positions_m[phase],
            case.depths_m[phase];
            resistivity_ohm_m = case.resistivity_ohm_m[phase, 1],
            relative_permittivity = case.insulation_relative_permittivity[phase, 3],
            relative_permeability = case.conductor_relative_permeability[phase, 1],
        ) for phase in 1:case.phase_count
    ]
    return cable_geometry_constants(
        conductors;
        phase_conductor_counts = ones(Int, case.phase_count),
        grounded_conductor_count = 0,
    )
end

function _cable_study_frequency_schedules(case::DeckCableConstantsCase)
    return CableFrequencyScanLoopSchedule[
        cable_frequency_scan_loop_schedule(
            card.start_frequency_hz,
            card.decade_count,
            card.points_per_decade;
            earth_resistivity_ohm_m = card.earth_resistivity_ohm_m,
            distance_m = card.distance_m,
            steady_state_frequency_hz = 60.0,
            card_output_flag = card.card_output_flag,
            transform_flag = card.transform_flag,
            modal_output_enabled = case.modal_output_flag != 0,
        ) for card in case.frequency_cards
    ]
end

function _cable_study_physical_checks(
    case::DeckCableConstantsCase,
    state::CablePipeSheathDerivedState,
    geometry::CableGeometryConstants,
    schedules::AbstractVector{CableFrequencyScanLoopSchedule},
)
    radii_ordered = all(
        phase -> all(diff(case.boundary_radii_m[phase, :]) .> 0.0),
        1:case.phase_count,
    )
    positive_materials = all(>(0.0), case.resistivity_ohm_m) &&
        all(>(0.0), case.conductor_relative_permeability) &&
        all(>(0.0), case.insulation_relative_permeability) &&
        all(>(0.0), case.insulation_relative_permittivity)
    symmetric_geometry = isapprox(
        geometry.direct_distance_m,
        transpose(geometry.direct_distance_m);
        atol = 1.0e-14,
        rtol = 1.0e-14,
    ) && isapprox(
        geometry.image_distance_m,
        transpose(geometry.image_distance_m);
        atol = 1.0e-14,
        rtol = 1.0e-14,
    )
    frequency_complete = !isempty(schedules) && all(schedule -> schedule.loop_executed, schedules)
    wave_speeds_finite = all(isfinite, state.layer_wave_speeds_m_per_s) &&
        all(>(0.0), state.layer_wave_speeds_m_per_s)
    return radii_ordered && positive_materials && symmetric_geometry &&
        frequency_complete && wave_speeds_finite && state.derived_state_executed
end

function run_cable_constants_study(parsed::DeckParseResult; case_index::Integer = 1)
    assert_deck_valid!(parsed)
    cases = deck_cable_constants_cases(parsed)
    isempty(cases) && throw(ArgumentError("deck contains no accepted CABLE CONSTANTS case"))
    index = Int(case_index)
    1 <= index <= length(cases) ||
        throw(ArgumentError("cable constants case_index must be between 1 and $(length(cases))"))
    case = cases[index]
    state = _cable_study_pipe_sheath_state(case)
    geometry = _cable_study_geometry(case)
    schedules = _cable_study_frequency_schedules(case)
    frequency_states = _cable_study_frequency_states(case, state, geometry, schedules)
    checks = _cable_study_physical_checks(case, state, geometry, schedules) &&
        !isempty(frequency_states) && all(row -> row.physical_checks_passed, frequency_states)
    return CableConstantsStudyResult(
        parsed.source,
        case,
        state,
        geometry,
        schedules,
        frequency_states,
        checks,
        (),
    )
end

function build_semlyen_frequency_dependent_line(
    result::CableConstantsStudyResult,
    line_length_m::Real,
    from_nodes::AbstractVector{<:Integer},
    to_nodes::AbstractVector{<:Integer},
    timestep_s::Real,
    propagation_fits::AbstractVector{LineStepResponseExponentialFitResult},
    admittance_fits::AbstractVector{LineStepResponseExponentialFitResult},
    characteristic_admittance_s::AbstractVector{<:Real};
    schedule_index::Integer = 1,
    phasor_frequency_hz::Real,
    propagation_relative_tolerance::Real = 0.10,
    admittance_relative_tolerance::Real = 0.10,
)
    result.physical_checks_passed || throw(ArgumentError(
        "cable-constants study must pass its physical checks before EMT line construction",
    ))
    index = Int(schedule_index)
    1 <= index <= length(result.frequency_schedules) || throw(ArgumentError(
        "schedule_index must select a generated cable frequency schedule",
    ))
    first_state = 1 + sum(
        schedule.frequency_count for schedule in result.frequency_schedules[1:(index - 1)];
        init = 0,
    )
    state_count = result.frequency_schedules[index].frequency_count
    last_state = first_state + state_count - 1
    last_state <= length(result.frequency_states) || throw(ArgumentError(
        "generated cable frequency-state table is incomplete for the selected schedule",
    ))
    return nested_cable_semlyen_frequency_dependent_line_from_fit(
        @view(result.frequency_states[first_state:last_state]),
        line_length_m,
        from_nodes,
        to_nodes,
        timestep_s,
        propagation_fits,
        admittance_fits,
        characteristic_admittance_s;
        phasor_frequency_hz = phasor_frequency_hz,
        propagation_relative_tolerance = propagation_relative_tolerance,
        admittance_relative_tolerance = admittance_relative_tolerance,
    )
end

end
