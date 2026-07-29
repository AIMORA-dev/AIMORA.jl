
using LinearAlgebra

using ..Branches

import ..Branches:
    EMTElement,
    stamp!,
    stamp_admittance_entry!,
    stamp_conductance!,
    stamp_history_current!,
    update!

export BergeronLine,
       LineFrequencyPoint,
       LineModalTransform,
       LineModeUnwindState,
       LineModalSolution,
       LineModalSolutionScan,
       CoupledLumpedSequenceImpedance,
       CoupledLumpedPhasePiSection,
       DistributedTransposedLineConstants,
       DistributedTransposedLineModalBranchState,
       DistributedTransposedLineSteadyStatePiEquivalent,
       DistributedTransposedLineHistoryState,
       DistributedTransposedLineInitialHistoryPhaseTerms,
       DistributedTransposedLineInitialHistoryModalComponents,
       DistributedTransposedLineInitialHistoryPhasors,
       DistributedTransposedLineInitialHistorySeed,
       DistributedTransposedLineModalTimestepUpdate,
       DistributedTransposedLinePhaseCurrentInjection,
       DistributedTransposedLineCompanionAdmittance,
       TerminalSurgeImpedanceAdmittance,
       CableGeometryConductor,
       CableGeometryConstants,
       CableElectrostaticAdmittance,
       CablePhaseElectrostaticAdmittance,
       CableConductorSeriesImpedance,
       CablePhaseLineConstants,
       CableFrequencyScanSeriesImpedance,
       CableFrequencyScanGeneratedLineConstants,
       CablePiSectionReportState,
       CablePipeSheathDerivedState,
       NestedCableFrequencyState,
       NestedCableTransientLineState,
       CableFrequencyScanLoopSchedule,
       CableFrequencyScanLineConstants,
       CableFrequencyScanRuntimeState,
       CableFrequencyScanRecursiveConvolutionState,
       CableFrequencyScanModalResponseSamples,
       SemlyenMartiFrequencyScan,
       SemlyenLineExponentialConvolutionCoefficients,
       SemlyenLineHarmonicHistoryUpdate,
       FrequencyDependentLineModalResponse,
       FrequencyDependentLineModalState,
       FrequencyDependentLineRuntimeState,
       LineFrequencySampleFitResult,
       LineRecursiveConvolutionFitResult,
       LineStepResponseExponentialFitResult,
       FrequencyDependentLineSampleRuntimeState,
       FrequencyDependentLineModalSampleRuntimeState,
       FrequencyDependentLineRecursiveConvolutionState,
       line_delay_steps,
       line_surge_admittance,
       line_history_currents,
       line_terminal_currents,
       line_terminal_voltages,
       line_traveling_waves,
       line_characteristic_impedance,
       line_propagation_constant,
       line_modal_transform!,
       line_modal_transform,
       line_phase_transform!,
       line_phase_transform,
       line_modal_eigen_order,
       line_mode_unwind!,
       line_mode_unwind,
       line_modal_solution,
       line_modal_solution_scan,
       line_modal_solution_scan!,
       coupled_lumped_sequence_impedance,
       coupled_lumped_matrix_impedance,
       coupled_lumped_phase_pi_section,
       distributed_transposed_line_constants,
       distributed_transposed_line_modal_branch_state,
       distributed_modal_line_branch_state,
       distributed_transposed_line_steady_state_pi_equivalent,
       distributed_transposed_line_history_state,
       distributed_transposed_line_initial_history_phase_matrix,
       distributed_transposed_line_initial_history_phase_terms,
       distributed_transposed_line_initial_history_modal_components,
       distributed_transposed_line_initial_history_phasors,
       distributed_transposed_line_initial_history_seed,
       distributed_transposed_line_history_current_injection!,
       distributed_transposed_line_modal_timestep_update,
       distributed_transposed_line_modal_timestep_update!,
       distributed_transposed_line_phase_current_injection,
       distributed_transposed_line_phase_current_injection!,
       distributed_transposed_line_companion_admittance,
       distributed_transposed_line_companion_admittance_stamp!,
       terminal_surge_impedance_admittance,
       cable_skin_effect_internal_impedance,
       cable_homogeneous_earth_return_impedance,
       cable_underground_homogeneous_earth_return_impedance,
       cable_homogeneous_earth_return_impedance_matrix,
       cable_bounded_skin_effect_internal_impedance,
       cable_bounded_earth_return_impedance_matrix,
       cable_geometry_constants,
       cable_pipe_sheath_derived_state,
       cable_geometry_electrostatic_admittance,
       cable_phase_electrostatic_admittance,
       cable_phase_series_impedance,
       nested_cable_frequency_state,
       nested_cable_transient_line_state,
       nested_cable_semlyen_frequency_dependent_line_from_fit,
       cable_frequency_scan_loop_schedule,
       cable_frequency_scan_series_impedance,
       cable_phase_line_constants,
       cable_frequency_scan_line_constants,
       cable_frequency_scan_line_constants_from_geometry,
       cable_frequency_scan_runtime_update,
       cable_frequency_scan_runtime_update!,
       cable_frequency_scan_runtime_update_from_geometry,
       cable_frequency_scan_recursive_convolution_update,
       cable_frequency_scan_recursive_convolution_update!,
       cable_frequency_scan_recursive_convolution_update_from_fit,
       cable_frequency_scan_modal_response_samples,
       cable_frequency_scan_recursive_convolution_fit,
       cable_frequency_scan_recursive_convolution_update_from_step_response_fit,
       cable_frequency_scan_recursive_convolution_update_from_scan_fit,
       cable_frequency_scan_recursive_convolution_update_from_geometry,
       semlyen_marti_frequency_scan,
       semlyen_marti_frequency_dependent_line_samples,
       semlyen_line_exponential_convolution_coefficients,
       semlyen_line_harmonic_history_update,
       frequency_dependent_line_point,
       frequency_dependent_line_sample_fit,
       line_step_response_exponential_fit,
       frequency_dependent_line_modal_response,
       frequency_dependent_line_modal_response!,
       frequency_dependent_line_runtime_update,
       frequency_dependent_line_runtime_update!,
       frequency_dependent_line_runtime_update_from_fit,
       frequency_dependent_line_runtime_update_from_fit!,
       frequency_dependent_line_sample_runtime_update,
       frequency_dependent_line_sample_runtime_update!,
       frequency_dependent_line_modal_sample_runtime_update,
       frequency_dependent_line_modal_sample_runtime_update!,
       frequency_dependent_line_recursive_convolution_fit,
       frequency_dependent_line_recursive_convolution_fit_from_step_response,
       frequency_dependent_line_recursive_convolution_update,
       frequency_dependent_line_recursive_convolution_update!,
       frequency_dependent_line_recursive_convolution_update_from_fit

const DISTRIBUTED_LINE_ENDPOINT_NORTON_SHARE = 0.5

mutable struct BergeronLine <: EMTElement
    a::Int
    b::Int
    zc::Float64
    travel_time_s::Float64
    dt_s::Float64
    attenuation::Float64
    delay_steps::Int
    from_wave_history::Vector{Float64}
    to_wave_history::Vector{Float64}
    write_index::Int
    h_from::Float64
    h_to::Float64
    v_from::Float64
    v_to::Float64
    i_from::Float64
    i_to::Float64
end

function line_delay_steps(travel_time_s::Real, dt_s::Real)::Int
    travel = Float64(travel_time_s)
    dt = Float64(dt_s)
    isfinite(travel) && travel > 0.0 || throw(ArgumentError("travel_time_s must be finite and positive"))
    isfinite(dt) && dt > 0.0 || throw(ArgumentError("dt_s must be finite and positive"))
    steps = Int(round(travel / dt))
    steps > 0 || throw(ArgumentError("travel_time_s must be at least one timestep"))
    abs(steps * dt - travel) <= max(1.0e-15, 16.0 * eps(Float64) * max(abs(travel), dt)) ||
        throw(ArgumentError("travel_time_s must be an integer multiple of dt_s for this Bergeron line"))
    return steps
end

function BergeronLine(
    a::Int,
    b::Int,
    zc::Real,
    travel_time_s::Real,
    dt_s::Real;
    attenuation::Real=1.0,
    initial_from_wave::Real=0.0,
    initial_to_wave::Real=0.0,
)
    a >= 0 && b >= 0 || throw(ArgumentError("line nodes must be nonnegative"))
    impedance = Float64(zc)
    isfinite(impedance) && impedance > 0.0 ||
        throw(ArgumentError("zc must be finite and positive"))
    atten = Float64(attenuation)
    isfinite(atten) && 0.0 <= atten <= 1.0 ||
        throw(ArgumentError("attenuation must be finite and between 0 and 1"))
    steps = line_delay_steps(travel_time_s, dt_s)
    from_initial = Float64(initial_from_wave)
    to_initial = Float64(initial_to_wave)
    isfinite(from_initial) && isfinite(to_initial) ||
        throw(ArgumentError("initial traveling waves must be finite"))
    return BergeronLine(
        a,
        b,
        impedance,
        Float64(travel_time_s),
        Float64(dt_s),
        atten,
        steps,
        fill(from_initial, steps),
        fill(to_initial, steps),
        1,
        -atten * to_initial,
        -atten * from_initial,
        0.0,
        0.0,
        0.0,
        0.0,
    )
end

line_surge_admittance(line::BergeronLine)::Float64 = inv(line.zc)

line_terminal_voltages(line::BergeronLine) = (from = line.v_from, to = line.v_to)

line_terminal_currents(line::BergeronLine) = (from = line.i_from, to = line.i_to)

line_history_currents(line::BergeronLine) = (from = line.h_from, to = line.h_to)

function line_traveling_waves(line::BergeronLine)
    last_index = line.write_index == 1 ? line.delay_steps : line.write_index - 1
    return (
        from = line.from_wave_history[last_index],
        to = line.to_wave_history[last_index],
    )
end

function stamp!(
    y::AbstractMatrix{Float64},
    rhs::AbstractVector{Float64},
    line::BergeronLine,
    _t::Float64,
    _dt::Float64,
)
    g = line_surge_admittance(line)
    line.h_from = -line.attenuation * line.to_wave_history[line.write_index]
    line.h_to = -line.attenuation * line.from_wave_history[line.write_index]
    stamp_conductance!(y, line.a, 0, g)
    stamp_conductance!(y, line.b, 0, g)
    stamp_history_current!(rhs, line.a, 0, line.h_from)
    stamp_history_current!(rhs, line.b, 0, line.h_to)
    return nothing
end

function update!(line::BergeronLine, voltages::AbstractVector{Float64}, _dt::Float64)
    g = line_surge_admittance(line)
    v_from = line.a == 0 ? 0.0 : voltages[line.a]
    v_to = line.b == 0 ? 0.0 : voltages[line.b]
    i_from = g * v_from + line.h_from
    i_to = g * v_to + line.h_to
    line.v_from = v_from
    line.v_to = v_to
    line.i_from = i_from
    line.i_to = i_to
    line.from_wave_history[line.write_index] = g * v_from + i_from
    line.to_wave_history[line.write_index] = g * v_to + i_to
    line.write_index = line.write_index == line.delay_steps ? 1 : line.write_index + 1
    return nothing
end

struct CoupledLumpedSequenceImpedance
    name::Symbol
    phase_indices::Vector{Int}
    from_nodes::Vector{Symbol}
    to_nodes::Vector{Symbol}
    from_node_indices::Vector{Int}
    to_node_indices::Vector{Int}
    line_numbers::Vector{Int}
    zero_sequence_resistance::Float64
    positive_sequence_resistance::Float64
    zero_sequence_inductance::Float64
    positive_sequence_inductance::Float64
    phase_resistance_matrix::Matrix{Float64}
    phase_inductance_matrix::Matrix{Float64}
    phase_capacitance_matrix::Matrix{Float64}
    stored_resistance_values::Vector{Float64}
    stored_inductance_values::Vector{Float64}
    stored_capacitance_values::Vector{Float64}
    input_kind::Symbol
    phase_count::Int
end

struct CoupledLumpedPhasePiSection
    name::Symbol
    phase_indices::Vector{Int}
    from_nodes::Vector{Symbol}
    to_nodes::Vector{Symbol}
    from_node_indices::Vector{Int}
    to_node_indices::Vector{Int}
    line_numbers::Vector{Int}
    phase_resistance_matrix::Matrix{Float64}
    phase_inductance_matrix::Matrix{Float64}
    phase_capacitance_matrix::Matrix{Float64}
    stored_resistance_values::Vector{Float64}
    stored_inductance_values::Vector{Float64}
    stored_capacitance_values::Vector{Float64}
    phase_count::Int
end

struct DistributedTransposedLineConstants
    name::Symbol
    phase_indices::Vector{Int}
    from_nodes::Vector{Symbol}
    to_nodes::Vector{Symbol}
    from_node_indices::Vector{Int}
    to_node_indices::Vector{Int}
    line_numbers::Vector{Int}
    line_length::Float64
    x_frequency_hz::Float64
    c_frequency_hz::Float64
    sequence_resistance_per_length::Vector{Float64}
    sequence_inductance_input_values::Vector{Float64}
    sequence_capacitance_input_values::Vector{Float64}
    sequence_inductance_h_per_length::Vector{Float64}
    sequence_capacitance_f_per_length::Vector{Float64}
    sequence_characteristic_impedances::Vector{Float64}
    sequence_total_resistances::Vector{Float64}
    sequence_propagation_times_s::Vector{Float64}
    ci_values::Vector{Float64}
    ck_values::Vector{Float64}
    cik_values::Vector{Float64}
    phase_resistance_matrix::Matrix{Float64}
    phase_inductance_matrix::Matrix{Float64}
    phase_capacitance_matrix::Matrix{Float64}
    stored_resistance_values::Vector{Float64}
    stored_inductance_values::Vector{Float64}
    stored_capacitance_values::Vector{Float64}
    phase_count::Int
end

struct DistributedTransposedLineModalBranchState
    name::Symbol
    phase_indices::Vector{Int}
    modal_sequence_indices::Vector{Int}
    continuation_copy_source_indices::Vector{Int}
    from_nodes::Vector{Symbol}
    to_nodes::Vector{Symbol}
    from_node_indices::Vector{Int}
    to_node_indices::Vector{Int}
    line_numbers::Vector{Int}
    group_length::Int
    modal_signed_characteristic_impedances::Vector{Float64}
    modal_total_resistances::Vector{Float64}
    modal_propagation_times_s::Vector{Float64}
    phase_to_modal_transform_matrix::Matrix{Float64}
    modal_to_phase_transform_matrix::Matrix{Float64}
    modal_admittance_denominator::Float64
    phase_count::Int
end

struct DistributedTransposedLineSteadyStatePiEquivalent
    name::Symbol
    phase_indices::Vector{Int}
    line_numbers::Vector{Int}
    modal_sequence_indices::Vector{Int}
    steady_state_frequency_hz::Float64
    angular_frequency_rad_s::Float64
    x_frequency_hz::Float64
    c_frequency_hz::Float64
    storage_start_index::Int
    storage_end_index::Int
    storage_row_indices::Vector{Int}
    phase_row_indices::Vector{Int}
    phase_column_indices::Vector{Int}
    phase_series_resistance_values::Vector{Float64}
    phase_series_reactance_values::Vector{Float64}
    phase_shunt_conductance_values::Vector{Float64}
    phase_shunt_susceptance_values::Vector{Float64}
    phase_series_impedance_matrix::Matrix{ComplexF64}
    phase_shunt_admittance_matrix::Matrix{ComplexF64}
    phase_count::Int
end

struct DistributedTransposedLineHistoryState
    name::Symbol
    phase_indices::Vector{Int}
    line_numbers::Vector{Int}
    modal_sequence_indices::Vector{Int}
    timestep_s::Float64
    steady_state_frequency_hz::Float64
    angular_step_rad::Float64
    history_storage_start_index::Int
    next_history_storage_index::Int
    storage_start_indices::Vector{Int}
    storage_end_indices::Vector{Int}
    storage_lengths::Vector{Int}
    history_sample_counts::Vector{Int}
    modal_history_interpolation_factors::Vector{Float64}
    history_read_indices::Vector{Int}
    outgoing_wave_from_real_values::Vector{Float64}
    outgoing_wave_from_imag_values::Vector{Float64}
    outgoing_wave_to_real_values::Vector{Float64}
    outgoing_wave_to_imag_values::Vector{Float64}
    modal_history_from_values::Vector{Vector{Float64}}
    modal_history_to_values::Vector{Vector{Float64}}
    packed_history_from_values::Vector{Float64}
    packed_history_to_values::Vector{Float64}
    initialized_from_steady_state::Bool
    phase_count::Int
end

struct DistributedTransposedLineInitialHistoryPhaseMatrix
    name::Symbol
    phase_indices::Vector{Int}
    line_numbers::Vector{Int}
    modal_sequence_indices::Vector{Int}
    modal_initial_history_impedance_values::Vector{Float64}
    modal_history_scale_values::Vector{Float64}
    phase_matrix_values::Matrix{Float64}
    phase_matrix_upper_values::Vector{Float64}
    phase_count::Int
end

struct DistributedTransposedLineInitialHistoryPhaseTerms
    name::Symbol
    phase_indices::Vector{Int}
    reduced_phase_matrix_upper_values::Vector{Float64}
    matrix_terminal_from_real_values::Vector{Float64}
    matrix_terminal_from_imag_values::Vector{Float64}
    matrix_terminal_to_real_values::Vector{Float64}
    matrix_terminal_to_imag_values::Vector{Float64}
    initial_terminal_from_real_values::Vector{Float64}
    initial_terminal_from_imag_values::Vector{Float64}
    initial_terminal_to_real_values::Vector{Float64}
    initial_terminal_to_imag_values::Vector{Float64}
    matrix_product_from_real_values::Vector{Float64}
    matrix_product_from_imag_values::Vector{Float64}
    matrix_product_to_real_values::Vector{Float64}
    matrix_product_to_imag_values::Vector{Float64}
    terminal_from_real_values::Vector{Float64}
    terminal_from_imag_values::Vector{Float64}
    terminal_to_real_values::Vector{Float64}
    terminal_to_imag_values::Vector{Float64}
    phase_count::Int
end

struct DistributedTransposedLineInitialHistoryModalComponents
    name::Symbol
    modal_indices::Vector{Int}
    term_modal_indices::Vector{Int}
    term_weight_values::Vector{Float64}
    term_from_real_values::Vector{Float64}
    term_from_imag_values::Vector{Float64}
    term_to_real_values::Vector{Float64}
    term_to_imag_values::Vector{Float64}
    term_from_real_products::Vector{Float64}
    term_from_imag_products::Vector{Float64}
    term_to_real_products::Vector{Float64}
    term_to_imag_products::Vector{Float64}
    modal_component_from_real_values::Vector{Float64}
    modal_component_from_imag_values::Vector{Float64}
    modal_component_to_real_values::Vector{Float64}
    modal_component_to_imag_values::Vector{Float64}
    phase_count::Int
end

struct DistributedTransposedLineInitialHistoryPhasors
    name::Symbol
    modal_indices::Vector{Int}
    modal_scale_values::Vector{Float64}
    modal_component_from_real_values::Vector{Float64}
    modal_component_from_imag_values::Vector{Float64}
    modal_component_to_real_values::Vector{Float64}
    modal_component_to_imag_values::Vector{Float64}
    modal_history_from_amplitudes::Vector{Float64}
    modal_history_from_phases::Vector{Float64}
    modal_history_to_amplitudes::Vector{Float64}
    modal_history_to_phases::Vector{Float64}
    phase_count::Int
end

struct DistributedTransposedLineInitialHistorySeed
    name::Symbol
    modal_indices::Vector{Int}
    history_start_indices::Vector{Int}
    modal_delay_counts::Vector{Int}
    sinusoidal_start_indices::Vector{Int}
    sinusoidal_end_indices::Vector{Int}
    interpolation_indices::Vector{Int}
    modal_history_interpolation_factors::Vector{Float64}
    modal_history_damping_values::Vector{Float64}
    modal_history_from_values::Vector{Vector{Float64}}
    modal_history_to_values::Vector{Vector{Float64}}
    modal_history_interpolation_from_values::Vector{Float64}
    modal_history_interpolation_to_values::Vector{Float64}
    packed_history_start_index::Int
    packed_history_end_index::Int
    packed_history_from_values::Vector{Float64}
    packed_history_to_values::Vector{Float64}
    phase_count::Int
end

struct DistributedTransposedLineModalTimestepUpdate
    name::Symbol
    phase_indices::Vector{Int}
    line_numbers::Vector{Int}
    modal_sequence_indices::Vector{Int}
    modal_terminal_voltage_from_values::Vector{Float64}
    modal_terminal_voltage_to_values::Vector{Float64}
    modal_current_from_values::Vector{Float64}
    modal_current_to_values::Vector{Float64}
    phase_current_injection_from_values::Vector{Float64}
    phase_current_injection_to_values::Vector{Float64}
    history_read_indices_before::Vector{Int}
    history_read_indices_after::Vector{Int}
    history_write_indices::Vector{Int}
    history_next_write_indices::Vector{Int}
    modal_history_interpolation_factors::Vector{Float64}
    phase_count::Int
end

struct DistributedTransposedLinePhaseCurrentInjection
    name::Symbol
    from_node_indices::Vector{Int}
    to_node_indices::Vector{Int}
    phase_current_injection_from_values::Vector{Float64}
    phase_current_injection_to_values::Vector{Float64}
    rhs_before_values::Vector{Float64}
    rhs_after_values::Vector{Float64}
    rhs_update_count::Int
    phase_count::Int
end

struct DistributedTransposedLineCompanionAdmittance <: EMTElement
    name::Symbol
    from_node_indices::Vector{Int}
    to_node_indices::Vector{Int}
    line_numbers::Vector{Int}
    modal_sequence_indices::Vector{Int}
    modal_companion_admittance_values::Vector{Float64}
    phase_companion_admittance_matrix::Matrix{Float64}
    phase_count::Int
end

struct TerminalSurgeImpedanceAdmittance <: EMTElement
    name::Symbol
    node_names::Vector{Symbol}
    node_indices::Vector{Int}
    line_numbers::Vector{Int}
    phase_impedance_matrix::Matrix{Float64}
    phase_admittance_matrix::Matrix{Float64}
    stamped_admittance_matrix::Matrix{Float64}
    phase_count::Int
end

_finite_impedance_value(name::AbstractString, value::Real) = begin
    numeric = Float64(value)
    isfinite(numeric) || throw(ArgumentError("$name must be finite"))
    numeric
end

function _positive_finite_impedance_value(name::AbstractString, value::Real)
    numeric = _finite_impedance_value(name, value)
    numeric > 0.0 || throw(ArgumentError("$name must be positive"))
    return numeric
end

function _sequence_phase_matrix(zero_sequence::Float64, positive_sequence::Float64)
    diagonal = (zero_sequence + 2.0 * positive_sequence) / 3.0
    mutual = (zero_sequence - positive_sequence) / 3.0
    return [
        diagonal mutual mutual
        mutual diagonal mutual
        mutual mutual diagonal
    ]
end

function _stored_symmetric_phase_values(matrix::AbstractMatrix{Float64})
    rows, columns = size(matrix)
    rows == columns || throw(ArgumentError("phase matrix must be square"))
    return Float64[matrix[row, column] for row in 1:rows for column in 1:row]
end

function _phase_triplet(values::AbstractVector{<:Integer})
    phases = Int.(values)
    phases == [1, 2, 3] ||
        throw(ArgumentError("coupled lumped sequence impedance requires phase indices [1, 2, 3]"))
    return phases
end

function _phase_index_vector(values::AbstractVector{<:Integer})
    phases = Int.(values)
    !isempty(phases) || throw(ArgumentError("phase indices must not be empty"))
    phases == collect(1:length(phases)) ||
        throw(ArgumentError("phase indices must be contiguous from 1"))
    return phases
end

function _symbol_triplet(name::AbstractString, values::AbstractVector{Symbol})
    length(values) == 3 || throw(ArgumentError("$name must have three phase entries"))
    return Symbol.(values)
end

function _symbol_phase_vector(name::AbstractString, values::AbstractVector{Symbol},
                              phase_count::Int)
    length(values) == phase_count ||
        throw(ArgumentError("$name must have $phase_count phase entries"))
    return Symbol.(values)
end

function _index_triplet(name::AbstractString, values::AbstractVector{<:Integer})
    length(values) == 3 || throw(ArgumentError("$name must have three phase entries"))
    indices = Int.(values)
    all(>=(0), indices) || throw(ArgumentError("$name must be nonnegative"))
    return indices
end

function _index_phase_vector(name::AbstractString, values::AbstractVector{<:Integer},
                             phase_count::Int)
    length(values) == phase_count ||
        throw(ArgumentError("$name must have $phase_count phase entries"))
    indices = Int.(values)
    all(>=(0), indices) || throw(ArgumentError("$name must be nonnegative"))
    return indices
end

function _line_number_triplet(values::AbstractVector{<:Integer})
    isempty(values) && return Int[]
    length(values) == 3 || throw(ArgumentError("line_numbers must be empty or have three phase entries"))
    return Int.(values)
end

function _line_number_vector(values::AbstractVector{<:Integer}, phase_count::Int)
    isempty(values) && return Int[]
    length(values) == phase_count ||
        throw(ArgumentError("line_numbers must be empty or have $phase_count phase entries"))
    return Int.(values)
end

function _finite_square_phase_matrix(matrix::AbstractMatrix{<:Real},
                                     phase_count::Int,
                                     field::AbstractString)
    size(matrix) == (phase_count, phase_count) ||
        throw(ArgumentError("$field matrix must be $phase_count by $phase_count"))
    values = Matrix{Float64}(matrix)
    all(isfinite, values) || throw(ArgumentError("$field matrix values must be finite"))
    return values
end

function coupled_lumped_sequence_impedance(
    name::Symbol,
    phase_indices::AbstractVector{<:Integer},
    from_nodes::AbstractVector{Symbol},
    to_nodes::AbstractVector{Symbol},
    from_node_indices::AbstractVector{<:Integer},
    to_node_indices::AbstractVector{<:Integer},
    zero_sequence_resistance::Real,
    positive_sequence_resistance::Real,
    zero_sequence_inductance::Real,
    positive_sequence_inductance::Real;
    line_numbers::AbstractVector{<:Integer}=Int[],
)
    phases = _phase_triplet(phase_indices)
    from_symbols = _symbol_triplet("from_nodes", from_nodes)
    to_symbols = _symbol_triplet("to_nodes", to_nodes)
    from_indices = _index_triplet("from_node_indices", from_node_indices)
    to_indices = _index_triplet("to_node_indices", to_node_indices)
    source_lines = _line_number_triplet(line_numbers)
    r0 = _finite_impedance_value("zero_sequence_resistance", zero_sequence_resistance)
    r1 = _finite_impedance_value("positive_sequence_resistance", positive_sequence_resistance)
    l0 = _finite_impedance_value("zero_sequence_inductance", zero_sequence_inductance)
    l1 = _finite_impedance_value("positive_sequence_inductance", positive_sequence_inductance)
    resistance_matrix = _sequence_phase_matrix(r0, r1)
    inductance_matrix = _sequence_phase_matrix(l0, l1)
    capacitance_matrix = zeros(Float64, 3, 3)
    return CoupledLumpedSequenceImpedance(
        name,
        phases,
        from_symbols,
        to_symbols,
        from_indices,
        to_indices,
        source_lines,
        r0,
        r1,
        l0,
        l1,
        resistance_matrix,
        inductance_matrix,
        capacitance_matrix,
        _stored_symmetric_phase_values(resistance_matrix),
        _stored_symmetric_phase_values(inductance_matrix),
        zeros(Float64, 6),
        :sequence_shorthand,
        3,
    )
end

function coupled_lumped_matrix_impedance(
    name::Symbol,
    phase_indices::AbstractVector{<:Integer},
    from_nodes::AbstractVector{Symbol},
    to_nodes::AbstractVector{Symbol},
    from_node_indices::AbstractVector{<:Integer},
    to_node_indices::AbstractVector{<:Integer},
    resistance_matrix::AbstractMatrix{<:Real},
    inductance_matrix::AbstractMatrix{<:Real};
    line_numbers::AbstractVector{<:Integer}=Int[],
)
    phase_count = length(phase_indices)
    phases = _phase_index_vector(phase_indices)
    from_symbols = _symbol_phase_vector("from_nodes", from_nodes, phase_count)
    to_symbols = _symbol_phase_vector("to_nodes", to_nodes, phase_count)
    from_indices = _index_phase_vector("from_node_indices", from_node_indices, phase_count)
    to_indices = _index_phase_vector("to_node_indices", to_node_indices, phase_count)
    source_lines = _line_number_vector(line_numbers, phase_count)
    resistance =
        _finite_square_phase_matrix(resistance_matrix, phase_count, "phase_resistance")
    inductance =
        _finite_square_phase_matrix(inductance_matrix, phase_count, "phase_inductance")
    tolerance = 64.0 * eps(Float64) * max(
        maximum(abs, resistance; init=0.0),
        maximum(abs, inductance; init=0.0),
        1.0,
    )
    maximum(abs, resistance - transpose(resistance); init=0.0) <= tolerance ||
        throw(ArgumentError("phase resistance matrix must be symmetric"))
    maximum(abs, inductance - transpose(inductance); init=0.0) <= tolerance ||
        throw(ArgumentError("phase inductance matrix must be symmetric"))
    resistance = 0.5 .* (resistance .+ transpose(resistance))
    inductance = 0.5 .* (inductance .+ transpose(inductance))
    zero_resistance = sum(resistance) / phase_count
    positive_resistance = phase_count == 1 ? zero_resistance :
        (tr(resistance) - zero_resistance) / (phase_count - 1)
    zero_inductance = sum(inductance) / phase_count
    positive_inductance = phase_count == 1 ? zero_inductance :
        (tr(inductance) - zero_inductance) / (phase_count - 1)
    stored_count = phase_count * (phase_count + 1) ÷ 2
    return CoupledLumpedSequenceImpedance(
        name,
        phases,
        from_symbols,
        to_symbols,
        from_indices,
        to_indices,
        source_lines,
        zero_resistance,
        positive_resistance,
        zero_inductance,
        positive_inductance,
        resistance,
        inductance,
        zeros(Float64, phase_count, phase_count),
        _stored_symmetric_phase_values(resistance),
        _stored_symmetric_phase_values(inductance),
        zeros(Float64, stored_count),
        :triangular_matrix,
        phase_count,
    )
end

function coupled_lumped_phase_pi_section(
    name::Symbol,
    phase_indices::AbstractVector{<:Integer},
    from_nodes::AbstractVector{Symbol},
    to_nodes::AbstractVector{Symbol},
    from_node_indices::AbstractVector{<:Integer},
    to_node_indices::AbstractVector{<:Integer},
    resistance_matrix::AbstractMatrix{<:Real},
    inductance_matrix::AbstractMatrix{<:Real},
    capacitance_matrix::AbstractMatrix{<:Real};
    line_numbers::AbstractVector{<:Integer}=Int[],
)
    phases = _phase_index_vector(phase_indices)
    phase_count = length(phases)
    from_symbols = _symbol_phase_vector("from_nodes", from_nodes, phase_count)
    to_symbols = _symbol_phase_vector("to_nodes", to_nodes, phase_count)
    from_indices = _index_phase_vector("from_node_indices", from_node_indices, phase_count)
    to_indices = _index_phase_vector("to_node_indices", to_node_indices, phase_count)
    source_lines = _line_number_vector(line_numbers, phase_count)
    resistance = _finite_square_phase_matrix(
        resistance_matrix,
        phase_count,
        "phase resistance",
    )
    inductance = _finite_square_phase_matrix(
        inductance_matrix,
        phase_count,
        "phase inductance",
    )
    capacitance = _finite_square_phase_matrix(
        capacitance_matrix,
        phase_count,
        "phase capacitance",
    )
    return CoupledLumpedPhasePiSection(
        name,
        phases,
        from_symbols,
        to_symbols,
        from_indices,
        to_indices,
        source_lines,
        resistance,
        inductance,
        capacitance,
        _stored_symmetric_phase_values(resistance),
        _stored_symmetric_phase_values(inductance),
        _stored_symmetric_phase_values(capacitance),
        phase_count,
    )
end

function terminal_surge_impedance_admittance(
    name::Symbol,
    node_names::AbstractVector{Symbol},
    node_indices::AbstractVector{<:Integer},
    phase_impedance_matrix::AbstractMatrix{<:Real};
    line_numbers::AbstractVector{<:Integer}=Int[],
    preexisting_admittance_matrix::Union{Nothing,AbstractMatrix{<:Real}}=nothing,
)
    phase_count = length(node_indices)
    phase_count > 0 || throw(ArgumentError("terminal surge impedance requires at least one terminal"))
    nodes = _symbol_phase_vector("node_names", node_names, phase_count)
    indices = _index_phase_vector("node_indices", node_indices, phase_count)
    lines = _line_number_vector(line_numbers, phase_count)
    impedance = _finite_square_phase_matrix(
        phase_impedance_matrix,
        phase_count,
        "terminal surge impedance",
    )
    symmetry_error = maximum(abs.(impedance - transpose(impedance)); init = 0.0)
    symmetry_error <= 1.0e-9 ||
        throw(ArgumentError("terminal surge impedance matrix must be symmetric"))
    admittance = Matrix{Float64}(inv(impedance))
    all(isfinite, admittance) ||
        throw(ArgumentError("terminal surge admittance matrix must be finite"))
    preexisting =
        preexisting_admittance_matrix === nothing ?
        zeros(Float64, phase_count, phase_count) :
        _finite_square_phase_matrix(
            preexisting_admittance_matrix,
            phase_count,
            "preexisting terminal admittance",
        )
    return TerminalSurgeImpedanceAdmittance(
        name,
        nodes,
        indices,
        lines,
        impedance,
        admittance,
        admittance - preexisting,
        phase_count,
    )
end

function _distributed_line_sequence_values(
    resistance_per_length::Float64,
    inductance_input::Float64,
    capacitance_input::Float64,
    line_length::Float64,
    x_frequency_hz::Float64,
    c_frequency_hz::Float64,
)
    inductance_h_per_length = inductance_input * 1.0e-3
    capacitance_f_per_length = capacitance_input / 1.0e6
    if x_frequency_hz > 0.0
        inductance_h_per_length *= 1000.0 / (2.0 * pi * x_frequency_hz)
    end
    if c_frequency_hz > 0.0
        capacitance_f_per_length /= 2.0 * pi * c_frequency_hz
    end
    inductance_h_per_length > 0.0 ||
        throw(ArgumentError("distributed line inductance must convert to a positive value"))
    capacitance_f_per_length > 0.0 ||
        throw(ArgumentError("distributed line capacitance must convert to a positive value"))
    characteristic_impedance = sqrt(inductance_h_per_length / capacitance_f_per_length)
    propagation_time_s = line_length * characteristic_impedance * capacitance_f_per_length
    total_resistance = resistance_per_length * line_length
    ci_value = resistance_per_length == 0.0 ? -characteristic_impedance : characteristic_impedance
    return (
        inductance_h_per_length = inductance_h_per_length,
        capacitance_f_per_length = capacitance_f_per_length,
        characteristic_impedance = characteristic_impedance,
        propagation_time_s = propagation_time_s,
        total_resistance = total_resistance,
        ci = ci_value,
        ck = total_resistance,
        cik = propagation_time_s,
    )
end

function distributed_transposed_line_constants(
    name::Symbol,
    phase_indices::AbstractVector{<:Integer},
    from_nodes::AbstractVector{Symbol},
    to_nodes::AbstractVector{Symbol},
    from_node_indices::AbstractVector{<:Integer},
    to_node_indices::AbstractVector{<:Integer},
    zero_sequence_resistance_per_length::Real,
    positive_sequence_resistance_per_length::Real,
    zero_sequence_inductance_input::Real,
    positive_sequence_inductance_input::Real,
    zero_sequence_capacitance_input::Real,
    positive_sequence_capacitance_input::Real,
    line_length::Real;
    x_frequency_hz::Real=0.0,
    c_frequency_hz::Real=0.0,
    line_numbers::AbstractVector{<:Integer}=Int[],
    third_sequence_resistance_per_length::Union{Nothing,Real}=nothing,
    third_sequence_inductance_input::Union{Nothing,Real}=nothing,
    third_sequence_capacitance_input::Union{Nothing,Real}=nothing,
)
    phases = _phase_triplet(phase_indices)
    from_symbols = _symbol_triplet("from_nodes", from_nodes)
    to_symbols = _symbol_triplet("to_nodes", to_nodes)
    from_indices = _index_triplet("from_node_indices", from_node_indices)
    to_indices = _index_triplet("to_node_indices", to_node_indices)
    source_lines = _line_number_triplet(line_numbers)
    length_value = _positive_finite_impedance_value("line_length", line_length)
    x_frequency = _finite_impedance_value("x_frequency_hz", x_frequency_hz)
    c_frequency = _finite_impedance_value("c_frequency_hz", c_frequency_hz)
    x_frequency >= 0.0 || throw(ArgumentError("x_frequency_hz must be nonnegative"))
    c_frequency >= 0.0 || throw(ArgumentError("c_frequency_hz must be nonnegative"))
    sequence_resistances = Float64[
        _finite_impedance_value("zero_sequence_resistance_per_length",
                                zero_sequence_resistance_per_length),
        _finite_impedance_value("positive_sequence_resistance_per_length",
                                positive_sequence_resistance_per_length),
    ]
    sequence_inductance_inputs = Float64[
        _positive_finite_impedance_value("zero_sequence_inductance_input",
                                         zero_sequence_inductance_input),
        _positive_finite_impedance_value("positive_sequence_inductance_input",
                                         positive_sequence_inductance_input),
    ]
    sequence_capacitance_inputs = Float64[
        _positive_finite_impedance_value("zero_sequence_capacitance_input",
                                         zero_sequence_capacitance_input),
        _positive_finite_impedance_value("positive_sequence_capacitance_input",
                                         positive_sequence_capacitance_input),
    ]
    explicit_third_sequence = any(
        value -> value !== nothing,
        (
            third_sequence_resistance_per_length,
            third_sequence_inductance_input,
            third_sequence_capacitance_input,
        ),
    )
    if explicit_third_sequence
        all(
            value -> value !== nothing,
            (
                third_sequence_resistance_per_length,
                third_sequence_inductance_input,
                third_sequence_capacitance_input,
            ),
        ) || throw(ArgumentError("explicit third distributed-line row requires resistance, inductance, and capacitance"))
        push!(
            sequence_resistances,
            _finite_impedance_value(
                "third_sequence_resistance_per_length",
                third_sequence_resistance_per_length,
            ),
        )
        push!(
            sequence_inductance_inputs,
            _positive_finite_impedance_value(
                "third_sequence_inductance_input",
                third_sequence_inductance_input,
            ),
        )
        push!(
            sequence_capacitance_inputs,
            _positive_finite_impedance_value(
                "third_sequence_capacitance_input",
                third_sequence_capacitance_input,
            ),
        )
    end
    sequence_values = [
        _distributed_line_sequence_values(
            sequence_resistances[index],
            sequence_inductance_inputs[index],
            sequence_capacitance_inputs[index],
            length_value,
            x_frequency,
            c_frequency,
        )
        for index in eachindex(sequence_resistances)
    ]
    inductances = Float64[value.inductance_h_per_length for value in sequence_values]
    capacitances = Float64[value.capacitance_f_per_length for value in sequence_values]
    characteristic_impedances =
        Float64[value.characteristic_impedance for value in sequence_values]
    propagation_times = Float64[value.propagation_time_s for value in sequence_values]
    total_resistances = Float64[value.total_resistance for value in sequence_values]
    ci_values = Float64[value.ci for value in sequence_values]
    ck_values = Float64[value.ck for value in sequence_values]
    cik_values = Float64[value.cik for value in sequence_values]
    resistance_matrix =
        _sequence_phase_matrix(sequence_resistances[1], sequence_resistances[2])
    inductance_matrix = _sequence_phase_matrix(inductances[1], inductances[2])
    capacitance_matrix = _sequence_phase_matrix(capacitances[1], capacitances[2])
    return DistributedTransposedLineConstants(
        name,
        phases,
        from_symbols,
        to_symbols,
        from_indices,
        to_indices,
        source_lines,
        length_value,
        x_frequency,
        c_frequency,
        sequence_resistances,
        sequence_inductance_inputs,
        sequence_capacitance_inputs,
        inductances,
        capacitances,
        characteristic_impedances,
        total_resistances,
        propagation_times,
        ci_values,
        ck_values,
        cik_values,
        resistance_matrix,
        inductance_matrix,
        capacitance_matrix,
        _stored_symmetric_phase_values(resistance_matrix),
        _stored_symmetric_phase_values(inductance_matrix),
        _stored_symmetric_phase_values(capacitance_matrix),
        3,
    )
end

function _finite_modal_transform_matrix(
    matrix::AbstractMatrix{<:Real},
    phase_count::Int,
    field::AbstractString,
)
    size(matrix) == (phase_count, phase_count) ||
        throw(ArgumentError("$field matrix must be $phase_count by $phase_count"))
    values = Matrix{Float64}(matrix)
    all(isfinite, values) || throw(ArgumentError("$field matrix entries must be finite"))
    return values
end

function distributed_modal_line_branch_state(
    name::Symbol,
    phase_indices::AbstractVector{<:Integer},
    from_nodes::AbstractVector{Symbol},
    to_nodes::AbstractVector{Symbol},
    from_node_indices::AbstractVector{<:Integer},
    to_node_indices::AbstractVector{<:Integer},
    modal_signed_characteristic_impedances::AbstractVector{<:Real},
    modal_total_resistances::AbstractVector{<:Real},
    modal_propagation_times_s::AbstractVector{<:Real},
    modal_to_phase_transform_matrix::AbstractMatrix{<:Real};
    phase_to_modal_transform_matrix::Union{Nothing,AbstractMatrix{<:Real}}=nothing,
    line_numbers::AbstractVector{<:Integer}=Int[],
    modal_sequence_indices::Union{Nothing,AbstractVector{<:Integer}}=nothing,
    continuation_copy_source_indices::Union{Nothing,AbstractVector{<:Integer}}=nothing,
    modal_admittance_denominator::Real=1.0,
)
    phase_count = length(phase_indices)
    phase_count > 0 || throw(ArgumentError("modal line state requires at least one phase"))
    phases = _copy_integer_modal_values("phase_indices", phase_indices, phase_count)
    phases == collect(1:phase_count) ||
        throw(ArgumentError("phase_indices must be contiguous from 1"))
    from_symbols = _symbol_phase_vector("from_nodes", from_nodes, phase_count)
    to_symbols = _symbol_phase_vector("to_nodes", to_nodes, phase_count)
    from_indices = _index_phase_vector("from_node_indices", from_node_indices, phase_count)
    to_indices = _index_phase_vector("to_node_indices", to_node_indices, phase_count)
    source_lines = _line_number_vector(line_numbers, phase_count)
    signed_impedances = _copy_real_modal_values(
        "modal_signed_characteristic_impedances",
        modal_signed_characteristic_impedances,
        phase_count,
    )
    all(!=(0.0), signed_impedances) ||
        throw(ArgumentError("modal characteristic impedances must be nonzero"))
    total_resistances = _copy_real_modal_values(
        "modal_total_resistances",
        modal_total_resistances,
        phase_count,
    )
    propagation_times = _copy_real_modal_values(
        "modal_propagation_times_s",
        modal_propagation_times_s,
        phase_count,
    )
    all(>(0.0), propagation_times) ||
        throw(ArgumentError("modal propagation times must be positive"))
    modal_to_phase = _finite_modal_transform_matrix(
        modal_to_phase_transform_matrix,
        phase_count,
        "modal_to_phase_transform_matrix",
    )
    phase_to_modal =
        phase_to_modal_transform_matrix === nothing ?
        Matrix{Float64}(inv(modal_to_phase)) :
        _finite_modal_transform_matrix(
            phase_to_modal_transform_matrix,
            phase_count,
            "phase_to_modal_transform_matrix",
        )
    inverse_error = maximum(abs.(phase_to_modal * modal_to_phase - I(phase_count)))
    inverse_error <= 1.0e-9 ||
        throw(ArgumentError("modal transform matrices must be inverses"))
    admittance_denominator =
        _positive_finite_impedance_value("modal_admittance_denominator",
                                         modal_admittance_denominator)
    sequences =
        modal_sequence_indices === nothing ?
        collect(1:phase_count) :
        _copy_integer_modal_values("modal_sequence_indices",
                                   modal_sequence_indices, phase_count)
    copies =
        continuation_copy_source_indices === nothing ?
        zeros(Int, phase_count) :
        Int.(continuation_copy_source_indices)
    length(copies) == phase_count ||
        throw(ArgumentError("continuation_copy_source_indices count must match phase count"))
    all(>=(0), copies) ||
        throw(ArgumentError("continuation_copy_source_indices entries must be nonnegative"))
    return DistributedTransposedLineModalBranchState(
        name,
        phases,
        sequences,
        copies,
        from_symbols,
        to_symbols,
        from_indices,
        to_indices,
        source_lines,
        phase_count,
        signed_impedances,
        total_resistances,
        propagation_times,
        phase_to_modal,
        modal_to_phase,
        admittance_denominator,
        phase_count,
    )
end

function distributed_transposed_line_modal_branch_state(
    constants::DistributedTransposedLineConstants;
    name::Symbol=Symbol(string(constants.name), "_modal_branch_state"),
)
    constants.phase_indices == [1, 2, 3] ||
        throw(ArgumentError("distributed transposed line modal state requires phase indices [1, 2, 3]"))
    length(constants.sequence_characteristic_impedances) in (2, 3) ||
        throw(ArgumentError("distributed transposed line modal state requires two or three sequence constants"))
    modal_sequence_indices =
        length(constants.sequence_characteristic_impedances) == 3 ? Int[1, 2, 3] : Int[1, 2, 2]
    continuation_copy_source_indices =
        length(constants.sequence_characteristic_impedances) == 3 ? Int[0, 0, 0] : Int[0, 0, 2]
    modal_to_phase_transform = _transposed_line_modal_transform()
    phase_to_modal_transform = constants.phase_count .* transpose(modal_to_phase_transform)
    return DistributedTransposedLineModalBranchState(
        name,
        copy(constants.phase_indices),
        modal_sequence_indices,
        continuation_copy_source_indices,
        copy(constants.from_nodes),
        copy(constants.to_nodes),
        copy(constants.from_node_indices),
        copy(constants.to_node_indices),
        copy(constants.line_numbers),
        length(modal_sequence_indices),
        Float64[constants.ci_values[index] for index in modal_sequence_indices],
        Float64[constants.ck_values[index] for index in modal_sequence_indices],
        Float64[constants.cik_values[index] for index in modal_sequence_indices],
        Matrix{Float64}(phase_to_modal_transform),
        Matrix{Float64}(modal_to_phase_transform),
        Float64(constants.phase_count),
        constants.phase_count,
    )
end

function _modal_complex_values(
    name::AbstractString,
    values::AbstractVector{<:Number},
    expected_count::Int,
)
    length(values) == expected_count ||
        throw(ArgumentError("$name must have $expected_count modal entries"))
    converted = ComplexF64.(values)
    all(value -> isfinite(real(value)) && isfinite(imag(value)), converted) ||
        throw(ArgumentError("$name entries must be finite"))
    return converted
end

function distributed_transposed_line_history_state(
    state::DistributedTransposedLineModalBranchState;
    timestep_s::Real,
    steady_state_frequency_hz::Real,
    history_storage_start_index::Integer=1,
    outgoing_wave_from::AbstractVector{<:Number}=zeros(ComplexF64, state.group_length),
    outgoing_wave_to::AbstractVector{<:Number}=zeros(ComplexF64, state.group_length),
    initialized_from_steady_state::Bool=true,
    name::Symbol=Symbol(string(state.name), "_history_state"),
)
    state.group_length > 0 ||
        throw(ArgumentError("distributed transposed line history state requires at least one mode"))
    timestep = _positive_finite_impedance_value("timestep_s", timestep_s)
    frequency = _positive_finite_impedance_value("steady_state_frequency_hz",
                                                 steady_state_frequency_hz)
    start_index = Int(history_storage_start_index)
    start_index > 0 || throw(ArgumentError("history_storage_start_index must be positive"))
    from_waves = _modal_complex_values("outgoing_wave_from", outgoing_wave_from, state.group_length)
    to_waves = _modal_complex_values("outgoing_wave_to", outgoing_wave_to, state.group_length)
    angular_step = 2.0 * pi * frequency * timestep

    starts = Int[]
    ends = Int[]
    storage_lengths = Int[]
    sample_counts = Int[]
    interpolation_factors = Float64[]
    from_histories = Vector{Float64}[]
    to_histories = Vector{Float64}[]
    next_index = start_index
    for mode_index in 1:state.group_length
        propagation_time = state.modal_propagation_times_s[mode_index]
        propagation_time > 0.0 ||
            throw(ArgumentError("modal propagation times must be positive"))
        storage_length = Int(floor(propagation_time / timestep + 2.0))
        storage_length >= 2 ||
            throw(ArgumentError("modal history storage length must be at least two"))
        sample_count = storage_length - 1
        delay_ratio = propagation_time / timestep
        interpolation_factor = delay_ratio - floor(delay_ratio)
        push!(starts, next_index)
        push!(ends, next_index + storage_length - 1)
        push!(storage_lengths, storage_length)
        push!(sample_counts, sample_count)
        push!(interpolation_factors, interpolation_factor)

        from_history = zeros(Float64, sample_count)
        to_history = zeros(Float64, sample_count)
        if initialized_from_steady_state
            from_wave = from_waves[mode_index]
            to_wave = to_waves[mode_index]
            for sample_index in 1:sample_count
                legacy_sample_index = sample_index + 1
                angle = (legacy_sample_index - storage_length - 1) * angular_step
                cosine = cos(angle)
                sine = sin(angle)
                from_history[sample_index] = real(from_wave) * cosine - imag(from_wave) * sine
                to_history[sample_index] = real(to_wave) * cosine - imag(to_wave) * sine
            end
        end
        push!(from_histories, from_history)
        push!(to_histories, to_history)
        next_index += storage_length + 1
    end

    return DistributedTransposedLineHistoryState(
        name,
        copy(state.phase_indices),
        copy(state.line_numbers),
        copy(state.modal_sequence_indices),
        timestep,
        frequency,
        angular_step,
        start_index,
        next_index,
        starts,
        ends,
        storage_lengths,
        sample_counts,
        interpolation_factors,
        ones(Int, state.group_length),
        real.(from_waves),
        imag.(from_waves),
        real.(to_waves),
        imag.(to_waves),
        from_histories,
        to_histories,
        Float64[],
        Float64[],
        initialized_from_steady_state,
        state.phase_count,
    )
end

function _distributed_history_interpolation_factors!(
    factors::Vector{Float64},
    state::DistributedTransposedLineModalBranchState,
    timestep::Float64,
)
    empty!(factors)
    for propagation_time in state.modal_propagation_times_s
        delay_ratio = propagation_time / timestep
        push!(factors, delay_ratio - floor(delay_ratio))
    end
    return factors
end

function _copy_real_modal_values(
    name::AbstractString,
    values::AbstractVector{<:Real},
    expected_count::Int,
)
    length(values) == expected_count ||
        throw(ArgumentError("$name must have $expected_count entries"))
    converted = Float64.(values)
    all(isfinite, converted) || throw(ArgumentError("$name entries must be finite"))
    return converted
end

function _copy_integer_modal_values(
    name::AbstractString,
    values::AbstractVector{<:Integer},
    expected_count::Int,
)
    length(values) == expected_count ||
        throw(ArgumentError("$name must have $expected_count entries"))
    converted = Int.(values)
    all(>(0), converted) || throw(ArgumentError("$name entries must be positive"))
    return converted
end

function _copy_packed_upper_symmetric_values(
    name::AbstractString,
    values::AbstractVector{<:Real},
    order::Int,
)
    expected_count = order * (order + 1) ÷ 2
    return _copy_real_modal_values(name, values, expected_count)
end

function _symmetric_packed_upper_product(
    matrix_upper::AbstractVector{Float64},
    vector::AbstractVector{Float64},
)
    order = length(vector)
    length(matrix_upper) == order * (order + 1) ÷ 2 ||
        throw(ArgumentError("matrix_upper has the wrong packed symmetric length"))
    result = zeros(Float64, order)
    matrix_index = 0
    for column in 1:order
        column_value = vector[column]
        for row in 1:column
            matrix_index += 1
            coefficient = matrix_upper[matrix_index]
            result[row] += coefficient * column_value
            row == column || (result[column] += coefficient * vector[row])
        end
    end
    return result
end

function distributed_transposed_line_initial_history_phase_matrix(
    state::DistributedTransposedLineModalBranchState;
    name::Symbol=Symbol(string(state.name), "_initial_history_phase_matrix"),
)
    state.phase_count == 3 && state.group_length == 3 ||
        throw(ArgumentError("distributed transposed line initial history phase matrix requires three phases and modes"))
    modal_initial_impedances = Float64[]
    modal_scale_values = Float64[]
    for index in 1:state.group_length
        signed_impedance = state.modal_signed_characteristic_impedances[index]
        characteristic_impedance = abs(signed_impedance)
        resistance = abs(state.modal_total_resistances[index])
        if signed_impedance < 0.0
            push!(modal_initial_impedances, signed_impedance)
            push!(modal_scale_values, inv(characteristic_impedance))
        else
            initial_impedance = characteristic_impedance - 0.25 * resistance
            shifted_impedance = characteristic_impedance + 0.25 * resistance
            initial_impedance > 0.0 ||
                throw(ArgumentError("distributed-line initial history impedance must be positive"))
            shifted_impedance > 0.0 ||
                throw(ArgumentError("distributed-line initial history scale impedance must be positive"))
            push!(modal_initial_impedances, initial_impedance)
            push!(modal_scale_values, inv(shifted_impedance))
        end
    end
    modal_transform = _transposed_line_modal_transform()
    phase_matrix = modal_transform *
        Diagonal(modal_initial_impedances) *
        transpose(modal_transform)
    return DistributedTransposedLineInitialHistoryPhaseMatrix(
        name,
        copy(state.phase_indices),
        copy(state.line_numbers),
        copy(state.modal_sequence_indices),
        modal_initial_impedances,
        modal_scale_values,
        Matrix{Float64}(phase_matrix),
        _stored_symmetric_phase_values(phase_matrix),
        state.phase_count,
    )
end

function distributed_transposed_line_initial_history_phase_terms(;
    reduced_phase_matrix_upper_values::AbstractVector{<:Real},
    matrix_terminal_from_real_values::AbstractVector{<:Real},
    matrix_terminal_from_imag_values::AbstractVector{<:Real},
    matrix_terminal_to_real_values::AbstractVector{<:Real},
    matrix_terminal_to_imag_values::AbstractVector{<:Real},
    initial_terminal_from_real_values::AbstractVector{<:Real},
    initial_terminal_from_imag_values::AbstractVector{<:Real},
    initial_terminal_to_real_values::AbstractVector{<:Real},
    initial_terminal_to_imag_values::AbstractVector{<:Real},
    phase_indices::AbstractVector{<:Integer} =
        collect(1:length(initial_terminal_from_real_values)),
    name::Symbol = :distributed_transposed_line_initial_history_phase_terms,
)
    phase_count = length(initial_terminal_from_real_values)
    phase_count > 0 ||
        throw(ArgumentError("distributed transposed line initial history phase terms require at least one phase"))
    phases = _copy_integer_modal_values("phase_indices", phase_indices, phase_count)
    length(unique(phases)) == phase_count ||
        throw(ArgumentError("phase_indices entries must be unique"))
    matrix_upper = _copy_packed_upper_symmetric_values(
        "reduced_phase_matrix_upper_values",
        reduced_phase_matrix_upper_values,
        phase_count,
    )
    from_real_source = _copy_real_modal_values(
        "matrix_terminal_from_real_values",
        matrix_terminal_from_real_values,
        phase_count,
    )
    from_imag_source = _copy_real_modal_values(
        "matrix_terminal_from_imag_values",
        matrix_terminal_from_imag_values,
        phase_count,
    )
    to_real_source = _copy_real_modal_values(
        "matrix_terminal_to_real_values",
        matrix_terminal_to_real_values,
        phase_count,
    )
    to_imag_source = _copy_real_modal_values(
        "matrix_terminal_to_imag_values",
        matrix_terminal_to_imag_values,
        phase_count,
    )
    initial_from_real = _copy_real_modal_values(
        "initial_terminal_from_real_values",
        initial_terminal_from_real_values,
        phase_count,
    )
    initial_from_imag = _copy_real_modal_values(
        "initial_terminal_from_imag_values",
        initial_terminal_from_imag_values,
        phase_count,
    )
    initial_to_real = _copy_real_modal_values(
        "initial_terminal_to_real_values",
        initial_terminal_to_real_values,
        phase_count,
    )
    initial_to_imag = _copy_real_modal_values(
        "initial_terminal_to_imag_values",
        initial_terminal_to_imag_values,
        phase_count,
    )

    from_real_product = _symmetric_packed_upper_product(matrix_upper, from_real_source)
    from_imag_product = _symmetric_packed_upper_product(matrix_upper, from_imag_source)
    to_real_product = _symmetric_packed_upper_product(matrix_upper, to_real_source)
    to_imag_product = _symmetric_packed_upper_product(matrix_upper, to_imag_source)

    return DistributedTransposedLineInitialHistoryPhaseTerms(
        name,
        phases,
        matrix_upper,
        from_real_source,
        from_imag_source,
        to_real_source,
        to_imag_source,
        initial_from_real,
        initial_from_imag,
        initial_to_real,
        initial_to_imag,
        from_real_product,
        from_imag_product,
        to_real_product,
        to_imag_product,
        initial_from_real .+ from_real_product,
        initial_from_imag .+ from_imag_product,
        initial_to_real .+ to_real_product,
        initial_to_imag .+ to_imag_product,
        phase_count,
    )
end

function distributed_transposed_line_initial_history_modal_components(;
    modal_indices::AbstractVector{<:Integer},
    term_modal_indices::AbstractVector{<:Integer},
    term_weight_values::AbstractVector{<:Real},
    term_from_real_values::AbstractVector{<:Real},
    term_from_imag_values::AbstractVector{<:Real},
    term_to_real_values::AbstractVector{<:Real},
    term_to_imag_values::AbstractVector{<:Real},
    name::Symbol = :distributed_transposed_line_initial_history_modal_components,
)
    mode_count = length(modal_indices)
    mode_count > 0 ||
        throw(ArgumentError("distributed transposed line initial history modal components require at least one mode"))
    term_count = length(term_modal_indices)
    term_count > 0 ||
        throw(ArgumentError("distributed transposed line initial history modal components require at least one term"))
    mode_indices = _copy_integer_modal_values("modal_indices", modal_indices, mode_count)
    output_index_by_mode = Dict(index => position for (position, index) in enumerate(mode_indices))
    length(output_index_by_mode) == mode_count ||
        throw(ArgumentError("modal_indices entries must be unique"))
    term_indices =
        _copy_integer_modal_values("term_modal_indices", term_modal_indices, term_count)
    weights = _copy_real_modal_values("term_weight_values", term_weight_values, term_count)
    from_real = _copy_real_modal_values(
        "term_from_real_values",
        term_from_real_values,
        term_count,
    )
    from_imag = _copy_real_modal_values(
        "term_from_imag_values",
        term_from_imag_values,
        term_count,
    )
    to_real = _copy_real_modal_values(
        "term_to_real_values",
        term_to_real_values,
        term_count,
    )
    to_imag = _copy_real_modal_values(
        "term_to_imag_values",
        term_to_imag_values,
        term_count,
    )

    from_real_products = weights .* from_real
    from_imag_products = weights .* from_imag
    to_real_products = weights .* to_real
    to_imag_products = weights .* to_imag
    component_from_real = zeros(Float64, mode_count)
    component_from_imag = zeros(Float64, mode_count)
    component_to_real = zeros(Float64, mode_count)
    component_to_imag = zeros(Float64, mode_count)
    seen_counts = zeros(Int, mode_count)
    for term_index in 1:term_count
        output_index = get(output_index_by_mode, term_indices[term_index], 0)
        output_index > 0 ||
            throw(ArgumentError("term_modal_indices entries must reference modal_indices"))
        component_from_real[output_index] += from_real_products[term_index]
        component_from_imag[output_index] += from_imag_products[term_index]
        component_to_real[output_index] += to_real_products[term_index]
        component_to_imag[output_index] += to_imag_products[term_index]
        seen_counts[output_index] += 1
    end
    all(>(0), seen_counts) ||
        throw(ArgumentError("each modal index must have at least one component term"))

    return DistributedTransposedLineInitialHistoryModalComponents(
        name,
        mode_indices,
        term_indices,
        weights,
        from_real,
        from_imag,
        to_real,
        to_imag,
        from_real_products,
        from_imag_products,
        to_real_products,
        to_imag_products,
        component_from_real,
        component_from_imag,
        component_to_real,
        component_to_imag,
        mode_count,
    )
end

function distributed_transposed_line_initial_history_phasors(;
    modal_scale_values::AbstractVector{<:Real},
    modal_component_from_real_values::AbstractVector{<:Real},
    modal_component_from_imag_values::AbstractVector{<:Real},
    modal_component_to_real_values::AbstractVector{<:Real},
    modal_component_to_imag_values::AbstractVector{<:Real},
    modal_indices::AbstractVector{<:Integer} = collect(1:length(modal_scale_values)),
    name::Symbol = :distributed_transposed_line_initial_history_phasors,
)
    mode_count = length(modal_scale_values)
    mode_count > 0 ||
        throw(ArgumentError("distributed transposed line initial history phasors require at least one mode"))
    scales = _copy_real_modal_values("modal_scale_values", modal_scale_values, mode_count)
    from_real = _copy_real_modal_values(
        "modal_component_from_real_values",
        modal_component_from_real_values,
        mode_count,
    )
    from_imag = _copy_real_modal_values(
        "modal_component_from_imag_values",
        modal_component_from_imag_values,
        mode_count,
    )
    to_real = _copy_real_modal_values(
        "modal_component_to_real_values",
        modal_component_to_real_values,
        mode_count,
    )
    to_imag = _copy_real_modal_values(
        "modal_component_to_imag_values",
        modal_component_to_imag_values,
        mode_count,
    )
    mode_indices = _copy_integer_modal_values("modal_indices", modal_indices, mode_count)

    from_amplitudes = Float64[]
    from_phases = Float64[]
    to_amplitudes = Float64[]
    to_phases = Float64[]
    for mode_index in 1:mode_count
        from_magnitude = hypot(from_real[mode_index], from_imag[mode_index])
        to_magnitude = hypot(to_real[mode_index], to_imag[mode_index])
        push!(from_amplitudes, from_magnitude * scales[mode_index])
        push!(
            from_phases,
            atan(from_imag[mode_index], from_magnitude == 0.0 ? 1.0 : from_real[mode_index]),
        )
        push!(to_amplitudes, to_magnitude * scales[mode_index])
        push!(
            to_phases,
            atan(to_imag[mode_index], to_magnitude == 0.0 ? 1.0 : to_real[mode_index]),
        )
    end

    return DistributedTransposedLineInitialHistoryPhasors(
        name,
        mode_indices,
        scales,
        from_real,
        from_imag,
        to_real,
        to_imag,
        from_amplitudes,
        from_phases,
        to_amplitudes,
        to_phases,
        mode_count,
    )
end

function distributed_transposed_line_initial_history_seed(;
    timestep_s::Real,
    steady_state_frequency_hz::Real,
    modal_delay_counts::AbstractVector{<:Integer},
    outgoing_wave_from_amplitudes::AbstractVector{<:Real},
    outgoing_wave_from_phases::AbstractVector{<:Real},
    outgoing_wave_to_amplitudes::AbstractVector{<:Real},
    outgoing_wave_to_phases::AbstractVector{<:Real},
    modal_history_interpolation_factors::AbstractVector{<:Real},
    modal_history_damping_values::AbstractVector{<:Real},
    modal_signed_characteristic_values::AbstractVector{<:Real} =
        ones(Float64, length(modal_delay_counts)),
    history_storage_start_index::Integer = 1,
    history_start_indices::Union{Nothing,AbstractVector{<:Integer}} = nothing,
    modal_indices::AbstractVector{<:Integer} = collect(1:length(modal_delay_counts)),
    name::Symbol = :distributed_transposed_line_initial_history_seed,
)
    mode_count = length(modal_delay_counts)
    mode_count > 0 ||
        throw(ArgumentError("distributed transposed line initial history seed requires at least one mode"))
    timestep = _positive_finite_impedance_value("timestep_s", timestep_s)
    frequency = _positive_finite_impedance_value("steady_state_frequency_hz",
                                                 steady_state_frequency_hz)
    angular_step = 2.0 * pi * frequency * timestep
    delay_counts = _copy_integer_modal_values("modal_delay_counts", modal_delay_counts, mode_count)
    from_amplitudes = _copy_real_modal_values(
        "outgoing_wave_from_amplitudes",
        outgoing_wave_from_amplitudes,
        mode_count,
    )
    from_phases = _copy_real_modal_values(
        "outgoing_wave_from_phases",
        outgoing_wave_from_phases,
        mode_count,
    )
    to_amplitudes = _copy_real_modal_values(
        "outgoing_wave_to_amplitudes",
        outgoing_wave_to_amplitudes,
        mode_count,
    )
    to_phases = _copy_real_modal_values(
        "outgoing_wave_to_phases",
        outgoing_wave_to_phases,
        mode_count,
    )
    interpolation_factors = _copy_real_modal_values(
        "modal_history_interpolation_factors",
        modal_history_interpolation_factors,
        mode_count,
    )
    all(value -> 0.0 <= value <= 1.0, interpolation_factors) ||
        throw(ArgumentError("modal_history_interpolation_factors must be between zero and one"))
    damping_values = _copy_real_modal_values(
        "modal_history_damping_values",
        modal_history_damping_values,
        mode_count,
    )
    characteristic_values = _copy_real_modal_values(
        "modal_signed_characteristic_values",
        modal_signed_characteristic_values,
        mode_count,
    )
    mode_indices = _copy_integer_modal_values("modal_indices", modal_indices, mode_count)
    start_index = Int(history_storage_start_index)
    start_index > 0 || throw(ArgumentError("history_storage_start_index must be positive"))
    starts = if history_start_indices === nothing
        generated = Int[]
        next_index = start_index
        for delay_count in delay_counts
            push!(generated, next_index)
            next_index += delay_count + 2
        end
        generated
    else
        _copy_integer_modal_values("history_start_indices", history_start_indices, mode_count)
    end
    all(diff(starts) .>= 0) ||
        throw(ArgumentError("history_start_indices must be ordered"))

    first_packed_index = minimum(starts)
    sinusoidal_starts = Int[]
    sinusoidal_ends = Int[]
    interpolation_indices = Int[]
    from_histories = Vector{Float64}[]
    to_histories = Vector{Float64}[]
    interpolation_from_values = Float64[]
    interpolation_to_values = Float64[]
    packed_end = first_packed_index

    for mode_index in 1:mode_count
        delay_count = delay_counts[mode_index]
        mode_start = starts[mode_index]
        fill_start = mode_start == first_packed_index ? mode_start : mode_start + 1
        fill_end = mode_start + delay_count + 1
        interpolation_index = fill_end + 1
        fill_start <= fill_end ||
            throw(ArgumentError("modal delay count leaves no sinusoidal history samples"))
        push!(sinusoidal_starts, fill_start)
        push!(sinusoidal_ends, fill_end)
        push!(interpolation_indices, interpolation_index)
        packed_end = max(packed_end, interpolation_index)

        from_history = Float64[]
        to_history = Float64[]
        for slot in fill_start:fill_end
            relative_slot = slot - mode_start
            angle = (relative_slot - delay_count - 1) * angular_step
            from_value =
                from_amplitudes[mode_index] * cos(angle + from_phases[mode_index])
            to_value =
                to_amplitudes[mode_index] * cos(angle + to_phases[mode_index])
            push!(from_history, from_value)
            push!(to_history, to_value)
        end
        push!(from_histories, from_history)
        push!(to_histories, to_history)

        damping = damping_values[mode_index]
        damping >= 0.0 ||
            throw(ArgumentError("negative modal_history_damping_values require frequency-dependent tail state"))
        first_read_slot = mode_start + 1
        first_read_offset = first_read_slot - fill_start + 1
        second_read_offset = first_read_offset + 1
        1 <= first_read_offset && second_read_offset <= length(from_history) ||
            throw(ArgumentError("initial distributed-line history seed requires two sinusoidal read samples"))
        factor = interpolation_factors[mode_index]
        from_terminal =
            from_history[first_read_offset] * factor +
            from_history[second_read_offset] * (1.0 - factor)
        to_terminal =
            to_history[first_read_offset] * factor +
            to_history[second_read_offset] * (1.0 - factor)
        if characteristic_values[mode_index] >= 0.0
            total_incident = from_terminal + to_terminal
            reflected_difference = (from_terminal - to_terminal) * damping
            push!(
                interpolation_from_values,
                0.5 * (total_incident + reflected_difference),
            )
            push!(
                interpolation_to_values,
                0.5 * (total_incident - reflected_difference),
            )
        else
            push!(interpolation_from_values, from_terminal)
            push!(interpolation_to_values, to_terminal)
        end
    end

    packed_from = zeros(Float64, packed_end - first_packed_index + 1)
    packed_to = zeros(Float64, packed_end - first_packed_index + 1)
    for mode_index in 1:mode_count
        for (offset, slot) in enumerate(
            sinusoidal_starts[mode_index]:sinusoidal_ends[mode_index],
        )
            packed_from[slot - first_packed_index + 1] =
                from_histories[mode_index][offset]
            packed_to[slot - first_packed_index + 1] =
                to_histories[mode_index][offset]
        end
        interpolation_offset =
            interpolation_indices[mode_index] - first_packed_index + 1
        packed_from[interpolation_offset] =
            interpolation_from_values[mode_index]
        packed_to[interpolation_offset] =
            interpolation_to_values[mode_index]
    end

    return DistributedTransposedLineInitialHistorySeed(
        name,
        mode_indices,
        starts,
        delay_counts,
        sinusoidal_starts,
        sinusoidal_ends,
        interpolation_indices,
        interpolation_factors,
        damping_values,
        from_histories,
        to_histories,
        interpolation_from_values,
        interpolation_to_values,
        first_packed_index,
        packed_end,
        packed_from,
        packed_to,
        mode_count,
    )
end

function distributed_transposed_line_history_state(
    state::DistributedTransposedLineModalBranchState,
    seed::DistributedTransposedLineInitialHistorySeed;
    timestep_s::Real,
    steady_state_frequency_hz::Real,
    initialized_from_steady_state::Bool=true,
    name::Symbol=Symbol(string(state.name), "_initial_history_state"),
)
    state.group_length == seed.phase_count ||
        throw(ArgumentError("distributed transposed line seed mode count must match modal branch state"))
    seed.modal_indices == collect(1:state.group_length) ||
        throw(ArgumentError("distributed transposed line seed modal indices must match branch state order"))
    timestep = _positive_finite_impedance_value("timestep_s", timestep_s)
    frequency = _positive_finite_impedance_value(
        "steady_state_frequency_hz",
        steady_state_frequency_hz,
    )
    history_from = [copy(values) for values in seed.modal_history_from_values]
    history_to = [copy(values) for values in seed.modal_history_to_values]
    sample_counts = Int[length(values) for values in history_from]
    length(history_to) == length(history_from) ||
        throw(ArgumentError("seed from/to modal history counts must match"))
    for index in eachindex(history_from)
        length(history_from[index]) == length(history_to[index]) ||
            throw(ArgumentError("seed from/to modal history lengths must match"))
        sample_counts[index] >= 2 ||
            throw(ArgumentError("seed modal history update requires at least two samples"))
    end
    storage_lengths = Int[
        seed.interpolation_indices[index] - seed.history_start_indices[index] + 1
        for index in 1:seed.phase_count
    ]
    return DistributedTransposedLineHistoryState(
        name,
        copy(state.phase_indices),
        copy(state.line_numbers),
        copy(state.modal_sequence_indices),
        timestep,
        frequency,
        2.0 * pi * frequency * timestep,
        seed.packed_history_start_index,
        seed.packed_history_end_index,
        copy(seed.history_start_indices),
        copy(seed.interpolation_indices),
        storage_lengths,
        sample_counts,
        copy(seed.modal_history_interpolation_factors),
        ones(Int, state.group_length),
        zeros(Float64, state.group_length),
        zeros(Float64, state.group_length),
        zeros(Float64, state.group_length),
        zeros(Float64, state.group_length),
        history_from,
        history_to,
        copy(seed.packed_history_from_values),
        copy(seed.packed_history_to_values),
        initialized_from_steady_state,
        state.phase_count,
    )
end
