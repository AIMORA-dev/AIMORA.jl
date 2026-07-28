
function frequency_dependent_line_modal_sample_runtime_update(
    yz_sample_matrices::AbstractVector,
    sample_rows::AbstractVector,
    target_frequency_hz::Real,
    line_length::Real,
    from_phase_voltage::AbstractVector,
    to_phase_voltage::AbstractVector;
    initial_frequency_hz::Real = target_frequency_hz,
    initial_ntol::Integer = 1,
    initial_nrp::Integer = 0,
    ntol::Integer = 2,
    nrp::Integer = 0,
)
    state = FrequencyDependentLineModalSampleRuntimeState(
        yz_sample_matrices,
        sample_rows,
        initial_frequency_hz,
        line_length;
        ntol = initial_ntol,
        nrp = initial_nrp,
    )
    return frequency_dependent_line_modal_sample_runtime_update!(
        state,
        target_frequency_hz,
        from_phase_voltage,
        to_phase_voltage;
        ntol = ntol,
        nrp = nrp,
    )
end

function _cable_frequency_scan_phase_zy_matrices(scan::CableFrequencyScanLineConstants)
    length(scan.frequencies_hz) == scan.frequency_count &&
        length(scan.frequency_constants) == scan.frequency_count &&
        length(scan.sample_rows) == scan.frequency_count ||
        throw(ArgumentError("cable frequency scan state is internally inconsistent"))
    matrices = Matrix{ComplexF64}[]
    sizehint!(matrices, scan.frequency_count)
    for (idx, constants) in pairs(scan.frequency_constants)
        abs(constants.frequency_hz - scan.frequencies_hz[idx]) <=
            _line_frequency_row_tolerance(scan.frequencies_hz[idx]) ||
            throw(ArgumentError("cable frequency scan constants must match scan frequencies"))
        constants.phase_count == scan.phase_count ||
            throw(ArgumentError("cable frequency scan phase count must be consistent"))
        push!(matrices, copy(constants.phase_zy_matrix_per_m2))
    end
    return matrices
end

function CableFrequencyScanRuntimeState(
    scan::CableFrequencyScanLineConstants,
    initial_frequency_hz::Real;
    ntol::Integer = 1,
    nrp::Integer = 0,
)
    modal_sample_state = FrequencyDependentLineModalSampleRuntimeState(
        _cable_frequency_scan_phase_zy_matrices(scan),
        scan.sample_rows,
        initial_frequency_hz,
        scan.line_length_m;
        ntol = ntol,
        nrp = nrp,
    )
    return CableFrequencyScanRuntimeState(scan, modal_sample_state, 0)
end

function cable_frequency_scan_runtime_update!(
    state::CableFrequencyScanRuntimeState,
    target_frequency_hz::Real,
    from_phase_voltage::AbstractVector,
    to_phase_voltage::AbstractVector;
    ntol::Integer = 2,
    nrp::Integer = 0,
)
    update_count_before = state.update_count
    base = frequency_dependent_line_modal_sample_runtime_update!(
        state.modal_sample_state,
        target_frequency_hz,
        from_phase_voltage,
        to_phase_voltage;
        ntol = ntol,
        nrp = nrp,
    )
    state.update_count += 1
    return merge(
        base,
        (
            source = :cable_frequency_scan_runtime_update,
            fortran_files = (:OVER47_FOR, :OVER10_FOR, :OVER11_FOR, :OVER12_FOR, :OVER13_FOR),
            fortran_routines = (:GUTS47, :ZYMX, :EIGEN, :OVER10, :SSEQIV, :OVER11, :OVER12, :OVER13),
            common_regions = (:LINEMODEL, :FDQLCL),
            cable_frequency_scan_frequencies_hz = copy(state.scan.frequencies_hz),
            cable_frequency_scan_row_count = state.scan.frequency_count,
            cable_phase_count = state.scan.phase_count,
            cable_line_length_m = state.scan.line_length_m,
            cable_frequency_scan_update_count = state.update_count,
            cable_runtime_update_count_mutated = state.update_count != update_count_before,
            bounded_cable_frequency_scan_runtime_executed = true,
            skin_effect_internal_impedance_executed =
                state.scan.skin_effect_internal_impedance_executed,
            earth_return_impedance_executed =
                state.scan.earth_return_impedance_executed,
            full_over47_series_impedance_generation_executed = false,
            full_over47_frequency_loop_executed = false,
            full_bpa_frequency_dependent_line_runtime_executed = false,
            deferred_calls = (
                :full_over47_series_impedance_generation,
                :full_over47_frequency_loop_side_effects,
                :full_dceign_lr_modal_eigenvector_bulk_package,
                :full_bpa_frequency_dependent_fitting,
                :recursive_convolution_line_runtime_oracle,
                :line_timestep_bulk_oracle,
            ),
        ),
    )
end

function cable_frequency_scan_runtime_update(
    scan::CableFrequencyScanLineConstants,
    target_frequency_hz::Real,
    from_phase_voltage::AbstractVector,
    to_phase_voltage::AbstractVector;
    initial_frequency_hz::Real = target_frequency_hz,
    initial_ntol::Integer = 1,
    initial_nrp::Integer = 0,
    ntol::Integer = 2,
    nrp::Integer = 0,
)
    state = CableFrequencyScanRuntimeState(
        scan,
        initial_frequency_hz;
        ntol = initial_ntol,
        nrp = initial_nrp,
    )
    return cable_frequency_scan_runtime_update!(
        state,
        target_frequency_hz,
        from_phase_voltage,
        to_phase_voltage;
        ntol = ntol,
        nrp = nrp,
    )
end

function CableFrequencyScanRecursiveConvolutionState(
    scan::CableFrequencyScanLineConstants,
    initial_frequency_hz::Real,
    pole_decay::AbstractMatrix,
    residue::AbstractMatrix;
    initial_ntol::Integer = 1,
    initial_nrp::Integer = 0,
    frequency_dependent_fitting_executed::Bool = false,
    frequency_loop_executed::Bool = false,
    pipe_sheath_side_effects_executed::Bool = false,
)
    recursive_state = FrequencyDependentLineRecursiveConvolutionState(
        _cable_frequency_scan_phase_zy_matrices(scan),
        scan.sample_rows,
        initial_frequency_hz,
        scan.line_length_m,
        pole_decay,
        residue;
        initial_ntol = initial_ntol,
        initial_nrp = initial_nrp,
        skin_effect_internal_impedance_executed =
            scan.skin_effect_internal_impedance_executed,
        earth_return_impedance_executed =
            scan.earth_return_impedance_executed,
        frequency_dependent_fitting_executed =
            frequency_dependent_fitting_executed,
        frequency_loop_executed = frequency_loop_executed,
        pipe_sheath_side_effects_executed = pipe_sheath_side_effects_executed,
    )
    return CableFrequencyScanRecursiveConvolutionState(scan, recursive_state, 0)
end

function CableFrequencyScanRecursiveConvolutionState(
    scan::CableFrequencyScanLineConstants,
    initial_frequency_hz::Real,
    fit::LineRecursiveConvolutionFitResult;
    initial_ntol::Integer = 1,
    initial_nrp::Integer = 0,
    frequency_dependent_fitting_executed::Bool = false,
    frequency_loop_executed::Bool = false,
    pipe_sheath_side_effects_executed::Bool = false,
)
    return CableFrequencyScanRecursiveConvolutionState(
        scan,
        initial_frequency_hz,
        fit.pole_decay,
        fit.residue;
        initial_ntol = initial_ntol,
        initial_nrp = initial_nrp,
        frequency_dependent_fitting_executed = frequency_dependent_fitting_executed,
        frequency_loop_executed = frequency_loop_executed,
        pipe_sheath_side_effects_executed = pipe_sheath_side_effects_executed,
    )
end

function cable_frequency_scan_recursive_convolution_update!(
    state::CableFrequencyScanRecursiveConvolutionState,
    target_frequency_hz::Real,
    from_phase_voltage::AbstractVector,
    to_phase_voltage::AbstractVector;
    ntol::Integer = 2,
    nrp::Integer = 0,
)
    update_count_before = state.update_count
    base = frequency_dependent_line_recursive_convolution_update!(
        state.recursive_state,
        target_frequency_hz,
        from_phase_voltage,
        to_phase_voltage;
        ntol = ntol,
        nrp = nrp,
    )
    state.update_count += 1
    return merge(
        base,
        (
            source = :cable_frequency_scan_recursive_convolution_update,
            fortran_files = (:OVER47_FOR, :OVER39_FOR, :OVER45_FOR, :OVER10_FOR, :OVER11_FOR, :OVER12_FOR, :OVER13_FOR),
            fortran_routines = (:GUTS47, :ZYMX, :EIGEN, :TDFIT, :RISE, :GUTS45, :OVER10, :SSEQIV, :OVER11, :OVER12, :OVER13),
            common_regions = (:LINEMODEL, :FDQLCL),
            cable_frequency_scan_frequencies_hz = copy(state.scan.frequencies_hz),
            cable_frequency_scan_row_count = state.scan.frequency_count,
            cable_phase_count = state.scan.phase_count,
            cable_line_length_m = state.scan.line_length_m,
            cable_recursive_convolution_update_count = state.update_count,
            cable_recursive_convolution_update_count_mutated =
                state.update_count != update_count_before,
            bounded_cable_frequency_scan_recursive_convolution_runtime_executed = true,
            skin_effect_internal_impedance_executed =
                state.scan.skin_effect_internal_impedance_executed,
            earth_return_impedance_executed =
                state.scan.earth_return_impedance_executed,
            full_over47_series_impedance_generation_executed = false,
            frequency_loop_executed =
                state.recursive_state.frequency_loop_executed,
            full_over47_frequency_loop_executed =
                state.recursive_state.frequency_loop_executed,
            pipe_sheath_side_effects_executed =
                state.recursive_state.pipe_sheath_side_effects_executed,
            full_over47_pipe_or_sheath_coupling_executed =
                state.recursive_state.pipe_sheath_side_effects_executed,
            full_bpa_frequency_dependent_fitting_executed = false,
            full_bpa_frequency_dependent_line_runtime_executed = false,
            deferred_calls = (
                :full_over47_series_impedance_generation,
                (state.recursive_state.pipe_sheath_side_effects_executed ? () : (:full_over47_pipe_or_sheath_coupling,))...,
                (state.recursive_state.frequency_loop_executed ? () : (:full_over47_frequency_loop_side_effects,))...,
                :full_bpa_frequency_dependent_fitting,
                :pole_zero_fit_oracle,
                :recursive_convolution_line_runtime_oracle,
                :line_timestep_bulk_oracle,
            ),
        ),
    )
end

function cable_frequency_scan_recursive_convolution_update(
    scan::CableFrequencyScanLineConstants,
    target_frequency_hz::Real,
    pole_decay::AbstractMatrix,
    residue::AbstractMatrix,
    from_phase_voltage::AbstractVector,
    to_phase_voltage::AbstractVector;
    initial_frequency_hz::Real = target_frequency_hz,
    initial_ntol::Integer = 1,
    initial_nrp::Integer = 0,
    ntol::Integer = 2,
    nrp::Integer = 0,
    frequency_dependent_fitting_executed::Bool = false,
    frequency_loop_executed::Bool = false,
    pipe_sheath_side_effects_executed::Bool = false,
)
    state = CableFrequencyScanRecursiveConvolutionState(
        scan,
        initial_frequency_hz,
        pole_decay,
        residue;
        initial_ntol = initial_ntol,
        initial_nrp = initial_nrp,
        frequency_dependent_fitting_executed = frequency_dependent_fitting_executed,
        frequency_loop_executed = frequency_loop_executed,
        pipe_sheath_side_effects_executed = pipe_sheath_side_effects_executed,
    )
    return cable_frequency_scan_recursive_convolution_update!(
        state,
        target_frequency_hz,
        from_phase_voltage,
        to_phase_voltage;
        ntol = ntol,
        nrp = nrp,
    )
end

function cable_frequency_scan_recursive_convolution_update_from_fit(
    scan::CableFrequencyScanLineConstants,
    fit::LineRecursiveConvolutionFitResult,
    target_frequency_hz::Real,
    from_phase_voltage::AbstractVector,
    to_phase_voltage::AbstractVector;
    initial_frequency_hz::Real = target_frequency_hz,
    initial_ntol::Integer = 1,
    initial_nrp::Integer = 0,
    ntol::Integer = 2,
    nrp::Integer = 0,
    frequency_dependent_fitting_executed::Bool = false,
    frequency_loop_executed::Bool = false,
    pipe_sheath_side_effects_executed::Bool = false,
)
    state = CableFrequencyScanRecursiveConvolutionState(
        scan,
        initial_frequency_hz,
        fit;
        initial_ntol = initial_ntol,
        initial_nrp = initial_nrp,
        frequency_dependent_fitting_executed = frequency_dependent_fitting_executed,
        frequency_loop_executed = frequency_loop_executed,
        pipe_sheath_side_effects_executed = pipe_sheath_side_effects_executed,
    )
    result = cable_frequency_scan_recursive_convolution_update!(
        state,
        target_frequency_hz,
        from_phase_voltage,
        to_phase_voltage;
        ntol = ntol,
        nrp = nrp,
    )
    return merge(
        result,
        (
            source = :cable_frequency_scan_recursive_convolution_update_from_fit,
            recursive_convolution_fit = fit,
            recursive_convolution_fit_frequency_hz = copy(fit.sample_frequencies_hz),
            recursive_convolution_fit_dt_s = fit.dt_s,
            recursive_convolution_fit_max_abs_error = fit.max_abs_error,
            bounded_recursive_convolution_fit_executed = true,
            bounded_cable_frequency_scan_recursive_convolution_runtime_executed = true,
            full_bpa_frequency_dependent_fitting_executed =
                frequency_dependent_fitting_executed,
            frequency_loop_executed = frequency_loop_executed,
            full_over47_frequency_loop_executed = frequency_loop_executed,
            pipe_sheath_side_effects_executed =
                pipe_sheath_side_effects_executed,
            full_over47_pipe_or_sheath_coupling_executed =
                pipe_sheath_side_effects_executed,
            full_bpa_frequency_dependent_line_runtime_executed = false,
            deferred_calls = (
                :full_over47_series_impedance_generation,
                (pipe_sheath_side_effects_executed ? () : (:full_over47_pipe_or_sheath_coupling,))...,
                (frequency_loop_executed ? () : (:full_over47_frequency_loop_side_effects,))...,
                :pole_zero_fit_oracle,
                :recursive_convolution_line_runtime_oracle,
                :line_timestep_bulk_oracle,
                (frequency_dependent_fitting_executed ? () : (:full_bpa_frequency_dependent_fitting,))...,
            ),
        ),
    )
end

function cable_frequency_scan_modal_response_samples(
    scan::CableFrequencyScanLineConstants;
    response_kind::Symbol = :propagated_modal_admittance,
)
    response_kind in (:propagated_modal_admittance, :modal_admittance, :propagation_factor) ||
        throw(ArgumentError("unsupported cable frequency-scan response kind"))
    length(scan.sample_rows) == scan.frequency_count ||
        throw(ArgumentError("cable frequency scan sample row count must match frequency count"))
    samples = Matrix{ComplexF64}(undef, scan.phase_count, scan.frequency_count)
    for frequency_index in 1:scan.frequency_count
        row = scan.sample_rows[frequency_index]
        length(row) == scan.phase_count ||
            throw(ArgumentError("cable frequency scan sample rows must match phase count"))
        for mode in 1:scan.phase_count
            point = row[mode]
            abs(point.frequency_hz - scan.frequencies_hz[frequency_index]) <=
                _line_frequency_row_tolerance(scan.frequencies_hz[frequency_index]) ||
                throw(ArgumentError("cable frequency scan response row frequency mismatch"))
            if response_kind == :propagation_factor
                samples[mode, frequency_index] = point.propagation_factor
            else
                abs(point.characteristic_impedance) > 0.0 ||
                    throw(ArgumentError("cable modal characteristic impedance must be nonzero"))
                modal_admittance = inv(point.characteristic_impedance)
                samples[mode, frequency_index] = response_kind == :modal_admittance ?
                    modal_admittance :
                    point.propagation_factor * modal_admittance
            end
        end
    end
    return CableFrequencyScanModalResponseSamples(
        copy(scan.frequencies_hz),
        samples,
        response_kind,
        scan.phase_count,
        scan.frequency_count,
    )
end

function cable_frequency_scan_recursive_convolution_update_from_step_response_fit(
    scan::CableFrequencyScanLineConstants,
    step_response_angular_frequencies_rad_s::AbstractVector,
    step_response_frequency_values::AbstractMatrix,
    dt_s::Real,
    target_frequency_hz::Real,
    from_phase_voltage::AbstractVector,
    to_phase_voltage::AbstractVector;
    response_kind::Symbol = :propagated_modal_admittance,
    final_values = ones(size(step_response_frequency_values, 1)),
    fit_span_s::Real,
    time_step_s::Real,
    time_start_s::Real = 0.0,
    fit_control_mode::Integer = 0,
    steady_state_angular_frequency_rad_s::Real = 2.0 * pi * 60.0,
    steady_state_shift = 1.0 + 0.0im,
    initial_frequency_hz::Real = target_frequency_hz,
    initial_ntol::Integer = 1,
    initial_nrp::Integer = 0,
    ntol::Integer = 2,
    nrp::Integer = 0,
    frequency_loop_executed::Bool = false,
    pipe_sheath_side_effects_executed::Bool = false,
)
    samples = cable_frequency_scan_modal_response_samples(
        scan;
        response_kind = response_kind,
    )
    fit_result = frequency_dependent_line_recursive_convolution_fit_from_step_response(
        samples.sample_frequencies_hz,
        samples.modal_response_samples,
        step_response_angular_frequencies_rad_s,
        step_response_frequency_values,
        dt_s;
        final_values = final_values,
        fit_span_s = fit_span_s,
        time_step_s = time_step_s,
        time_start_s = time_start_s,
        fit_control_mode = fit_control_mode,
        steady_state_angular_frequency_rad_s = steady_state_angular_frequency_rad_s,
        steady_state_shift = steady_state_shift,
    )
    result = cable_frequency_scan_recursive_convolution_update_from_fit(
        scan,
        fit_result.recursive_fit,
        target_frequency_hz,
        from_phase_voltage,
        to_phase_voltage;
        initial_frequency_hz = initial_frequency_hz,
        initial_ntol = initial_ntol,
        initial_nrp = initial_nrp,
        ntol = ntol,
        nrp = nrp,
        frequency_dependent_fitting_executed = true,
        frequency_loop_executed = frequency_loop_executed,
        pipe_sheath_side_effects_executed = pipe_sheath_side_effects_executed,
    )
    return merge(
        result,
        (
            source = :cable_frequency_scan_recursive_convolution_update_from_step_response_fit,
            cable_modal_response_samples = samples,
            cable_modal_response_kind = response_kind,
            step_response_fits = fit_result.step_response_fits,
            step_response_pole_decay = fit_result.pole_decay,
            tdfit_step_response_fit_executed = true,
            bounded_cable_frequency_scan_response_fit_executed = true,
            full_bpa_frequency_dependent_fitting_executed = true,
            deferred_calls = Tuple(
                call for call in result.deferred_calls
                if call != :full_bpa_frequency_dependent_fitting
            ),
        ),
    )
end

function cable_frequency_scan_recursive_convolution_fit(
    scan::CableFrequencyScanLineConstants,
    pole_decay::AbstractMatrix,
    dt_s::Real;
    response_kind::Symbol = :propagated_modal_admittance,
)
    samples = cable_frequency_scan_modal_response_samples(
        scan;
        response_kind = response_kind,
    )
    return frequency_dependent_line_recursive_convolution_fit(
        samples.sample_frequencies_hz,
        samples.modal_response_samples,
        pole_decay,
        dt_s,
    )
end

function cable_frequency_scan_recursive_convolution_update_from_scan_fit(
    scan::CableFrequencyScanLineConstants,
    pole_decay::AbstractMatrix,
    dt_s::Real,
    target_frequency_hz::Real,
    from_phase_voltage::AbstractVector,
    to_phase_voltage::AbstractVector;
    response_kind::Symbol = :propagated_modal_admittance,
    initial_frequency_hz::Real = target_frequency_hz,
    initial_ntol::Integer = 1,
    initial_nrp::Integer = 0,
    ntol::Integer = 2,
    nrp::Integer = 0,
    frequency_loop_executed::Bool = false,
    pipe_sheath_side_effects_executed::Bool = false,
)
    samples = cable_frequency_scan_modal_response_samples(
        scan;
        response_kind = response_kind,
    )
    fit = frequency_dependent_line_recursive_convolution_fit(
        samples.sample_frequencies_hz,
        samples.modal_response_samples,
        pole_decay,
        dt_s,
    )
    result = cable_frequency_scan_recursive_convolution_update_from_fit(
        scan,
        fit,
        target_frequency_hz,
        from_phase_voltage,
        to_phase_voltage;
        initial_frequency_hz = initial_frequency_hz,
        initial_ntol = initial_ntol,
        initial_nrp = initial_nrp,
        ntol = ntol,
        nrp = nrp,
        frequency_loop_executed = frequency_loop_executed,
        pipe_sheath_side_effects_executed = pipe_sheath_side_effects_executed,
    )
    return merge(
        result,
        (
            source = :cable_frequency_scan_recursive_convolution_update_from_scan_fit,
            cable_modal_response_samples = samples,
            cable_modal_response_kind = response_kind,
            cable_recursive_convolution_fit = fit,
            bounded_cable_frequency_scan_response_fit_executed = true,
            bounded_cable_frequency_scan_recursive_convolution_runtime_executed = true,
            frequency_loop_executed = frequency_loop_executed,
            full_over47_frequency_loop_executed = frequency_loop_executed,
            pipe_sheath_side_effects_executed =
                pipe_sheath_side_effects_executed,
            full_over47_pipe_or_sheath_coupling_executed =
                pipe_sheath_side_effects_executed,
            full_bpa_frequency_dependent_fitting_executed = false,
            full_bpa_frequency_dependent_line_runtime_executed = false,
            deferred_calls = (
                :full_over47_series_impedance_generation,
                (pipe_sheath_side_effects_executed ? () : (:full_over47_pipe_or_sheath_coupling,))...,
                (frequency_loop_executed ? () : (:full_over47_frequency_loop_side_effects,))...,
                :full_bpa_frequency_dependent_fitting,
                :pole_zero_fit_oracle,
                :recursive_convolution_line_runtime_oracle,
                :line_timestep_bulk_oracle,
            ),
        ),
    )
end

function frequency_dependent_line_recursive_convolution_update!(
    state::FrequencyDependentLineRecursiveConvolutionState,
    target_frequency_hz::Real,
    from_phase_voltage::AbstractVector,
    to_phase_voltage::AbstractVector;
    ntol::Integer = 2,
    nrp::Integer = 0,
)
    base = frequency_dependent_line_modal_sample_runtime_update!(
        state.modal_sample_state,
        target_frequency_hz,
        from_phase_voltage,
        to_phase_voltage;
        ntol = ntol,
        nrp = nrp,
    )
    mode_count, term_count = size(state.modal_history_state)
    length(base.modal_response.modal_voltage) == mode_count ||
        throw(ArgumentError("line recursive convolution modal-voltage count must match state"))
    state.previous_modal_history_state .= state.modal_history_state
    state.previous_convolution_modal_current .= state.convolution_modal_current
    state.previous_convolution_phase_current .= state.convolution_phase_current
    for mode in 1:mode_count
        modal_voltage = base.modal_response.modal_voltage[mode]
        for term in 1:term_count
            state.modal_history_state[mode, term] =
                state.pole_decay[mode, term] * state.modal_history_state[mode, term] +
                state.residue[mode, term] * modal_voltage
        end
        state.convolution_modal_current[mode] =
            sum(view(state.modal_history_state, mode, :); init = 0.0 + 0.0im)
    end
    line_phase_transform!(
        state.convolution_phase_current,
        base.modal_solution.transform,
        state.convolution_modal_current,
    )
    state.sending_phase_current .= base.sending_phase_current .+
        state.convolution_phase_current
    state.receiving_phase_current .= base.receiving_phase_current .-
        state.convolution_phase_current
    state.update_count += 1
    return merge(
        base,
        (
            source = :frequency_dependent_line_recursive_convolution_update,
            fortran_files = (:OVER39_FOR, :OVER45_FOR, :OVER10_FOR, :OVER11_FOR, :OVER12_FOR, :OVER13_FOR),
            fortran_routines = (:TDFIT, :RISE, :GUTS45, :OVER10, :SSEQIV, :OVER11, :OVER12, :OVER13),
            pole_decay = copy(state.pole_decay),
            residue = copy(state.residue),
            previous_modal_history_state =
                copy(state.previous_modal_history_state),
            modal_history_state = copy(state.modal_history_state),
            previous_convolution_modal_current =
                copy(state.previous_convolution_modal_current),
            convolution_modal_current = copy(state.convolution_modal_current),
            previous_convolution_phase_current =
                copy(state.previous_convolution_phase_current),
            convolution_phase_current = copy(state.convolution_phase_current),
            base_sending_phase_current = base.sending_phase_current,
            base_receiving_phase_current = base.receiving_phase_current,
            sending_phase_current = copy(state.sending_phase_current),
            receiving_phase_current = copy(state.receiving_phase_current),
            recursive_convolution_update_count = state.update_count,
            recursive_convolution_term_count = term_count,
            bounded_recursive_convolution_runtime_executed = true,
            skin_effect_internal_impedance_executed =
                state.skin_effect_internal_impedance_executed,
            earth_return_impedance_executed =
                state.earth_return_impedance_executed,
            frequency_dependent_fitting_executed =
                state.frequency_dependent_fitting_executed,
            full_bpa_frequency_dependent_fitting_executed =
                state.frequency_dependent_fitting_executed,
            frequency_loop_executed = state.frequency_loop_executed,
            full_over47_frequency_loop_executed = state.frequency_loop_executed,
            pipe_sheath_side_effects_executed =
                state.pipe_sheath_side_effects_executed,
            full_over47_pipe_or_sheath_coupling_executed =
                state.pipe_sheath_side_effects_executed,
            full_frequency_dependent_line_runtime_executed = false,
            deferred_calls = (
                (state.frequency_dependent_fitting_executed ? () : (:full_bpa_frequency_dependent_fitting,))...,
                (state.pipe_sheath_side_effects_executed ? () : (:full_over47_pipe_or_sheath_coupling,))...,
                (state.frequency_loop_executed ? () : (:full_over47_frequency_loop_side_effects,))...,
                :pole_zero_fit_oracle,
                :recursive_convolution_line_runtime_oracle,
                :line_timestep_bulk_oracle,
            ),
        ),
    )
end

function frequency_dependent_line_recursive_convolution_update(
    yz_sample_matrices::AbstractVector,
    sample_rows::AbstractVector,
    target_frequency_hz::Real,
    line_length::Real,
    pole_decay::AbstractMatrix,
    residue::AbstractMatrix,
    from_phase_voltage::AbstractVector,
    to_phase_voltage::AbstractVector;
    initial_frequency_hz::Real = target_frequency_hz,
    initial_ntol::Integer = 1,
    initial_nrp::Integer = 0,
    ntol::Integer = 2,
    nrp::Integer = 0,
    frequency_loop_executed::Bool = false,
    pipe_sheath_side_effects_executed::Bool = false,
)
    state = FrequencyDependentLineRecursiveConvolutionState(
        yz_sample_matrices,
        sample_rows,
        initial_frequency_hz,
        line_length,
        pole_decay,
        residue;
        initial_ntol = initial_ntol,
        initial_nrp = initial_nrp,
        frequency_loop_executed = frequency_loop_executed,
        pipe_sheath_side_effects_executed = pipe_sheath_side_effects_executed,
    )
    return frequency_dependent_line_recursive_convolution_update!(
        state,
        target_frequency_hz,
        from_phase_voltage,
        to_phase_voltage;
        ntol = ntol,
        nrp = nrp,
    )
end

function frequency_dependent_line_recursive_convolution_update_from_fit(
    yz_sample_matrices::AbstractVector,
    sample_rows::AbstractVector,
    fit::LineRecursiveConvolutionFitResult,
    target_frequency_hz::Real,
    line_length::Real,
    from_phase_voltage::AbstractVector,
    to_phase_voltage::AbstractVector;
    initial_frequency_hz::Real = target_frequency_hz,
    initial_ntol::Integer = 1,
    initial_nrp::Integer = 0,
    ntol::Integer = 2,
    nrp::Integer = 0,
    frequency_loop_executed::Bool = false,
    pipe_sheath_side_effects_executed::Bool = false,
)
    state = FrequencyDependentLineRecursiveConvolutionState(
        yz_sample_matrices,
        sample_rows,
        initial_frequency_hz,
        line_length,
        fit;
        initial_ntol = initial_ntol,
        initial_nrp = initial_nrp,
        frequency_loop_executed = frequency_loop_executed,
        pipe_sheath_side_effects_executed = pipe_sheath_side_effects_executed,
    )
    result = frequency_dependent_line_recursive_convolution_update!(
        state,
        target_frequency_hz,
        from_phase_voltage,
        to_phase_voltage;
        ntol = ntol,
        nrp = nrp,
    )
    return merge(
        result,
        (
            source = :frequency_dependent_line_recursive_convolution_update_from_fit,
            recursive_convolution_fit = fit,
            recursive_convolution_fit_frequency_hz = copy(fit.sample_frequencies_hz),
            recursive_convolution_fit_dt_s = fit.dt_s,
            recursive_convolution_fit_max_abs_error = fit.max_abs_error,
            bounded_recursive_convolution_fit_executed = true,
            full_bpa_frequency_dependent_fitting_executed = false,
            frequency_loop_executed = frequency_loop_executed,
            full_over47_frequency_loop_executed = frequency_loop_executed,
            pipe_sheath_side_effects_executed =
                pipe_sheath_side_effects_executed,
            full_over47_pipe_or_sheath_coupling_executed =
                pipe_sheath_side_effects_executed,
            full_frequency_dependent_line_runtime_executed = false,
            deferred_calls = (
                :full_bpa_frequency_dependent_fitting,
                (pipe_sheath_side_effects_executed ? () : (:full_over47_pipe_or_sheath_coupling,))...,
                (frequency_loop_executed ? () : (:full_over47_frequency_loop_side_effects,))...,
                :pole_zero_fit_oracle,
                :recursive_convolution_line_runtime_oracle,
                :line_timestep_bulk_oracle,
            ),
        ),
    )
end

const LINE_VACUUM_PERMITTIVITY_F_PER_M = 8.8541878128e-12
const LINE_VACUUM_PERMEABILITY_H_PER_M = 4.0e-7 * pi
const LINE_SKIN_EFFECT_I_APPROX_LIMIT = 3.75
const LINE_SKIN_EFFECT_K_APPROX_LIMIT = 2.0
const LINE_SKIN_EFFECT_EXPONENT_LIMIT = 43.0
const LINE_SKIN_EFFECT_DEFAULT_EPSILN = 1.0e-8
const LINE_EARTH_RETURN_EULER_SCALE = 1.781072417990
const LINE_EARTH_RETURN_SERIES_LIMIT = 5.0
const LINE_EARTH_RETURN_SERIES_TERMS = 21
const LINE_EARTH_RETURN_SERIES_CONVERGENCE = 1.0e-12
const CABLE_FREQUENCY_SCAN_LOG_DECADE = 2.302585093
const CABLE_FREQUENCY_SCAN_DISTANCE_RECORD_SCALE = 1.0e-3
const CABLE_PIPE_INFINITY = 1.0e20
const LINE_DEFAULT_LIGHT_SPEED_M_PER_S = 2.997925e8

function _checked_line_real_vector(values::AbstractVector, expected_count::Int, label::AbstractString)
    length(values) == expected_count ||
        throw(ArgumentError("$label count must be $expected_count"))
    checked = Float64.(values)
    all(isfinite, checked) ||
        throw(ArgumentError("$label entries must be finite"))
    return checked
end

function _checked_line_real_matrix(
    values::AbstractMatrix,
    expected_rows::Int,
    expected_columns::Int,
    label::AbstractString,
)
    size(values) == (expected_rows, expected_columns) ||
        throw(ArgumentError("$label size must be $(expected_rows)x$(expected_columns)"))
    checked = Matrix{Float64}(values)
    all(isfinite, checked) ||
        throw(ArgumentError("$label entries must be finite"))
    return checked
end

function _cable_surface_position_kind(isyst::Int)
    isyst < 0 && return :underground
    isyst == 0 && return :earth_surface
    return :overhead
end

function _cable_kind(itypec::Int)
    itypec == 1 && return :overhead_line
    itypec == 2 && return :non_pipe_cable
    itypec == 3 && return :pipe_type_cable
    return :cable
end

function _cable_selected_grounded_count(
    requested_grounded_count::Int,
    conductor_count::Int,
    active_phase_count::Int,
    pipe_count::Int,
    cable_kind_code::Int,
    grounded_flags::AbstractVector{<:Integer},
)
    if requested_grounded_count > 3
        length(grounded_flags) >= conductor_count ||
            throw(ArgumentError("grounded_flags must cover conductor_count"))
        for conductor in conductor_count:-1:1
            if Int(grounded_flags[conductor]) <= 0
                return conductor, false
            end
        end
        return 0, true
    elseif requested_grounded_count == 0
        return conductor_count, false
    elseif requested_grounded_count <= 1
        return cable_kind_code == 2 ? conductor_count : conductor_count - 1, false
    elseif requested_grounded_count <= 2
        return active_phase_count, false
    end
    return pipe_count, false
end

function _cable_pi_section_report_state(
    cable_kind_code::Int,
    pi_section_count::Int,
    crossbonded::Int,
    total_length_m::Float64,
    section_length_m::Float64,
    sheath_grounding_resistance_ohm::Float64,
)
    punch_requested = pi_section_count != 0 || crossbonded != 0
    section_count = abs(pi_section_count)
    modeling_kind =
        !punch_requested ? :none :
        cable_kind_code == 1 ? :overhead_line :
        pi_section_count < 0 ? :discrete_cable :
        :homogeneous_cable
    return CablePiSectionReportState(
        punch_requested,
        total_length_m,
        section_length_m,
        section_count,
        crossbonded != 0,
        modeling_kind,
        sheath_grounding_resistance_ohm,
        punch_requested,
    )
end

function _cable_crossbond_layout_validation_error(
    cable_kind_code::Int,
    phase_count::Int,
    conductor_count::Int,
    requested_grounded_count::Int,
    crossbonded::Int,
    pi_section_count::Int,
    layer_counts::AbstractVector{Int},
)
    (cable_kind_code == 1 || crossbonded == 0 || pi_section_count < 0) && return :none
    phase_count == 3 || return :phase_count_must_be_three
    length(layer_counts) >= 3 || return :phase_layer_count_missing
    for phase in 1:3
        (layer_counts[phase] != 1 && layer_counts[phase] <= 3) ||
            return :invalid_phase_layer_count
    end
    conductor_count == 6 + requested_grounded_count ||
        return :crossbond_reduction_count_mismatch
    layer_counts[1] <= layer_counts[2] <= layer_counts[3] ||
        return :crossbond_layer_order
    return :none
end

function _cable_pi_model_validation_error(
    phase_count::Int,
    reduction_count::Int,
    pi_section_count::Int,
    crossbonded::Int,
    separate_insulation_sections::Bool,
    layer_counts::AbstractVector{Int},
)
    pi_section_count == 0 && return :none
    if pi_section_count >= 0 && crossbonded == 0
        return separate_insulation_sections ? :separate_insulation_without_crossbond : :none
    end

    if pi_section_count < 0 && crossbonded != 0
        reduction_count in (6, 7) || return :discrete_crossbond_reduction_count
    elseif pi_section_count < 0
        reduction_count > 7 && return :discrete_reduction_count_too_large
        if reduction_count >= 6
            nothing
        elseif reduction_count == 2
            phase_count == 1 || return :two_conductor_discrete_requires_single_phase
            return :none
        elseif reduction_count == 3
            separate_insulation_sections ||
                return :three_conductor_discrete_requires_separate_insulation
            phase_count == 1 || return :three_conductor_discrete_requires_single_phase
            return :none
        else
            return :invalid_discrete_reduction_count
        end
    end

    phase_count == 3 || return :phase_count_must_be_three
    length(layer_counts) >= 3 || return :phase_layer_count_missing
    for phase in 1:3
        (layer_counts[phase] != 1 && layer_counts[phase] <= 3) ||
            return :invalid_phase_layer_count
    end
    layer_counts[1] <= layer_counts[2] <= layer_counts[3] ||
        return :crossbond_layer_order
    reduction_count == 6 && return :none
    (pi_section_count > 0 && crossbonded != 0 && reduction_count == 4) &&
        return :none
    layer_counts[2] == 2 || return :middle_phase_requires_two_conductors
    return :none
end

function cable_pipe_sheath_derived_state(;
    cable_kind_code::Integer,
    surface_position_code::Integer,
    conductor_count::Integer,
    active_phase_count::Integer,
    pipe_count::Integer,
    active_phase_count_without_pipe::Integer,
    requested_grounded_count::Integer,
    grounded_flags::AbstractVector{<:Integer},
    layer_counts::AbstractVector{<:Integer},
    boundary_radii_m::AbstractMatrix,
    resistivity_ohm_m::AbstractMatrix,
    relative_permeability::AbstractMatrix,
    relative_permittivity::AbstractMatrix,
    pipe_radii_m::AbstractVector,
    pipe_resistivity_ohm_m::Real,
    pipe_relative_permeability::Real,
    pipe_inner_insulator_relative_permittivity::Real,
    pipe_outer_insulator_relative_permittivity::Real,
    conductor_depths_m::AbstractVector,
    conductor_distances_m::AbstractVector,
    conductor_pipe_center_distances_m::AbstractVector,
    conductor_angles_rad::AbstractVector,
    pi_section_count::Integer = 0,
    crossbonded::Integer = 0,
    separate_insulation_sections::Bool = false,
    total_length_m::Real = 0.0,
    section_length_m::Real = 0.0,
    sheath_grounding_resistance_ohm::Real = 0.0,
    light_speed_m_per_s::Real = LINE_DEFAULT_LIGHT_SPEED_M_PER_S,
)
    cable_kind = Int(cable_kind_code)
    surface_position = Int(surface_position_code)
    conductors = Int(conductor_count)
    phases = Int(active_phase_count)
    pipes = Int(pipe_count)
    phases_without_pipe = Int(active_phase_count_without_pipe)
    requested_grounded = Int(requested_grounded_count)
    conductors > 0 || throw(ArgumentError("conductor_count must be positive"))
    phases > 0 || throw(ArgumentError("active_phase_count must be positive"))
    pipes >= 0 || throw(ArgumentError("pipe_count must be nonnegative"))
    phases_without_pipe >= 0 ||
        throw(ArgumentError("active_phase_count_without_pipe must be nonnegative"))
    layer_count_values = Int.(collect(layer_counts))
    length(layer_count_values) == phases ||
        throw(ArgumentError("layer_counts must cover active_phase_count"))
    all(count -> 1 <= count <= 3, layer_count_values) ||
        throw(ArgumentError("layer_counts entries must be between 1 and 3"))
    radii = _checked_line_real_matrix(boundary_radii_m, phases, 7, "boundary_radii_m")
    resistivity = _checked_line_real_matrix(resistivity_ohm_m, phases, 3, "resistivity_ohm_m")
    permeability = _checked_line_real_matrix(relative_permeability, phases, 3, "relative_permeability")
    permittivity = _checked_line_real_matrix(relative_permittivity, phases, 3, "relative_permittivity")
    pipe_radii = _checked_line_real_vector(pipe_radii_m, 3, "pipe_radii_m")
    depths = _checked_line_real_vector(conductor_depths_m, phases, "conductor_depths_m")
    distances = _checked_line_real_vector(conductor_distances_m, phases, "conductor_distances_m")
    pipe_center_distances =
        _checked_line_real_vector(conductor_pipe_center_distances_m, phases, "conductor_pipe_center_distances_m")
    angles = _checked_line_real_vector(conductor_angles_rad, phases, "conductor_angles_rad")
    pipe_resistivity = _checked_line_positive(pipe_resistivity_ohm_m, "pipe_resistivity_ohm_m")
    pipe_permeability = _checked_line_positive(pipe_relative_permeability, "pipe_relative_permeability")
    pipe_inner_permittivity =
        _checked_line_positive(pipe_inner_insulator_relative_permittivity, "pipe_inner_insulator_relative_permittivity")
    raw_outer_permittivity = Float64(pipe_outer_insulator_relative_permittivity)
    isfinite(raw_outer_permittivity) && raw_outer_permittivity >= 0.0 ||
        throw(ArgumentError("pipe_outer_insulator_relative_permittivity must be finite and nonnegative"))
    outer_defaulted = raw_outer_permittivity == 0.0
    pipe_outer_permittivity = outer_defaulted ? 1.0 : raw_outer_permittivity
    light_speed = _checked_line_positive(light_speed_m_per_s, "light_speed_m_per_s")
    total_length = Float64(total_length_m)
    section_length = Float64(section_length_m)
    sheath_resistance = Float64(sheath_grounding_resistance_ohm)
    all(isfinite, (total_length, section_length, sheath_resistance)) ||
        throw(ArgumentError("pi-section report scalar values must be finite"))

    selected_grounded, no_ungrounded_stop = _cable_selected_grounded_count(
        requested_grounded,
        conductors,
        phases_without_pipe,
        pipes,
        cable_kind,
        grounded_flags,
    )

    pipe_branch_executed = cable_kind != 2
    pipe_return_included = pipe_branch_executed && pipes != 0
    pipe_radii_mutated = copy(pipe_radii)
    if pipe_branch_executed && !pipe_return_included
        pipe_radii_mutated[2] = CABLE_PIPE_INFINITY
        pipe_radii_mutated[3] = CABLE_PIPE_INFINITY
    end

    layer_wave_speeds = zeros(Float64, phases, 3)
    for phase in 1:phases, layer in 1:3
        layer_wave_speeds[phase, layer] = light_speed / sqrt(permittivity[phase, layer])
    end

    u0 = LINE_VACUUM_PERMEABILITY_H_PER_M
    core_inner = similar(depths)
    core_outer = similar(depths)
    sheath_inner = similar(depths)
    sheath_outer = similar(depths)
    armor_inner = similar(depths)
    armor_outer = similar(depths)
    core_log = similar(depths)
    sheath_log = similar(depths)
    armor_log = similar(depths)
    for phase in 1:phases
        core_inner[phase] = radii[phase, 1] * sqrt(u0 / resistivity[phase, 1] * permeability[phase, 1])
        core_outer[phase] = radii[phase, 2] * sqrt(u0 / resistivity[phase, 1] * permeability[phase, 1])
        core_log[phase] = log(radii[phase, 3] / radii[phase, 2])
        sheath_inner[phase] = radii[phase, 3] * sqrt(u0 / resistivity[phase, 2] * permeability[phase, 2])
        sheath_outer[phase] = radii[phase, 4] * sqrt(u0 / resistivity[phase, 2] * permeability[phase, 2])
        sheath_log[phase] = log(radii[phase, 5] / radii[phase, 4])
        armor_inner[phase] = radii[phase, 5] * sqrt(u0 / resistivity[phase, 3] * permeability[phase, 3])
        armor_outer[phase] = radii[phase, 6] * sqrt(u0 / resistivity[phase, 3] * permeability[phase, 3])
        armor_log[phase] = log(radii[phase, 7] / radii[phase, 6])
    end
    all(isfinite, vcat(core_inner, core_outer, sheath_inner, sheath_outer, armor_inner, armor_outer, core_log, sheath_log, armor_log)) ||
        throw(ArgumentError("derived cable pipe/sheath values must be finite"))

    pipe_limit = NaN
    pipe_radius_valid = true
    if cable_kind == 3
        pipe_limit = (1.0 + 2.0 / sqrt(3.0)) * radii[1, 7]
        phases == 2 && (pipe_limit = 2.0 * radii[1, 7])
        phases == 1 && (pipe_limit = radii[1, 7])
        pipe_radius_valid = pipe_radii_mutated[1] > pipe_limit
    end

    layer_count_order_valid = true
    for phase in 2:phases
        if layer_count_values[phase - 1] < layer_count_values[phase]
            layer_count_order_valid = false
            break
        end
    end

    matrix_width =
        cable_kind == 3 ? pipes : phases
    admittance_count = conductors
    pi_reduction_count =
        Int(pi_section_count) >= 0 && Int(crossbonded) != 0 ? 4 : selected_grounded
    crossbond_validation_error = _cable_crossbond_layout_validation_error(
        cable_kind,
        phases,
        conductors,
        requested_grounded,
        Int(crossbonded),
        Int(pi_section_count),
        layer_count_values,
    )
    pi_validation_error = _cable_pi_model_validation_error(
        phases,
        pi_reduction_count,
        Int(pi_section_count),
        Int(crossbonded),
        separate_insulation_sections,
        layer_count_values,
    )
    pi_report = _cable_pi_section_report_state(
        cable_kind,
        Int(pi_section_count),
        Int(crossbonded),
        total_length,
        section_length,
        sheath_resistance,
    )

    return CablePipeSheathDerivedState(
        _cable_kind(cable_kind),
        _cable_surface_position_kind(surface_position),
        phases,
        conductors,
        pipes,
        selected_grounded,
        matrix_width,
        admittance_count,
        pi_reduction_count,
        vec(radii[:, 7]),
        depths,
        distances,
        pipe_center_distances,
        angles,
        pipe_return_included,
        pipe_radii_mutated,
        pipe_resistivity,
        pipe_permeability,
        pipe_inner_permittivity,
        pipe_outer_permittivity,
        outer_defaulted,
        layer_wave_speeds,
        core_inner,
        core_outer,
        sheath_inner,
        sheath_outer,
        armor_inner,
        armor_outer,
        core_log,
        sheath_log,
        armor_log,
        pipe_limit,
        pipe_radius_valid,
        layer_count_order_valid,
        no_ungrounded_stop,
        crossbond_validation_error == :none,
        crossbond_validation_error,
        pi_validation_error == :none,
        pi_validation_error,
        pi_report,
        true,
    )
end

function _cable_frequency_scan_retained_value(
    value,
    fallback::Float64,
    label::AbstractString,
)
    if value === nothing
        return fallback > 0.0 ? fallback : NaN
    end
    checked = Float64(value)
    checked == -Inf && return fallback > 0.0 ? fallback : NaN
    isfinite(checked) || throw(ArgumentError("$label must be finite"))
    return checked
end

function cable_frequency_scan_loop_schedule(
    start_frequency_hz::Real,
    decade_count::Integer,
    points_per_decade::Integer;
    earth_resistivity_ohm_m::Real = 0.0,
    distance_m::Real = 0.0,
    steady_state_frequency_hz::Real = 60.0,
    retained_distance_m = nothing,
    retained_earth_resistivity_ohm_m = nothing,
    alteration_mode::Integer = 0,
    last_output_family::Integer = 0,
    card_output_flag::Integer = 0,
    transform_flag::Integer = 0,
    modal_output_enabled::Bool = false,
    matching_root_count::Integer = 0,
    retry_flag::Integer = 0,
    history_capacity::Integer = typemax(Int),
)
    raw_start = Float64(start_frequency_hz)
    steady = Float64(steady_state_frequency_hz)
    isfinite(raw_start) && raw_start >= 0.0 ||
        throw(ArgumentError("start_frequency_hz must be finite and nonnegative"))
    isfinite(steady) && steady > 0.0 ||
        throw(ArgumentError("steady_state_frequency_hz must be finite and positive"))
    start = raw_start == 0.0 ? steady : raw_start
    decades = Int(decade_count)
    decades >= 0 || throw(ArgumentError("decade_count must be nonnegative"))
    points = Int(points_per_decade)
    points >= 0 || throw(ArgumentError("points_per_decade must be nonnegative"))
    points == 0 && (points = 1)
    earth = Float64(earth_resistivity_ohm_m)
    isfinite(earth) || throw(ArgumentError("earth_resistivity_ohm_m must be finite"))
    distance = Float64(distance_m)
    isfinite(distance) || throw(ArgumentError("distance_m must be finite"))
    alteration = Int(alteration_mode)
    retained_distance = _cable_frequency_scan_retained_value(
        retained_distance_m,
        distance,
        "retained_distance_m",
    )
    retained_earth = _cable_frequency_scan_retained_value(
        retained_earth_resistivity_ohm_m,
        earth,
        "retained_earth_resistivity_ohm_m",
    )
    if alteration == 2 && isfinite(retained_earth)
        earth = retained_earth
    end

    ratio = exp(CABLE_FREQUENCY_SCAN_LOG_DECADE / points)
    frequencies = Float64[]
    base_frequencies = Float64[]
    decade_indices = Int[]
    interval_indices = Int[]
    output_indices = Int[]
    iprint = 0
    iii = 0
    kkk = 1
    fdecad = start
    final_distance = distance
    voltage_control_values = zeros(Float64, 6)

    function push_frequency!(frequency::Float64, base::Float64, decade_index::Int, interval_index::Int)
        push!(frequencies, frequency)
        push!(base_frequencies, base)
        push!(decade_indices, decade_index)
        push!(interval_indices, interval_index)
        push!(output_indices, iprint)
        return nothing
    end

    if decades > 0
        iprint += 1
        voltage_control_values .= (
            start,
            start * 10.0^decades,
            ratio,
            log((start * 10.0^decades) / start) / log(ratio) + 1.5,
            Float64(decades),
            Float64(points),
        )
        final_distance = 0.0
        push_frequency!(start, start, 0, 0)
        while iii != decades
            pkkk = kkk
            frequency = fdecad * exp(Float64(pkkk) * CABLE_FREQUENCY_SCAN_LOG_DECADE / points)
            kkk += 1
            if kkk > points
                iii += 1
                fdecad *= 10.0
                iprint += 1
                kkk = 1
                push_frequency!(fdecad, fdecad / 10.0, iii - 1, points)
            else
                iprint += 1
                push_frequency!(frequency, fdecad, iii, pkkk)
            end
        end
    else
        iprint += 1
        push_frequency!(start, start, 0, 0)
    end

    history_limit = Int(history_capacity)
    history_limit >= 0 || throw(ArgumentError("history_capacity must be nonnegative"))
    resistivity_history_values = Float64[]
    resistivity_history_frequencies = Float64[]
    if alteration == 1 && decades <= 0
        iprint <= history_limit ||
            throw(ArgumentError("frequency scan resistivity history capacity was exceeded"))
        push!(resistivity_history_values, earth)
        push!(resistivity_history_frequencies, start)
    end

    retry = false
    final_retry_flag = Int(retry_flag)
    matching_roots = Int(matching_root_count)
    matching_roots >= 0 || throw(ArgumentError("matching_root_count must be nonnegative"))
    if length(frequencies) > 1 && final_retry_flag <= 0
        ratio_found = matching_roots / Float64(iprint - 1)
        if ratio_found >= 0.75
            retry = true
            final_retry_flag = 1
            output_indices .= output_indices .+ 1
            iprint += 1
        end
    end

    unit_distance_record_count =
        decades <= 0 && Int(last_output_family) == 39 ? 1 : 0
    unit_distance_km =
        unit_distance_record_count == 1 ?
        distance * CABLE_FREQUENCY_SCAN_DISTANCE_RECORD_SCALE :
        NaN

    return CableFrequencyScanLoopSchedule(
        start,
        steady,
        decades,
        points,
        length(frequencies),
        frequencies,
        base_frequencies,
        decade_indices,
        interval_indices,
        output_indices,
        ratio,
        Float64(earth_resistivity_ohm_m),
        earth,
        fill(earth, length(frequencies)),
        distance,
        final_distance,
        fill(final_distance, length(frequencies)),
        retained_distance,
        retained_earth,
        voltage_control_values,
        resistivity_history_values,
        resistivity_history_frequencies,
        unit_distance_record_count,
        unit_distance_km,
        Int(card_output_flag),
        Int(transform_flag),
        alteration,
        modal_output_enabled || unit_distance_record_count == 1,
        retry,
        final_retry_flag,
        isempty(output_indices) ? 0 : output_indices[end],
        true,
    )
end

function _checked_line_positive(value::Real, label::AbstractString)
    checked = Float64(value)
    isfinite(checked) && checked > 0.0 ||
        throw(ArgumentError("$label must be finite and positive"))
    return checked
end

function _checked_line_nonnegative(value::Real, label::AbstractString)
    checked = Float64(value)
    isfinite(checked) && checked >= 0.0 ||
        throw(ArgumentError("$label must be finite and nonnegative"))
    return checked
end

function CableGeometryConductor(
    radius_m::Real,
    horizontal_position_m::Real,
    depth_m::Real;
    resistivity_ohm_m::Real = 1.724e-8,
    relative_permittivity::Real = 1.0,
    relative_permeability::Real = 1.0,
)
    radius = _checked_line_positive(radius_m, "cable conductor radius_m")
    horizontal = Float64(horizontal_position_m)
    isfinite(horizontal) ||
        throw(ArgumentError("cable conductor horizontal_position_m must be finite"))
    depth = _checked_line_positive(depth_m, "cable conductor depth_m")
    2.0 * depth > radius ||
        throw(ArgumentError("cable conductor image distance must exceed its radius"))
    resistivity = _checked_line_positive(resistivity_ohm_m, "cable conductor resistivity_ohm_m")
    epsilon_r = _checked_line_positive(relative_permittivity, "cable conductor relative_permittivity")
    mu_r = _checked_line_positive(relative_permeability, "cable conductor relative_permeability")
    return CableGeometryConductor(
        radius,
        horizontal,
        depth,
        resistivity,
        epsilon_r,
        mu_r,
    )
end

function _checked_cable_phase_counts(
    phase_conductor_counts,
    grounded_conductor_count::Integer,
    conductor_count::Int,
)
    grounded = Int(grounded_conductor_count)
    grounded >= 0 ||
        throw(ArgumentError("grounded_conductor_count must be nonnegative"))
    counts = phase_conductor_counts === nothing ?
        fill(1, conductor_count - grounded) : collect(Int.(phase_conductor_counts))
    !isempty(counts) || throw(ArgumentError("phase_conductor_counts must be nonempty"))
    all(>(0), counts) ||
        throw(ArgumentError("phase_conductor_counts entries must be positive"))
    sum(counts) + grounded == conductor_count ||
        throw(ArgumentError("phase_conductor_counts plus grounded_conductor_count must match conductor count"))
    return counts, grounded
end

function cable_geometry_constants(
    conductors::AbstractVector{CableGeometryConductor};
    phase_conductor_counts = nothing,
    grounded_conductor_count::Integer = 0,
)
    conductor_count = length(conductors)
    conductor_count > 0 ||
        throw(ArgumentError("cable geometry requires at least one conductor"))
    counts, grounded = _checked_cable_phase_counts(
        phase_conductor_counts,
        grounded_conductor_count,
        conductor_count,
    )
    radius = [conductor.radius_m for conductor in conductors]
    horizontal = [conductor.horizontal_position_m for conductor in conductors]
    depth = [conductor.depth_m for conductor in conductors]
    direct = Matrix{Float64}(undef, conductor_count, conductor_count)
    image = Matrix{Float64}(undef, conductor_count, conductor_count)
    angle = Matrix{Float64}(undef, conductor_count, conductor_count)
    potential_log = Matrix{Float64}(undef, conductor_count, conductor_count)
    for row in 1:conductor_count, col in 1:conductor_count
        dx = horizontal[row] - horizontal[col]
        dy = depth[row] - depth[col]
        direct_distance = row == col ? radius[row] : hypot(dx, dy)
        direct_distance > 0.0 ||
            throw(ArgumentError("cable conductor pair distance must be positive"))
        image_distance = hypot(dx, depth[row] + depth[col])
        image_distance > direct_distance ||
            throw(ArgumentError("cable image distance must exceed direct conductor distance"))
        direct[row, col] = direct_distance
        image[row, col] = image_distance
        angle[row, col] = row == col ? 0.0 :
            atan(abs(dx) / (depth[row] + depth[col]))
        potential_log[row, col] = log(image_distance / direct_distance)
    end
    resistivity = [conductor.resistivity_ohm_m for conductor in conductors]
    epsilon_r = [conductor.relative_permittivity for conductor in conductors]
    mu_r = [conductor.relative_permeability for conductor in conductors]
    return CableGeometryConstants(
        collect(conductors),
        counts,
        grounded,
        radius,
        horizontal,
        depth,
        direct,
        image,
        angle,
        potential_log,
        resistivity,
        inv.(resistivity),
        epsilon_r,
        LINE_VACUUM_PERMITTIVITY_F_PER_M .* epsilon_r,
        mu_r,
        LINE_VACUUM_PERMEABILITY_H_PER_M .* mu_r,
        conductor_count,
        length(counts),
    )
end

function _checked_cable_electrostatic_permittivity(
    constants::CableGeometryConstants,
    relative_permittivity,
)
    if relative_permittivity === nothing
        epsilon_r = constants.relative_permittivity[1]
        tolerance = max(1.0e-12, 64.0 * eps(Float64) * max(1.0, abs(epsilon_r)))
        all(value -> abs(value - epsilon_r) <= tolerance, constants.relative_permittivity) ||
            throw(ArgumentError("relative_permittivity must be supplied for nonhomogeneous cable dielectrics"))
        return epsilon_r
    end
    return _checked_line_positive(relative_permittivity, "cable electrostatic relative_permittivity")
end

function cable_geometry_electrostatic_admittance(
    constants::CableGeometryConstants,
    frequency_hz::Real;
    relative_permittivity = nothing,
)
    frequency = Float64(frequency_hz)
    isfinite(frequency) && frequency >= 0.0 ||
        throw(ArgumentError("frequency_hz must be finite and nonnegative"))
    epsilon_r = _checked_cable_electrostatic_permittivity(constants, relative_permittivity)
    scale = inv(2.0 * pi * LINE_VACUUM_PERMITTIVITY_F_PER_M * epsilon_r)
    potential = scale .* constants.potential_log_matrix
    all(isfinite, potential) ||
        throw(ArgumentError("cable potential coefficient matrix entries must be finite"))
    capacitance = inv(Symmetric(potential))
    all(isfinite, capacitance) ||
        throw(ArgumentError("cable capacitance matrix entries must be finite"))
    omega = 2.0 * pi * frequency
    shunt = ComplexF64.(0.0, omega .* capacitance)
    return CableElectrostaticAdmittance(
        frequency,
        omega,
        epsilon_r,
        Matrix{Float64}(potential),
        Matrix{Float64}(capacitance),
        shunt,
        constants.conductor_count,
        constants.phase_count,
    )
end

function _cable_phase_conductor_matrix(constants::CableGeometryConstants)
    incidence = zeros(Float64, constants.phase_count, constants.conductor_count)
    conductor = 1
    for phase in 1:constants.phase_count
        for _ in 1:constants.phase_conductor_counts[phase]
            incidence[phase, conductor] = 1.0
            conductor += 1
        end
    end
    return incidence
end

function _cable_phase_average_matrix(constants::CableGeometryConstants)
    average = _cable_phase_conductor_matrix(constants)
    for phase in 1:constants.phase_count
        average[phase, :] ./= constants.phase_conductor_counts[phase]
    end
    return average
end

function _line_skin_effect_bessel_i0_i1_k0_k1(x::ComplexF64, scaled::Bool)
    xa = abs(x)
    if xa <= LINE_SKIN_EFFECT_I_APPROX_LIMIT
        y = x / LINE_SKIN_EFFECT_I_APPROX_LIMIT
        y1 = y * y
        y2 = y1 * y1
        y3 = y2 * y1
        y4 = y3 * y1
        y5 = y4 * y1
        y6 = y5 * y1
        i0 = 1.0 + 3.5156229 * y1 + 3.0899424 * y2 + 1.2067492 * y3 +
            0.2659732 * y4 + 0.0360768 * y5 + 0.0045813 * y6
        i1 = x * (0.5 + 0.87890594 * y1 + 0.51498869 * y2 +
            0.15084934 * y3 + 0.02658733 * y4 + 0.00301532 * y5 +
            0.00032411 * y6)
    else
        y = x / LINE_SKIN_EFFECT_I_APPROX_LIMIT
        y0 = sqrt(x)
        scaled || (y0 *= exp(-x))
        y1 = inv(y)
        y2 = y1 / y
        y3 = y2 / y
        y4 = y3 / y
        y5 = y4 / y
        y6 = y5 / y
        y7 = y6 / y
        y8 = y7 / y
        i0 = (0.39894228 + 0.01328592 * y1 + 0.00225319 * y2 -
            0.00157565 * y3 + 0.00916281 * y4 - 0.02057706 * y5 +
            0.02635537 * y6 - 0.01647633 * y7 + 0.00392377 * y8) / y0
        i1 = (0.39894228 - 0.03988024 * y1 - 0.00362018 * y2 +
            0.00163801 * y3 - 0.01031555 * y4 + 0.02282967 * y5 -
            0.02895312 * y6 + 0.01787654 * y7 - 0.00420059 * y8) / y0
    end

    if xa <= LINE_SKIN_EFFECT_K_APPROX_LIMIT
        y = x / LINE_SKIN_EFFECT_K_APPROX_LIMIT
        y0 = log(y)
        y1 = y * y
        y2 = y1 * y1
        y3 = y2 * y1
        y4 = y3 * y1
        y5 = y4 * y1
        y6 = y5 * y1
        k0 = -y0 * i0 - 0.57721566 + 0.42278420 * y1 + 0.23069756 * y2 +
            0.03488590 * y3 + 0.00262698 * y4 + 0.00010750 * y5 +
            0.00000740 * y6
        k1 = y0 * i1 + (1.0 + 0.15443144 * y1 - 0.67278579 * y2 -
            0.18156897 * y3 - 0.01919402 * y4 - 0.00110404 * y5 -
            0.00004686 * y6) / x
    else
        y = ComplexF64(LINE_SKIN_EFFECT_K_APPROX_LIMIT) / x
        y0 = sqrt(x)
        scaled || (y0 *= exp(x))
        y1 = y * y
        y2 = y1 * y
        y3 = y2 * y
        y4 = y3 * y
        y5 = y4 * y
        k0 = (1.25331414 - 0.07832358 * y + 0.02189568 * y1 -
            0.01062446 * y2 + 0.00587872 * y3 - 0.00251540 * y4 +
            0.00053208 * y5) / y0
        k1 = (1.25331414 + 0.23498619 * y - 0.03655620 * y1 +
            0.01504268 * y2 - 0.00780353 * y3 + 0.00325614 * y4 -
            0.00068245 * y5) / y0
    end
    return ComplexF64(i0), ComplexF64(i1), ComplexF64(k0), ComplexF64(k1)
end

function cable_skin_effect_internal_impedance(
    inner_diffusion_factor::Real,
    outer_diffusion_factor::Real,
    relative_permeability::Real,
    frequency_hz::Real;
    epsiln::Real = LINE_SKIN_EFFECT_DEFAULT_EPSILN,
    exponent_limit::Real = LINE_SKIN_EFFECT_EXPONENT_LIMIT,
)
    inner = _checked_line_nonnegative(inner_diffusion_factor, "skin-effect inner diffusion factor")
    outer = _checked_line_positive(outer_diffusion_factor, "skin-effect outer diffusion factor")
    inner <= outer ||
        throw(ArgumentError("skin-effect inner diffusion factor must not exceed the outer factor"))
    permeability_ratio =
        _checked_line_positive(relative_permeability, "skin-effect relative_permeability")
    frequency = _checked_line_positive(frequency_hz, "skin-effect frequency_hz")
    tolerance = _checked_line_positive(epsiln, "skin-effect epsiln")
    overflow_limit =
        _checked_line_positive(exponent_limit, "skin-effect exponent_limit")

    omega = 2.0 * pi * frequency
    cjw = ComplexF64(0.0, omega)
    sjw = sqrt(cjw)
    u2p = LINE_VACUUM_PERMEABILITY_H_PER_M / (2.0 * pi)
    x2 = ComplexF64(outer) * sjw
    small_inner_threshold = 0.5 * tolerance / 1000.0 * 1.0e-3 / 100.0
    if inner < small_inner_threshold
        i0, i1, _, _ = _line_skin_effect_bessel_i0_i1_k0_k1(x2, true)
        return ComplexF64(cjw * u2p / x2 * i0 / i1 * permeability_ratio)
    end

    x1 = ComplexF64(inner) * sjw
    use_scaled = abs(x2) >= 10.0 && abs(x1) > LINE_SKIN_EFFECT_I_APPROX_LIMIT
    _, inner_i1, _, inner_k1 =
        _line_skin_effect_bessel_i0_i1_k0_k1(x1, use_scaled)
    outer_i0, outer_i1, outer_k0, outer_k1 =
        _line_skin_effect_bessel_i0_i1_k0_k1(x2, use_scaled)
    if use_scaled
        exponent = x2 - x1
        if abs(exponent) > overflow_limit
            return ComplexF64(cjw * u2p / x2 * outer_i0 / outer_i1 * permeability_ratio)
        end
        scale = exp(exponent)
        numerator = outer_i0 * inner_k1 * scale + outer_k0 * inner_i1 / scale
        denominator = outer_i1 * inner_k1 * scale - outer_k1 * inner_i1 / scale
        return ComplexF64(cjw * u2p / x2 * numerator / denominator * permeability_ratio)
    end
    numerator = outer_i0 * inner_k1 + outer_k0 * inner_i1
    denominator = outer_i1 * inner_k1 - outer_k1 * inner_i1
    return ComplexF64(cjw * u2p / x2 * numerator / denominator * permeability_ratio)
end

function cable_bounded_skin_effect_internal_impedance(
    radius_m::Real,
    resistivity_ohm_m::Real,
    permeability_h_per_m::Real,
    frequency_hz::Real,
)
    radius = _checked_line_positive(radius_m, "skin-effect conductor radius_m")
    resistivity = _checked_line_positive(resistivity_ohm_m, "skin-effect conductor resistivity_ohm_m")
    permeability = _checked_line_positive(permeability_h_per_m, "skin-effect conductor permeability_h_per_m")
    frequency = _checked_line_positive(frequency_hz, "skin-effect frequency_hz")
    omega = 2.0 * pi * frequency
    dc_resistance = resistivity / (pi * radius^2)
    low_frequency_reactance = omega * permeability / (8.0 * pi)
    skin_depth = sqrt(2.0 * resistivity / (omega * permeability))
    skin_ratio = radius / skin_depth
    surface_impedance = resistivity / (2.0 * pi * radius * skin_depth)
    blend = skin_ratio^4 / (1.0 + skin_ratio^4)
    low_frequency_resistance = dc_resistance * (1.0 + skin_ratio^4 / 48.0)
    low_frequency_internal_reactance = low_frequency_reactance / (1.0 + skin_ratio^4 / 48.0)
    resistance = (1.0 - blend) * low_frequency_resistance + blend * surface_impedance
    reactance = (1.0 - blend) * low_frequency_internal_reactance + blend * surface_impedance
    return ComplexF64(resistance, reactance)
end

function _line_homogeneous_earth_return_terms(e_value::Float64, angle_rad::Float64)
    if e_value <= LINE_EARTH_RETURN_SERIES_LIMIT
        r2 = e_value^4
        r1 = r2 / 16.0
        sn = e_value^2 / 8.0
        bn = r1 / 12.0
        cn = e_value / 3.0
        dn = 5.0 / 4.0
        en = e_value^3 / 45.0
        fn = 5.0 / 3.0
        a1 = 0.0
        a2 = 0.0
        a3 = 0.0
        a4 = 0.0
        b1 = 0.0
        b2 = 0.0
        b3 = 0.0
        b4 = 0.0
        evennn = 0.0
        for iteration in 1:LINE_EARTH_RETURN_SERIES_TERMS
            t = iteration - 1.0
            t1 = 2.0 * t
            t2 = 4.0 * t
            cs1 = cos((t2 + 2.0) * angle_rad)
            ss1 = sin((t2 + 2.0) * angle_rad)
            cs2 = cos((t2 + 4.0) * angle_rad)
            ss2 = sin((t2 + 4.0) * angle_rad)
            cs3 = cos((t2 + 1.0) * angle_rad)
            cs4 = cos((t2 + 3.0) * angle_rad)
            if iteration == 1
                a1 = sn * cs1
                a2 = sn * ss1
                a3 = bn * cs2
                a4 = bn * ss2
                b1 = cn * cs3
                b2 = dn * a1
                b3 = en * cs4
                b4 = fn * a3
                evennn = a1 + a2 + a3 + a4 + b1 + b2 + b3 + b4
                continue
            end
            t3 = -t1 * (t1 + 1.0)^2 * (t1 + 2.0)
            t4 = -(t1 + 1.0) * (t1 + 2.0)^2 * (t1 + 3.0)
            t5 = -(t2 - 1.0) * (t2 + 1.0)^2 * (t2 + 3.0)
            t6 = inv(t2) + inv(t1 + 1.0) + inv(t1 + 2.0) - inv(t2 + 4.0)
            t7 = -(t2 + 1.0) * (t2 + 3.0)^2 * (t2 + 5.0)
            t8 = inv(t2 + 2.0) + inv(t1 + 2.0) + inv(t1 + 3.0) -
                inv(t2 + 6.0)
            sn *= r1 / t3
            bn *= r1 / t4
            cn *= r2 / t5
            dn += t6
            en *= r2 / t7
            fn += t8
            a1 += sn * cs1
            a2 += sn * ss1
            a3 += bn * cs2
            a4 += bn * ss2
            b1 += cn * cs3
            b2 += dn * sn * cs1
            b3 += en * cs4
            b4 += fn * bn * cs2
            previous = evennn
            evennn = a1 + a2 + a3 + a4 + b1 + b2 + b3 + b4
            if evennn != 0.0 &&
               (1.0 - previous / evennn)^2 < LINE_EARTH_RETURN_SERIES_CONVERGENCE
                break
            end
        end
        euler_scale = 2.0 / LINE_EARTH_RETURN_EULER_SCALE
        sqrt_two = sqrt(2.0)
        p1 = pi * (1.0 - a3) / 4.0 + a1 * log(euler_scale / e_value) +
            angle_rad * a2 + b2 + sqrt_two * (b3 - b1)
        q1 = 0.5 + (1.0 - a3) * log(euler_scale / e_value) -
            angle_rad * a4 - pi * a1 / 4.0 - b4 + sqrt_two * (b1 + b3)
        return p1, q1
    end
    sqrt_two = sqrt(2.0)
    cs1 = sqrt_two * cos(angle_rad)
    cs2 = 2.0 * cos(2.0 * angle_rad)
    cs3 = sqrt_two * cos(3.0 * angle_rad)
    cs4 = 3.0 * sqrt_two * cos(5.0 * angle_rad)
    p1 = (cs1 + (cs4 / e_value^3 + cs3 / e_value - cs2) / e_value) / e_value
    q1 = (cs1 + (cs4 / e_value^2 - cs3) / e_value^2) / e_value
    return p1, q1
end

function cable_homogeneous_earth_return_impedance(
    image_distance_m::Real,
    image_angle_rad::Real,
    earth_resistivity_ohm_m::Real,
    frequency_hz::Real,
)
    image_distance = _checked_line_positive(image_distance_m, "earth-return image_distance_m")
    angle = Float64(image_angle_rad)
    isfinite(angle) ||
        throw(ArgumentError("earth-return image_angle_rad must be finite"))
    earth_resistivity = _checked_line_positive(earth_resistivity_ohm_m, "earth_resistivity_ohm_m")
    frequency = _checked_line_positive(frequency_hz, "earth-return frequency_hz")
    omega = 2.0 * pi * frequency
    e_value = image_distance * sqrt(LINE_VACUUM_PERMEABILITY_H_PER_M / earth_resistivity) *
        sqrt(omega)
    p1, q1 = _line_homogeneous_earth_return_terms(e_value, angle)
    scale = omega * (LINE_VACUUM_PERMEABILITY_H_PER_M / (2.0 * pi))
    return ComplexF64(scale * p1, scale * q1)
end
