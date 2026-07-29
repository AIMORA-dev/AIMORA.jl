const DECK_CASE_SEQUENCE_SUMMARY_SCHEMA = "aimora.deck_case_sequence_summary.v1"

"""
The physical frequency points implied by an FS or LMFS request.

Linear schedules expose a physical step in hertz. Logarithmic schedules expose
the dimensionless geometric ratio without retaining overloaded compatibility
state.
"""
struct DeckFrequencyScanSchedule
    request_line_no::Int
    request_kind::Symbol
    minimum_frequency_hz::Float64
    maximum_frequency_hz::Float64
    points_per_decade::Int
    linear_step_hz::Union{Nothing,Float64}
    geometric_ratio::Float64
    frequencies_hz::Vector{Float64}
end

struct DeckNetworkFrequencyScanPoint
    frequency_hz::Float64
    node_voltage_phasors::Vector{ComplexF64}
    nodal_current_phasors::Vector{ComplexF64}
    source_injection_phasors::Vector{ComplexF64}
    output_voltage_phasors::Vector{ComplexF64}
    residual_max_abs::Float64
    relative_residual_max_abs::Float64
    admittance_symmetry_max_abs_error::Float64
    minimum_dissipative_eigenvalue_s::Float64
    physical_checks_passed::Bool
end

struct DeckNetworkFrequencyScanStudy
    request_line_no::Int
    node_names::Vector{Symbol}
    output_node_indices::Vector{Int}
    frequencies_hz::Vector{Float64}
    points::Vector{DeckNetworkFrequencyScanPoint}
    physical_checks_passed::Bool
end

struct DeckImpulseResponseFitControl
    fit_span_s::Float64
    time_step_s::Float64

    function DeckImpulseResponseFitControl(
        fit_span_s::Real,
        time_step_s::Real,
    )
        span = Float64(fit_span_s)
        step = Float64(time_step_s)
        isfinite(span) && span > 0.0 ||
            throw(ArgumentError("impulse-response fit span must be finite and positive"))
        isfinite(step) && step > 0.0 && step < span ||
            throw(ArgumentError("impulse-response timestep must be positive and shorter than the fit span"))
        return new(span, step)
    end
end

struct DeckLineFrequencyScanStudy
    line_name::Symbol
    mode_index::Int
    frequencies_hz::Vector{Float64}
    characteristic_impedance_ohm::Vector{ComplexF64}
    propagation_without_delay::Vector{ComplexF64}
    physical_checks_passed::Bool
end

struct DeckLineImpulseResponseStudy
    line_name::Symbol
    mode_index::Int
    frequencies_hz::Vector{Float64}
    propagation_fit::LineStepResponseExponentialFitResult
    characteristic_admittance_fit::Union{Nothing,LineStepResponseExponentialFitResult}
    characteristic_admittance_is_constant::Bool
    propagation_dc_value::Float64
    characteristic_admittance_dc_siemens::Float64
    physical_checks_passed::Bool
end

struct DeckAuxiliaryStudyRun
    schedules::Vector{DeckFrequencyScanSchedule}
    network_frequency_scans::Vector{DeckNetworkFrequencyScanStudy}
    line_frequency_scans::Vector{DeckLineFrequencyScanStudy}
    impulse_responses::Vector{DeckLineImpulseResponseStudy}
    line_constants::Union{Nothing,LineConstantsStudyResult}
    handled_requests::Vector{Symbol}
    compatibility_exclusions::Vector{Symbol}
    deferred_requests::Vector{Symbol}
end

struct DeckCaseRun
    case_index::Int
    start_line_no::Int
    end_line_no::Int
    execution_kind::Symbol
    trace::Union{Nothing,DeckEMTTrace}
    auxiliary::DeckAuxiliaryStudyRun
end

struct DeckCaseSequenceRun
    source::String
    cases::Vector{DeckCaseRun}
    aborted_case_count::Int
    discarded_card_count::Int
    run_terminated::Bool
end

function deck_frequency_scan_schedule(row::DeckParser.DeckStudyOptionRequestRow)
    row.request_kind in (:frequency_scan, :line_model_frequency_scan) ||
        throw(ArgumentError("request $(row.request_kind) is not a frequency scan"))
    length(row.numeric_values) == 3 && length(row.integer_values) == 1 ||
        throw(ArgumentError("frequency scan requires FMIN, DELF, FMAX, and N8"))
    fmin, requested_delta, fmax = row.numeric_values
    points_per_decade = only(row.integer_values)
    isfinite(fmin) && fmin > 0.0 ||
        throw(ArgumentError("minimum scan frequency must be finite and positive"))
    isfinite(fmax) && fmax > fmin ||
        throw(ArgumentError("maximum scan frequency must exceed the minimum"))
    points_per_decade >= 0 ||
        throw(ArgumentError("frequency-scan points per decade must be nonnegative"))

    frequencies = Float64[]
    if points_per_decade == 0
        isfinite(requested_delta) && requested_delta > 0.0 ||
            throw(ArgumentError("linear scan delta must be finite and positive"))
        step_count = floor(
            Int,
            (fmax - fmin) / requested_delta +
            32.0 * eps(Float64) * max(abs(fmax / requested_delta), 1.0),
        )
        frequencies = [fmin + index * requested_delta for index in 0:step_count]
        linear_step = requested_delta
        ratio = 1.0
    else
        ratio = exp(log(10.0) / points_per_decade)
        step_count = floor(
            Int,
            log(fmax / fmin) / log(ratio) +
            32.0 * eps(Float64) * max(abs(log(fmax / fmin) / log(ratio)), 1.0),
        )
        frequencies = [fmin * ratio^index for index in 0:step_count]
        linear_step = nothing
    end
    if last(frequencies) > fmax &&
       last(frequencies) <= fmax * (1.0 + 64.0 * eps(Float64))
        frequencies[end] = fmax
    end
    all(isfinite, frequencies) && all(>(0.0), frequencies) ||
        throw(ArgumentError("frequency scan generated invalid samples"))
    issorted(frequencies) && all(diff(frequencies) .> 0.0) ||
        throw(ArgumentError("frequency scan must be strictly increasing"))
    last(frequencies) <= fmax * (1.0 + 64.0 * eps(Float64)) ||
        throw(ArgumentError("frequency scan exceeds its upper bound"))
    return DeckFrequencyScanSchedule(
        row.line_no,
        row.request_kind,
        fmin,
        fmax,
        points_per_decade,
        linear_step,
        ratio,
        frequencies,
    )
end

function _impulse_response_schedule(schedules::Vector{DeckFrequencyScanSchedule})
    candidates = filter(
        schedule ->
            schedule.points_per_decade > 0 && length(schedule.frequencies_hz) >= 4,
        schedules,
    )
    isempty(candidates) && return nothing
    return first(candidates)
end

function _network_frequency_scan_partition(
    parsed::DeckParser.DeckParseResult,
    frequency_hz::Float64,
)
    base = DeckParser.deck_steady_state_frequency_partition(parsed)
    source_frequencies = copy(base.source_frequencies_hz)
    for source_index in base.active_source_row_indices
        source_frequencies[source_index] = frequency_hz
    end
    return DeckParser.DeckSteadyStateFrequencyPartition(
        copy(base.node_source_row_indices),
        fill(frequency_hz, length(base.node_frequencies_hz)),
        copy(base.active_source_row_indices),
        source_frequencies,
        copy(base.source_successor_indices),
        copy.(base.source_groups),
        fill(frequency_hz, length(base.subnetwork_frequencies_hz)),
        copy.(base.subnetwork_node_indices),
        copy(base.inactive_node_indices),
        base.scalar_branch_count,
        base.initially_closed_switch_count,
        copy(base.unsupported_topology_kinds),
    )
end

function _network_frequency_scan_point(
    parsed::DeckParser.DeckParseResult,
    frequency_hz::Float64,
    output_node_indices::Vector{Int},
)
    node_count = maximum(values(parsed.node_map); init = 0)
    node_count > 0 ||
        throw(ArgumentError("whole-network frequency scan requires at least one network node"))
    partition = _network_frequency_scan_partition(parsed, frequency_hz)
    admittance, rhs, switch_representatives =
        _deck_steady_state_nodal_equations(
            parsed,
            node_count;
            frequency_partition = partition,
        )
    voltage = _solve_grouped_steady_state_admittance(
        admittance,
        rhs,
        switch_representatives,
    )
    nodal_current = admittance * voltage
    residual = nodal_current - rhs
    residual_max = maximum(abs, residual; init = 0.0)
    residual_scale = max(
        norm(rhs, Inf),
        norm(admittance, Inf) * norm(voltage, Inf),
        1.0,
    )
    relative_residual = residual_max / residual_scale
    symmetry_error = maximum(
        abs,
        admittance - transpose(admittance);
        init = 0.0,
    )
    dissipative = Hermitian(0.5 .* (admittance + adjoint(admittance)))
    minimum_dissipative_eigenvalue =
        minimum(eigvals(dissipative); init = 0.0)
    output_voltage = ComplexF64[
        node == 0 ? 0.0 + 0.0im : voltage[node]
        for node in output_node_indices
    ]
    finite_values =
        all(isfinite, real.(voltage)) &&
        all(isfinite, imag.(voltage)) &&
        all(isfinite, real.(nodal_current)) &&
        all(isfinite, imag.(nodal_current))
    physical_checks =
        finite_values &&
        relative_residual <= 1.0e-10 &&
        symmetry_error <= 1.0e-10
    physical_checks || throw(ArgumentError(
        "whole-network frequency scan failed its nodal residual or reciprocity check at $frequency_hz Hz",
    ))
    return DeckNetworkFrequencyScanPoint(
        frequency_hz,
        voltage,
        nodal_current,
        rhs,
        output_voltage,
        residual_max,
        relative_residual,
        symmetry_error,
        minimum_dissipative_eigenvalue,
        physical_checks,
    )
end

function _network_frequency_scan_study(
    parsed::DeckParser.DeckParseResult,
    schedule::DeckFrequencyScanSchedule,
)
    schedule.request_kind == :frequency_scan ||
        throw(ArgumentError("whole-network scan requires a FREQUENCY SCAN schedule"))
    output_node_indices = DeckParser.deck_over16_output_node_indices(parsed)
    points = DeckNetworkFrequencyScanPoint[
        _network_frequency_scan_point(
            parsed,
            frequency_hz,
            output_node_indices,
        )
        for frequency_hz in schedule.frequencies_hz
    ]
    checks =
        length(points) == length(schedule.frequencies_hz) &&
        all(point -> point.physical_checks_passed, points)
    checks ||
        throw(ArgumentError("whole-network frequency scan did not execute every requested point"))
    return DeckNetworkFrequencyScanStudy(
        schedule.request_line_no,
        ordered_node_names(parsed.node_map),
        output_node_indices,
        copy(schedule.frequencies_hz),
        points,
        checks,
    )
end

function _line_impulse_response_study(group, mode, schedule, fit_control)
    frequencies_hz = schedule.frequencies_hz
    angular_frequencies = 2.0 * pi .* vcat(0.0, frequencies_hz)
    propagation = mode.conversion.propagation_without_delay
    impedance = mode.conversion.characteristic_impedance
    propagation_values = Float64[
        real(pole_residue_transfer_value(propagation, im * omega))
        for omega in angular_frequencies
    ]
    admittance_values = Float64[
        real(inv(pole_residue_transfer_value(impedance, im * omega)))
        for omega in angular_frequencies
    ]
    propagation_dc = first(propagation_values)
    admittance_dc = first(admittance_values)
    propagation_dc > 0.0 && admittance_dc > 0.0 ||
        throw(ArgumentError("line impulse response requires positive DC transfer values"))
    propagation_fit = line_step_response_exponential_fit(
        angular_frequencies,
        propagation_values;
        final_value = propagation_dc,
        fit_span_s=fit_control.fit_span_s,
        time_step_s=fit_control.time_step_s,
    )
    admittance_is_constant = all(
        value -> isapprox(value, admittance_dc; rtol=1.0e-12, atol=1.0e-15),
        admittance_values,
    )
    admittance_fit = admittance_is_constant ?
        nothing :
        line_step_response_exponential_fit(
            angular_frequencies,
            admittance_values;
            final_value = admittance_dc,
            fit_span_s=fit_control.fit_span_s,
            time_step_s=fit_control.time_step_s,
            fit_control_mode = -1,
        )
    physical_checks_passed =
        mode.conversion.physical_checks_passed &&
        propagation_fit.fit_executed &&
        (admittance_is_constant || admittance_fit.fit_executed) &&
        isfinite(propagation_fit.normalized_square_error) &&
        (admittance_is_constant || isfinite(admittance_fit.normalized_square_error))
    physical_checks_passed ||
        throw(ArgumentError("line impulse-response physical checks did not pass"))
    return DeckLineImpulseResponseStudy(
        group.name,
        mode.mode_index,
        copy(frequencies_hz),
        propagation_fit,
        admittance_fit,
        admittance_is_constant,
        propagation_dc,
        admittance_dc,
        physical_checks_passed,
    )
end

function _line_frequency_scan_study(group, mode, schedule)
    frequencies_hz = schedule.frequencies_hz
    impedance = ComplexF64[
        pole_residue_transfer_value(
            mode.conversion.characteristic_impedance,
            2.0im * pi * frequency_hz,
        )
        for frequency_hz in frequencies_hz
    ]
    propagation = ComplexF64[
        pole_residue_transfer_value(
            mode.conversion.propagation_without_delay,
            2.0im * pi * frequency_hz,
        )
        for frequency_hz in frequencies_hz
    ]
    physical_checks_passed =
        mode.conversion.physical_checks_passed &&
        all(value -> isfinite(real(value)) && isfinite(imag(value)), impedance) &&
        all(value -> real(value) > 0.0, impedance) &&
        all(value -> isfinite(real(value)) && isfinite(imag(value)), propagation) &&
        all(value -> abs(value) <= 1.0 + 1.0e-10, propagation)
    physical_checks_passed ||
        throw(ArgumentError("line-model frequency-scan physical checks did not pass"))
    return DeckLineFrequencyScanStudy(
        group.name,
        mode.mode_index,
        copy(frequencies_hz),
        impedance,
        propagation,
        physical_checks_passed,
    )
end

function run_deck_auxiliary_studies(
    parsed::DeckParser.DeckParseResult;
    impulse_response_fit_control::Union{
        Nothing,
        DeckImpulseResponseFitControl,
    }=nothing,
)
    assert_deck_valid!(parsed)
    rows = DeckParser.deck_study_option_request_rows(parsed)
    simulation_control_rows =
        DeckParser.deck_simulation_control_request_rows(parsed)
    schedules = DeckFrequencyScanSchedule[
        deck_frequency_scan_schedule(row)
        for row in rows
        if row.request_kind in (:frequency_scan, :line_model_frequency_scan)
    ]
    requested_kinds = Set(getfield.(rows, :request_kind))
    simulation_control_kinds =
        Set(getfield.(simulation_control_rows, :request_kind))
    groups = DeckParser.deck_rational_frequency_line_groups(parsed)
    network_frequency_scans = DeckNetworkFrequencyScanStudy[]
    line_frequency_scans = DeckLineFrequencyScanStudy[]
    impulse_responses = DeckLineImpulseResponseStudy[]
    handled_requests = Symbol[]
    compatibility_exclusions = Symbol[]
    deferred_requests = Symbol[]
    line_constants = if !isempty(
        DeckParser.deck_line_constants_physical_conductors(parsed),
    ) || !isempty(DeckParser.deck_line_constants_frequency_cards(parsed))
        run_line_constants_study(parsed)
    else
        nothing
    end

    line_model_schedules = filter(
        schedule -> schedule.request_kind == :line_model_frequency_scan,
        schedules,
    )
    if !isempty(line_model_schedules)
        if isempty(groups)
            push!(deferred_requests, :line_model_frequency_scan_model)
        else
            for schedule in line_model_schedules, group in groups, mode in group.modes
                push!(
                    line_frequency_scans,
                    _line_frequency_scan_study(group, mode, schedule),
                )
            end
        end
    end

    for schedule in schedules
        schedule.request_kind == :frequency_scan || continue
        push!(
            network_frequency_scans,
            _network_frequency_scan_study(parsed, schedule),
        )
        push!(handled_requests, :whole_network_frequency_scan)
    end

    if :hauer_impulse_response_setup in requested_kinds
        schedule = _impulse_response_schedule(schedules)
        if schedule === nothing
            push!(deferred_requests, :impulse_response_logarithmic_scan)
        elseif isempty(groups)
            push!(deferred_requests, :impulse_response_line_model)
        elseif impulse_response_fit_control === nothing
            push!(deferred_requests, :impulse_response_fit_control)
        else
            for group in groups, mode in group.modes
                push!(
                    impulse_responses,
                    _line_impulse_response_study(
                        group,
                        mode,
                        schedule,
                        impulse_response_fit_control,
                    ),
                )
            end
        end
    end
    # MIDOV1 only reset scratch/output units in the monolithic executable.
    # Each Julia study already owns isolated output streams and typed archives,
    # so the request is satisfied without emulating global file units.
    :file_request in requested_kinds &&
        push!(handled_requests, :isolated_output_file_lifecycle)
    # CONVERT ZNO rewrote pre-M37 card layouts and was never a simulation
    # model. Keep it visible as an explicit retired input-conversion utility,
    # not as a silently deferred production side effect.
    :zinc_oxide_format_conversion in requested_kinds &&
        push!(compatibility_exclusions, :pre_m37_zinc_oxide_card_converter)
    # CHANGE SWITCH (special request 37) invoked OVER41/CRDCHG only to
    # rewrite pre-M37 type-91/92/93 input cards into newer formats. It
    # never participated in a simulation. Keep the accepted request visible
    # as a retired input-conversion utility instead of treating CRDCHG as
    # transformer-parameter behavior or a deferred runtime effect.
    :switch_pseudononlinear_conversion in simulation_control_kinds &&
        push!(
            compatibility_exclusions,
            :pre_m37_switch_pseudononlinear_card_converter,
        )
    return DeckAuxiliaryStudyRun(
        schedules,
        network_frequency_scans,
        line_frequency_scans,
        impulse_responses,
        line_constants,
        unique(handled_requests),
        unique(compatibility_exclusions),
        unique(deferred_requests),
    )
end

function _case_execution_kind(
    parsed::DeckParser.DeckParseResult,
    auxiliary::DeckAuxiliaryStudyRun,
)
    has_auxiliary_request = !isempty(
        DeckParser.deck_study_option_request_rows(parsed),
    )
    if auxiliary.line_constants !== nothing
        return :line_constants
    elseif !isempty(parsed.elements)
        return has_auxiliary_request ? :emt_with_auxiliary : :emt
    elseif has_auxiliary_request
        return :auxiliary
    elseif isempty(parsed.node_map)
        return :metadata_only
    end
    throw(ArgumentError(
        "case $(parsed.source) has model nodes but no ordinary-network or auxiliary execution owner",
    ))
end

function run_deck_case_sequence_emt(
    sequence::DeckParser.DeckCaseSequence;
    dt_s::Float64=20e-6,
    t_end_s::Float64=0.0,
    time_horizon::Symbol=:arguments,
    output_schedule::Symbol=:all_steps,
    impulse_response_fit_control::Union{
        Nothing,
        DeckImpulseResponseFitControl,
    }=nothing,
)
    runs = DeckCaseRun[]
    for data_case in sequence.cases
        parsed = data_case.parsed
        assert_deck_valid!(parsed)
        auxiliary = run_deck_auxiliary_studies(
            parsed;
            impulse_response_fit_control,
        )
        execution_kind = _case_execution_kind(parsed, auxiliary)
        trace = execution_kind in (:emt, :emt_with_auxiliary) ?
            run_deck_emt(
                parsed;
                dt_s,
                t_end_s,
                time_horizon,
                output_schedule,
            ) :
            nothing
        push!(
            runs,
            DeckCaseRun(
                data_case.case_index,
                data_case.start_line_no,
                data_case.end_line_no,
                execution_kind,
                trace,
                auxiliary,
            ),
        )
    end
    return DeckCaseSequenceRun(
        sequence.source,
        runs,
        sequence.aborted_case_count,
        sequence.discarded_card_count,
        sequence.run_terminated,
    )
end

function deck_case_sequence_result(run::DeckCaseSequenceRun; elapsed_s::Float64=0.0)
    executed_count = count(data_case -> data_case.trace !== nothing, run.cases)
    schedule_count = sum(
        data_case -> length(data_case.auxiliary.schedules),
        run.cases;
        init=0,
    )
    line_frequency_scan_count = sum(
        data_case -> length(data_case.auxiliary.line_frequency_scans),
        run.cases;
        init=0,
    )
    network_frequency_scan_count = sum(
        data_case -> length(data_case.auxiliary.network_frequency_scans),
        run.cases;
        init=0,
    )
    impulse_response_count = sum(
        data_case -> length(data_case.auxiliary.impulse_responses),
        run.cases;
        init=0,
    )
    line_constants_count = count(
        data_case -> data_case.auxiliary.line_constants !== nothing,
        run.cases,
    )
    deferred = unique(vcat(
        [data_case.auxiliary.deferred_requests for data_case in run.cases]...,
        Symbol[],
    ))
    exclusions = unique(vcat(
        [data_case.auxiliary.compatibility_exclusions for data_case in run.cases]...,
        Symbol[],
    ))
    warnings = StudyWarning[
        study_warning(
            request,
            "The deck sequence retained an explicit deferred auxiliary request: $(request).",
        )
        for request in deferred
    ]
    append!(
        warnings,
        StudyWarning[
            study_warning(
                exclusion,
                "The request is an explicitly retired input-conversion utility, not a production simulation effect: $(exclusion).",
            )
            for exclusion in exclusions
        ],
    )
    return study_result(
        :emt;
        status=isempty(warnings) ? :ok : :warning,
        quantities=[
            result_quantity(:case_count, length(run.cases); unit="count"),
            result_quantity(:executed_emt_case_count, executed_count; unit="count"),
            result_quantity(:aborted_case_count, run.aborted_case_count; unit="count"),
            result_quantity(:frequency_scan_count, schedule_count; unit="count"),
            result_quantity(
                :whole_network_frequency_scan_count,
                network_frequency_scan_count;
                unit="count",
            ),
            result_quantity(
                :line_model_frequency_scan_count,
                line_frequency_scan_count;
                unit="count",
            ),
            result_quantity(
                :line_impulse_response_count,
                impulse_response_count;
                unit="count",
            ),
            result_quantity(
                :line_constants_case_count,
                line_constants_count;
                unit="count",
            ),
            result_quantity(:elapsed_s, elapsed_s; unit="s"),
        ],
        assumptions=[
            study_assumption(
                :case_state_ownership,
                "isolated";
                description="Every accepted data case owns a separate parsed model and runtime state.",
            ),
            study_assumption(
                :external_reference_in_runtime,
                false;
                description="The external reference executable is comparison evidence only.",
            ),
        ],
        warnings,
        metadata=Dict{Symbol,Any}(
            :source => run.source,
            :run_terminated => run.run_terminated,
            :discarded_card_count => run.discarded_card_count,
            :deferred_requests => deferred,
            :compatibility_exclusions => exclusions,
        ),
    )
end

_toml_text(value::AbstractString) = sprint(show, String(value))

function write_deck_case_sequence_summary(
    path::AbstractString,
    run::DeckCaseSequenceRun;
    elapsed_s::Float64=0.0,
)
    result = deck_case_sequence_result(run; elapsed_s)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "schema = ", _toml_text(DECK_CASE_SEQUENCE_SUMMARY_SCHEMA))
        println(io, "source = ", _toml_text(run.source))
        println(io, "status = ", _toml_text(String(result.status)))
        println(io, "case_count = ", length(run.cases))
        println(io, "aborted_case_count = ", run.aborted_case_count)
        println(io, "discarded_card_count = ", run.discarded_card_count)
        println(io, "run_terminated = ", run.run_terminated)
        println(io, "elapsed_s = ", elapsed_s)
        for data_case in run.cases
            println(io)
            println(io, "[[cases]]")
            println(io, "case_index = ", data_case.case_index)
            println(io, "start_line_no = ", data_case.start_line_no)
            println(io, "end_line_no = ", data_case.end_line_no)
            println(
                io,
                "execution_kind = ",
                _toml_text(String(data_case.execution_kind)),
            )
            println(io, "emt_executed = ", data_case.trace !== nothing)
            println(
                io,
                "frequency_scan_count = ",
                length(data_case.auxiliary.schedules),
            )
            println(
                io,
                "whole_network_frequency_scan_count = ",
                length(data_case.auxiliary.network_frequency_scans),
            )
            println(
                io,
                "whole_network_frequency_sample_count = ",
                sum(
                    scan -> length(scan.points),
                    data_case.auxiliary.network_frequency_scans;
                    init=0,
                ),
            )
            println(
                io,
                "line_model_frequency_scan_count = ",
                length(data_case.auxiliary.line_frequency_scans),
            )
            println(
                io,
                "frequency_sample_count = ",
                sum(
                    scan -> length(scan.frequencies_hz),
                    data_case.auxiliary.line_frequency_scans;
                    init=0,
                ),
            )
            println(
                io,
                "frequency_schedules_hz = [",
                join(
                    [
                        "[" * join(string.(schedule.frequencies_hz), ", ") * "]"
                        for schedule in data_case.auxiliary.schedules
                    ],
                    ", ",
                ),
                "]",
            )
            println(
                io,
                "impulse_response_count = ",
                length(data_case.auxiliary.impulse_responses),
            )
            println(
                io,
                "impulse_response_fit_errors = [",
                join(
                    string.(
                        getfield.(
                            getfield.(
                                data_case.auxiliary.impulse_responses,
                                :propagation_fit,
                            ),
                            :normalized_square_error,
                        ),
                    ),
                    ", ",
                ),
                "]",
            )
            println(
                io,
                "handled_requests = [",
                join(
                    _toml_text.(String.(data_case.auxiliary.handled_requests)),
                    ", ",
                ),
                "]",
            )
            println(
                io,
                "compatibility_exclusions = [",
                join(
                    _toml_text.(
                        String.(data_case.auxiliary.compatibility_exclusions),
                    ),
                    ", ",
                ),
                "]",
            )
            println(
                io,
                "deferred_requests = [",
                join(
                    _toml_text.(String.(data_case.auxiliary.deferred_requests)),
                    ", ",
                ),
                "]",
            )
            if data_case.trace !== nothing
                println(
                    io,
                    "final_voltage_pu = [",
                    join(string.(data_case.trace.voltage_pu[:, end]), ", "),
                    "]",
                )
            end
        end
    end
    return path
end
