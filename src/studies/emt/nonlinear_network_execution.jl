"""One prelocalized fixed-step discontinuity whose accepted interval uses a declared companion treatment."""
struct NonlinearEMTDiscontinuity
    time_s::Float64
    reason::Symbol
    event_name::Symbol
    priority::Int
    root_iteration_count::Int
    root_bracket_width_s::Float64

    function NonlinearEMTDiscontinuity(
        time_s::Real,
        reason::Symbol;
        event_name::Symbol=reason,
        priority::Integer=0,
        root_iteration_count::Integer=0,
        root_bracket_width_s::Real=0.0,
    )
        time = Float64(time_s)
        isfinite(time) && time > 0.0 || throw(ArgumentError(
            "nonlinear EMT discontinuity time must be finite and positive",
        ))
        reason in (:localized_event, :topology_change) || throw(ArgumentError(
            "nonlinear EMT discontinuity reason must be :localized_event or :topology_change",
        ))
        isempty(String(event_name)) && throw(ArgumentError(
            "nonlinear EMT discontinuity event name must not be empty",
        ))
        root_iterations = Int(root_iteration_count)
        root_iterations >= 0 || throw(ArgumentError(
            "nonlinear EMT event root-iteration count must be nonnegative",
        ))
        bracket_width = Float64(root_bracket_width_s)
        isfinite(bracket_width) && bracket_width >= 0.0 || throw(ArgumentError(
            "nonlinear EMT event root-bracket width must be finite and nonnegative",
        ))
        return new(
            time,
            reason,
            event_name,
            Int(priority),
            root_iterations,
            bracket_width,
        )
    end
end

"""Adapt one already accepted hybrid-event occurrence without replaying its transition."""
function NonlinearEMTDiscontinuity(occurrence::HybridEventOccurrence)
    return NonlinearEMTDiscontinuity(
        occurrence.time_s,
        occurrence.topology_invalidating ? :topology_change : :localized_event;
        event_name=occurrence.name,
        priority=occurrence.priority,
        root_iteration_count=occurrence.root_iteration_count,
        root_bracket_width_s=occurrence.root_bracket_width_s,
    )
end

_nonlinear_event_order_key(discontinuity::NonlinearEMTDiscontinuity) =
    (discontinuity.priority, String(discontinuity.event_name))

_nonlinear_task_order_key(occurrence::SampledTaskOccurrence) =
    (occurrence.priority, String(occurrence.name), occurrence.execution_index)

"""Exact fixed-step calendar for a nonlinear EMT network and its prelocalized discontinuities."""
struct NonlinearEMTStudySchedule
    step_s::Float64
    final_time_s::Float64
    step_count::Int
    discontinuity_reason_by_step::Dict{Int,Symbol}
    discontinuities_by_step::Dict{Int,Vector{NonlinearEMTDiscontinuity}}
    accepted_tasks_by_step::Dict{Int,Vector{SampledTaskOccurrence}}

    function NonlinearEMTStudySchedule(
        step_s::Real,
        final_time_s::Real;
        discontinuities=NonlinearEMTDiscontinuity[],
        accepted_tasks=SampledTaskOccurrence[],
    )
        step = Float64(step_s)
        final_time = Float64(final_time_s)
        isfinite(step) && step > 0.0 || throw(ArgumentError(
            "nonlinear EMT schedule step must be finite and positive",
        ))
        isfinite(final_time) && final_time >= step || throw(ArgumentError(
            "nonlinear EMT final time must be finite and at least one step",
        ))
        count = round(Int, final_time / step)
        abs(count * step - final_time) <= 16.0 * eps(Float64) * final_time ||
            throw(ArgumentError(
                "nonlinear EMT final time must lie on the fixed-step calendar",
            ))
        reason_by_step = Dict{Int,Symbol}()
        discontinuities_by_step = Dict{Int,Vector{NonlinearEMTDiscontinuity}}()
        for discontinuity in discontinuities
            discontinuity isa NonlinearEMTDiscontinuity || throw(ArgumentError(
                "nonlinear EMT schedule requires typed discontinuities",
            ))
            discontinuity_step = round(Int, discontinuity.time_s / step)
            1 <= discontinuity_step <= count || throw(ArgumentError(
                "nonlinear EMT discontinuity lies outside the scheduled horizon",
            ))
            abs(discontinuity_step * step - discontinuity.time_s) <=
                16.0 * eps(Float64) * max(discontinuity.time_s, step) ||
                throw(ArgumentError(
                    "nonlinear EMT discontinuity must be prelocalized to the fixed-step calendar",
                ))
            push!(
                get!(
                    discontinuities_by_step,
                    discontinuity_step,
                    NonlinearEMTDiscontinuity[],
                ),
                discontinuity,
            )
        end
        for (discontinuity_step, boundary_events) in discontinuities_by_step
            event_names = getproperty.(boundary_events, :event_name)
            length(unique(event_names)) == length(event_names) || throw(ArgumentError(
                "nonlinear EMT schedule repeats an event name at one boundary",
            ))
            order_keys = _nonlinear_event_order_key.(boundary_events)
            issorted(order_keys) || throw(ArgumentError(
                "nonlinear EMT simultaneous events must retain accepted priority/name order",
            ))
            reason_by_step[discontinuity_step] = any(
                event -> event.reason === :topology_change,
                boundary_events,
            ) ? :topology_change : :localized_event
        end
        tasks_by_step = Dict{Int,Vector{SampledTaskOccurrence}}()
        for occurrence in accepted_tasks
            occurrence isa SampledTaskOccurrence || throw(ArgumentError(
                "nonlinear EMT schedule requires typed accepted sampled-task occurrences",
            ))
            isfinite(occurrence.time_s) && occurrence.time_s > 0.0 || throw(ArgumentError(
                "nonlinear EMT accepted sampled-task time must be finite and positive",
            ))
            occurrence.tick >= 0 || throw(ArgumentError(
                "nonlinear EMT accepted sampled-task tick must be nonnegative",
            ))
            occurrence.execution_index > 0 || throw(ArgumentError(
                "nonlinear EMT accepted sampled-task execution index must be positive",
            ))
            task_step = round(Int, occurrence.time_s / step)
            1 <= task_step <= count || throw(ArgumentError(
                "nonlinear EMT accepted sampled task lies outside the scheduled horizon",
            ))
            abs(task_step * step - occurrence.time_s) <=
                16.0 * eps(Float64) * max(occurrence.time_s, step) ||
                throw(ArgumentError(
                    "nonlinear EMT accepted sampled task must lie on the fixed-step calendar",
                ))
            push!(
                get!(tasks_by_step, task_step, SampledTaskOccurrence[]),
                occurrence,
            )
        end
        for boundary_tasks in values(tasks_by_step)
            task_names = getproperty.(boundary_tasks, :name)
            length(unique(task_names)) == length(task_names) || throw(ArgumentError(
                "nonlinear EMT schedule repeats a sampled-task name at one boundary",
            ))
            order_keys = _nonlinear_task_order_key.(boundary_tasks)
            issorted(order_keys) || throw(ArgumentError(
                "nonlinear EMT simultaneous sampled tasks must retain accepted priority/name order",
            ))
        end
        return new(
            step,
            final_time,
            count,
            reason_by_step,
            discontinuities_by_step,
            tasks_by_step,
        )
    end
end

"""Accepted nonlinear EMT voltage, ideal-constraint current, and numerical diagnostics on an exact schedule."""
struct NonlinearEMTStudyTrace
    time_s::Vector{Float64}
    voltage_v::Matrix{Float64}
    constraint_current_a::Matrix{Float64}
    diagnostics::Vector{NonlinearSolveDiagnostics}
    accepted_discontinuities::Vector{NonlinearEMTDiscontinuity}
    accepted_tasks::Vector{SampledTaskOccurrence}
end

"""Execute a typed nonlinear EMT schedule without duplicating accepted steps, histories, events, or outputs."""
function evaluate_nonlinear_emt_network!(
    system::NonlinearNodalSystem,
    schedule::NonlinearEMTStudySchedule;
    chatter_decisions::AbstractDict{Int,<:NonlinearChatterDecision}=
        Dict{Int,NonlinearChatterDecision}(),
)
    all(step -> 1 <= step <= schedule.step_count, keys(chatter_decisions)) ||
        throw(ArgumentError("nonlinear chatter decision step lies outside the schedule"))
    node_count = nonlinear_linear_system(system).node_count
    constraint_count = length(system.ideal_constraints)
    time_s = collect(0:schedule.step_count) .* schedule.step_s
    voltage_v = zeros(Float64, node_count, schedule.step_count + 1)
    constraint_current_a = zeros(Float64, constraint_count, schedule.step_count + 1)
    voltage_v[:, 1] .= nonlinear_linear_system(system).v
    constraint_current_a[:, 1] .= system.accepted_state.constraint_current_a
    diagnostics = NonlinearSolveDiagnostics[]
    accepted_discontinuities = NonlinearEMTDiscontinuity[]
    accepted_tasks = SampledTaskOccurrence[]
    sizehint!(diagnostics, schedule.step_count)
    for step in 1:schedule.step_count
        discontinuity_reason = get(
            schedule.discontinuity_reason_by_step,
            step,
            :none,
        )
        discontinuity_treatment = discontinuity_reason === :none ?
            :none : :two_backward_euler_half_steps
        chatter_decision = get(chatter_decisions, step, nothing)
        boundary_events = get(
            schedule.discontinuities_by_step,
            step,
            NonlinearEMTDiscontinuity[],
        )
        boundary_tasks = get(
            schedule.accepted_tasks_by_step,
            step,
            SampledTaskOccurrence[],
        )
        if chatter_decision !== nothing && chatter_decision.critical_damping_allowed &&
           (!isempty(boundary_events) || !isempty(boundary_tasks))
            throw(ArgumentError(
                "nonlinear chatter damping cannot cross an accepted event or sampled-task boundary",
            ))
        end
        result = advance_nonlinear_step!(
            system,
            time_s[step + 1],
            schedule.step_s;
            discontinuity_treatment,
            discontinuity_reason,
            chatter_decision,
        )
        result.accepted || throw(result.failure)
        voltage_v[:, step + 1] .= result.voltage_v
        constraint_count > 0 &&
            (constraint_current_a[:, step + 1] .= result.constraint_current_a)
        push!(diagnostics, result.diagnostics)
        append!(accepted_discontinuities, boundary_events)
        append!(accepted_tasks, boundary_tasks)
    end
    return NonlinearEMTStudyTrace(
        time_s,
        voltage_v,
        constraint_current_a,
        diagnostics,
        accepted_discontinuities,
        accepted_tasks,
    )
end
