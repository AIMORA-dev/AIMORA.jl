struct DCSimulatorInitializationResult
    primary_rows::Vector{Int}
    voltage_differences::Vector{Float64}
    control_values::Vector{Float64}
end

function initialize_dc_simulator_sources!(
    plan::DeckOVER16BoundaryPlan,
    initial_node_voltages::AbstractVector{<:Real},
)
    primary_rows = findall(==(:dc_simulator_primary), plan.source_layout_kinds)
    voltage_differences = Float64[]
    control_values = Float64[]
    for primary_row in primary_rows
        successor_row = primary_row + 1
        balancing_row = primary_row + 2
        balancing_row <= plan.source_row_count || throw(ArgumentError(
            "DC-simulator source topology requires controller, successor, and balancing rows",
        ))
        plan.source_layout_kinds[successor_row] == :dc_simulator_successor ||
            throw(ArgumentError("DC-simulator controller row is missing its typed successor"))
        plan.source_layout_kinds[balancing_row] == :dc_simulator_balancing_source ||
            throw(ArgumentError("DC-simulator controller row is missing its balancing source"))

        primary_node = abs(plan.source_node_values[primary_row])
        successor_node = abs(plan.source_node_values[successor_row])
        maximum((primary_node, successor_node)) <= length(initial_node_voltages) ||
            throw(ArgumentError("DC-simulator nodes exceed the initialized voltage state"))
        voltage_difference =
            Float64(initial_node_voltages[primary_node]) -
            Float64(initial_node_voltages[successor_node])

        old_successor_stop = plan.source_tstop_values[successor_row]
        plan.source_tstop_values[successor_row] = 0.0
        intermediate_control =
            plan.source_crest_values[balancing_row] -
            voltage_difference * plan.source_tstop_values[primary_row]
        selected_control = plan.source_iform_values[successor_row] > 1 ?
            old_successor_stop : intermediate_control
        voltage_scaled_control =
            voltage_difference * plan.source_tstop_values[primary_row] -
            plan.source_crest_values[primary_row]
        updated_control =
            voltage_scaled_control / plan.source_time1_values[primary_row] +
            selected_control * plan.source_sfreq_values[primary_row] /
            plan.source_crest_values[successor_row]
        lower_state =
            plan.source_time1_values[successor_row] *
            plan.source_tstop_values[primary_row] -
            plan.source_crest_values[primary_row]
        upper_state =
            plan.source_tstart_values[successor_row] *
            plan.source_tstop_values[primary_row] -
            plan.source_crest_values[primary_row]

        plan.source_crest_values[primary_row] += 2.0 * updated_control
        plan.source_sfreq_values[primary_row] =
            1.0 - 2.0 * plan.source_sfreq_values[primary_row] /
            plan.source_crest_values[successor_row]
        plan.source_time1_values[primary_row] =
            1.0 - 2.0 / plan.source_time1_values[primary_row]
        plan.source_crest_values[successor_row] = 2.0 * updated_control
        plan.source_sfreq_values[successor_row] = voltage_scaled_control
        plan.source_time1_values[successor_row] = lower_state
        plan.source_time2_values[successor_row] = upper_state
        plan.source_tstart_values[successor_row] =
            plan.source_sfreq_values[primary_row] * intermediate_control +
            plan.source_time1_values[primary_row] * voltage_scaled_control

        push!(voltage_differences, voltage_difference)
        push!(control_values, updated_control)
    end
    return DCSimulatorInitializationResult(
        primary_rows,
        voltage_differences,
        control_values,
    )
end
