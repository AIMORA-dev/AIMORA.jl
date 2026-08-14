export TransformerApparatusRuntime,
       TransformerApparatusRuntimeSnapshot,
       transformer_apparatus_runtime,
       transformer_apparatus_runtime_snapshot,
       restore_transformer_apparatus_runtime_snapshot!,
       transformer_apparatus_runtime_diagnostics

abstract type AbstractTransformerAcceptedRuntimeState end

mutable struct TransformerTerminalMatrixRuntimeState <: AbstractTransformerAcceptedRuntimeState
    coil_current_a::Vector{Float64}
    capacitor_current_a::Vector{Float64}
    coil_voltage_v::Vector{Float64}
    winding_loss_state::Vector{Float64}
    winding_loss_voltage_v::Vector{Float64}
    terminal_voltage_v::Vector{Float64}
    terminal_current_a::Vector{Float64}
    previous_terminal_power_w::Float64
    terminal_power_w::Float64
    supplied_energy_j::Float64
    winding_loss_energy_j::Float64
    frequency_dependent_winding_loss_energy_j::Float64
    dielectric_loss_energy_j::Float64
    stored_magnetic_energy_j::Float64
    stored_frequency_dependent_winding_energy_j::Float64
    stored_electric_energy_j::Float64
    maximum_kcl_residual_a::Float64
    accepted_time_s::Float64
    accepted_step_count::Int
end

mutable struct TransformerWidebandRuntimeState <: AbstractTransformerAcceptedRuntimeState
    rational_state::Vector{Float64}
    terminal_voltage_v::Vector{Float64}
    terminal_current_a::Vector{Float64}
    previous_terminal_power_w::Float64
    terminal_power_w::Float64
    supplied_energy_j::Float64
    minimum_supplied_energy_j::Float64
    dissipated_energy_j::Float64
    stored_energy_j::Float64
    maximum_energy_residual_j::Float64
    maximum_state_magnitude::Float64
    accepted_time_s::Float64
    accepted_step_count::Int
end

mutable struct TransformerMagneticRuntimeState <: AbstractTransformerAcceptedRuntimeState
    coil_current_a::Vector{Float64}
    capacitor_current_a::Vector{Float64}
    coil_voltage_v::Vector{Float64}
    branch_flux_wb::Vector{Float64}
    branch_mmf_drop_at::Vector{Float64}
    magnetic_node_potential_at::Vector{Float64}
    tellinen_state::Vector{Union{Nothing,TellinenMagneticState}}
    winding_loss_state::Vector{Float64}
    winding_loss_voltage_v::Vector{Float64}
    terminal_voltage_v::Vector{Float64}
    terminal_current_a::Vector{Float64}
    previous_terminal_power_w::Float64
    terminal_power_w::Float64
    supplied_energy_j::Float64
    winding_loss_energy_j::Float64
    frequency_dependent_winding_loss_energy_j::Float64
    dielectric_loss_energy_j::Float64
    hysteresis_loss_energy_j::Float64
    classical_eddy_core_loss_energy_j::Float64
    excess_core_loss_energy_j::Float64
    stored_leakage_energy_j::Float64
    stored_magnetic_energy_j::Float64
    stored_frequency_dependent_winding_energy_j::Float64
    stored_electric_energy_j::Float64
    last_classical_dynamic_mmf_at::Vector{Float64}
    last_excess_dynamic_mmf_at::Vector{Float64}
    maximum_magnetic_continuity_residual_wb::Float64
    maximum_magnetic_constitutive_residual_at::Float64
    maximum_local_nonlinear_iterations::Int
    accepted_time_s::Float64
    accepted_step_count::Int
end

mutable struct TransformerNetworkRuntimeState <: AbstractTransformerAcceptedRuntimeState
    represented_node_voltage_v::Vector{Float64}
    terminal_voltage_v::Vector{Float64}
    terminal_current_a::Vector{Float64}
    branch_current_a::Vector{Float64}
    branch_voltage_v::Vector{Float64}
    capacitor_current_a::Vector{Float64}
    previous_terminal_power_w::Float64
    terminal_power_w::Float64
    supplied_energy_j::Float64
    winding_loss_energy_j::Float64
    dielectric_loss_energy_j::Float64
    stored_magnetic_energy_j::Float64
    stored_electric_energy_j::Float64
    maximum_internal_kcl_residual_a::Float64
    accepted_time_s::Float64
    accepted_step_count::Int
end

mutable struct _TransformerTerminalMatrixCandidate
    conductance_s::Matrix{Float64}
    companion_impedance_ohm::Matrix{Float64}
    winding_loss_impedance_ohm::Matrix{Float64}
    history_current_a::Vector{Float64}
    history_voltage_v::Vector{Float64}
    coil_voltage_v::Vector{Float64}
    coil_current_a::Vector{Float64}
    coil_current_workspace_a::Vector{Float64}
    capacitor_current_a::Vector{Float64}
    terminal_current_a::Vector{Float64}
    terminal_jacobian_s::Matrix{Float64}
    coil_jacobian_s::Matrix{Float64}
    terminal_coil_jacobian_s::Matrix{Float64}
    winding_loss_state_transition::Matrix{Float64}
    winding_loss_endpoint_input::Matrix{Float64}
    winding_loss_state::Vector{Float64}
    winding_loss_voltage_v::Vector{Float64}
    prepared_step_s::Float64
    prepared_companion_method::Symbol
end

mutable struct _TransformerWidebandCandidate
    state_transition::Matrix{Float64}
    endpoint_input::Matrix{Float64}
    companion_admittance_s::Matrix{Float64}
    history_current_a::Vector{Float64}
    state::Vector{Float64}
    terminal_current_a::Vector{Float64}
end

mutable struct _TransformerMagneticCandidate
    resistance_ohm::Matrix{Float64}
    leakage_inductance_h::Matrix{Float64}
    capacitance_f::Matrix{Float64}
    conductance_s::Matrix{Float64}
    electrical_impedance_ohm::Matrix{Float64}
    electrical_history_voltage_v::Vector{Float64}
    coil_voltage_v::Vector{Float64}
    coil_current_a::Vector{Float64}
    capacitor_current_a::Vector{Float64}
    branch_flux_wb::Vector{Float64}
    branch_mmf_drop_at::Vector{Float64}
    branch_classical_dynamic_mmf_at::Vector{Float64}
    branch_excess_dynamic_mmf_at::Vector{Float64}
    branch_classical_loss_power_w::Vector{Float64}
    branch_excess_loss_power_w::Vector{Float64}
    branch_differential_reluctance_at_per_wb::Vector{Float64}
    magnetic_node_potential_at::Vector{Float64}
    tellinen_state::Vector{Union{Nothing,TellinenMagneticState}}
    terminal_current_a::Vector{Float64}
    terminal_jacobian_s::Matrix{Float64}
    winding_loss_state_transition::Matrix{Float64}
    winding_loss_endpoint_input::Matrix{Float64}
    winding_loss_state::Vector{Float64}
    winding_loss_voltage_v::Vector{Float64}
    local_iterations::Int
    electrical_residual_v::Float64
    magnetic_constitutive_residual_at::Float64
    magnetic_continuity_residual_wb::Float64
end

mutable struct _TransformerNetworkBranchGroup
    incidence::Matrix{Float64}
    resistance_ohm::Matrix{Float64}
    inductance_h::Matrix{Float64}
    accepted_current_a::Vector{Float64}
    accepted_voltage_v::Vector{Float64}
    conductance_s::Matrix{Float64}
    history_current_a::Vector{Float64}
    current_offset::UnitRange{Int}
end

mutable struct _TransformerNetworkCandidate
    node_order::Vector{Symbol}
    external_node_indices::Vector{Int}
    internal_node_indices::Vector{Int}
    branch_groups::Vector{_TransformerNetworkBranchGroup}
    capacitance_f::Matrix{Float64}
    conductance_s::Matrix{Float64}
    nodal_admittance_s::Matrix{Float64}
    nodal_history_current_a::Vector{Float64}
    companion_admittance_s::Matrix{Float64}
    history_current_a::Vector{Float64}
    represented_node_voltage_v::Vector{Float64}
    branch_current_a::Vector{Float64}
    branch_voltage_v::Vector{Float64}
    capacitor_current_a::Vector{Float64}
end

function _transformer_total_stored_energy(state::TransformerTerminalMatrixRuntimeState)
    return state.stored_magnetic_energy_j +
        state.stored_frequency_dependent_winding_energy_j +
        state.stored_electric_energy_j
end

function _transformer_total_stored_energy(state::TransformerMagneticRuntimeState)
    return state.stored_leakage_energy_j + state.stored_magnetic_energy_j +
        state.stored_frequency_dependent_winding_energy_j +
        state.stored_electric_energy_j
end

_transformer_total_stored_energy(state::TransformerWidebandRuntimeState) =
    state.stored_energy_j

function _transformer_total_stored_energy(state::TransformerNetworkRuntimeState)
    return state.stored_magnetic_energy_j + state.stored_electric_energy_j
end

function _transformer_total_dissipated_energy(state::TransformerTerminalMatrixRuntimeState)
    return state.winding_loss_energy_j +
        state.frequency_dependent_winding_loss_energy_j +
        state.dielectric_loss_energy_j
end

function _transformer_total_dissipated_energy(state::TransformerMagneticRuntimeState)
    return state.winding_loss_energy_j +
        state.frequency_dependent_winding_loss_energy_j +
        state.dielectric_loss_energy_j + state.hysteresis_loss_energy_j +
        state.classical_eddy_core_loss_energy_j + state.excess_core_loss_energy_j
end

_transformer_total_dissipated_energy(state::TransformerWidebandRuntimeState) =
    state.dissipated_energy_j

function _transformer_total_dissipated_energy(state::TransformerNetworkRuntimeState)
    return state.winding_loss_energy_j + state.dielectric_loss_energy_j
end

"""Public physical device state integrated by the existing private nodal/Newton transaction."""
mutable struct TransformerApparatusRuntime{P<:TransformerApparatusPreparation} <:
               AbstractNonlinearCurrentDevice
    preparation::P
    terminal_nodes::Vector{Int}
    accepted_state::AbstractTransformerAcceptedRuntimeState
    candidate::Union{
        _TransformerTerminalMatrixCandidate,
        _TransformerWidebandCandidate,
        _TransformerMagneticCandidate,
        _TransformerNetworkCandidate,
    }
    candidate_time_s::Float64
    candidate_step_s::Float64
    companion_method::Symbol
    prepared::Bool
    preparation_count::Int
    trial_evaluation_count::Int
    accepted_evaluation_count::Int
    rejected_trial_count::Int
    initial_stored_energy_j::Float64
    event_state::TransformerApparatusEventState
    candidate_apparatus_terminal_voltage_v::Vector{Float64}
    candidate_network_terminal_current_a::Vector{Float64}
    candidate_network_terminal_jacobian_s::Matrix{Float64}
end

function _transformer_terminal_runtime_state(preparation::TransformerApparatusPreparation)
    coil_count = length(preparation.coil_order)
    terminal_count = length(preparation.terminal_order)
    model = preparation.specification.model
    winding_loss_state_count = model isa HybridTransformerModel ?
        size(model.winding_loss_state_matrix_per_s, 1) : 0
    matrices = _transformer_terminal_matrices(preparation)
    connection = preparation.specification.connection
    coil_voltage = transpose(connection.incidence) *
        preparation.initial_node_voltage_v
    terminal_voltage_discrete_derivative = _transformer_initial_discrete_derivative(
        preparation.initial_node_voltage_v,
        preparation.initial_node_voltage_derivative_v_per_s,
        preparation.initialization_mode,
        preparation.specification.settings,
    )
    coil_voltage_derivative = transpose(connection.incidence) *
        terminal_voltage_discrete_derivative
    capacitor_current = matrices.capacitance_f * coil_voltage_derivative
    winding_loss_state = zeros(winding_loss_state_count)
    winding_loss_voltage = zeros(coil_count)
    if model isa HybridTransformerModel
        winding_loss = _transformer_initial_winding_loss_state(
            model,
            preparation.initialization_mode,
            preparation.initial_coil_current_a,
            preparation.initial_coil_current_derivative_a_per_s,
            preparation.specification.settings,
        )
        winding_loss_state .= winding_loss.state
        winding_loss_voltage .= winding_loss.voltage
    end
    total_coil_current = preparation.initial_coil_current_a .+
        capacitor_current .+ matrices.conductance_s * coil_voltage
    terminal_current = connection.incidence * total_coil_current
    terminal_power = dot(preparation.initial_node_voltage_v, terminal_current)
    stored_magnetic_energy = 0.5 * dot(
        preparation.initial_coil_current_a,
        matrices.inductance_h * preparation.initial_coil_current_a,
    )
    stored_winding_loss_energy = model isa HybridTransformerModel ?
        0.5 * dot(
            winding_loss_state,
            model.winding_loss_storage_matrix_j * winding_loss_state,
        ) : 0.0
    stored_electric_energy =
        0.5 * dot(coil_voltage, matrices.capacitance_f * coil_voltage)
    state = TransformerTerminalMatrixRuntimeState(
        copy(preparation.initial_coil_current_a),
        capacitor_current,
        coil_voltage,
        winding_loss_state,
        winding_loss_voltage,
        copy(preparation.initial_node_voltage_v),
        terminal_current,
        terminal_power,
        terminal_power,
        0.0,
        0.0,
        0.0,
        0.0,
        stored_magnetic_energy,
        stored_winding_loss_energy,
        stored_electric_energy,
        0.0,
        preparation.initial_time_s,
        0,
    )
    return state
end

function _transformer_terminal_matrices(preparation::TransformerApparatusPreparation)
    model = preparation.specification.model
    if model isa LowFrequencyTransformerModel || model isa BCTRANTransformerModel
        return model.matrices
    elseif model isa HybridTransformerModel
        preparation.effective_linear_inductance_h === nothing && _transformer_refusal(
            :nonlinear_magnetic_runtime_required,
            :runtime_construction,
            preparation.specification.id,
            preparation.specification.tier,
            "hybrid nonlinear magnetic material requires the coupled magnetic runtime",
        )
        return TransformerTerminalMatrices(
            model.leakage.resistance_ohm,
            something(preparation.effective_linear_inductance_h);
            capacitance_f=model.leakage.capacitance_f,
            conductance_s=model.leakage.conductance_s,
        )
    elseif model isa MagneticEquivalentCircuitModel
        preparation.effective_linear_inductance_h === nothing && _transformer_refusal(
            :nonlinear_magnetic_runtime_required,
            :runtime_construction,
            preparation.specification.id,
            preparation.specification.tier,
            "MEC nonlinear magnetic material requires the coupled magnetic runtime",
        )
        return TransformerTerminalMatrices(
            model.winding_resistance_ohm,
            something(preparation.effective_linear_inductance_h);
            capacitance_f=model.terminal_capacitance_f,
            conductance_s=model.terminal_conductance_s,
        )
    end
    _transformer_refusal(
        :incompatible_runtime_model,
        :runtime_construction,
        preparation.specification.id,
        preparation.specification.tier,
        "selected transformer tier is not a coupled terminal-matrix runtime",
    )
end

function _transformer_terminal_candidate(
    preparation::TransformerApparatusPreparation,
)
    terminal_count = length(preparation.terminal_order)
    coil_count = length(preparation.coil_order)
    model = preparation.specification.model
    winding_loss_state_count = model isa HybridTransformerModel ?
        size(model.winding_loss_state_matrix_per_s, 1) : 0
    return _TransformerTerminalMatrixCandidate(
        zeros(coil_count, coil_count),
        zeros(coil_count, coil_count),
        zeros(coil_count, coil_count),
        zeros(coil_count),
        zeros(coil_count),
        zeros(coil_count),
        zeros(coil_count),
        zeros(coil_count),
        zeros(coil_count),
        zeros(terminal_count),
        zeros(terminal_count, terminal_count),
        zeros(coil_count, coil_count),
        zeros(terminal_count, coil_count),
        zeros(winding_loss_state_count, winding_loss_state_count),
        zeros(winding_loss_state_count, coil_count),
        zeros(winding_loss_state_count),
        zeros(coil_count),
        NaN,
        :unprepared,
    )
end

function _transformer_magnetic_model_parts(
    preparation::TransformerApparatusPreparation,
)
    model = preparation.specification.model
    if model isa HybridTransformerModel
        return (
            resistance_ohm=model.leakage.resistance_ohm,
            leakage_inductance_h=model.leakage.inductance_h,
            capacitance_f=model.leakage.capacitance_f,
            conductance_s=model.leakage.conductance_s,
            magnetic_graph=model.magnetic_graph,
            winding_loss_model=model,
        )
    end
    model = model::MagneticEquivalentCircuitModel
    return (
        resistance_ohm=model.winding_resistance_ohm,
        leakage_inductance_h=model.leakage_inductance_h,
        capacitance_f=model.terminal_capacitance_f,
        conductance_s=model.terminal_conductance_s,
        magnetic_graph=model.magnetic_graph,
        winding_loss_model=nothing,
    )
end

function _transformer_piecewise_integral(grid_x, grid_y, start_x::Float64, end_x::Float64)
    start_x == end_x && return 0.0
    lower = min(start_x, end_x)
    upper = max(start_x, end_x)
    first(grid_x) <= lower <= upper <= last(grid_x) || throw(DomainError(
        (start_x, end_x),
        "transformer magnetic energy integral is outside its material domain",
    ))
    points = Float64[lower]
    append!(points, value for value in grid_x if lower < value < upper)
    push!(points, upper)
    integral = 0.0
    for index in 1:(length(points) - 1)
        left_value, _ = _linear_interpolation_and_slope(
            points[index],
            grid_x,
            grid_y,
        )
        right_value, _ = _linear_interpolation_and_slope(
            points[index + 1],
            grid_x,
            grid_y,
        )
        integral += 0.5 * (left_value + right_value) *
            (points[index + 1] - points[index])
    end
    return start_x <= end_x ? integral : -integral
end

function _transformer_magnetic_energy_density(
    material::LinearTransformerMagneticMaterial,
    flux_density_t::Float64,
)
    field_strength = magnetic_material_field(material, flux_density_t)
    return 0.5 * flux_density_t * field_strength
end

function _transformer_magnetic_energy_density(
    material::PiecewiseLinearTransformerMagneticMaterial,
    flux_density_t::Float64,
)
    magnitude = abs(flux_density_t)
    return _transformer_piecewise_integral(
        material.flux_density_t,
        material.field_strength_a_per_m,
        0.0,
        magnitude,
    )
end

function _transformer_magnetic_energy_density(
    material::TellinenTransformerMagneticMaterial,
    flux_density_t::Float64,
)
    anhysteretic_flux_density = 0.5 .* (
        material.lower_branch.flux_density_t .+
        material.upper_branch.flux_density_t
    )
    return _transformer_piecewise_integral(
        anhysteretic_flux_density,
        material.lower_branch.field_strength_a_per_m,
        0.0,
        flux_density_t,
    )
end

function _transformer_branch_magnetic_trial(
    graph::TransformerMagneticGraph,
    branch_index::Int,
    branch_flux_wb::Float64,
    accepted_tellinen_state::Union{Nothing,TellinenMagneticState},
)
    branch = graph.branches[branch_index]
    material = graph.materials[branch.material_index]
    flux_density = branch_flux_wb / branch.cross_section_m2
    air_gap_flux_density = flux_density / branch.air_gap_effective_area_factor
    core_length = branch.length_m - branch.air_gap_length_m
    trial_tellinen_state = nothing
    if material isa TellinenTransformerMagneticMaterial
        accepted_tellinen_state === nothing && throw(ArgumentError(
            "Tellinen transformer branch is missing its accepted hysteresis state",
        ))
        tellinen_trial = if flux_density == accepted_tellinen_state.flux_density_t
            (
                state=accepted_tellinen_state,
                differential_reluctivity_m_per_h=inv(
                    _tellinen_flux_derivative(
                        material,
                        accepted_tellinen_state.field_strength_a_per_m,
                        accepted_tellinen_state.flux_density_t,
                        accepted_tellinen_state.direction,
                    ),
                ),
            )
        else
            tellinen_trial_from_flux_density(
                material,
                accepted_tellinen_state,
                flux_density,
            )
        end
        field_strength = tellinen_trial.state.field_strength_a_per_m
        differential_reluctivity = tellinen_trial.differential_reluctivity_m_per_h
        trial_tellinen_state = tellinen_trial.state
    else
        field_strength = magnetic_material_field(material, flux_density)
        differential_reluctivity = magnetic_material_differential_reluctivity(
            material,
            flux_density,
        )
    end
    mmf_drop = core_length * field_strength +
        branch.air_gap_length_m * air_gap_flux_density /
        TRANSFORMER_VACUUM_PERMEABILITY_H_PER_M
    differential_mmf_per_flux = (
        core_length * differential_reluctivity / branch.cross_section_m2 +
        branch.air_gap_length_m / (
            TRANSFORMER_VACUUM_PERMEABILITY_H_PER_M *
            branch.cross_section_m2 * branch.air_gap_effective_area_factor
        )
    )
    isfinite(mmf_drop) && isfinite(differential_mmf_per_flux) &&
        differential_mmf_per_flux > 0.0 || throw(ArgumentError(
            "transformer magnetic branch has a nonpositive differential reluctance",
        ))
    core_energy = _transformer_magnetic_energy_density(material, flux_density) *
        core_length * branch.cross_section_m2
    gap_energy = 0.5 * branch.air_gap_length_m *
        branch.cross_section_m2 * branch.air_gap_effective_area_factor *
        air_gap_flux_density^2 / TRANSFORMER_VACUUM_PERMEABILITY_H_PER_M
    return (
        flux_density_t=flux_density,
        mmf_drop_at=mmf_drop,
        differential_reluctance_at_per_wb=differential_mmf_per_flux,
        tellinen_state=trial_tellinen_state,
        stored_energy_j=core_energy + gap_energy,
    )
end

function _transformer_magnetic_residual!(
    branch_mmf_drop_at,
    branch_differential_reluctance_at_per_wb,
    trial_tellinen_state,
    graph::TransformerMagneticGraph,
    branch_flux_wb,
    magnetic_node_potential_at,
    coil_current_a,
    accepted_tellinen_state,
    ;
    previous_branch_flux_wb=nothing,
    previous_classical_dynamic_mmf_at=nothing,
    previous_excess_dynamic_mmf_at=nothing,
    step_s=nothing,
    branch_classical_dynamic_mmf_at=nothing,
    branch_excess_dynamic_mmf_at=nothing,
    branch_classical_loss_power_w=nothing,
    branch_excess_loss_power_w=nothing,
    companion_method::Symbol=:trapezoidal,
)
    for branch_index in eachindex(graph.branches)
        trial = _transformer_branch_magnetic_trial(
            graph,
            branch_index,
            branch_flux_wb[branch_index],
            accepted_tellinen_state[branch_index],
        )
        branch_mmf_drop_at[branch_index] = trial.mmf_drop_at
        branch_differential_reluctance_at_per_wb[branch_index] =
            trial.differential_reluctance_at_per_wb
        trial_tellinen_state[branch_index] = trial.tellinen_state
        classical_mmf = 0.0
        excess_mmf = 0.0
        classical_power = 0.0
        excess_power = 0.0
        loss = graph.dynamic_core_loss[branch_index]
        if loss !== nothing && previous_branch_flux_wb !== nothing && step_s !== nothing
            branch = graph.branches[branch_index]
            flux_density_rate = (
                branch_flux_wb[branch_index] - previous_branch_flux_wb[branch_index]
            ) / (branch.cross_section_m2 * step_s)
            dynamic = _transformer_dynamic_core_loss_field(loss, flux_density_rate)
            core_length = branch.length_m - branch.air_gap_length_m
            core_volume = core_length * branch.cross_section_m2
            previous_classical_mmf = previous_classical_dynamic_mmf_at === nothing ?
                0.0 : previous_classical_dynamic_mmf_at[branch_index]
            previous_excess_mmf = previous_excess_dynamic_mmf_at === nothing ?
                0.0 : previous_excess_dynamic_mmf_at[branch_index]
            integration_factor = companion_method === :backward_euler ? 1.0 : 2.0
            classical_mmf = integration_factor * core_length *
                dynamic.classical_field_a_per_m -
                (companion_method === :backward_euler ? 0.0 : previous_classical_mmf)
            excess_mmf = integration_factor * core_length *
                dynamic.excess_field_a_per_m -
                (companion_method === :backward_euler ? 0.0 : previous_excess_mmf)
            branch_mmf_drop_at[branch_index] += classical_mmf + excess_mmf
            branch_differential_reluctance_at_per_wb[branch_index] +=
                integration_factor * core_length *
                dynamic.differential_field_a_s_per_m_t /
                (branch.cross_section_m2 * step_s)
            classical_power = core_volume * dynamic.classical_loss_density_w_per_m3
            excess_power = core_volume * dynamic.excess_loss_density_w_per_m3
        end
        branch_classical_dynamic_mmf_at === nothing ||
            (branch_classical_dynamic_mmf_at[branch_index] = classical_mmf)
        branch_excess_dynamic_mmf_at === nothing ||
            (branch_excess_dynamic_mmf_at[branch_index] = excess_mmf)
        branch_classical_loss_power_w === nothing ||
            (branch_classical_loss_power_w[branch_index] = classical_power)
        branch_excess_loss_power_w === nothing ||
            (branch_excess_loss_power_w[branch_index] = excess_power)
    end
    constitutive_residual = branch_mmf_drop_at .+
        transpose(graph.incidence) * magnetic_node_potential_at .-
        graph.winding_turns * coil_current_a
    continuity_residual = graph.incidence * branch_flux_wb
    return constitutive_residual, continuity_residual
end

function _transformer_initial_tellinen_states(
    graph::TransformerMagneticGraph,
    branch_flux_wb=zeros(length(graph.branches)),
    branch_flux_derivative_wb_per_s=zeros(length(graph.branches)),
)
    result = Vector{Union{Nothing,TellinenMagneticState}}(undef, length(graph.branches))
    for branch_index in eachindex(graph.branches)
        branch = graph.branches[branch_index]
        material = graph.materials[branch.material_index]
        result[branch_index] = material isa TellinenTransformerMagneticMaterial ?
            tellinen_state(
                material;
                field_strength_a_per_m=0.0,
                flux_density_t=branch_flux_wb[branch_index] /
                    branch.cross_section_m2,
                direction=branch_flux_derivative_wb_per_s[branch_index] < 0.0 ?
                    -1 : 1,
            ) : nothing
    end
    return result
end

function _transformer_initial_dynamic_core_loss(
    graph::TransformerMagneticGraph,
    branch_flux_derivative_wb_per_s,
)
    classical_dynamic_mmf = zeros(length(graph.branches))
    excess_dynamic_mmf = zeros(length(graph.branches))
    for branch_index in eachindex(graph.branches)
        loss = graph.dynamic_core_loss[branch_index]
        loss === nothing && continue
        branch = graph.branches[branch_index]
        flux_density_rate = branch_flux_derivative_wb_per_s[branch_index] /
            branch.cross_section_m2
        dynamic = _transformer_dynamic_core_loss_field(loss, flux_density_rate)
        core_length = branch.length_m - branch.air_gap_length_m
        classical_dynamic_mmf[branch_index] =
            core_length * dynamic.classical_field_a_per_m
        excess_dynamic_mmf[branch_index] =
            core_length * dynamic.excess_field_a_per_m
    end
    return classical_dynamic_mmf, excess_dynamic_mmf
end

function _solve_transformer_prescribed_magnetic_state!(
    branch_mmf_drop_at,
    branch_differential_reluctance_at_per_wb,
    magnetic_node_potential_at,
    trial_tellinen_state,
    graph::TransformerMagneticGraph,
    branch_flux_wb,
    coil_current_a,
    accepted_tellinen_state,
    additional_branch_mmf_at,
    settings::TransformerRuntimeSettings,
)
    _transformer_magnetic_residual!(
        branch_mmf_drop_at,
        branch_differential_reluctance_at_per_wb,
        trial_tellinen_state,
        graph,
        branch_flux_wb,
        magnetic_node_potential_at,
        coil_current_a,
        accepted_tellinen_state,
    )
    branch_mmf_drop_at .+= additional_branch_mmf_at
    right_hand_side = graph.winding_turns * coil_current_a .- branch_mmf_drop_at
    magnetic_node_potential_at .=
        (graph.incidence * transpose(graph.incidence)) \
        (graph.incidence * right_hand_side)
    constitutive_residual = branch_mmf_drop_at .+
        transpose(graph.incidence) * magnetic_node_potential_at .-
        graph.winding_turns * coil_current_a
    continuity_residual = graph.incidence * branch_flux_wb
    constitutive_scale = max(
        maximum(abs, graph.winding_turns * coil_current_a; init=0.0),
        maximum(abs, branch_mmf_drop_at; init=0.0),
        1.0,
    )
    flux_scale = max(maximum(abs, branch_flux_wb; init=0.0), 1.0e-12)
    maximum_constitutive_residual =
        maximum(abs, constitutive_residual; init=0.0)
    maximum_continuity_residual = maximum(abs, continuity_residual; init=0.0)
    maximum_constitutive_residual <=
        settings.nonlinear_residual_relative_tolerance * constitutive_scale ||
        throw(DomainError(
            maximum_constitutive_residual,
            "prescribed transformer residual flux is incompatible with its material law and winding currents",
        ))
    maximum_continuity_residual <=
        settings.magnetic_continuity_absolute_tolerance_wb +
        settings.nonlinear_residual_relative_tolerance * flux_scale ||
        throw(DomainError(
            maximum_continuity_residual,
            "prescribed transformer residual flux violates magnetic continuity",
        ))
    return maximum_constitutive_residual, maximum_continuity_residual
end

function _solve_transformer_initial_magnetic_state!(
    branch_flux_wb,
    branch_mmf_drop_at,
    branch_differential_reluctance_at_per_wb,
    magnetic_node_potential_at,
    trial_tellinen_state,
    graph::TransformerMagneticGraph,
    coil_current_a,
    accepted_tellinen_state,
    settings::TransformerRuntimeSettings,
)
    branch_count = length(graph.branches)
    node_count = length(graph.node_order)
    for iteration in 1:settings.maximum_local_nonlinear_iterations
        constitutive_residual, continuity_residual = _transformer_magnetic_residual!(
            branch_mmf_drop_at,
            branch_differential_reluctance_at_per_wb,
            trial_tellinen_state,
            graph,
            branch_flux_wb,
            magnetic_node_potential_at,
            coil_current_a,
            accepted_tellinen_state,
        )
        constitutive_scale = max(
            maximum(abs, graph.winding_turns * coil_current_a; init=0.0),
            maximum(abs, branch_mmf_drop_at; init=0.0),
            1.0,
        )
        flux_scale = max(maximum(abs, branch_flux_wb; init=0.0), 1.0e-12)
        if maximum(abs, constitutive_residual; init=0.0) <=
               settings.nonlinear_residual_relative_tolerance * constitutive_scale &&
           maximum(abs, continuity_residual; init=0.0) <=
               settings.magnetic_continuity_absolute_tolerance_wb +
               settings.nonlinear_residual_relative_tolerance * flux_scale
            return iteration
        end
        system = [
            Diagonal(branch_differential_reluctance_at_per_wb) transpose(graph.incidence)
            graph.incidence zeros(Float64, node_count, node_count)
        ]
        correction = -(system \ vcat(constitutive_residual, continuity_residual))
        branch_flux_wb .+= correction[1:branch_count]
        magnetic_node_potential_at .+= correction[(branch_count + 1):end]
    end
    throw(ArgumentError("transformer initial nonlinear magnetic state did not converge"))
end

function _transformer_magnetic_runtime_state(
    preparation::TransformerApparatusPreparation,
)
    parts = _transformer_magnetic_model_parts(preparation)
    graph = parts.magnetic_graph
    coil_count = length(preparation.coil_order)
    terminal_count = length(preparation.terminal_order)
    branch_count = length(graph.branches)
    node_count = length(graph.node_order)
    branch_flux = preparation.initial_branch_flux_wb === nothing ?
        zeros(branch_count) : copy(preparation.initial_branch_flux_wb)
    supplied_branch_flux_derivative =
        preparation.initial_branch_flux_derivative_wb_per_s === nothing ?
        zeros(branch_count) :
        copy(preparation.initial_branch_flux_derivative_wb_per_s)
    branch_flux_derivative = _transformer_initial_discrete_derivative(
        branch_flux,
        supplied_branch_flux_derivative,
        preparation.initialization_mode,
        preparation.specification.settings,
    )
    branch_mmf = zeros(branch_count)
    branch_differential = zeros(branch_count)
    magnetic_potential = zeros(node_count)
    accepted_tellinen = _transformer_initial_tellinen_states(
        graph,
        branch_flux,
        branch_flux_derivative,
    )
    trial_tellinen = copy(accepted_tellinen)
    classical_dynamic_mmf, excess_dynamic_mmf =
        _transformer_initial_dynamic_core_loss(graph, branch_flux_derivative)
    initial_iteration_count = 1
    maximum_constitutive_residual = 0.0
    maximum_continuity_residual = 0.0
    if preparation.initial_branch_flux_wb === nothing
        initial_iteration_count = _solve_transformer_initial_magnetic_state!(
            branch_flux,
            branch_mmf,
            branch_differential,
            magnetic_potential,
            trial_tellinen,
            graph,
            preparation.initial_coil_current_a,
            accepted_tellinen,
            preparation.specification.settings,
        )
        accepted_tellinen .= trial_tellinen
    end
    try
        maximum_constitutive_residual, maximum_continuity_residual =
            _solve_transformer_prescribed_magnetic_state!(
                branch_mmf,
                branch_differential,
                magnetic_potential,
                trial_tellinen,
                graph,
                branch_flux,
                preparation.initial_coil_current_a,
                accepted_tellinen,
                classical_dynamic_mmf .+ excess_dynamic_mmf,
                preparation.specification.settings,
            )
    catch error
        _transformer_refusal(
            :inconsistent_initial_magnetic_state,
            :initialize,
            preparation.specification.id,
            preparation.specification.tier,
            sprint(showerror, error),
        )
    end
    accepted_tellinen .= trial_tellinen
    loss_model = parts.winding_loss_model
    winding_loss_state_count = loss_model === nothing ? 0 :
        size(loss_model.winding_loss_state_matrix_per_s, 1)
    winding_loss_state = zeros(winding_loss_state_count)
    winding_loss_voltage = zeros(coil_count)
    if loss_model !== nothing
        winding_loss = _transformer_initial_winding_loss_state(
            loss_model,
            preparation.initialization_mode,
            preparation.initial_coil_current_a,
            preparation.initial_coil_current_derivative_a_per_s,
            preparation.specification.settings,
        )
        winding_loss_state .= winding_loss.state
        winding_loss_voltage .= winding_loss.voltage
    end
    connection = preparation.specification.connection
    coil_voltage = transpose(connection.incidence) *
        preparation.initial_node_voltage_v
    terminal_voltage_discrete_derivative = _transformer_initial_discrete_derivative(
        preparation.initial_node_voltage_v,
        preparation.initial_node_voltage_derivative_v_per_s,
        preparation.initialization_mode,
        preparation.specification.settings,
    )
    coil_voltage_derivative = transpose(connection.incidence) *
        terminal_voltage_discrete_derivative
    capacitor_current = parts.capacitance_f * coil_voltage_derivative
    total_coil_current = preparation.initial_coil_current_a .+
        capacitor_current .+ parts.conductance_s * coil_voltage
    terminal_current = connection.incidence * total_coil_current
    terminal_power = dot(preparation.initial_node_voltage_v, terminal_current)
    stored_leakage_energy = 0.5 * dot(
        preparation.initial_coil_current_a,
        parts.leakage_inductance_h * preparation.initial_coil_current_a,
    )
    stored_electric_energy =
        0.5 * dot(coil_voltage, parts.capacitance_f * coil_voltage)
    stored_winding_loss_energy = loss_model === nothing ? 0.0 :
        0.5 * dot(
            winding_loss_state,
            loss_model.winding_loss_storage_matrix_j * winding_loss_state,
        )
    state = TransformerMagneticRuntimeState(
        copy(preparation.initial_coil_current_a),
        capacitor_current,
        coil_voltage,
        branch_flux,
        branch_mmf,
        magnetic_potential,
        accepted_tellinen,
        winding_loss_state,
        winding_loss_voltage,
        copy(preparation.initial_node_voltage_v),
        terminal_current,
        terminal_power,
        terminal_power,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        stored_leakage_energy,
        0.0,
        stored_winding_loss_energy,
        stored_electric_energy,
        classical_dynamic_mmf,
        excess_dynamic_mmf,
        maximum_continuity_residual,
        maximum_constitutive_residual,
        initial_iteration_count,
        preparation.initial_time_s,
        0,
    )
    stored_magnetic_energy = 0.0
    for branch_index in eachindex(graph.branches)
        stored_magnetic_energy += _transformer_branch_magnetic_trial(
            graph,
            branch_index,
            state.branch_flux_wb[branch_index],
            state.tellinen_state[branch_index],
        ).stored_energy_j
    end
    state.stored_magnetic_energy_j = stored_magnetic_energy
    return state
end

function _transformer_magnetic_candidate(
    preparation::TransformerApparatusPreparation,
)
    parts = _transformer_magnetic_model_parts(preparation)
    graph = parts.magnetic_graph
    coil_count = length(preparation.coil_order)
    terminal_count = length(preparation.terminal_order)
    branch_count = length(graph.branches)
    node_count = length(graph.node_order)
    loss_model = parts.winding_loss_model
    winding_loss_state_count = loss_model === nothing ? 0 :
        size(loss_model.winding_loss_state_matrix_per_s, 1)
    return _TransformerMagneticCandidate(
        copy(parts.resistance_ohm),
        copy(parts.leakage_inductance_h),
        copy(parts.capacitance_f),
        copy(parts.conductance_s),
        zeros(coil_count, coil_count),
        zeros(coil_count),
        zeros(coil_count),
        zeros(coil_count),
        zeros(coil_count),
        zeros(branch_count),
        zeros(branch_count),
        zeros(branch_count),
        zeros(branch_count),
        zeros(branch_count),
        zeros(branch_count),
        zeros(branch_count),
        zeros(node_count),
        Vector{Union{Nothing,TellinenMagneticState}}(nothing, branch_count),
        zeros(terminal_count),
        zeros(terminal_count, terminal_count),
        zeros(winding_loss_state_count, winding_loss_state_count),
        zeros(winding_loss_state_count, coil_count),
        zeros(winding_loss_state_count),
        zeros(coil_count),
        0,
        Inf,
        Inf,
        Inf,
    )
end

function _transformer_wideband_runtime_state(preparation::TransformerApparatusPreparation)
    model = preparation.specification.model::WidebandTransformerModel
    rational_state = if preparation.initialization_mode ===
                        SinusoidalTransformerOperatingPoint
        _transformer_sinusoidal_linear_state(
            model.state_matrix_per_s,
            model.input_matrix,
            preparation.initial_node_voltage_v,
            preparation.initial_node_voltage_derivative_v_per_s,
            preparation.specification.settings,
        ).state
    else
        zeros(size(model.state_matrix_per_s, 1))
    end
    terminal_current = model.output_matrix_s_per_s * rational_state .+
        model.direct_admittance_s * preparation.initial_node_voltage_v
    terminal_power = dot(preparation.initial_node_voltage_v, terminal_current)
    stored_energy =
        0.5 * dot(rational_state, model.storage_matrix_j * rational_state)
    return TransformerWidebandRuntimeState(
        rational_state,
        copy(preparation.initial_node_voltage_v),
        terminal_current,
        terminal_power,
        terminal_power,
        0.0,
        0.0,
        0.0,
        stored_energy,
        0.0,
        maximum(abs, rational_state; init=0.0),
        preparation.initial_time_s,
        0,
    )
end

function _transformer_wideband_candidate(preparation::TransformerApparatusPreparation)
    model = preparation.specification.model::WidebandTransformerModel
    state_count = size(model.state_matrix_per_s, 1)
    terminal_count = length(preparation.terminal_order)
    return _TransformerWidebandCandidate(
        zeros(state_count, state_count),
        zeros(state_count, terminal_count),
        zeros(terminal_count, terminal_count),
        zeros(terminal_count),
        zeros(state_count),
        zeros(terminal_count),
    )
end

function _single_ladder_branch_group(
    branch::TransformerLadderBranch,
    node_count::Int,
    current_index::Int,
)
    incidence = zeros(Float64, node_count, 1)
    branch.from_node_index == 0 || (incidence[branch.from_node_index, 1] = 1.0)
    branch.to_node_index == 0 || (incidence[branch.to_node_index, 1] = -1.0)
    return _TransformerNetworkBranchGroup(
        incidence,
        reshape([branch.resistance_ohm], 1, 1),
        reshape([branch.inductance_h], 1, 1),
        zeros(1),
        zeros(1),
        zeros(1, 1),
        zeros(1),
        current_index:current_index,
    )
end

function _grey_box_candidate(preparation::TransformerApparatusPreparation)
    model = preparation.specification.model::GreyBoxTransformerModel
    node_count = length(model.node_order)
    branch_groups = _TransformerNetworkBranchGroup[
        _single_ladder_branch_group(branch, node_count, index)
        for (index, branch) in pairs(model.branches)
    ]
    internal_indices = setdiff(collect(1:node_count), model.terminal_node_indices)
    return _TransformerNetworkCandidate(
        copy(model.node_order),
        copy(model.terminal_node_indices),
        internal_indices,
        branch_groups,
        copy(model.capacitance_f),
        copy(model.conductance_s),
        zeros(node_count, node_count),
        zeros(node_count),
        zeros(length(model.terminal_node_indices), length(model.terminal_node_indices)),
        zeros(length(model.terminal_node_indices)),
        zeros(node_count),
        zeros(length(model.branches)),
        zeros(length(model.branches)),
        zeros(node_count),
    )
end

function _white_box_candidate(preparation::TransformerApparatusPreparation)
    model = preparation.specification.model::WhiteBoxTransformerModel
    winding_count = length(model.winding_order)
    maximum_section_count = maximum(model.section_count_per_winding)
    connection = preparation.specification.connection
    length(connection.coil_order) == winding_count || _transformer_refusal(
        :winding_coil_count_mismatch,
        :runtime_construction,
        preparation.specification.id,
        preparation.specification.tier,
        "white-box runtime requires one declared coil per winding in this coupled section slice",
    )
    terminal_count = length(connection.node_order)
    terminal_count == 2 * winding_count || _transformer_refusal(
        :white_box_terminal_topology,
        :runtime_construction,
        preparation.specification.id,
        preparation.specification.tier,
        "white-box section runtime requires explicit start and end terminals for every winding",
    )
    internal_count = sum(model.section_count_per_winding .- 1)
    node_count = terminal_count + internal_count
    node_order = vcat(
        connection.node_order,
        Symbol[
            Symbol(winding, :_section_, section)
            for (winding, section_count) in
                zip(model.winding_order, model.section_count_per_winding)
            for section in 1:(section_count - 1)
        ],
    )
    internal_start = Vector{Int}(undef, winding_count)
    internal_cursor = terminal_count + 1
    for winding in 1:winding_count
        internal_start[winding] = internal_cursor
        internal_cursor += model.section_count_per_winding[winding] - 1
    end
    function winding_section_node(winding::Int, boundary::Int)
        section_count = model.section_count_per_winding[winding]
        0 <= boundary <= section_count || throw(BoundsError(
            0:section_count,
            boundary,
        ))
        boundary == 0 && return 2 * winding - 1
        boundary == section_count && return 2 * winding
        return internal_start[winding] + boundary - 1
    end
    branch_groups = _TransformerNetworkBranchGroup[]
    current_cursor = 1
    capacitance = zeros(Float64, node_count, node_count)
    conductance = zeros(Float64, node_count, node_count)
    for section in 1:maximum_section_count
        active_windings = findall(
            section_count -> section <= section_count,
            model.section_count_per_winding,
        )
        active_count = length(active_windings)
        incidence = zeros(Float64, node_count, active_count)
        for (active_index, winding) in pairs(active_windings)
            incidence[winding_section_node(winding, section - 1), active_index] = 1.0
            incidence[winding_section_node(winding, section), active_index] = -1.0
        end
        length_m = model.section_length_m[section]
        resistance = model.series_resistance_ohm_per_m[section][
            active_windings,
            active_windings,
        ] .* length_m
        inductance = model.series_inductance_h_per_m[section][
            active_windings,
            active_windings,
        ] .* length_m
        group = _TransformerNetworkBranchGroup(
            incidence,
            resistance,
            inductance,
            zeros(active_count),
            zeros(active_count),
            zeros(active_count, active_count),
            zeros(active_count),
            current_cursor:(current_cursor + active_count - 1),
        )
        push!(branch_groups, group)
        current_cursor += active_count
        shunt_incidence = zeros(Float64, node_count, active_count)
        for (active_index, winding) in pairs(active_windings)
            shunt_incidence[winding_section_node(winding, section), active_index] = 1.0
        end
        capacitance .+= shunt_incidence *
            (model.shunt_capacitance_f_per_m[section][
                active_windings,
                active_windings,
            ] .* length_m) *
            transpose(shunt_incidence)
        conductance .+= shunt_incidence *
            (model.shunt_conductance_s_per_m[section][
                active_windings,
                active_windings,
            ] .* length_m) *
            transpose(shunt_incidence)
    end
    external_indices = collect(1:terminal_count)
    internal_indices = collect((terminal_count + 1):node_count)
    branch_value_count = current_cursor - 1
    return _TransformerNetworkCandidate(
        node_order,
        external_indices,
        internal_indices,
        branch_groups,
        capacitance,
        conductance,
        zeros(node_count, node_count),
        zeros(node_count),
        zeros(terminal_count, terminal_count),
        zeros(terminal_count),
        zeros(node_count),
        zeros(branch_value_count),
        zeros(branch_value_count),
        zeros(node_count),
    )
end

function _transformer_sinusoidal_network_initial_state(
    preparation::TransformerApparatusPreparation,
    candidate::_TransformerNetworkCandidate,
)
    physical_frequency, discrete_frequency =
        _transformer_sinusoidal_frequencies(preparation.specification.settings)
    external_voltage_phasor = _transformer_sinusoidal_phasor(
        preparation.initial_node_voltage_v,
        preparation.initial_node_voltage_derivative_v_per_s,
        physical_frequency,
    )
    node_count = length(candidate.node_order)
    nodal_admittance = ComplexF64.(candidate.conductance_s) .+
        im * discrete_frequency .* ComplexF64.(candidate.capacitance_f)
    branch_current_phasor = zeros(ComplexF64, length(candidate.branch_current_a))
    branch_voltage_phasor = zeros(ComplexF64, length(candidate.branch_voltage_v))
    branch_admittances = Matrix{ComplexF64}[]
    for group in candidate.branch_groups
        branch_impedance = ComplexF64.(group.resistance_ohm) .+
            im * discrete_frequency .* ComplexF64.(group.inductance_h)
        branch_admittance = inv(branch_impedance)
        push!(branch_admittances, branch_admittance)
        complex_incidence = ComplexF64.(group.incidence)
        nodal_admittance .+= complex_incidence * branch_admittance *
            transpose(complex_incidence)
    end
    represented_voltage_phasor = zeros(ComplexF64, node_count)
    external = candidate.external_node_indices
    internal = candidate.internal_node_indices
    represented_voltage_phasor[external] .= external_voltage_phasor
    if !isempty(internal)
        represented_voltage_phasor[internal] .=
            -(nodal_admittance[internal, internal] \
              (
                  nodal_admittance[internal, external] *
                  external_voltage_phasor
              ))
    end
    for (group, branch_admittance) in
        zip(candidate.branch_groups, branch_admittances)
        branch_voltage = transpose(ComplexF64.(group.incidence)) *
            represented_voltage_phasor
        branch_current = branch_admittance * branch_voltage
        branch_voltage_phasor[group.current_offset] .= branch_voltage
        branch_current_phasor[group.current_offset] .= branch_current
    end
    nodal_current_phasor = nodal_admittance * represented_voltage_phasor
    return (
        represented_node_voltage_v=real.(represented_voltage_phasor),
        represented_node_voltage_derivative_v_per_s=
            real.(im * discrete_frequency .* represented_voltage_phasor),
        terminal_current_a=real.(nodal_current_phasor[external]),
        branch_current_a=real.(branch_current_phasor),
        branch_voltage_v=real.(branch_voltage_phasor),
        internal_kcl_residual_a=maximum(
            abs,
            nodal_current_phasor[internal];
            init=0.0,
        ),
    )
end

function _transformer_network_runtime_state(
    preparation::TransformerApparatusPreparation,
    candidate::_TransformerNetworkCandidate,
)
    node_count = length(candidate.node_order)
    if preparation.initialization_mode === SinusoidalTransformerOperatingPoint
        initial = _transformer_sinusoidal_network_initial_state(
            preparation,
            candidate,
        )
    else
        represented_voltage = zeros(node_count)
        represented_voltage[candidate.external_node_indices] .=
            preparation.initial_node_voltage_v
        represented_voltage_derivative = zeros(node_count)
        represented_voltage_derivative[candidate.external_node_indices] .=
            _transformer_initial_discrete_derivative(
                preparation.initial_node_voltage_v,
                preparation.initial_node_voltage_derivative_v_per_s,
                preparation.initialization_mode,
                preparation.specification.settings,
            )
        branch_current = zeros(length(candidate.branch_current_a))
        branch_voltage = zeros(length(candidate.branch_voltage_v))
        nodal_current = candidate.conductance_s * represented_voltage .+
            candidate.capacitance_f * represented_voltage_derivative
        for group in candidate.branch_groups
            voltage = transpose(group.incidence) * represented_voltage
            branch_voltage[group.current_offset] .= voltage
            nodal_current .+= group.incidence * branch_current[group.current_offset]
        end
        initial = (
            represented_node_voltage_v=represented_voltage,
            represented_node_voltage_derivative_v_per_s=
                represented_voltage_derivative,
            terminal_current_a=
                nodal_current[candidate.external_node_indices],
            branch_current_a=branch_current,
            branch_voltage_v=branch_voltage,
            internal_kcl_residual_a=maximum(
                abs,
                nodal_current[candidate.internal_node_indices];
                init=0.0,
            ),
        )
    end
    capacitor_current = candidate.capacitance_f *
        initial.represented_node_voltage_derivative_v_per_s
    for group in candidate.branch_groups
        group.accepted_current_a .= initial.branch_current_a[group.current_offset]
        group.accepted_voltage_v .= initial.branch_voltage_v[group.current_offset]
    end
    stored_magnetic_energy = sum(
        0.5 * dot(
            group.accepted_current_a,
            group.inductance_h * group.accepted_current_a,
        ) for group in candidate.branch_groups
    )
    stored_electric_energy = 0.5 * dot(
        initial.represented_node_voltage_v,
        candidate.capacitance_f * initial.represented_node_voltage_v,
    )
    terminal_power = dot(
        preparation.initial_node_voltage_v,
        initial.terminal_current_a,
    )
    return TransformerNetworkRuntimeState(
        initial.represented_node_voltage_v,
        copy(preparation.initial_node_voltage_v),
        initial.terminal_current_a,
        initial.branch_current_a,
        initial.branch_voltage_v,
        capacitor_current,
        terminal_power,
        terminal_power,
        0.0,
        0.0,
        0.0,
        stored_magnetic_energy,
        stored_electric_energy,
        initial.internal_kcl_residual_a,
        preparation.initial_time_s,
        0,
    )
end

function transformer_apparatus_runtime(
    preparation::TransformerApparatusPreparation,
    terminal_nodes,
)
    nodes = Int.(terminal_nodes)
    length(nodes) == length(preparation.terminal_order) || throw(DimensionMismatch(
        "transformer runtime terminal nodes must match the preparation terminal order",
    ))
    all(>(0), nodes) && length(unique(nodes)) == length(nodes) ||
        throw(ArgumentError(
            "transformer runtime terminal nodes must be unique positive network indices",
        ))
    tier = preparation.specification.tier
    accepted_state, candidate = if tier in (
        HybridTransformerTier,
        MagneticEquivalentCircuitTier,
    )
        _transformer_magnetic_runtime_state(preparation),
        _transformer_magnetic_candidate(preparation)
    elseif tier in (
        LowFrequencyTerminalTier,
        BCTRANTerminalTier,
    )
        _transformer_terminal_runtime_state(preparation),
        _transformer_terminal_candidate(preparation)
    elseif tier === WidebandBlackBoxTier
        _transformer_wideband_runtime_state(preparation),
        _transformer_wideband_candidate(preparation)
    elseif tier === GreyBoxLadderTier
        network_candidate = _grey_box_candidate(preparation)
        _transformer_network_runtime_state(preparation, network_candidate), network_candidate
    else
        network_candidate = _white_box_candidate(preparation)
        _transformer_network_runtime_state(preparation, network_candidate), network_candidate
    end
    return TransformerApparatusRuntime(
        preparation,
        nodes,
        accepted_state,
        candidate,
        preparation.initial_time_s,
        preparation.specification.settings.timestep_s,
        :trapezoidal,
        false,
        0,
        0,
        0,
        0,
        _transformer_total_stored_energy(accepted_state),
        TransformerApparatusEventState(length(nodes)),
        zeros(length(nodes)),
        zeros(length(nodes)),
        zeros(length(nodes), length(nodes)),
    )
end

nonlinear_terminal_nodes(runtime::TransformerApparatusRuntime) = runtime.terminal_nodes
nonlinear_device_formulation(::TransformerApparatusRuntime) = PhysicalConstitutiveCurrent
nonlinear_device_provenance(runtime::TransformerApparatusRuntime) =
    first(runtime.preparation.specification.sources).provenance

function _prepare_transformer_terminal_candidate!(
    runtime::TransformerApparatusRuntime,
    candidate::_TransformerTerminalMatrixCandidate,
    state::TransformerTerminalMatrixRuntimeState,
    step_s::Float64,
)
    matrices = _transformer_terminal_matrices(runtime.preparation)
    model = runtime.preparation.specification.model
    backward_euler = runtime.companion_method === :backward_euler
    alpha = (backward_euler ? 1.0 : 2.0) / step_s
    companion_parameters_changed =
        candidate.prepared_companion_method !== runtime.companion_method ||
        candidate.prepared_step_s != step_s
    if companion_parameters_changed
        fill!(candidate.winding_loss_impedance_ohm, 0.0)
        if model isa HybridTransformerModel
            state_count = size(model.winding_loss_state_matrix_per_s, 1)
            if state_count > 0
                identity_state = Matrix{Float64}(I, state_count, state_count)
                state_matrix = backward_euler ?
                    identity_state .- step_s .* model.winding_loss_state_matrix_per_s :
                    alpha .* identity_state .- model.winding_loss_state_matrix_per_s
                condition_number = cond(state_matrix)
                isfinite(condition_number) && condition_number <= 1.0e14 ||
                    _transformer_refusal(
                        :ill_conditioned_winding_loss_map,
                        :prepare_step,
                        runtime.preparation.specification.id,
                        runtime.preparation.specification.tier,
                        "transformer winding-loss bilinear state map is ill-conditioned";
                        diagnostics=(condition_number=condition_number,),
                    )
                factorization = lu(state_matrix)
                if backward_euler
                    candidate.winding_loss_state_transition .=
                        factorization \ identity_state
                    candidate.winding_loss_endpoint_input .= factorization \
                        (step_s .* model.winding_loss_input_matrix)
                else
                    candidate.winding_loss_state_transition .= factorization \
                        (alpha .* identity_state .+ model.winding_loss_state_matrix_per_s)
                    candidate.winding_loss_endpoint_input .= factorization \
                        model.winding_loss_input_matrix
                end
            end
            candidate.winding_loss_impedance_ohm .=
                model.winding_loss_direct_ohm .+
                model.winding_loss_output_matrix_ohm_per_s *
                candidate.winding_loss_endpoint_input
        end
        candidate.companion_impedance_ohm .= matrices.resistance_ohm .+
            candidate.winding_loss_impedance_ohm .+
            alpha .* matrices.inductance_h
        condition_number = cond(candidate.companion_impedance_ohm)
        isfinite(condition_number) && condition_number <= 1.0e14 ||
            _transformer_refusal(
                :ill_conditioned_companion,
                :prepare_step,
                runtime.preparation.specification.id,
                runtime.preparation.specification.tier,
                "transformer terminal companion impedance is singular or ill-conditioned";
                diagnostics=(condition_number=condition_number,),
            )
        candidate.conductance_s .= inv(Symmetric(candidate.companion_impedance_ohm))
        candidate.prepared_step_s = step_s
        candidate.prepared_companion_method = runtime.companion_method
    end
    history_voltage = candidate.history_voltage_v
    mul!(history_voltage, matrices.inductance_h, state.coil_current_a, alpha, 0.0)
    if !backward_euler
        mul!(history_voltage, matrices.resistance_ohm, state.coil_current_a, -1.0, 1.0)
        history_voltage .+= state.coil_voltage_v
    end
    if model isa HybridTransformerModel
        if backward_euler
            mul!(
                candidate.winding_loss_state,
                candidate.winding_loss_state_transition,
                state.winding_loss_state,
            )
        else
            history_voltage .-= state.winding_loss_voltage_v
            mul!(
                candidate.winding_loss_state,
                candidate.winding_loss_state_transition,
                state.winding_loss_state,
            )
            mul!(
                candidate.winding_loss_state,
                candidate.winding_loss_endpoint_input,
                state.coil_current_a,
                1.0,
                1.0,
            )
        end
        mul!(
            candidate.coil_voltage_v,
            model.winding_loss_output_matrix_ohm_per_s,
            candidate.winding_loss_state,
        )
        history_voltage .-= candidate.coil_voltage_v
    end
    mul!(candidate.history_current_a, candidate.conductance_s, history_voltage)
    fill!(candidate.coil_voltage_v, 0.0)
    fill!(candidate.coil_current_a, 0.0)
    fill!(candidate.capacitor_current_a, 0.0)
    fill!(candidate.terminal_current_a, 0.0)
    fill!(candidate.terminal_jacobian_s, 0.0)
    fill!(candidate.winding_loss_state, 0.0)
    fill!(candidate.winding_loss_voltage_v, 0.0)
    return candidate
end

function _prepare_transformer_magnetic_candidate!(
    runtime::TransformerApparatusRuntime,
    candidate::_TransformerMagneticCandidate,
    state::TransformerMagneticRuntimeState,
    step_s::Float64,
)
    graph = _transformer_magnetic_model_parts(runtime.preparation).magnetic_graph
    backward_euler = runtime.companion_method === :backward_euler
    alpha = (backward_euler ? 1.0 : 2.0) / step_s
    candidate.electrical_impedance_ohm .= candidate.resistance_ohm .+
        alpha .* candidate.leakage_inductance_h
    candidate.electrical_history_voltage_v .= if backward_euler
        alpha .* candidate.leakage_inductance_h * state.coil_current_a .+
        alpha .* (transpose(graph.winding_turns) * state.branch_flux_wb)
    else
        state.coil_voltage_v .+
        (alpha .* candidate.leakage_inductance_h .- candidate.resistance_ohm) *
            state.coil_current_a .+
        alpha .* (transpose(graph.winding_turns) * state.branch_flux_wb)
    end
    model = runtime.preparation.specification.model
    if model isa HybridTransformerModel
        state_count = size(model.winding_loss_state_matrix_per_s, 1)
        if state_count > 0
            identity_state = Matrix{Float64}(I, state_count, state_count)
            state_matrix = backward_euler ?
                identity_state .- step_s .* model.winding_loss_state_matrix_per_s :
                alpha .* identity_state .- model.winding_loss_state_matrix_per_s
            condition_number = cond(state_matrix)
            isfinite(condition_number) && condition_number <= 1.0e14 ||
                _transformer_refusal(
                    :ill_conditioned_winding_loss_map,
                    :prepare_step,
                    runtime.preparation.specification.id,
                    runtime.preparation.specification.tier,
                    "transformer nonlinear magnetic winding-loss map is ill-conditioned";
                    diagnostics=(condition_number=condition_number,),
                )
            factorization = lu(state_matrix)
            if backward_euler
                candidate.winding_loss_state_transition .=
                    factorization \ identity_state
                candidate.winding_loss_endpoint_input .= factorization \
                    (step_s .* model.winding_loss_input_matrix)
            else
                candidate.winding_loss_state_transition .= factorization \
                    (alpha .* identity_state .+ model.winding_loss_state_matrix_per_s)
                candidate.winding_loss_endpoint_input .= factorization \
                    model.winding_loss_input_matrix
            end
        end
        candidate.electrical_impedance_ohm .+= model.winding_loss_direct_ohm .+
            model.winding_loss_output_matrix_ohm_per_s *
            candidate.winding_loss_endpoint_input
        if backward_euler
            candidate.electrical_history_voltage_v .-=
                model.winding_loss_output_matrix_ohm_per_s *
                (candidate.winding_loss_state_transition * state.winding_loss_state)
        else
            candidate.electrical_history_voltage_v .-= state.winding_loss_voltage_v
            candidate.electrical_history_voltage_v .-=
                model.winding_loss_output_matrix_ohm_per_s * (
                    candidate.winding_loss_state_transition * state.winding_loss_state .+
                    candidate.winding_loss_endpoint_input * state.coil_current_a
                )
        end
    end
    candidate.coil_current_a .= state.coil_current_a
    candidate.branch_flux_wb .= state.branch_flux_wb
    candidate.branch_mmf_drop_at .= state.branch_mmf_drop_at
    candidate.magnetic_node_potential_at .= state.magnetic_node_potential_at
    candidate.tellinen_state .= state.tellinen_state
    fill!(candidate.coil_voltage_v, 0.0)
    fill!(candidate.capacitor_current_a, 0.0)
    fill!(candidate.terminal_current_a, 0.0)
    fill!(candidate.terminal_jacobian_s, 0.0)
    fill!(candidate.winding_loss_state, 0.0)
    fill!(candidate.winding_loss_voltage_v, 0.0)
    candidate.local_iterations = 0
    candidate.electrical_residual_v = Inf
    candidate.magnetic_constitutive_residual_at = Inf
    candidate.magnetic_continuity_residual_wb = Inf
    return candidate
end

function _prepare_transformer_wideband_candidate!(
    runtime::TransformerApparatusRuntime,
    candidate::_TransformerWidebandCandidate,
    state::TransformerWidebandRuntimeState,
    step_s::Float64,
)
    model = runtime.preparation.specification.model::WidebandTransformerModel
    state_count = size(model.state_matrix_per_s, 1)
    backward_euler = runtime.companion_method === :backward_euler
    alpha = (backward_euler ? 1.0 : 2.0) / step_s
    identity_state = Matrix{Float64}(I, state_count, state_count)
    bilinear_matrix = alpha .* identity_state .- model.state_matrix_per_s
    condition_number = cond(bilinear_matrix)
    isfinite(condition_number) && condition_number <= 1.0e14 ||
        _transformer_refusal(
            :ill_conditioned_bilinear_map,
            :prepare_step,
            runtime.preparation.specification.id,
            runtime.preparation.specification.tier,
            "wideband transformer bilinear state map is ill-conditioned";
            diagnostics=(condition_number=condition_number,),
        )
    factorization = lu(bilinear_matrix)
    candidate.state_transition .= backward_euler ?
        factorization \ (alpha .* identity_state) :
        factorization \ (alpha .* identity_state .+ model.state_matrix_per_s)
    candidate.endpoint_input .= factorization \ model.input_matrix
    candidate.companion_admittance_s .= model.direct_admittance_s .+
        model.output_matrix_s_per_s * candidate.endpoint_input
    scale = max(maximum(abs, candidate.companion_admittance_s; init=0.0), 1.0)
    maximum(
        abs,
        candidate.companion_admittance_s - transpose(candidate.companion_admittance_s);
        init=0.0,
    ) <= 1.0e-8 * scale || _transformer_refusal(
        :nonreciprocal_discrete_companion,
        :prepare_step,
        runtime.preparation.specification.id,
        runtime.preparation.specification.tier,
        "wideband transformer discrete companion is not reciprocal",
    )
    minimum_eigenvalue = minimum(
        eigvals(Symmetric(0.5 .* (
            candidate.companion_admittance_s .+
            transpose(candidate.companion_admittance_s)
        )));
        init=Inf,
    )
    minimum_eigenvalue >= -1.0e-10 * scale || _transformer_refusal(
        :active_discrete_companion,
        :prepare_step,
        runtime.preparation.specification.id,
        runtime.preparation.specification.tier,
        "wideband transformer discrete companion violates passivity";
        diagnostics=(minimum_eigenvalue_s=minimum_eigenvalue,),
    )
    candidate.history_current_a .= model.output_matrix_s_per_s * (
        candidate.state_transition * state.rational_state .+
        (backward_euler ? zeros(length(state.rational_state)) :
         candidate.endpoint_input * state.terminal_voltage_v)
    )
    fill!(candidate.state, 0.0)
    fill!(candidate.terminal_current_a, 0.0)
    return candidate
end

function _prepare_transformer_network_candidate!(
    runtime::TransformerApparatusRuntime,
    candidate::_TransformerNetworkCandidate,
    state::TransformerNetworkRuntimeState,
    step_s::Float64,
)
    backward_euler = runtime.companion_method === :backward_euler
    alpha = (backward_euler ? 1.0 : 2.0) / step_s
    fill!(candidate.nodal_admittance_s, 0.0)
    fill!(candidate.nodal_history_current_a, 0.0)
    for group in candidate.branch_groups
        impedance = group.resistance_ohm .+
            alpha .* group.inductance_h
        condition_number = cond(impedance)
        isfinite(condition_number) && condition_number <= 1.0e14 ||
            _transformer_refusal(
                :ill_conditioned_ladder_branch,
                :prepare_step,
                runtime.preparation.specification.id,
                runtime.preparation.specification.tier,
                "transformer internal coupled branch companion is ill-conditioned";
                diagnostics=(condition_number=condition_number,),
            )
        group.conductance_s .= inv(Symmetric(impedance))
        history_voltage = backward_euler ?
            alpha .* group.inductance_h * group.accepted_current_a :
            group.accepted_voltage_v .+
                (alpha .* group.inductance_h .- group.resistance_ohm) *
                group.accepted_current_a
        group.history_current_a .= group.conductance_s * history_voltage
        candidate.nodal_admittance_s .+=
            group.incidence * group.conductance_s * transpose(group.incidence)
        candidate.nodal_history_current_a .+=
            group.incidence * group.history_current_a
    end
    capacitor_conductance = alpha .* candidate.capacitance_f
    candidate.nodal_admittance_s .+= candidate.conductance_s .+ capacitor_conductance
    candidate.nodal_history_current_a .+=
        -(backward_euler ? zeros(length(state.capacitor_current_a)) :
          state.capacitor_current_a) .-
        capacitor_conductance * state.represented_node_voltage_v
    for fault_id in sort!(
        collect(keys(runtime.event_state.active_internal_faults));
        by=String,
    )
        from_index, to_index, conductance_s =
            runtime.event_state.active_internal_faults[fault_id]
        candidate.nodal_admittance_s[from_index, from_index] += conductance_s
        candidate.nodal_admittance_s[to_index, to_index] += conductance_s
        candidate.nodal_admittance_s[from_index, to_index] -= conductance_s
        candidate.nodal_admittance_s[to_index, from_index] -= conductance_s
    end
    external = candidate.external_node_indices
    internal = candidate.internal_node_indices
    if isempty(internal)
        candidate.companion_admittance_s .=
            candidate.nodal_admittance_s[external, external]
        candidate.history_current_a .= candidate.nodal_history_current_a[external]
    else
        internal_matrix = candidate.nodal_admittance_s[internal, internal]
        condition_number = cond(internal_matrix)
        isfinite(condition_number) && condition_number <= 1.0e14 ||
            _transformer_refusal(
                :ill_conditioned_internal_network,
                :prepare_step,
                runtime.preparation.specification.id,
                runtime.preparation.specification.tier,
                "transformer represented internal network is singular or ill-conditioned";
                diagnostics=(condition_number=condition_number,),
            )
        internal_factor = lu(internal_matrix)
        internal_to_external = internal_factor \
            candidate.nodal_admittance_s[internal, external]
        internal_history = internal_factor \ candidate.nodal_history_current_a[internal]
        candidate.companion_admittance_s .=
            candidate.nodal_admittance_s[external, external] .-
            candidate.nodal_admittance_s[external, internal] * internal_to_external
        candidate.history_current_a .=
            candidate.nodal_history_current_a[external] .-
            candidate.nodal_admittance_s[external, internal] * internal_history
    end
    fill!(candidate.represented_node_voltage_v, 0.0)
    fill!(candidate.branch_current_a, 0.0)
    fill!(candidate.branch_voltage_v, 0.0)
    fill!(candidate.capacitor_current_a, 0.0)
    return candidate
end

function prepare_nonlinear_device_step!(
    runtime::TransformerApparatusRuntime,
    time_s::Float64,
    step_s::Float64,
    companion_method::Symbol,
)
    isfinite(time_s) && isfinite(step_s) && step_s > 0.0 || throw(ArgumentError(
        "transformer candidate time and step must be finite with a positive step",
    ))
    accepted_time = runtime.accepted_state.accepted_time_s
    expected_time = accepted_time + step_s
    abs(time_s - expected_time) <=
        64.0 * eps(Float64) * max(abs(time_s), abs(expected_time), 1.0) ||
        _transformer_refusal(
            :nonforward_candidate_time,
            :prepare_step,
            runtime.preparation.specification.id,
            runtime.preparation.specification.tier,
            "transformer candidate time must be exactly one configured step after the accepted state";
            diagnostics=(
                accepted_time_s=accepted_time,
                expected_time_s=expected_time,
                candidate_time_s=time_s,
            ),
        )
    normalized_companion_method = companion_method in (
        :trapezoidal,
        :TrapezoidalCompanion,
    ) ? :trapezoidal : companion_method in (
        :backward_euler,
        :BackwardEulerCompanion,
    ) ? :backward_euler : companion_method
    normalized_companion_method in (:trapezoidal, :backward_euler) || _transformer_refusal(
        :unsupported_companion_method,
        :prepare_step,
        runtime.preparation.specification.id,
        runtime.preparation.specification.tier,
        "transformer apparatus admits trapezoidal execution and event-localization backward Euler only";
        diagnostics=(requested_method=companion_method,),
    )
    configured_step = runtime.preparation.specification.settings.timestep_s
    step_tolerance = 64.0 * eps(Float64) * max(step_s, configured_step)
    step_valid = normalized_companion_method === :trapezoidal ?
        abs(step_s - configured_step) <= step_tolerance :
        step_s <= configured_step + step_tolerance
    step_valid || _transformer_refusal(
            :timestep_identity_mismatch,
            :prepare_step,
            runtime.preparation.specification.id,
            runtime.preparation.specification.tier,
            "candidate timestep is incompatible with the transformer preparation";
            diagnostics=(configured_step_s=configured_step, candidate_step_s=step_s),
        )
    previous_companion_method = runtime.companion_method
    runtime.companion_method = normalized_companion_method
    try
        if runtime.candidate isa _TransformerTerminalMatrixCandidate
            _prepare_transformer_terminal_candidate!(
                runtime,
                runtime.candidate,
                runtime.accepted_state::TransformerTerminalMatrixRuntimeState,
                step_s,
            )
        elseif runtime.candidate isa _TransformerMagneticCandidate
            _prepare_transformer_magnetic_candidate!(
                runtime,
                runtime.candidate,
                runtime.accepted_state::TransformerMagneticRuntimeState,
                step_s,
            )
        elseif runtime.candidate isa _TransformerWidebandCandidate
            _prepare_transformer_wideband_candidate!(
                runtime,
                runtime.candidate,
                runtime.accepted_state::TransformerWidebandRuntimeState,
                step_s,
            )
        else
            _prepare_transformer_network_candidate!(
                runtime,
                runtime.candidate,
                runtime.accepted_state::TransformerNetworkRuntimeState,
                step_s,
            )
        end
    catch error
        runtime.companion_method = previous_companion_method
        rethrow(error)
    end
    runtime.candidate_time_s = time_s
    runtime.candidate_step_s = step_s
    runtime.prepared = true
    runtime.preparation_count += 1
    return nothing
end

function _transformer_terminal_trial!(
    runtime::TransformerApparatusRuntime,
    terminal_voltage_v::AbstractVector{Float64},
)
    candidate = runtime.candidate::_TransformerTerminalMatrixCandidate
    state = runtime.accepted_state::TransformerTerminalMatrixRuntimeState
    matrices = _transformer_terminal_matrices(runtime.preparation)
    connection = runtime.preparation.specification.connection
    transformer_coil_voltages!(candidate.coil_voltage_v, connection, terminal_voltage_v)
    mul!(candidate.coil_current_a, candidate.conductance_s, candidate.coil_voltage_v)
    candidate.coil_current_a .+= candidate.history_current_a
    model = runtime.preparation.specification.model
    if model isa HybridTransformerModel
        candidate.winding_loss_state .= if runtime.companion_method === :backward_euler
            candidate.winding_loss_state_transition * state.winding_loss_state .+
            candidate.winding_loss_endpoint_input * candidate.coil_current_a
        else
            candidate.winding_loss_state_transition * state.winding_loss_state .+
            candidate.winding_loss_endpoint_input * (
                candidate.coil_current_a .+ state.coil_current_a
            )
        end
        candidate.winding_loss_voltage_v .=
            model.winding_loss_output_matrix_ohm_per_s *
            candidate.winding_loss_state .+
            model.winding_loss_direct_ohm * candidate.coil_current_a
    end
    alpha = (runtime.companion_method === :backward_euler ? 1.0 : 2.0) /
        runtime.candidate_step_s
    candidate.coil_current_workspace_a .=
        candidate.coil_voltage_v .- state.coil_voltage_v
    mul!(
        candidate.capacitor_current_a,
        matrices.capacitance_f,
        candidate.coil_current_workspace_a,
        alpha,
        0.0,
    )
    runtime.companion_method === :backward_euler ||
        (candidate.capacitor_current_a .-= state.capacitor_current_a)
    candidate.coil_current_workspace_a .=
        candidate.coil_current_a .+ candidate.capacitor_current_a
    mul!(
        candidate.coil_current_workspace_a,
        matrices.conductance_s,
        candidate.coil_voltage_v,
        1.0,
        1.0,
    )
    transformer_terminal_currents!(
        candidate.terminal_current_a,
        connection,
        candidate.coil_current_workspace_a,
    )
    candidate.coil_jacobian_s .= candidate.conductance_s .+
        alpha .* matrices.capacitance_f .+
        matrices.conductance_s
    mul!(
        candidate.terminal_coil_jacobian_s,
        connection.incidence,
        candidate.coil_jacobian_s,
    )
    mul!(
        candidate.terminal_jacobian_s,
        candidate.terminal_coil_jacobian_s,
        transpose(connection.incidence),
    )
    return candidate.terminal_current_a, candidate.terminal_jacobian_s
end

function _transformer_magnetic_trial!(
    runtime::TransformerApparatusRuntime,
    terminal_voltage_v::AbstractVector{Float64},
)
    candidate = runtime.candidate::_TransformerMagneticCandidate
    state = runtime.accepted_state::TransformerMagneticRuntimeState
    parts = _transformer_magnetic_model_parts(runtime.preparation)
    graph = parts.magnetic_graph
    connection = runtime.preparation.specification.connection
    settings = runtime.preparation.specification.settings
    transformer_coil_voltages!(
        candidate.coil_voltage_v,
        connection,
        terminal_voltage_v,
    )
    coil_count = length(candidate.coil_current_a)
    branch_count = length(candidate.branch_flux_wb)
    magnetic_node_count = length(candidate.magnetic_node_potential_at)
    alpha = (runtime.companion_method === :backward_euler ? 1.0 : 2.0) /
        runtime.candidate_step_s
    local_system = zeros(
        Float64,
        coil_count + branch_count + magnetic_node_count,
        coil_count + branch_count + magnetic_node_count,
    )
    converged = false
    for iteration in 1:settings.maximum_local_nonlinear_iterations
        model = runtime.preparation.specification.model
        if model isa HybridTransformerModel
            candidate.winding_loss_state .= if runtime.companion_method === :backward_euler
                candidate.winding_loss_state_transition * state.winding_loss_state .+
                candidate.winding_loss_endpoint_input * candidate.coil_current_a
            else
                candidate.winding_loss_state_transition * state.winding_loss_state .+
                candidate.winding_loss_endpoint_input * (
                    candidate.coil_current_a .+ state.coil_current_a
                )
            end
            candidate.winding_loss_voltage_v .=
                model.winding_loss_output_matrix_ohm_per_s *
                candidate.winding_loss_state .+
                model.winding_loss_direct_ohm * candidate.coil_current_a
        end
        constitutive_residual, continuity_residual = _transformer_magnetic_residual!(
            candidate.branch_mmf_drop_at,
            candidate.branch_differential_reluctance_at_per_wb,
            candidate.tellinen_state,
            graph,
            candidate.branch_flux_wb,
            candidate.magnetic_node_potential_at,
            candidate.coil_current_a,
            state.tellinen_state,
            previous_branch_flux_wb=state.branch_flux_wb,
            previous_classical_dynamic_mmf_at=
                state.last_classical_dynamic_mmf_at,
            previous_excess_dynamic_mmf_at=state.last_excess_dynamic_mmf_at,
            step_s=runtime.candidate_step_s,
            branch_classical_dynamic_mmf_at=
                candidate.branch_classical_dynamic_mmf_at,
            branch_excess_dynamic_mmf_at=candidate.branch_excess_dynamic_mmf_at,
            branch_classical_loss_power_w=candidate.branch_classical_loss_power_w,
            branch_excess_loss_power_w=candidate.branch_excess_loss_power_w,
            companion_method=runtime.companion_method,
        )
        electrical_residual = candidate.electrical_impedance_ohm *
            candidate.coil_current_a .+
            alpha .* (transpose(graph.winding_turns) * candidate.branch_flux_wb) .-
            candidate.coil_voltage_v .-
            candidate.electrical_history_voltage_v
        electrical_scale = max(
            maximum(abs, candidate.coil_voltage_v; init=0.0),
            maximum(abs, candidate.electrical_history_voltage_v; init=0.0),
            1.0,
        )
        magnetic_scale = max(
            maximum(abs, graph.winding_turns * candidate.coil_current_a; init=0.0),
            maximum(abs, candidate.branch_mmf_drop_at; init=0.0),
            1.0,
        )
        flux_scale = max(maximum(abs, candidate.branch_flux_wb; init=0.0), 1.0e-12)
        candidate.electrical_residual_v = maximum(abs, electrical_residual; init=0.0)
        candidate.magnetic_constitutive_residual_at =
            maximum(abs, constitutive_residual; init=0.0)
        candidate.magnetic_continuity_residual_wb =
            maximum(abs, continuity_residual; init=0.0)
        candidate.local_iterations = iteration
        converged =
            candidate.electrical_residual_v <=
                settings.nonlinear_residual_relative_tolerance * electrical_scale &&
            candidate.magnetic_constitutive_residual_at <=
                settings.nonlinear_residual_relative_tolerance * magnetic_scale &&
            candidate.magnetic_continuity_residual_wb <=
                settings.magnetic_continuity_absolute_tolerance_wb +
                settings.nonlinear_residual_relative_tolerance * flux_scale
        local_system .= 0.0
        local_system[1:coil_count, 1:coil_count] .=
            candidate.electrical_impedance_ohm
        local_system[
            1:coil_count,
            (coil_count + 1):(coil_count + branch_count),
        ] .= alpha .* transpose(graph.winding_turns)
        local_system[
            (coil_count + 1):(coil_count + branch_count),
            1:coil_count,
        ] .= -graph.winding_turns
        local_system[
            (coil_count + 1):(coil_count + branch_count),
            (coil_count + 1):(coil_count + branch_count),
        ] .= Diagonal(candidate.branch_differential_reluctance_at_per_wb)
        local_system[
            (coil_count + 1):(coil_count + branch_count),
            (coil_count + branch_count + 1):end,
        ] .= transpose(graph.incidence)
        local_system[
            (coil_count + branch_count + 1):end,
            (coil_count + 1):(coil_count + branch_count),
        ] .= graph.incidence
        converged && break
        correction = -(local_system \ vcat(
            electrical_residual,
            constitutive_residual,
            continuity_residual,
        ))
        all(isfinite, correction) || throw(ArgumentError(
            "transformer nonlinear magnetic correction became nonfinite",
        ))
        candidate.coil_current_a .+= correction[1:coil_count]
        candidate.branch_flux_wb .+= correction[
            (coil_count + 1):(coil_count + branch_count)
        ]
        candidate.magnetic_node_potential_at .+= correction[
            (coil_count + branch_count + 1):end
        ]
    end
    converged || _transformer_refusal(
        :nonlinear_magnetic_nonconvergence,
        :trial,
        runtime.preparation.specification.id,
        runtime.preparation.specification.tier,
        "transformer coupled electrical and magnetic endpoint did not converge";
        diagnostics=(
            iterations=candidate.local_iterations,
            electrical_residual_v=candidate.electrical_residual_v,
            magnetic_constitutive_residual_at=
                candidate.magnetic_constitutive_residual_at,
            magnetic_continuity_residual_wb=
                candidate.magnetic_continuity_residual_wb,
        ),
    )
    voltage_right_hand_side = zeros(
        Float64,
        coil_count + branch_count + magnetic_node_count,
        coil_count,
    )
    voltage_right_hand_side[1:coil_count, :] .= Matrix{Float64}(
        I,
        coil_count,
        coil_count,
    )
    coil_current_jacobian = (local_system \ voltage_right_hand_side)[
        1:coil_count,
        :,
    ]
    candidate.capacitor_current_a .= alpha .* candidate.capacitance_f *
        (candidate.coil_voltage_v .- state.coil_voltage_v)
    runtime.companion_method === :backward_euler ||
        (candidate.capacitor_current_a .-= state.capacitor_current_a)
    total_coil_current = candidate.coil_current_a .+
        candidate.capacitor_current_a .+
        candidate.conductance_s * candidate.coil_voltage_v
    transformer_terminal_currents!(
        candidate.terminal_current_a,
        connection,
        total_coil_current,
    )
    coil_terminal_jacobian = coil_current_jacobian .+
        alpha .* candidate.capacitance_f .+ candidate.conductance_s
    candidate.terminal_jacobian_s .= connection.incidence *
        coil_terminal_jacobian * transpose(connection.incidence)
    return candidate.terminal_current_a, candidate.terminal_jacobian_s
end

function _transformer_wideband_trial!(
    runtime::TransformerApparatusRuntime,
    terminal_voltage_v::AbstractVector{Float64},
)
    candidate = runtime.candidate::_TransformerWidebandCandidate
    state = runtime.accepted_state::TransformerWidebandRuntimeState
    model = runtime.preparation.specification.model::WidebandTransformerModel
    model = runtime.preparation.specification.model::WidebandTransformerModel
    candidate.state .= candidate.state_transition * state.rational_state .+
        candidate.endpoint_input * (
            runtime.companion_method === :backward_euler ?
                terminal_voltage_v : terminal_voltage_v .+ state.terminal_voltage_v
        )
    candidate.terminal_current_a .=
        model.output_matrix_s_per_s * candidate.state .+
        model.direct_admittance_s * terminal_voltage_v
    return candidate.terminal_current_a, candidate.companion_admittance_s
end

function _transformer_network_trial!(
    runtime::TransformerApparatusRuntime,
    terminal_voltage_v::AbstractVector{Float64},
)
    candidate = runtime.candidate::_TransformerNetworkCandidate
    state = runtime.accepted_state::TransformerNetworkRuntimeState
    external = candidate.external_node_indices
    internal = candidate.internal_node_indices
    candidate.represented_node_voltage_v[external] .= terminal_voltage_v
    if !isempty(internal)
        internal_matrix = candidate.nodal_admittance_s[internal, internal]
        right_hand_side = -(
            candidate.nodal_admittance_s[internal, external] * terminal_voltage_v .+
            candidate.nodal_history_current_a[internal]
        )
        candidate.represented_node_voltage_v[internal] .= internal_matrix \ right_hand_side
    end
    terminal_current = candidate.companion_admittance_s * terminal_voltage_v .+
        candidate.history_current_a
    current_cursor = 1
    for group in candidate.branch_groups
        branch_voltage = transpose(group.incidence) * candidate.represented_node_voltage_v
        branch_current = group.conductance_s * branch_voltage .+ group.history_current_a
        candidate.branch_voltage_v[group.current_offset] .= branch_voltage
        candidate.branch_current_a[group.current_offset] .= branch_current
        current_cursor += length(branch_current)
    end
    alpha = (runtime.companion_method === :backward_euler ? 1.0 : 2.0) /
        runtime.candidate_step_s
    candidate.capacitor_current_a .= alpha .* candidate.capacitance_f *
        (candidate.represented_node_voltage_v .- state.represented_node_voltage_v)
    runtime.companion_method === :backward_euler ||
        (candidate.capacitor_current_a .-= state.capacitor_current_a)
    return terminal_current, candidate.companion_admittance_s
end

function _transformer_raw_trial!(
    runtime::TransformerApparatusRuntime,
    terminal_voltage_v::AbstractVector{Float64},
)
    return if runtime.candidate isa _TransformerTerminalMatrixCandidate
        _transformer_terminal_trial!(runtime, terminal_voltage_v)
    elseif runtime.candidate isa _TransformerMagneticCandidate
        _transformer_magnetic_trial!(runtime, terminal_voltage_v)
    elseif runtime.candidate isa _TransformerWidebandCandidate
        _transformer_wideband_trial!(runtime, terminal_voltage_v)
    else
        _transformer_network_trial!(runtime, terminal_voltage_v)
    end
end

function nonlinear_current_jacobian!(
    terminal_current_a::AbstractVector{Float64},
    terminal_jacobian_s::AbstractMatrix{Float64},
    runtime::TransformerApparatusRuntime,
    terminal_voltage_v::AbstractVector{Float64},
    time_s::Float64,
)
    runtime.prepared || throw(ArgumentError(
        "transformer apparatus step must be prepared before trial evaluation",
    ))
    length(terminal_voltage_v) == length(runtime.terminal_nodes) ||
        throw(DimensionMismatch("transformer trial terminal voltage count is incompatible"))
    length(terminal_current_a) >= length(runtime.terminal_nodes) &&
        size(terminal_jacobian_s, 1) >= length(runtime.terminal_nodes) &&
        size(terminal_jacobian_s, 2) >= length(runtime.terminal_nodes) ||
        throw(DimensionMismatch("transformer nonlinear workspaces are too small"))
    all(isfinite, terminal_voltage_v) && isfinite(time_s) || throw(ArgumentError(
        "transformer trial voltage and time must be finite",
    ))
    current, jacobian = _transformer_event_trial!(runtime, terminal_voltage_v)
    terminal_count = length(runtime.terminal_nodes)
    terminal_current_a[1:terminal_count] .= current
    terminal_jacobian_s[1:terminal_count, 1:terminal_count] .= jacobian
    runtime.trial_evaluation_count += 1
    return nothing
end

function _transformer_endpoint_inner_product(
    previous_left,
    accepted_left,
    previous_right,
    accepted_right,
    backward_euler::Bool,
)
    if backward_euler
        return dot(accepted_left, accepted_right)
    end
    result = 0.0
    @inbounds for index in eachindex(
        previous_left,
        accepted_left,
        previous_right,
        accepted_right,
    )
        result += 0.25 *
            (previous_left[index] + accepted_left[index]) *
            (previous_right[index] + accepted_right[index])
    end
    return result
end

function _transformer_quadratic_form!(workspace, matrix, values)
    mul!(workspace, matrix, values)
    return dot(values, workspace)
end

function _transformer_endpoint_quadratic_form!(
    average_workspace,
    product_workspace,
    matrix,
    previous_values,
    accepted_values,
    backward_euler::Bool,
)
    if backward_euler
        average_workspace .= accepted_values
    else
        average_workspace .= 0.5 .* (previous_values .+ accepted_values)
    end
    return _transformer_quadratic_form!(
        product_workspace,
        matrix,
        average_workspace,
    )
end

function _accept_transformer_energy!(
    state,
    previous_voltage,
    previous_current,
    voltage,
    current,
    step_s,
    settings,
    companion_method::Symbol=:trapezoidal,
)
    power = dot(voltage, current)
    energy_increment = step_s * _transformer_endpoint_inner_product(
        previous_voltage,
        voltage,
        previous_current,
        current,
        companion_method === :backward_euler,
    )
    state.supplied_energy_j += energy_increment
    state.previous_terminal_power_w = power
    state.terminal_power_w = power
    state.accepted_step_count += 1
    state.accepted_time_s += step_s
    all(isfinite, (
        state.supplied_energy_j,
        state.terminal_power_w,
        state.accepted_time_s,
    )) || throw(ArgumentError("transformer accepted energy state became nonfinite"))
    return nothing
end

function _accept_transformer_linear_terminal_state!(
    runtime::TransformerApparatusRuntime,
    candidate::_TransformerTerminalMatrixCandidate,
    state::TransformerTerminalMatrixRuntimeState,
    matrices::TransformerTerminalMatrices,
    voltage,
    current,
)
    backward_euler = runtime.companion_method === :backward_euler
    magnetic_energy = 0.5 * _transformer_quadratic_form!(
        candidate.history_voltage_v,
        matrices.inductance_h,
        candidate.coil_current_a,
    )
    electric_energy = 0.5 * _transformer_quadratic_form!(
        candidate.history_voltage_v,
        matrices.capacitance_f,
        candidate.coil_voltage_v,
    )
    winding_loss_increment = runtime.candidate_step_s *
        _transformer_endpoint_quadratic_form!(
            candidate.history_voltage_v,
            candidate.coil_current_workspace_a,
            matrices.resistance_ohm,
            state.coil_current_a,
            candidate.coil_current_a,
            backward_euler,
        )
    dielectric_loss_increment = runtime.candidate_step_s *
        _transformer_endpoint_quadratic_form!(
            candidate.history_voltage_v,
            candidate.coil_current_workspace_a,
            matrices.conductance_s,
            state.coil_voltage_v,
            candidate.coil_voltage_v,
            backward_euler,
        )
    _accept_transformer_energy!(
        state,
        state.terminal_voltage_v,
        state.terminal_current_a,
        voltage,
        current,
        runtime.candidate_step_s,
        runtime.preparation.specification.settings,
        runtime.companion_method,
    )
    state.coil_current_a .= candidate.coil_current_a
    state.capacitor_current_a .= candidate.capacitor_current_a
    state.coil_voltage_v .= candidate.coil_voltage_v
    state.terminal_voltage_v .= voltage
    state.terminal_current_a .= current
    state.stored_magnetic_energy_j = magnetic_energy
    state.stored_electric_energy_j = electric_energy
    state.winding_loss_energy_j += winding_loss_increment
    state.dielectric_loss_energy_j += dielectric_loss_increment
    return nothing
end

function _accept_transformer_terminal_state!(
    runtime::TransformerApparatusRuntime,
    voltage,
)
    candidate = runtime.candidate::_TransformerTerminalMatrixCandidate
    state = runtime.accepted_state::TransformerTerminalMatrixRuntimeState
    matrices = _transformer_terminal_matrices(runtime.preparation)
    current, _ = _transformer_terminal_trial!(runtime, voltage)
    model = runtime.preparation.specification.model
    if !(model isa HybridTransformerModel)
        return _accept_transformer_linear_terminal_state!(
            runtime,
            candidate,
            state,
            matrices,
            voltage,
            current,
        )
    end
    previous_coil_current = copy(state.coil_current_a)
    previous_coil_voltage = copy(state.coil_voltage_v)
    previous_winding_loss_state = copy(state.winding_loss_state)
    previous_winding_loss_voltage = copy(state.winding_loss_voltage_v)
    previous_terminal_voltage = copy(state.terminal_voltage_v)
    previous_terminal_current = copy(state.terminal_current_a)
    backward_euler = runtime.companion_method === :backward_euler
    state.coil_current_a .= candidate.coil_current_a
    state.capacitor_current_a .= candidate.capacitor_current_a
    state.coil_voltage_v .= candidate.coil_voltage_v
    state.winding_loss_state .= candidate.winding_loss_state
    state.winding_loss_voltage_v .= candidate.winding_loss_voltage_v
    state.terminal_voltage_v .= voltage
    state.terminal_current_a .= current
    state.stored_magnetic_energy_j =
        0.5 * dot(state.coil_current_a, matrices.inductance_h * state.coil_current_a)
    state.stored_electric_energy_j =
        0.5 * dot(state.coil_voltage_v, matrices.capacitance_f * state.coil_voltage_v)
    step_s = runtime.candidate_step_s
    state.stored_frequency_dependent_winding_energy_j = 0.5 * dot(
        state.winding_loss_state,
        model.winding_loss_storage_matrix_j * state.winding_loss_state,
    )
    average_winding_loss_state = backward_euler ?
        state.winding_loss_state :
        0.5 .* (previous_winding_loss_state .+ state.winding_loss_state)
    average_winding_loss_current = backward_euler ?
        state.coil_current_a :
        0.5 .* (previous_coil_current .+ state.coil_current_a)
    state.frequency_dependent_winding_loss_energy_j += step_s * (
        0.5 * dot(
            average_winding_loss_state,
            model.winding_loss_dissipation_matrix_j_per_s *
            average_winding_loss_state,
        ) +
        dot(
            average_winding_loss_current,
            model.winding_loss_direct_ohm * average_winding_loss_current,
        )
    )
    average_winding_loss_voltage = backward_euler ?
        state.winding_loss_voltage_v :
        0.5 .* (previous_winding_loss_voltage .+ state.winding_loss_voltage_v)
    winding_loss_supplied_energy = step_s * dot(
        average_winding_loss_current,
        average_winding_loss_voltage,
    )
    winding_loss_balance = winding_loss_supplied_energy -
        (
            state.stored_frequency_dependent_winding_energy_j -
            0.5 * dot(
                previous_winding_loss_state,
                model.winding_loss_storage_matrix_j *
                previous_winding_loss_state,
            )
        ) -
        step_s * (
            0.5 * dot(
                average_winding_loss_state,
                model.winding_loss_dissipation_matrix_j_per_s *
                average_winding_loss_state,
            ) +
            dot(
                average_winding_loss_current,
                model.winding_loss_direct_ohm *
                average_winding_loss_current,
            )
        )
    backward_euler || abs(winding_loss_balance) <=
        runtime.preparation.specification.settings.energy_absolute_tolerance_j +
        512.0 * eps(Float64) * max(abs(winding_loss_supplied_energy), 1.0) ||
        throw(ArgumentError(
            "transformer winding-loss storage and dissipation balance failed",
        ))
    average_coil_current = backward_euler ? state.coil_current_a :
        0.5 .* (previous_coil_current .+ state.coil_current_a)
    average_coil_voltage = backward_euler ? state.coil_voltage_v :
        0.5 .* (previous_coil_voltage .+ state.coil_voltage_v)
    state.winding_loss_energy_j += step_s * dot(
        average_coil_current,
        matrices.resistance_ohm * average_coil_current,
    )
    state.dielectric_loss_energy_j += step_s * dot(
        average_coil_voltage,
        matrices.conductance_s * average_coil_voltage,
    )
    _accept_transformer_energy!(
        state,
        previous_terminal_voltage,
        previous_terminal_current,
        voltage,
        current,
        step_s,
        runtime.preparation.specification.settings,
        runtime.companion_method,
    )
    return nothing
end

function _accept_transformer_wideband_state!(
    runtime::TransformerApparatusRuntime,
    voltage,
)
    candidate = runtime.candidate::_TransformerWidebandCandidate
    state = runtime.accepted_state::TransformerWidebandRuntimeState
    model = runtime.preparation.specification.model::WidebandTransformerModel
    current, _ = _transformer_wideband_trial!(runtime, voltage)
    previous_rational_state = copy(state.rational_state)
    previous_terminal_voltage = copy(state.terminal_voltage_v)
    previous_terminal_current = copy(state.terminal_current_a)
    backward_euler = runtime.companion_method === :backward_euler
    state.rational_state .= candidate.state
    state.terminal_voltage_v .= voltage
    state.terminal_current_a .= current
    state.maximum_state_magnitude = max(
        state.maximum_state_magnitude,
        maximum(abs, state.rational_state; init=0.0),
    )
    previous_stored_energy = state.stored_energy_j
    state.stored_energy_j = 0.5 * dot(
        state.rational_state,
        model.storage_matrix_j * state.rational_state,
    )
    average_rational_state = backward_euler ? state.rational_state :
        0.5 .* (previous_rational_state .+ state.rational_state)
    average_terminal_voltage = backward_euler ? state.terminal_voltage_v :
        0.5 .* (previous_terminal_voltage .+ state.terminal_voltage_v)
    loss_increment = runtime.candidate_step_s * (
        0.5 * dot(
            average_rational_state,
            model.dissipation_matrix_j_per_s * average_rational_state,
        ) +
        dot(
            average_terminal_voltage,
            model.direct_admittance_s * average_terminal_voltage,
        )
    )
    state.dissipated_energy_j += loss_increment
    previous_supplied_energy = state.supplied_energy_j
    _accept_transformer_energy!(
        state,
        previous_terminal_voltage,
        previous_terminal_current,
        voltage,
        current,
        runtime.candidate_step_s,
        runtime.preparation.specification.settings,
        runtime.companion_method,
    )
    energy_residual = (state.supplied_energy_j - previous_supplied_energy) -
        (state.stored_energy_j - previous_stored_energy) - loss_increment
    state.maximum_energy_residual_j = max(
        state.maximum_energy_residual_j,
        abs(energy_residual),
    )
    backward_euler || abs(energy_residual) <=
        runtime.preparation.specification.settings.energy_absolute_tolerance_j +
        512.0 * eps(Float64) * max(abs(state.supplied_energy_j), 1.0) ||
        throw(ArgumentError("wideband transformer discrete energy balance failed"))
    state.minimum_supplied_energy_j = min(
        state.minimum_supplied_energy_j,
        state.supplied_energy_j,
    )
    return nothing
end

function _accept_transformer_magnetic_state!(
    runtime::TransformerApparatusRuntime,
    voltage,
)
    candidate = runtime.candidate::_TransformerMagneticCandidate
    state = runtime.accepted_state::TransformerMagneticRuntimeState
    parts = _transformer_magnetic_model_parts(runtime.preparation)
    graph = parts.magnetic_graph
    current, _ = _transformer_magnetic_trial!(runtime, voltage)
    previous_coil_current = copy(state.coil_current_a)
    previous_coil_voltage = copy(state.coil_voltage_v)
    previous_branch_flux = copy(state.branch_flux_wb)
    previous_tellinen_state = copy(state.tellinen_state)
    previous_winding_loss_state = copy(state.winding_loss_state)
    previous_winding_loss_voltage = copy(state.winding_loss_voltage_v)
    previous_terminal_voltage = copy(state.terminal_voltage_v)
    previous_terminal_current = copy(state.terminal_current_a)
    backward_euler = runtime.companion_method === :backward_euler
    state.coil_current_a .= candidate.coil_current_a
    state.capacitor_current_a .= candidate.capacitor_current_a
    state.coil_voltage_v .= candidate.coil_voltage_v
    state.branch_flux_wb .= candidate.branch_flux_wb
    state.branch_mmf_drop_at .= candidate.branch_mmf_drop_at
    state.magnetic_node_potential_at .= candidate.magnetic_node_potential_at
    state.tellinen_state .= candidate.tellinen_state
    state.winding_loss_state .= candidate.winding_loss_state
    state.winding_loss_voltage_v .= candidate.winding_loss_voltage_v
    state.terminal_voltage_v .= voltage
    state.terminal_current_a .= current
    state.stored_leakage_energy_j = 0.5 * dot(
        state.coil_current_a,
        candidate.leakage_inductance_h * state.coil_current_a,
    )
    state.stored_electric_energy_j = 0.5 * dot(
        state.coil_voltage_v,
        candidate.capacitance_f * state.coil_voltage_v,
    )
    state.stored_magnetic_energy_j = 0.0
    hysteresis_increment = 0.0
    for branch_index in eachindex(graph.branches)
        trial = _transformer_branch_magnetic_trial(
            graph,
            branch_index,
            state.branch_flux_wb[branch_index],
            state.tellinen_state[branch_index],
        )
        state.stored_magnetic_energy_j += trial.stored_energy_j
        material = graph.materials[graph.branches[branch_index].material_index]
        material isa TellinenTransformerMagneticMaterial || continue
        previous_material_state = previous_tellinen_state[branch_index]
        accepted_material_state = state.tellinen_state[branch_index]
        previous_material_state === nothing && throw(ArgumentError(
            "transformer accepted Tellinen state is missing",
        ))
        accepted_material_state === nothing && throw(ArgumentError(
            "transformer candidate Tellinen state is missing",
        ))
        branch = graph.branches[branch_index]
        core_volume = (branch.length_m - branch.air_gap_length_m) *
            branch.cross_section_m2
        integration_field_strength = backward_euler ?
            accepted_material_state.field_strength_a_per_m :
            0.5 * (
                previous_material_state.field_strength_a_per_m +
                accepted_material_state.field_strength_a_per_m
            )
        magnetic_work = integration_field_strength * (
            accepted_material_state.flux_density_t -
            previous_material_state.flux_density_t
        ) * core_volume
        stored_change = (
            _transformer_magnetic_energy_density(
                material,
                accepted_material_state.flux_density_t,
            ) -
            _transformer_magnetic_energy_density(
                material,
                previous_material_state.flux_density_t,
            )
        ) * core_volume
        hysteresis_increment += magnetic_work - stored_change
    end
    state.hysteresis_loss_energy_j += hysteresis_increment
    step_s = runtime.candidate_step_s
    state.last_classical_dynamic_mmf_at .=
        candidate.branch_classical_dynamic_mmf_at
    state.last_excess_dynamic_mmf_at .= candidate.branch_excess_dynamic_mmf_at
    classical_dynamic_loss_increment = step_s *
        sum(candidate.branch_classical_loss_power_w)
    excess_dynamic_loss_increment = step_s * sum(candidate.branch_excess_loss_power_w)
    loss_roundoff = 512.0 * eps(Float64) * max(
        abs(classical_dynamic_loss_increment),
        abs(excess_dynamic_loss_increment),
        1.0,
    )
    classical_dynamic_loss_increment >= -loss_roundoff || throw(ArgumentError(
        "transformer classical eddy core loss became active",
    ))
    excess_dynamic_loss_increment >= -loss_roundoff || throw(ArgumentError(
        "transformer excess core loss became active",
    ))
    state.classical_eddy_core_loss_energy_j +=
        max(0.0, classical_dynamic_loss_increment)
    state.excess_core_loss_energy_j += max(0.0, excess_dynamic_loss_increment)
    average_coil_current = backward_euler ? state.coil_current_a :
        0.5 .* (previous_coil_current .+ state.coil_current_a)
    average_coil_voltage = backward_euler ? state.coil_voltage_v :
        0.5 .* (previous_coil_voltage .+ state.coil_voltage_v)
    state.winding_loss_energy_j += step_s * dot(
        average_coil_current,
        candidate.resistance_ohm * average_coil_current,
    )
    state.dielectric_loss_energy_j += step_s * dot(
        average_coil_voltage,
        candidate.conductance_s * average_coil_voltage,
    )
    model = runtime.preparation.specification.model
    if model isa HybridTransformerModel
        state.stored_frequency_dependent_winding_energy_j = 0.5 * dot(
            state.winding_loss_state,
            model.winding_loss_storage_matrix_j * state.winding_loss_state,
        )
        average_winding_loss_state = backward_euler ? state.winding_loss_state :
            0.5 .* (previous_winding_loss_state .+ state.winding_loss_state)
        state.frequency_dependent_winding_loss_energy_j += step_s * (
            0.5 * dot(
                average_winding_loss_state,
                model.winding_loss_dissipation_matrix_j_per_s *
                average_winding_loss_state,
            ) +
            dot(
                average_coil_current,
                model.winding_loss_direct_ohm * average_coil_current,
            )
        )
        average_winding_loss_voltage = backward_euler ? state.winding_loss_voltage_v :
            0.5 .* (previous_winding_loss_voltage .+ state.winding_loss_voltage_v)
        winding_loss_supplied_energy = step_s * dot(
            average_coil_current,
            average_winding_loss_voltage,
        )
        previous_stored_loss_energy = 0.5 * dot(
            previous_winding_loss_state,
            model.winding_loss_storage_matrix_j * previous_winding_loss_state,
        )
        winding_loss_dissipation = step_s * (
            0.5 * dot(
                average_winding_loss_state,
                model.winding_loss_dissipation_matrix_j_per_s *
                average_winding_loss_state,
            ) +
            dot(
                average_coil_current,
                model.winding_loss_direct_ohm * average_coil_current,
            )
        )
        winding_loss_balance = winding_loss_supplied_energy -
            (
                state.stored_frequency_dependent_winding_energy_j -
                previous_stored_loss_energy
            ) - winding_loss_dissipation
        backward_euler || abs(winding_loss_balance) <=
            runtime.preparation.specification.settings.energy_absolute_tolerance_j +
            512.0 * eps(Float64) * max(abs(winding_loss_supplied_energy), 1.0) ||
            throw(ArgumentError(
                "transformer nonlinear magnetic winding-loss energy balance failed",
            ))
    end
    state.maximum_magnetic_continuity_residual_wb = max(
        state.maximum_magnetic_continuity_residual_wb,
        candidate.magnetic_continuity_residual_wb,
    )
    state.maximum_magnetic_constitutive_residual_at = max(
        state.maximum_magnetic_constitutive_residual_at,
        candidate.magnetic_constitutive_residual_at,
    )
    state.maximum_local_nonlinear_iterations = max(
        state.maximum_local_nonlinear_iterations,
        candidate.local_iterations,
    )
    _accept_transformer_energy!(
        state,
        previous_terminal_voltage,
        previous_terminal_current,
        voltage,
        current,
        step_s,
        runtime.preparation.specification.settings,
        runtime.companion_method,
    )
    all(isfinite, previous_branch_flux) || throw(ArgumentError(
        "transformer previous magnetic flux state became nonfinite",
    ))
    return nothing
end

function _accept_transformer_network_state!(
    runtime::TransformerApparatusRuntime,
    voltage,
)
    candidate = runtime.candidate::_TransformerNetworkCandidate
    state = runtime.accepted_state::TransformerNetworkRuntimeState
    current, _ = _transformer_network_trial!(runtime, voltage)
    previous_node_voltage = copy(state.represented_node_voltage_v)
    previous_terminal_voltage = copy(state.terminal_voltage_v)
    accepted_time = state.accepted_time_s
    event_time_tolerance = 64.0 * eps(Float64) * max(abs(accepted_time), 1.0)
    apparatus_current_changed_at_boundary = any(
        occurrence ->
            _transformer_event_requires_apparatus_current_reconstruction(
                occurrence.kind,
            ) &&
            abs(occurrence.time_s - accepted_time) <= event_time_tolerance,
        runtime.event_state.occurrences,
    )
    previous_terminal_current = apparatus_current_changed_at_boundary ?
        candidate.companion_admittance_s * previous_terminal_voltage .+
            candidate.history_current_a :
        copy(state.terminal_current_a)
    backward_euler = runtime.companion_method === :backward_euler
    state.represented_node_voltage_v .= candidate.represented_node_voltage_v
    state.terminal_voltage_v .= voltage
    state.terminal_current_a .= current
    state.branch_current_a .= candidate.branch_current_a
    state.branch_voltage_v .= candidate.branch_voltage_v
    state.capacitor_current_a .= candidate.capacitor_current_a
    internal_residual = candidate.nodal_admittance_s *
        state.represented_node_voltage_v .+ candidate.nodal_history_current_a
    state.maximum_internal_kcl_residual_a = max(
        state.maximum_internal_kcl_residual_a,
        maximum(
            abs,
            internal_residual[candidate.internal_node_indices];
            init=0.0,
        ),
    )
    state.stored_magnetic_energy_j = 0.0
    average_node_voltage = backward_euler ? state.represented_node_voltage_v :
        0.5 .* (previous_node_voltage .+ state.represented_node_voltage_v)
    state.dielectric_loss_energy_j += runtime.candidate_step_s * dot(
        average_node_voltage,
        candidate.conductance_s * average_node_voltage,
    )
    for group in candidate.branch_groups
        previous_branch_current = copy(group.accepted_current_a)
        group.accepted_current_a .= state.branch_current_a[group.current_offset]
        group.accepted_voltage_v .= state.branch_voltage_v[group.current_offset]
        state.stored_magnetic_energy_j += 0.5 * dot(
            group.accepted_current_a,
            group.inductance_h * group.accepted_current_a,
        )
        average_branch_current = backward_euler ? group.accepted_current_a :
            0.5 .* (previous_branch_current .+ group.accepted_current_a)
        state.winding_loss_energy_j += runtime.candidate_step_s * dot(
            average_branch_current,
            group.resistance_ohm * average_branch_current,
        )
    end
    state.stored_electric_energy_j = 0.5 * dot(
        state.represented_node_voltage_v,
        candidate.capacitance_f * state.represented_node_voltage_v,
    )
    for fault_id in sort!(
        collect(keys(runtime.event_state.active_internal_faults));
        by=String,
    )
        from_index, to_index, conductance_s =
            runtime.event_state.active_internal_faults[fault_id]
        previous_fault_voltage = previous_node_voltage[from_index] -
            previous_node_voltage[to_index]
        accepted_fault_voltage = state.represented_node_voltage_v[from_index] -
            state.represented_node_voltage_v[to_index]
        integration_fault_voltage = backward_euler ? accepted_fault_voltage :
            0.5 * (previous_fault_voltage + accepted_fault_voltage)
        internal_fault_energy_increment = runtime.candidate_step_s *
            conductance_s * integration_fault_voltage^2
        runtime.event_state.event_energy_j += internal_fault_energy_increment
        runtime.event_state.internal_fault_energy_j += internal_fault_energy_increment
    end
    _accept_transformer_energy!(
        state,
        previous_terminal_voltage,
        previous_terminal_current,
        voltage,
        current,
        runtime.candidate_step_s,
        runtime.preparation.specification.settings,
        runtime.companion_method,
    )
    return nothing
end

function accept_nonlinear_device_state!(
    runtime::TransformerApparatusRuntime,
    terminal_voltage_v::AbstractVector{Float64},
    _terminal_current_a::AbstractVector{Float64},
    _terminal_jacobian_s::AbstractMatrix{Float64},
    time_s::Float64,
)
    runtime.prepared || throw(ArgumentError(
        "transformer apparatus step must be prepared before state acceptance",
    ))
    abs(time_s - runtime.candidate_time_s) <=
        64.0 * eps(Float64) * max(abs(time_s), abs(runtime.candidate_time_s), 1.0) ||
        throw(ArgumentError("transformer accepted time does not match the prepared candidate"))
    network_current, _ = _transformer_event_trial!(runtime, terminal_voltage_v)
    apparatus_voltage = runtime.candidate_apparatus_terminal_voltage_v
    if runtime.candidate isa _TransformerTerminalMatrixCandidate
        _accept_transformer_terminal_state!(runtime, apparatus_voltage)
    elseif runtime.candidate isa _TransformerMagneticCandidate
        _accept_transformer_magnetic_state!(runtime, apparatus_voltage)
    elseif runtime.candidate isa _TransformerWidebandCandidate
        _accept_transformer_wideband_state!(runtime, apparatus_voltage)
    else
        _accept_transformer_network_state!(runtime, apparatus_voltage)
    end
    _accept_transformer_event_state!(runtime, terminal_voltage_v, network_current)
    runtime.accepted_evaluation_count += 1
    return nothing
end

function finish_nonlinear_device_step!(runtime::TransformerApparatusRuntime)
    runtime.prepared || throw(ArgumentError(
        "transformer apparatus has no prepared step to finish",
    ))
    runtime.prepared = false
    return nothing
end

struct TransformerApparatusRuntimeSnapshot
    schema_version::Int
    preparation_signature_sha256::String
    state::NamedTuple
    deterministic_signature_sha256::String
end

function _transformer_runtime_state_tuple(state::TransformerTerminalMatrixRuntimeState)
    return (
        kind=:terminal_matrix,
        coil_current_a=copy(state.coil_current_a),
        capacitor_current_a=copy(state.capacitor_current_a),
        coil_voltage_v=copy(state.coil_voltage_v),
        winding_loss_state=copy(state.winding_loss_state),
        winding_loss_voltage_v=copy(state.winding_loss_voltage_v),
        terminal_voltage_v=copy(state.terminal_voltage_v),
        terminal_current_a=copy(state.terminal_current_a),
        previous_terminal_power_w=state.previous_terminal_power_w,
        terminal_power_w=state.terminal_power_w,
        supplied_energy_j=state.supplied_energy_j,
        winding_loss_energy_j=state.winding_loss_energy_j,
        frequency_dependent_winding_loss_energy_j=
            state.frequency_dependent_winding_loss_energy_j,
        dielectric_loss_energy_j=state.dielectric_loss_energy_j,
        stored_magnetic_energy_j=state.stored_magnetic_energy_j,
        stored_frequency_dependent_winding_energy_j=
            state.stored_frequency_dependent_winding_energy_j,
        stored_electric_energy_j=state.stored_electric_energy_j,
        maximum_kcl_residual_a=state.maximum_kcl_residual_a,
        accepted_time_s=state.accepted_time_s,
        accepted_step_count=state.accepted_step_count,
    )
end

function _transformer_runtime_state_tuple(state::TransformerMagneticRuntimeState)
    return (
        kind=:nonlinear_magnetic,
        coil_current_a=copy(state.coil_current_a),
        capacitor_current_a=copy(state.capacitor_current_a),
        coil_voltage_v=copy(state.coil_voltage_v),
        branch_flux_wb=copy(state.branch_flux_wb),
        branch_mmf_drop_at=copy(state.branch_mmf_drop_at),
        magnetic_node_potential_at=copy(state.magnetic_node_potential_at),
        tellinen_state=copy(state.tellinen_state),
        winding_loss_state=copy(state.winding_loss_state),
        winding_loss_voltage_v=copy(state.winding_loss_voltage_v),
        terminal_voltage_v=copy(state.terminal_voltage_v),
        terminal_current_a=copy(state.terminal_current_a),
        previous_terminal_power_w=state.previous_terminal_power_w,
        terminal_power_w=state.terminal_power_w,
        supplied_energy_j=state.supplied_energy_j,
        winding_loss_energy_j=state.winding_loss_energy_j,
        frequency_dependent_winding_loss_energy_j=
            state.frequency_dependent_winding_loss_energy_j,
        dielectric_loss_energy_j=state.dielectric_loss_energy_j,
        hysteresis_loss_energy_j=state.hysteresis_loss_energy_j,
        classical_eddy_core_loss_energy_j=state.classical_eddy_core_loss_energy_j,
        excess_core_loss_energy_j=state.excess_core_loss_energy_j,
        stored_leakage_energy_j=state.stored_leakage_energy_j,
        stored_magnetic_energy_j=state.stored_magnetic_energy_j,
        stored_frequency_dependent_winding_energy_j=
            state.stored_frequency_dependent_winding_energy_j,
        stored_electric_energy_j=state.stored_electric_energy_j,
        last_classical_dynamic_mmf_at=copy(state.last_classical_dynamic_mmf_at),
        last_excess_dynamic_mmf_at=copy(state.last_excess_dynamic_mmf_at),
        maximum_magnetic_continuity_residual_wb=
            state.maximum_magnetic_continuity_residual_wb,
        maximum_magnetic_constitutive_residual_at=
            state.maximum_magnetic_constitutive_residual_at,
        maximum_local_nonlinear_iterations=state.maximum_local_nonlinear_iterations,
        accepted_time_s=state.accepted_time_s,
        accepted_step_count=state.accepted_step_count,
    )
end

function _transformer_runtime_state_tuple(state::TransformerWidebandRuntimeState)
    return (
        kind=:wideband,
        rational_state=copy(state.rational_state),
        terminal_voltage_v=copy(state.terminal_voltage_v),
        terminal_current_a=copy(state.terminal_current_a),
        previous_terminal_power_w=state.previous_terminal_power_w,
        terminal_power_w=state.terminal_power_w,
        supplied_energy_j=state.supplied_energy_j,
        minimum_supplied_energy_j=state.minimum_supplied_energy_j,
        dissipated_energy_j=state.dissipated_energy_j,
        stored_energy_j=state.stored_energy_j,
        maximum_energy_residual_j=state.maximum_energy_residual_j,
        maximum_state_magnitude=state.maximum_state_magnitude,
        accepted_time_s=state.accepted_time_s,
        accepted_step_count=state.accepted_step_count,
    )
end

function _transformer_runtime_state_tuple(state::TransformerNetworkRuntimeState)
    return (
        kind=:represented_network,
        represented_node_voltage_v=copy(state.represented_node_voltage_v),
        terminal_voltage_v=copy(state.terminal_voltage_v),
        terminal_current_a=copy(state.terminal_current_a),
        branch_current_a=copy(state.branch_current_a),
        branch_voltage_v=copy(state.branch_voltage_v),
        capacitor_current_a=copy(state.capacitor_current_a),
        previous_terminal_power_w=state.previous_terminal_power_w,
        terminal_power_w=state.terminal_power_w,
        supplied_energy_j=state.supplied_energy_j,
        winding_loss_energy_j=state.winding_loss_energy_j,
        dielectric_loss_energy_j=state.dielectric_loss_energy_j,
        stored_magnetic_energy_j=state.stored_magnetic_energy_j,
        stored_electric_energy_j=state.stored_electric_energy_j,
        maximum_internal_kcl_residual_a=state.maximum_internal_kcl_residual_a,
        accepted_time_s=state.accepted_time_s,
        accepted_step_count=state.accepted_step_count,
    )
end

function _write_transformer_snapshot_value(io, value)
    if value isa AbstractArray
        println(io, ndims(value))
        println(io, join(size(value), ','))
        for element in value
            _write_transformer_snapshot_value(io, element)
        end
    elseif value isa NamedTuple
        for name in keys(value)
            println(io, name)
            _write_transformer_snapshot_value(io, getfield(value, name))
        end
    elseif value isa Tuple
        println(io, length(value))
        for element in value
            _write_transformer_snapshot_value(io, element)
        end
    elseif value isa AbstractDict
        ordered_keys = sort!(collect(keys(value)); by=String)
        println(io, length(ordered_keys))
        for key in ordered_keys
            _write_transformer_snapshot_value(io, key)
            _write_transformer_snapshot_value(io, value[key])
        end
    elseif value isa TransformerApparatusEventCommand ||
           value isa TransformerApparatusEventOccurrence
        for name in fieldnames(typeof(value))
            println(io, name)
            _write_transformer_snapshot_value(io, getfield(value, name))
        end
    elseif value isa TellinenMagneticState
        for name in fieldnames(TellinenMagneticState)
            println(io, name)
            _write_transformer_snapshot_value(io, getfield(value, name))
        end
    elseif value isa AbstractFloat
        println(io, bitstring(Float64(value)))
    elseif value === nothing
        println(io, "nothing")
    else
        println(io, value)
    end
    return nothing
end

function _transformer_event_snapshot_state(runtime::TransformerApparatusRuntime)
    state = runtime.event_state
    return (
        terminal_nodes=copy(runtime.terminal_nodes),
        terminal_closed=copy(state.terminal_closed),
        terminal_transform=copy(state.terminal_transform),
        active_terminal_faults=deepcopy(state.active_terminal_faults),
        active_winding_faults=deepcopy(state.active_winding_faults),
        active_grounding_paths=deepcopy(state.active_grounding_paths),
        active_internal_faults=deepcopy(state.active_internal_faults),
        terminal_fault_admittance_s=copy(state.terminal_fault_admittance_s),
        scheduled_commands=copy(state.scheduled_commands),
        applied_event_ids=copy(state.applied_event_ids),
        occurrences=copy(state.occurrences),
        topology_transition_count=state.topology_transition_count,
        tap_change_count=state.tap_change_count,
        phase_shift_change_count=state.phase_shift_change_count,
        event_energy_j=state.event_energy_j,
        external_fault_energy_j=state.external_fault_energy_j,
        internal_fault_energy_j=state.internal_fault_energy_j,
        numerical_dissipation_energy_j=state.numerical_dissipation_energy_j,
        maximum_energy_balance_residual_j=state.maximum_energy_balance_residual_j,
        previous_network_voltage_v=copy(state.previous_network_voltage_v),
        previous_fault_current_a=copy(state.previous_fault_current_a),
        preparation_count=runtime.preparation_count,
        trial_evaluation_count=runtime.trial_evaluation_count,
        accepted_evaluation_count=runtime.accepted_evaluation_count,
        rejected_trial_count=runtime.rejected_trial_count,
        initial_stored_energy_j=runtime.initial_stored_energy_j,
        candidate_time_s=runtime.candidate_time_s,
        candidate_step_s=runtime.candidate_step_s,
        companion_method=runtime.companion_method,
    )
end

function _transformer_snapshot_signature(preparation_signature, state)
    io = IOBuffer()
    println(io, preparation_signature)
    for name in keys(state)
        println(io, name)
        value = getfield(state, name)
        _write_transformer_snapshot_value(io, value)
    end
    return bytes2hex(sha256(take!(io)))
end

function transformer_apparatus_runtime_snapshot(runtime::TransformerApparatusRuntime)
    runtime.prepared && throw(ArgumentError(
        "transformer snapshot cannot be taken during an active trial step",
    ))
    state = merge(
        _transformer_runtime_state_tuple(runtime.accepted_state),
        (event_state=_transformer_event_snapshot_state(runtime),),
    )
    signature = _transformer_snapshot_signature(
        runtime.preparation.preparation_signature_sha256,
        state,
    )
    return TransformerApparatusRuntimeSnapshot(
        1,
        runtime.preparation.preparation_signature_sha256,
        state,
        signature,
    )
end

function _restore_transformer_event_snapshot!(runtime, snapshot)
    hasproperty(snapshot, :event_state) || return runtime
    source = snapshot.event_state
    state = runtime.event_state
    source.terminal_nodes == runtime.terminal_nodes || throw(ArgumentError(
        "transformer snapshot network-terminal identity is incompatible",
    ))
    length(source.terminal_closed) == length(state.terminal_closed) ||
        throw(DimensionMismatch("transformer snapshot event terminal count is incompatible"))
    size(source.terminal_transform) == size(state.terminal_transform) ||
        throw(DimensionMismatch("transformer snapshot event transform is incompatible"))
    size(source.terminal_fault_admittance_s) ==
        size(state.terminal_fault_admittance_s) || throw(DimensionMismatch(
            "transformer snapshot event fault admittance is incompatible",
        ))
    state.terminal_closed .= source.terminal_closed
    state.terminal_transform .= source.terminal_transform
    empty!(state.active_terminal_faults)
    merge!(state.active_terminal_faults, deepcopy(source.active_terminal_faults))
    empty!(state.active_winding_faults)
    merge!(state.active_winding_faults, deepcopy(source.active_winding_faults))
    empty!(state.active_grounding_paths)
    merge!(state.active_grounding_paths, deepcopy(source.active_grounding_paths))
    empty!(state.active_internal_faults)
    merge!(state.active_internal_faults, deepcopy(source.active_internal_faults))
    state.terminal_fault_admittance_s .= source.terminal_fault_admittance_s
    empty!(state.scheduled_commands)
    append!(state.scheduled_commands, deepcopy(source.scheduled_commands))
    empty!(state.applied_event_ids)
    append!(state.applied_event_ids, source.applied_event_ids)
    empty!(state.occurrences)
    append!(state.occurrences, deepcopy(source.occurrences))
    state.topology_transition_count = source.topology_transition_count
    state.tap_change_count = source.tap_change_count
    state.phase_shift_change_count = source.phase_shift_change_count
    state.event_energy_j = source.event_energy_j
    state.external_fault_energy_j = source.external_fault_energy_j
    state.internal_fault_energy_j = source.internal_fault_energy_j
    state.numerical_dissipation_energy_j = source.numerical_dissipation_energy_j
    state.maximum_energy_balance_residual_j = source.maximum_energy_balance_residual_j
    state.previous_network_voltage_v .= source.previous_network_voltage_v
    state.previous_fault_current_a .= source.previous_fault_current_a
    runtime.preparation_count = source.preparation_count
    runtime.trial_evaluation_count = source.trial_evaluation_count
    runtime.accepted_evaluation_count = source.accepted_evaluation_count
    runtime.rejected_trial_count = source.rejected_trial_count
    runtime.initial_stored_energy_j == source.initial_stored_energy_j ||
        throw(ArgumentError("transformer snapshot initial energy identity is incompatible"))
    runtime.candidate_time_s = source.candidate_time_s
    runtime.candidate_step_s = source.candidate_step_s
    runtime.companion_method = source.companion_method
    return runtime
end

function _restore_transformer_state!(
    state::TransformerTerminalMatrixRuntimeState,
    snapshot::NamedTuple,
)
    snapshot.kind === :terminal_matrix || throw(ArgumentError(
        "transformer snapshot state family is incompatible",
    ))
    for name in (
        :coil_current_a,
        :capacitor_current_a,
        :coil_voltage_v,
        :winding_loss_state,
        :winding_loss_voltage_v,
        :terminal_voltage_v,
        :terminal_current_a,
    )
        target = getfield(state, name)
        source = getfield(snapshot, name)
        length(target) == length(source) || throw(DimensionMismatch(
            "transformer snapshot vector size is incompatible",
        ))
        target .= source
    end
    for name in (
        :previous_terminal_power_w,
        :terminal_power_w,
        :supplied_energy_j,
        :winding_loss_energy_j,
        :frequency_dependent_winding_loss_energy_j,
        :dielectric_loss_energy_j,
        :stored_magnetic_energy_j,
        :stored_frequency_dependent_winding_energy_j,
        :stored_electric_energy_j,
        :maximum_kcl_residual_a,
        :accepted_time_s,
        :accepted_step_count,
    )
        setfield!(state, name, getfield(snapshot, name))
    end
end

function _restore_transformer_state!(
    state::TransformerWidebandRuntimeState,
    snapshot::NamedTuple,
)
    snapshot.kind === :wideband || throw(ArgumentError(
        "transformer snapshot state family is incompatible",
    ))
    for name in (:rational_state, :terminal_voltage_v, :terminal_current_a)
        target = getfield(state, name)
        source = getfield(snapshot, name)
        length(target) == length(source) || throw(DimensionMismatch(
            "transformer snapshot vector size is incompatible",
        ))
        target .= source
    end
    for name in (
        :previous_terminal_power_w,
        :terminal_power_w,
        :supplied_energy_j,
        :minimum_supplied_energy_j,
        :dissipated_energy_j,
        :stored_energy_j,
        :maximum_energy_residual_j,
        :maximum_state_magnitude,
        :accepted_time_s,
        :accepted_step_count,
    )
        setfield!(state, name, getfield(snapshot, name))
    end
end

function _restore_transformer_state!(
    state::TransformerMagneticRuntimeState,
    snapshot::NamedTuple,
)
    snapshot.kind === :nonlinear_magnetic || throw(ArgumentError(
        "transformer snapshot state family is incompatible",
    ))
    for name in (
        :coil_current_a,
        :capacitor_current_a,
        :coil_voltage_v,
        :branch_flux_wb,
        :branch_mmf_drop_at,
        :magnetic_node_potential_at,
        :tellinen_state,
        :winding_loss_state,
        :winding_loss_voltage_v,
        :terminal_voltage_v,
        :terminal_current_a,
        :last_classical_dynamic_mmf_at,
        :last_excess_dynamic_mmf_at,
    )
        target = getfield(state, name)
        source = getfield(snapshot, name)
        length(target) == length(source) || throw(DimensionMismatch(
            "transformer nonlinear magnetic snapshot vector size is incompatible",
        ))
        target .= source
    end
    for name in (
        :previous_terminal_power_w,
        :terminal_power_w,
        :supplied_energy_j,
        :winding_loss_energy_j,
        :frequency_dependent_winding_loss_energy_j,
        :dielectric_loss_energy_j,
        :hysteresis_loss_energy_j,
        :classical_eddy_core_loss_energy_j,
        :excess_core_loss_energy_j,
        :stored_leakage_energy_j,
        :stored_magnetic_energy_j,
        :stored_frequency_dependent_winding_energy_j,
        :stored_electric_energy_j,
        :maximum_magnetic_continuity_residual_wb,
        :maximum_magnetic_constitutive_residual_at,
        :maximum_local_nonlinear_iterations,
        :accepted_time_s,
        :accepted_step_count,
    )
        setfield!(state, name, getfield(snapshot, name))
    end
end

function _restore_transformer_state!(
    state::TransformerNetworkRuntimeState,
    snapshot::NamedTuple,
)
    snapshot.kind === :represented_network || throw(ArgumentError(
        "transformer snapshot state family is incompatible",
    ))
    for name in (
        :represented_node_voltage_v,
        :terminal_voltage_v,
        :terminal_current_a,
        :branch_current_a,
        :branch_voltage_v,
        :capacitor_current_a,
    )
        target = getfield(state, name)
        source = getfield(snapshot, name)
        length(target) == length(source) || throw(DimensionMismatch(
            "transformer snapshot vector size is incompatible",
        ))
        target .= source
    end
    for name in (
        :previous_terminal_power_w,
        :terminal_power_w,
        :supplied_energy_j,
        :winding_loss_energy_j,
        :dielectric_loss_energy_j,
        :stored_magnetic_energy_j,
        :stored_electric_energy_j,
        :maximum_internal_kcl_residual_a,
        :accepted_time_s,
        :accepted_step_count,
    )
        setfield!(state, name, getfield(snapshot, name))
    end
end

function restore_transformer_apparatus_runtime_snapshot!(
    runtime::TransformerApparatusRuntime,
    snapshot::TransformerApparatusRuntimeSnapshot,
)
    runtime.prepared && throw(ArgumentError(
        "transformer snapshot cannot be restored during an active trial step",
    ))
    snapshot.schema_version == 1 || throw(ArgumentError(
        "transformer snapshot schema version is unsupported",
    ))
    snapshot.preparation_signature_sha256 ==
        runtime.preparation.preparation_signature_sha256 || throw(ArgumentError(
            "transformer snapshot preparation identity is stale",
        ))
    expected_signature = _transformer_snapshot_signature(
        snapshot.preparation_signature_sha256,
        snapshot.state,
    )
    snapshot.deterministic_signature_sha256 == expected_signature ||
        throw(ArgumentError("transformer snapshot integrity signature does not match"))
    _restore_transformer_state!(runtime.accepted_state, snapshot.state)
    _restore_transformer_event_snapshot!(runtime, snapshot.state)
    if runtime.candidate isa _TransformerNetworkCandidate
        state = runtime.accepted_state::TransformerNetworkRuntimeState
        for group in runtime.candidate.branch_groups
            group.accepted_current_a .= state.branch_current_a[group.current_offset]
            group.accepted_voltage_v .= state.branch_voltage_v[group.current_offset]
        end
    end
    return runtime
end

function transformer_apparatus_runtime_diagnostics(runtime::TransformerApparatusRuntime)
    state = runtime.accepted_state
    common = (
        apparatus=runtime.preparation.specification.id,
        tier=_TRANSFORMER_TIER_IDS[runtime.preparation.specification.tier],
        initialization_mode=_TRANSFORMER_INITIALIZATION_MODE_IDS[
            runtime.preparation.initialization_mode
        ],
        initial_time_s=runtime.preparation.initial_time_s,
        residual_flux_projection_correction_wb=
            runtime.preparation.residual_flux_projection_correction_wb,
        operating_point_electrical_residual_v=
            runtime.preparation.operating_point_electrical_residual_v,
        preparation_signature_sha256=
            runtime.preparation.preparation_signature_sha256,
        accepted_time_s=state.accepted_time_s,
        accepted_step_count=state.accepted_step_count,
        terminal_power_w=state.terminal_power_w,
        supplied_energy_j=state.supplied_energy_j,
        preparation_count=runtime.preparation_count,
        trial_evaluation_count=runtime.trial_evaluation_count,
        accepted_evaluation_count=runtime.accepted_evaluation_count,
        prepared=runtime.prepared,
        terminal_closed=copy(runtime.event_state.terminal_closed),
        active_terminal_fault_count=length(runtime.event_state.active_terminal_faults),
        active_winding_fault_count=length(runtime.event_state.active_winding_faults),
        active_grounding_path_count=length(runtime.event_state.active_grounding_paths),
        active_internal_fault_count=length(runtime.event_state.active_internal_faults),
        event_occurrence_count=length(runtime.event_state.occurrences),
        topology_transition_count=runtime.event_state.topology_transition_count,
        tap_change_count=runtime.event_state.tap_change_count,
        phase_shift_change_count=runtime.event_state.phase_shift_change_count,
        event_energy_j=runtime.event_state.event_energy_j,
        external_fault_energy_j=runtime.event_state.external_fault_energy_j,
        internal_fault_energy_j=runtime.event_state.internal_fault_energy_j,
        numerical_dissipation_energy_j=
            runtime.event_state.numerical_dissipation_energy_j,
        maximum_energy_balance_residual_j=
            runtime.event_state.maximum_energy_balance_residual_j,
        total_supplied_energy_j=state.supplied_energy_j +
            runtime.event_state.external_fault_energy_j,
    )
    physical = if state isa TransformerTerminalMatrixRuntimeState
        (
            winding_loss_energy_j=state.winding_loss_energy_j,
            frequency_dependent_winding_loss_energy_j=
                state.frequency_dependent_winding_loss_energy_j,
            dielectric_loss_energy_j=state.dielectric_loss_energy_j,
            stored_magnetic_energy_j=state.stored_magnetic_energy_j,
            stored_frequency_dependent_winding_energy_j=
                state.stored_frequency_dependent_winding_energy_j,
            stored_electric_energy_j=state.stored_electric_energy_j,
            maximum_kcl_residual_a=state.maximum_kcl_residual_a,
        )
    elseif state isa TransformerMagneticRuntimeState
        (
            winding_loss_energy_j=state.winding_loss_energy_j,
            frequency_dependent_winding_loss_energy_j=
                state.frequency_dependent_winding_loss_energy_j,
            dielectric_loss_energy_j=state.dielectric_loss_energy_j,
            hysteresis_loss_energy_j=state.hysteresis_loss_energy_j,
            classical_eddy_core_loss_energy_j=
                state.classical_eddy_core_loss_energy_j,
            excess_core_loss_energy_j=state.excess_core_loss_energy_j,
            stored_leakage_energy_j=state.stored_leakage_energy_j,
            stored_magnetic_energy_j=state.stored_magnetic_energy_j,
            stored_frequency_dependent_winding_energy_j=
                state.stored_frequency_dependent_winding_energy_j,
            stored_electric_energy_j=state.stored_electric_energy_j,
            maximum_magnetic_continuity_residual_wb=
                state.maximum_magnetic_continuity_residual_wb,
            maximum_magnetic_constitutive_residual_at=
                state.maximum_magnetic_constitutive_residual_at,
            maximum_local_nonlinear_iterations=
                state.maximum_local_nonlinear_iterations,
        )
    elseif state isa TransformerWidebandRuntimeState
        (
            dissipated_energy_j=state.dissipated_energy_j,
            stored_energy_j=state.stored_energy_j,
            minimum_supplied_energy_j=state.minimum_supplied_energy_j,
            maximum_energy_residual_j=state.maximum_energy_residual_j,
            maximum_state_magnitude=state.maximum_state_magnitude,
        )
    else
        state = state::TransformerNetworkRuntimeState
        (
            winding_loss_energy_j=state.winding_loss_energy_j,
            dielectric_loss_energy_j=state.dielectric_loss_energy_j,
            stored_magnetic_energy_j=state.stored_magnetic_energy_j,
            stored_electric_energy_j=state.stored_electric_energy_j,
            maximum_internal_kcl_residual_a=
                state.maximum_internal_kcl_residual_a,
        )
    end
    return merge(common, physical)
end
