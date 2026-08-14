export ModernMachineFamily,
       WoundFieldSynchronousMachine,
       CageInductionMachine,
       WoundRotorInductionMachine,
       PermanentMagnetSynchronousMachine,
       DoublyFedInductionMachine,
       SynchronousCondenserMachine,
       ModernMachineOperatingMode,
       MachineMotorMode,
       MachineGeneratorMode,
       MachineCondenserMode,
       ModernMachineInitializationMode,
       DeenergizedMachineInitialization,
       SpecifiedMachineInitialization,
       SinusoidalMachineOperatingPoint,
       MachineRotorBranch,
       MachineElectricalParameters,
       MachineMagneticCoenergyLaw,
       MachineShaftMass,
       MachineShaftCoupling,
       MachineControlParameters,
       MachineRuntimeSettings,
       ModernMachineSpecification,
       ModernMachinePreparation,
       ModernMachineReadiness,
       ModernMachineRefusal,
       modern_machine_contract,
       prepare_modern_machine,
       modern_machine_readiness,
       modern_machine_signature

@enum ModernMachineFamily begin
    WoundFieldSynchronousMachine
    CageInductionMachine
    WoundRotorInductionMachine
    PermanentMagnetSynchronousMachine
    DoublyFedInductionMachine
    SynchronousCondenserMachine
end

@enum ModernMachineOperatingMode begin
    MachineMotorMode
    MachineGeneratorMode
    MachineCondenserMode
end

@enum ModernMachineInitializationMode begin
    DeenergizedMachineInitialization
    SpecifiedMachineInitialization
    SinusoidalMachineOperatingPoint
end

const _MODERN_MACHINE_FAMILY_IDS = Dict(
    WoundFieldSynchronousMachine => :wound_field_synchronous,
    CageInductionMachine => :cage_induction,
    WoundRotorInductionMachine => :wound_rotor_induction,
    PermanentMagnetSynchronousMachine => :permanent_magnet_synchronous,
    DoublyFedInductionMachine => :doubly_fed_induction,
    SynchronousCondenserMachine => :synchronous_condenser,
)

const _MODERN_MACHINE_MODE_IDS = Dict(
    MachineMotorMode => :motor,
    MachineGeneratorMode => :generator,
    MachineCondenserMode => :condenser,
)

struct ModernMachineRefusal <: Exception
    code::Symbol
    operation::Symbol
    machine::Symbol
    family::ModernMachineFamily
    message::String
    diagnostics::NamedTuple
end

Base.showerror(io::IO, refusal::ModernMachineRefusal) = print(
    io,
    String(refusal.code),
    " during ",
    String(refusal.operation),
    " for modern machine ",
    String(refusal.machine),
    " (",
    String(_MODERN_MACHINE_FAMILY_IDS[refusal.family]),
    "): ",
    refusal.message,
)

function _modern_machine_refusal(
    code::Symbol,
    operation::Symbol,
    specification,
    message::AbstractString;
    diagnostics=NamedTuple(),
)
    throw(ModernMachineRefusal(
        code,
        operation,
        specification.id,
        specification.family,
        String(message),
        diagnostics,
    ))
end

struct MachineRotorBranch
    id::Symbol
    resistance_ohm::Float64
    leakage_inductance_h::Float64
    terminal_exposed::Bool

    function MachineRotorBranch(
        id::Symbol;
        resistance_ohm::Real,
        leakage_inductance_h::Real,
        terminal_exposed::Bool=false,
    )
        isempty(String(id)) && throw(ArgumentError("machine rotor-branch id must not be empty"))
        resistance = Float64(resistance_ohm)
        leakage = Float64(leakage_inductance_h)
        isfinite(resistance) && resistance > 0.0 || throw(ArgumentError(
            "machine rotor-branch resistance must be finite and positive",
        ))
        isfinite(leakage) && leakage > 0.0 || throw(ArgumentError(
            "machine rotor-branch leakage inductance must be finite and positive",
        ))
        return new(id, resistance, leakage, terminal_exposed)
    end
end

struct MachineElectricalParameters
    stator_resistance_ohm::Float64
    zero_sequence_inductance_h::Float64
    stator_d_leakage_inductance_h::Float64
    stator_q_leakage_inductance_h::Float64
    d_axis_magnetizing_inductance_h::Float64
    q_axis_magnetizing_inductance_h::Float64
    field_resistance_ohm::Float64
    field_leakage_inductance_h::Float64
    d_damper_resistance_ohm::Float64
    d_damper_leakage_inductance_h::Float64
    q_damper_resistance_ohm::Float64
    q_damper_leakage_inductance_h::Float64
    permanent_magnet_flux_wb::Float64
    iron_loss_conductance_s::Float64

    function MachineElectricalParameters(;
        stator_resistance_ohm::Real,
        zero_sequence_inductance_h::Real,
        stator_d_leakage_inductance_h::Real,
        stator_q_leakage_inductance_h::Real,
        d_axis_magnetizing_inductance_h::Real,
        q_axis_magnetizing_inductance_h::Real,
        field_resistance_ohm::Real=0.0,
        field_leakage_inductance_h::Real=0.0,
        d_damper_resistance_ohm::Real=0.0,
        d_damper_leakage_inductance_h::Real=0.0,
        q_damper_resistance_ohm::Real=0.0,
        q_damper_leakage_inductance_h::Real=0.0,
        permanent_magnet_flux_wb::Real=0.0,
        iron_loss_conductance_s::Real=0.0,
    )
        values = Float64.((
            stator_resistance_ohm,
            zero_sequence_inductance_h,
            stator_d_leakage_inductance_h,
            stator_q_leakage_inductance_h,
            d_axis_magnetizing_inductance_h,
            q_axis_magnetizing_inductance_h,
            field_resistance_ohm,
            field_leakage_inductance_h,
            d_damper_resistance_ohm,
            d_damper_leakage_inductance_h,
            q_damper_resistance_ohm,
            q_damper_leakage_inductance_h,
            permanent_magnet_flux_wb,
            iron_loss_conductance_s,
        ))
        all(isfinite, values) || throw(ArgumentError(
            "machine electrical parameters must be finite",
        ))
        all(>(0.0), values[1:6]) || throw(ArgumentError(
            "machine stator resistance and base inductances must be positive",
        ))
        all(>=(0.0), values[7:end]) || throw(ArgumentError(
            "machine optional electrical parameters must be nonnegative",
        ))
        return new(values...)
    end
end

"""Convex magnetic coenergy correction with symmetric d/q cross derivatives."""
struct MachineMagneticCoenergyLaw
    radial_coefficient_per_wb2_h::Float64
    cross_coefficient_per_wb2_h::Float64
    maximum_flux_wb::Float64

    function MachineMagneticCoenergyLaw(;
        radial_coefficient_per_wb2_h::Real=0.0,
        cross_coefficient_per_wb2_h::Real=0.0,
        maximum_flux_wb::Real=Inf,
    )
        radial = Float64(radial_coefficient_per_wb2_h)
        cross = Float64(cross_coefficient_per_wb2_h)
        maximum_flux = Float64(maximum_flux_wb)
        isfinite(radial) && radial >= 0.0 || throw(ArgumentError(
            "machine radial saturation coefficient must be finite and nonnegative",
        ))
        isfinite(cross) && cross >= 0.0 || throw(ArgumentError(
            "machine cross-saturation coefficient must be finite and nonnegative",
        ))
        cross <= 2.0 * radial || throw(ArgumentError(
            "machine cross-saturation coefficient must preserve global coenergy convexity",
        ))
        (isfinite(maximum_flux) && maximum_flux > 0.0) || maximum_flux == Inf ||
            throw(ArgumentError("machine maximum flux must be positive or infinite"))
        return new(radial, cross, maximum_flux)
    end
end

struct MachineShaftMass
    id::Symbol
    inertia_kg_m2::Float64
    damping_nm_s_per_rad::Float64
    initial_angle_rad::Float64
    initial_speed_rad_s::Float64

    function MachineShaftMass(
        id::Symbol;
        inertia_kg_m2::Real,
        damping_nm_s_per_rad::Real=0.0,
        initial_angle_rad::Real=0.0,
        initial_speed_rad_s::Real=0.0,
    )
        inertia = Float64(inertia_kg_m2)
        damping = Float64(damping_nm_s_per_rad)
        angle = Float64(initial_angle_rad)
        speed = Float64(initial_speed_rad_s)
        isempty(String(id)) && throw(ArgumentError("machine shaft-mass id must not be empty"))
        isfinite(inertia) && inertia > 0.0 || throw(ArgumentError(
            "machine shaft inertia must be finite and positive",
        ))
        isfinite(damping) && damping >= 0.0 || throw(ArgumentError(
            "machine shaft damping must be finite and nonnegative",
        ))
        isfinite(angle) && isfinite(speed) || throw(ArgumentError(
            "machine initial shaft angle and speed must be finite",
        ))
        return new(id, inertia, damping, angle, speed)
    end
end

struct MachineShaftCoupling
    id::Symbol
    left_mass::Symbol
    right_mass::Symbol
    stiffness_nm_per_rad::Float64
    damping_nm_s_per_rad::Float64
    initial_twist_rad::Float64

    function MachineShaftCoupling(
        id::Symbol,
        left_mass::Symbol,
        right_mass::Symbol;
        stiffness_nm_per_rad::Real,
        damping_nm_s_per_rad::Real=0.0,
        initial_twist_rad::Real=0.0,
    )
        isempty(String(id)) && throw(ArgumentError("machine shaft-coupling id must not be empty"))
        left_mass != right_mass || throw(ArgumentError(
            "machine shaft coupling must connect distinct masses",
        ))
        stiffness = Float64(stiffness_nm_per_rad)
        damping = Float64(damping_nm_s_per_rad)
        twist = Float64(initial_twist_rad)
        isfinite(stiffness) && stiffness > 0.0 || throw(ArgumentError(
            "machine shaft stiffness must be finite and positive",
        ))
        isfinite(damping) && damping >= 0.0 || throw(ArgumentError(
            "machine shaft coupling damping must be finite and nonnegative",
        ))
        isfinite(twist) || throw(ArgumentError(
            "machine shaft initial twist must be finite",
        ))
        return new(id, left_mass, right_mass, stiffness, damping, twist)
    end
end

struct MachineControlParameters
    enabled::Bool
    task_period_s::Float64
    task_phase_s::Float64
    voltage_reference_v::Float64
    excitation_gain::Float64
    excitation_time_constant_s::Float64
    field_voltage_min_v::Float64
    field_voltage_max_v::Float64
    speed_reference_rad_s::Float64
    governor_droop_rad_s_per_nm::Float64
    governor_time_constant_s::Float64
    torque_min_nm::Float64
    torque_max_nm::Float64
    stabilizer_gain::Float64
    stabilizer_washout_s::Float64
    stabilizer_lead_s::Float64
    stabilizer_lag_s::Float64

    function MachineControlParameters(;
        enabled::Bool=false,
        task_period_s::Real=1.0e-3,
        task_phase_s::Real=0.0,
        voltage_reference_v::Real=1.0,
        excitation_gain::Real=0.0,
        excitation_time_constant_s::Real=0.05,
        field_voltage_min_v::Real=0.0,
        field_voltage_max_v::Real=1.0,
        speed_reference_rad_s::Real=0.0,
        governor_droop_rad_s_per_nm::Real=0.0,
        governor_time_constant_s::Real=0.1,
        torque_min_nm::Real=-Inf,
        torque_max_nm::Real=Inf,
        stabilizer_gain::Real=0.0,
        stabilizer_washout_s::Real=1.0,
        stabilizer_lead_s::Real=0.05,
        stabilizer_lag_s::Real=0.1,
    )
        period = Float64(task_period_s)
        phase = Float64(task_phase_s)
        values = Float64.((
            voltage_reference_v,
            excitation_gain,
            excitation_time_constant_s,
            field_voltage_min_v,
            field_voltage_max_v,
            speed_reference_rad_s,
            governor_droop_rad_s_per_nm,
            governor_time_constant_s,
            torque_min_nm,
            torque_max_nm,
            stabilizer_gain,
            stabilizer_washout_s,
            stabilizer_lead_s,
            stabilizer_lag_s,
        ))
        isfinite(period) && period > 0.0 || throw(ArgumentError(
            "machine control task period must be finite and positive",
        ))
        isfinite(phase) && 0.0 <= phase < period || throw(ArgumentError(
            "machine control task phase must lie in [0, task period)",
        ))
        all(isfinite, values[[1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14]]) ||
            throw(ArgumentError("machine control parameters must be finite"))
        values[3] > 0.0 && values[8] > 0.0 && values[12] > 0.0 &&
            values[13] >= 0.0 && values[14] > 0.0 || throw(ArgumentError(
            "machine control time constants must be positive and lead nonnegative",
        ))
        values[4] <= values[5] || throw(ArgumentError(
            "machine field-voltage limits are reversed",
        ))
        values[9] <= values[10] || throw(ArgumentError(
            "machine torque limits are reversed",
        ))
        return new(enabled, period, phase, values...)
    end
end

struct MachineRuntimeSettings
    timestep_s::Float64
    nonlinear_tolerance::Float64
    maximum_nonlinear_iterations::Int
    event_time_tolerance_s::Float64
    energy_tolerance_j::Float64

    function MachineRuntimeSettings(;
        timestep_s::Real,
        nonlinear_tolerance::Real=1.0e-10,
        maximum_nonlinear_iterations::Integer=16,
        event_time_tolerance_s::Real=1.0e-13,
        energy_tolerance_j::Real=1.0e-8,
    )
        values = Float64.((
            timestep_s,
            nonlinear_tolerance,
            event_time_tolerance_s,
            energy_tolerance_j,
        ))
        all(value -> isfinite(value) && value > 0.0, values) || throw(ArgumentError(
            "machine runtime tolerances and timestep must be finite and positive",
        ))
        iterations = Int(maximum_nonlinear_iterations)
        iterations >= 2 || throw(ArgumentError(
            "machine nonlinear iteration limit must be at least two",
        ))
        return new(values[1], values[2], iterations, values[3], values[4])
    end
end

struct ModernMachineSpecification
    id::Symbol
    family::ModernMachineFamily
    operating_mode::ModernMachineOperatingMode
    pole_pairs::Int
    phase_order::NTuple{3,Symbol}
    terminal_order::NTuple{4,Symbol}
    electrical::MachineElectricalParameters
    rotor_branches::Vector{MachineRotorBranch}
    saturation::MachineMagneticCoenergyLaw
    shaft_masses::Vector{MachineShaftMass}
    shaft_couplings::Vector{MachineShaftCoupling}
    electromagnetic_mass::Symbol
    controls::MachineControlParameters
    settings::MachineRuntimeSettings
    initialization_mode::ModernMachineInitializationMode
    initial_phase_voltage_v::NTuple{3,Float64}
    initial_phase_current_a::NTuple{3,Float64}
    initial_field_voltage_v::Float64
    initial_rotor_voltage_dq_v::NTuple{2,Float64}
    initial_mechanical_torque_nm::Float64
    provenance::ParameterProvenance
    uncertainty::String
    validity_domain::String

    function ModernMachineSpecification(
        id::Symbol,
        family::ModernMachineFamily;
        operating_mode::ModernMachineOperatingMode=MachineMotorMode,
        pole_pairs::Integer,
        phase_order=(:a, :b, :c),
        terminal_order=(:a, :b, :c, :neutral),
        electrical::MachineElectricalParameters,
        rotor_branches=MachineRotorBranch[],
        saturation::MachineMagneticCoenergyLaw=MachineMagneticCoenergyLaw(),
        shaft_masses,
        shaft_couplings=MachineShaftCoupling[],
        electromagnetic_mass::Symbol=first(shaft_masses).id,
        controls::MachineControlParameters=MachineControlParameters(),
        settings::MachineRuntimeSettings,
        initialization_mode::ModernMachineInitializationMode=
            DeenergizedMachineInitialization,
        initial_phase_voltage_v=(0.0, 0.0, 0.0),
        initial_phase_current_a=(0.0, 0.0, 0.0),
        initial_field_voltage_v::Real=0.0,
        initial_rotor_voltage_dq_v=(0.0, 0.0),
        initial_mechanical_torque_nm::Real=0.0,
        provenance::ParameterProvenance=ParameterProvenance(
            "caller-supplied generic modern machine data",
            "SI",
            "explicit typed values; no inferred per-unit or family conversion",
            "unknown; caller must replace for quantitative claims",
            "declared generic machine family and runtime settings",
            PhysicalModelParameter,
        ),
        uncertainty::AbstractString="unknown",
        validity_domain::AbstractString="generic phase-domain fixed-step EMT machine",
    )
        isempty(String(id)) && throw(ArgumentError("modern machine id must not be empty"))
        pairs = Int(pole_pairs)
        pairs > 0 || throw(ArgumentError("modern machine pole-pair count must be positive"))
        phases = Tuple(Symbol.(phase_order))
        terminals = Tuple(Symbol.(terminal_order))
        length(phases) == 3 && length(unique(phases)) == 3 || throw(ArgumentError(
            "modern machine requires three unique ordered phases",
        ))
        length(terminals) == 4 && length(unique(terminals)) == 4 || throw(ArgumentError(
            "modern machine requires three phase terminals and one unique neutral",
        ))
        masses = MachineShaftMass[shaft_masses...]
        couplings = MachineShaftCoupling[shaft_couplings...]
        branches = MachineRotorBranch[rotor_branches...]
        isempty(masses) && throw(ArgumentError("modern machine requires at least one shaft mass"))
        length(masses) <= 16 || throw(ArgumentError(
            "modern machine supports at most sixteen shaft masses",
        ))
        mass_ids = getfield.(masses, :id)
        length(unique(mass_ids)) == length(mass_ids) || throw(ArgumentError(
            "modern machine shaft-mass ids must be unique",
        ))
        electromagnetic_mass in mass_ids || throw(ArgumentError(
            "modern machine electromagnetic mass is not present",
        ))
        branch_ids = getfield.(branches, :id)
        length(unique(branch_ids)) == length(branch_ids) || throw(ArgumentError(
            "modern machine rotor-branch ids must be unique",
        ))
        phase_voltage = Tuple(Float64.(initial_phase_voltage_v))
        phase_current = Tuple(Float64.(initial_phase_current_a))
        rotor_voltage = Tuple(Float64.(initial_rotor_voltage_dq_v))
        length(phase_voltage) == 3 && length(phase_current) == 3 &&
            length(rotor_voltage) == 2 || throw(DimensionMismatch(
                "modern machine initial phase and rotor values have invalid dimensions",
            ))
        initial_values = (
            phase_voltage...,
            phase_current...,
            rotor_voltage...,
            Float64(initial_field_voltage_v),
            Float64(initial_mechanical_torque_nm),
        )
        all(isfinite, initial_values) || throw(ArgumentError(
            "modern machine initial values must be finite",
        ))
        isempty(strip(uncertainty)) && throw(ArgumentError(
            "modern machine uncertainty must be explicit",
        ))
        isempty(strip(validity_domain)) && throw(ArgumentError(
            "modern machine validity domain must be explicit",
        ))
        return new(
            id,
            family,
            operating_mode,
            pairs,
            phases,
            terminals,
            electrical,
            branches,
            saturation,
            masses,
            couplings,
            electromagnetic_mass,
            controls,
            settings,
            initialization_mode,
            phase_voltage,
            phase_current,
            Float64(initial_field_voltage_v),
            rotor_voltage,
            Float64(initial_mechanical_torque_nm),
            provenance,
            String(uncertainty),
            String(validity_domain),
        )
    end
end

struct ModernMachineElectricalLayout
    zero_index::Int
    d_indices::UnitRange{Int}
    q_indices::UnitRange{Int}
    stator_d_index::Int
    stator_q_index::Int
    field_index::Union{Nothing,Int}
    d_damper_index::Union{Nothing,Int}
    q_damper_index::Union{Nothing,Int}
    rotor_d_indices::Vector{Int}
    rotor_q_indices::Vector{Int}
end

struct ModernMachinePreparation
    specification::ModernMachineSpecification
    layout::ModernMachineElectricalLayout
    inductance_h::Matrix{Float64}
    inverse_inductance_per_h::Matrix{Float64}
    resistance_ohm::Vector{Float64}
    permanent_flux_offset_wb::Vector{Float64}
    terminal_voltage_map::Matrix{Float64}
    terminal_current_map::Matrix{Float64}
    initial_flux_wb::Vector{Float64}
    deterministic_signature_sha256::String
end

struct ModernMachineReadiness
    ready::Bool
    code::Symbol
    family::Symbol
    state_count::Int
    rotor_branch_count::Int
    shaft_mass_count::Int
    control_enabled::Bool
    deterministic_signature_sha256::String
    limitations::Vector{String}
end

function modern_machine_contract()
    validity = ModelValidityDomain(
        :modern_phase_domain_machine;
        description="Three ordered stator phases and explicit neutral with declared rotor, field, shaft, saturation, control, and fixed-step domains.",
        bounds=(
            NumericDomainBound(:phase_count; unit="count", lower=3, upper=3),
            NumericDomainBound(:shaft_mass_count; unit="count", lower=1, upper=16),
            NumericDomainBound(:deep_bar_branch_count; unit="count", lower=0, upper=8),
        ),
        unsupported_phenomena=(
            :finite_element_air_gap,
            :magnetic_hysteresis,
            :thermal_aging,
            :vendor_controller_equivalence,
            :protection_certification,
        ),
    )
    inventory = DynamicStateInventory(
        differential=(
            :electrical_flux_linkage,
            :shaft_angle,
            :shaft_speed,
            :control_state,
        ),
        algebraic=(
            :terminal_voltage,
            :terminal_current,
            :winding_current,
            :electromagnetic_torque,
        ),
        discrete=(
            :family_identity,
            :topology_identity,
            :limiter_mode,
            :event_counter,
        ),
        delayed_history=(
            :trapezoidal_companion_history,
            :sampled_control_input,
            :energy_accumulator,
        ),
        scheduler=(:control_task_calendar, :machine_event_calendar),
        random=(),
    )
    return ScientificModelContract(
        :modern_emt_machine,
        :phase_domain_electromechanical_machine;
        owner="AIMORA.ModernMachines",
        maturity=:implemented,
        fidelity=FieldCoupledDetailed,
        validity_domain=validity,
        state_inventory=inventory,
        inputs=(
            ContractQuantity(:phase_terminal_voltage; unit="V", orientation="phase-to-neutral"),
            ContractQuantity(:field_voltage; unit="V", orientation="into field winding"),
            ContractQuantity(:rotor_port_voltage; unit="V", orientation="into rotor port"),
            ContractQuantity(:mechanical_torque; unit="N*m", orientation="positive into shaft"),
        ),
        outputs=(
            ContractQuantity(:phase_terminal_current; unit="A", orientation="positive into machine"),
            ContractQuantity(:electromagnetic_torque; unit="N*m", orientation="positive rotor direction"),
            ContractQuantity(:stored_energy; unit="J"),
            ContractQuantity(:energy_residual; unit="J"),
        ),
        assumptions=(
            "power-invariant zero/d/q transformation at the accepted rotor angle",
            "fixed-step trapezoidal electrical state with exact analytic terminal tangent",
            "convex coenergy saturation and passive winding, rotor, shaft, and damping data",
        ),
        mutation_order=(
            :validate_identity,
            :capture_transaction,
            :apply_due_events,
            :sample_controls,
            :solve_electrical_trial,
            :solve_shaft_trial,
            :verify_energy_and_domain,
            :accept_once,
        ),
    )
end

function _machine_layout(specification::ModernMachineSpecification)
    family = specification.family
    synchronous = family in (
        WoundFieldSynchronousMachine,
        SynchronousCondenserMachine,
    )
    permanent = family === PermanentMagnetSynchronousMachine
    induction = family in (
        CageInductionMachine,
        WoundRotorInductionMachine,
        DoublyFedInductionMachine,
    )
    d_extra = synchronous ? 2 : induction ? length(specification.rotor_branches) : 0
    q_extra = synchronous ? 1 : induction ? length(specification.rotor_branches) : 0
    permanent && (d_extra = 0; q_extra = 0)
    d_indices = 2:(2 + d_extra)
    q_start = last(d_indices) + 1
    q_indices = q_start:(q_start + q_extra)
    field_index = synchronous ? first(d_indices) + 1 : nothing
    d_damper_index = synchronous ? first(d_indices) + 2 : nothing
    q_damper_index = synchronous ? first(q_indices) + 1 : nothing
    rotor_d = induction ? collect((first(d_indices) + 1):last(d_indices)) : Int[]
    rotor_q = induction ? collect((first(q_indices) + 1):last(q_indices)) : Int[]
    return ModernMachineElectricalLayout(
        1,
        d_indices,
        q_indices,
        first(d_indices),
        first(q_indices),
        field_index,
        d_damper_index,
        q_damper_index,
        rotor_d,
        rotor_q,
    )
end

function _validate_machine_family!(specification::ModernMachineSpecification)
    family = specification.family
    branches = specification.rotor_branches
    electrical = specification.electrical
    if family in (WoundFieldSynchronousMachine, SynchronousCondenserMachine)
        isempty(branches) || _modern_machine_refusal(
            :unexpected_rotor_branch,
            :prepare,
            specification,
            "synchronous families use explicit field and damper circuits, not induction rotor branches",
        )
        all(>(0.0), (
            electrical.field_resistance_ohm,
            electrical.field_leakage_inductance_h,
            electrical.d_damper_resistance_ohm,
            electrical.d_damper_leakage_inductance_h,
            electrical.q_damper_resistance_ohm,
            electrical.q_damper_leakage_inductance_h,
        )) || _modern_machine_refusal(
            :missing_field_or_damper_data,
            :prepare,
            specification,
            "wound-field synchronous families require positive field and d/q damper data",
        )
    elseif family === PermanentMagnetSynchronousMachine
        isempty(branches) || _modern_machine_refusal(
            :unexpected_rotor_branch,
            :prepare,
            specification,
            "permanent-magnet family does not admit induction rotor branches",
        )
        electrical.permanent_magnet_flux_wb > 0.0 || _modern_machine_refusal(
            :missing_permanent_flux,
            :prepare,
            specification,
            "permanent-magnet family requires positive permanent flux",
        )
    else
        1 <= length(branches) <= 8 || _modern_machine_refusal(
            :invalid_rotor_branch_count,
            :prepare,
            specification,
            "induction families require one through eight passive rotor branches",
        )
        exposed = count(getfield.(branches, :terminal_exposed))
        if family === CageInductionMachine
            exposed == 0 || _modern_machine_refusal(
                :exposed_cage_terminal,
                :prepare,
                specification,
                "cage induction rotor branches cannot expose terminals",
            )
        else
            exposed == 1 || _modern_machine_refusal(
                :missing_rotor_port,
                :prepare,
                specification,
                "wound-rotor and doubly-fed machines require exactly one exposed rotor branch",
            )
        end
    end
    specification.operating_mode === MachineCondenserMode &&
        family !== SynchronousCondenserMachine && _modern_machine_refusal(
            :invalid_condenser_mode,
            :prepare,
            specification,
            "condenser mode is reserved for the synchronous-condenser family",
        )
    family === SynchronousCondenserMachine &&
        specification.operating_mode !== MachineCondenserMode && _modern_machine_refusal(
            :missing_condenser_mode,
            :prepare,
            specification,
            "synchronous-condenser family requires condenser operating mode",
        )
    return nothing
end

function _validate_machine_shaft!(specification::ModernMachineSpecification)
    ids = getfield.(specification.shaft_masses, :id)
    edges = Set{Tuple{Symbol,Symbol}}()
    for coupling in specification.shaft_couplings
        coupling.left_mass in ids && coupling.right_mass in ids || _modern_machine_refusal(
            :unknown_shaft_mass,
            :prepare,
            specification,
            "shaft coupling references an unknown mass",
        )
        edge = coupling.left_mass < coupling.right_mass ?
            (coupling.left_mass, coupling.right_mass) :
            (coupling.right_mass, coupling.left_mass)
        edge in edges && _modern_machine_refusal(
            :duplicate_shaft_coupling,
            :prepare,
            specification,
            "shaft graph contains a duplicate coupling",
        )
        push!(edges, edge)
    end
    length(ids) == 1 && return nothing
    adjacency = Dict(id => Symbol[] for id in ids)
    for coupling in specification.shaft_couplings
        push!(adjacency[coupling.left_mass], coupling.right_mass)
        push!(adjacency[coupling.right_mass], coupling.left_mass)
    end
    visited = Set{Symbol}([first(ids)])
    frontier = Symbol[first(ids)]
    while !isempty(frontier)
        id = pop!(frontier)
        for neighbor in adjacency[id]
            neighbor in visited && continue
            push!(visited, neighbor)
            push!(frontier, neighbor)
        end
    end
    length(visited) == length(ids) || _modern_machine_refusal(
        :disconnected_shaft,
        :prepare,
        specification,
        "every shaft mass must belong to one connected graph",
    )
    return nothing
end

function _machine_axis_inductance(leakages::Vector{Float64}, magnetizing::Float64)
    matrix = Diagonal(leakages) + magnetizing .* ones(length(leakages), length(leakages))
    minimum(eigvals(Symmetric(Matrix(matrix)))) > 0.0 || throw(ArgumentError(
        "machine axis inductance matrix must be positive definite",
    ))
    return Matrix(matrix)
end

function _machine_matrices(
    specification::ModernMachineSpecification,
    layout::ModernMachineElectricalLayout,
)
    electrical = specification.electrical
    state_count = last(layout.q_indices)
    inductance = zeros(state_count, state_count)
    resistance = zeros(state_count)
    inductance[layout.zero_index, layout.zero_index] =
        electrical.zero_sequence_inductance_h
    resistance[layout.zero_index] = electrical.stator_resistance_ohm
    d_leakages = Float64[electrical.stator_d_leakage_inductance_h]
    q_leakages = Float64[electrical.stator_q_leakage_inductance_h]
    d_resistance = Float64[electrical.stator_resistance_ohm]
    q_resistance = Float64[electrical.stator_resistance_ohm]
    if specification.family in (
        WoundFieldSynchronousMachine,
        SynchronousCondenserMachine,
    )
        append!(d_leakages, (
            electrical.field_leakage_inductance_h,
            electrical.d_damper_leakage_inductance_h,
        ))
        push!(q_leakages, electrical.q_damper_leakage_inductance_h)
        append!(d_resistance, (
            electrical.field_resistance_ohm,
            electrical.d_damper_resistance_ohm,
        ))
        push!(q_resistance, electrical.q_damper_resistance_ohm)
    elseif specification.family in (
        CageInductionMachine,
        WoundRotorInductionMachine,
        DoublyFedInductionMachine,
    )
        append!(d_leakages, getfield.(specification.rotor_branches, :leakage_inductance_h))
        append!(q_leakages, getfield.(specification.rotor_branches, :leakage_inductance_h))
        append!(d_resistance, getfield.(specification.rotor_branches, :resistance_ohm))
        append!(q_resistance, getfield.(specification.rotor_branches, :resistance_ohm))
    end
    inductance[layout.d_indices, layout.d_indices] .= _machine_axis_inductance(
        d_leakages,
        electrical.d_axis_magnetizing_inductance_h,
    )
    inductance[layout.q_indices, layout.q_indices] .= _machine_axis_inductance(
        q_leakages,
        electrical.q_axis_magnetizing_inductance_h,
    )
    resistance[layout.d_indices] .= d_resistance
    resistance[layout.q_indices] .= q_resistance
    offset = zeros(state_count)
    if specification.family === PermanentMagnetSynchronousMachine
        offset[layout.stator_d_index] = electrical.permanent_magnet_flux_wb
    end
    return inductance, resistance, offset
end

function _machine_terminal_maps(layout::ModernMachineElectricalLayout)
    phase_to_terminal = [
        1.0 0.0 0.0 -1.0
        0.0 1.0 0.0 -1.0
        0.0 0.0 1.0 -1.0
    ]
    return phase_to_terminal, transpose(phase_to_terminal)
end

function _machine_signature(
    specification::ModernMachineSpecification,
    inductance,
    resistance,
    offset,
)
    io = IOBuffer()
    println(io, specification.id)
    println(io, _MODERN_MACHINE_FAMILY_IDS[specification.family])
    println(io, _MODERN_MACHINE_MODE_IDS[specification.operating_mode])
    println(io, specification.pole_pairs)
    println(io, join(specification.phase_order, ','))
    println(io, join(specification.terminal_order, ','))
    for value in (inductance..., resistance..., offset...)
        println(io, bitstring(Float64(value)))
    end
    for branch in specification.rotor_branches
        println(io, branch.id, ':', bitstring(branch.resistance_ohm), ':',
            bitstring(branch.leakage_inductance_h), ':', branch.terminal_exposed)
    end
    for mass in specification.shaft_masses
        println(io, mass.id, ':', bitstring(mass.inertia_kg_m2), ':',
            bitstring(mass.damping_nm_s_per_rad))
    end
    println(io, specification.provenance.source)
    println(io, specification.provenance.units)
    println(io, specification.uncertainty)
    println(io, specification.validity_domain)
    return bytes2hex(sha256(take!(io)))
end

function _machine_initial_flux(
    specification::ModernMachineSpecification,
    layout::ModernMachineElectricalLayout,
    inductance::Matrix{Float64},
    offset::Vector{Float64},
)
    current = zeros(size(inductance, 1))
    angle = first(specification.shaft_masses).initial_angle_rad * specification.pole_pairs
    transform = machine_phase_transform(angle)
    current[[layout.zero_index, layout.stator_d_index, layout.stator_q_index]] .=
        transform * collect(specification.initial_phase_current_a)
    if layout.field_index !== nothing
        electrical = specification.electrical
        current[layout.field_index] = specification.initial_field_voltage_v /
            electrical.field_resistance_ohm
    end
    return inductance * current + offset
end

function prepare_modern_machine(specification::ModernMachineSpecification)
    _validate_machine_family!(specification)
    _validate_machine_shaft!(specification)
    layout = _machine_layout(specification)
    inductance, resistance, offset = _machine_matrices(specification, layout)
    terminal_voltage_map, terminal_current_map = _machine_terminal_maps(layout)
    initial_flux = _machine_initial_flux(
        specification,
        layout,
        inductance,
        offset,
    )
    signature = _machine_signature(
        specification,
        inductance,
        resistance,
        offset,
    )
    return ModernMachinePreparation(
        specification,
        layout,
        inductance,
        inv(Symmetric(inductance)),
        resistance,
        offset,
        terminal_voltage_map,
        terminal_current_map,
        initial_flux,
        signature,
    )
end

modern_machine_signature(preparation::ModernMachinePreparation) =
    preparation.deterministic_signature_sha256

function modern_machine_readiness(specification::ModernMachineSpecification)
    try
        preparation = prepare_modern_machine(specification)
        return ModernMachineReadiness(
            true,
            :ready,
            _MODERN_MACHINE_FAMILY_IDS[specification.family],
            length(preparation.initial_flux_wb),
            length(specification.rotor_branches),
            length(specification.shaft_masses),
            specification.controls.enabled,
            preparation.deterministic_signature_sha256,
            String[
                "fixed-step phase-domain generic machine; no vendor or certification claim",
                "thermal, hysteresis, spatial harmonic, protection, and arbitrary converter behavior are unsupported",
            ],
        )
    catch error
        if error isa ModernMachineRefusal || error isa ArgumentError ||
            error isa DimensionMismatch
            return ModernMachineReadiness(
                false,
                error isa ModernMachineRefusal ? error.code : :invalid_input,
                _MODERN_MACHINE_FAMILY_IDS[specification.family],
                0,
                length(specification.rotor_branches),
                length(specification.shaft_masses),
                specification.controls.enabled,
                "",
                String[sprint(showerror, error)],
            )
        end
        rethrow()
    end
end
