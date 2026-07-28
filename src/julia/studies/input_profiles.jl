module StudyInputProfiles

using ..StudyInputs:
    InputSpec,
    StudyInputProfile,
    required_keys,
    optional_keys,
    missing_required_inputs,
    validate_study_inputs_issues,
    validate_study_inputs

export InputSpec,
       StudyInputProfile,
       input_profile,
       available_input_profiles,
       required_keys,
       optional_keys,
       missing_required_inputs,
       validate_study_inputs_issues,
       validate_study_inputs

spec(key::Symbol, label::AbstractString, description::AbstractString; unit = nothing) =
    InputSpec(key; label = label, unit = unit, description = description)

const PROFILES = Dict{Symbol,StudyInputProfile}(
    :power_flow => StudyInputProfile(
        :power_flow,
        [
            spec(:buses, "Buses", "Bus names, voltage bases, and bus types."),
            spec(:branches, "Branches", "Lines, cables, and transformers needed to connect the network."),
            spec(:sources, "Sources", "Slack/grid sources or generator voltage setpoints."),
            spec(:loads, "Loads", "Active and reactive power demand.", unit = "MW/Mvar or pu"),
        ],
        [
            spec(:controls, "Controls", "Tap changers, switched shunts, and voltage control settings."),
            spec(:limits, "Limits", "Voltage, loading, and reactive limits."),
            spec(:scenarios, "Scenarios", "Load/generation cases and operating variants."),
        ];
        notes = "Power flow needs the electrical network and operating injections, not EMT timestep data.",
    ),

    :short_circuit => StudyInputProfile(
        :short_circuit,
        [
            spec(:buses, "Buses", "Fault locations and voltage bases."),
            spec(:source_equivalents, "Source equivalents", "Utility or generator short-circuit equivalents."),
            spec(:impedances, "Impedances", "Equipment positive/negative/zero sequence or phase impedances."),
            spec(:fault_definition, "Fault definition", "Fault type, location, grounding, and prefault assumptions."),
        ],
        [
            spec(:standards, "Standards", "IEC, ANSI, or project-specific calculation method."),
            spec(:motor_contribution, "Motor contribution", "Motor and converter fault contribution assumptions."),
            spec(:protective_devices, "Protective devices", "Devices to report duty against."),
        ];
        notes = "Short-circuit studies should not require load-flow controls unless the selected method needs them.",
    ),

    :protection => StudyInputProfile(
        :protection,
        [
            spec(:protective_devices, "Protective devices", "Relays, fuses, breakers, reclosers, settings, and curves."),
            spec(:fault_duties, "Fault duties", "Minimum and maximum fault currents at protected locations.", unit = "A or kA"),
            spec(:equipment_limits, "Equipment limits", "Protected cable, transformer, bus, motor, or generator damage limits."),
        ],
        [
            spec(:full_network, "Full network", "Network model for deriving fault currents instead of entering them."),
            spec(:load_currents, "Load currents", "Normal and emergency load currents.", unit = "A"),
            spec(:ct_vt, "CT/VT data", "Instrument transformer ratios, classes, and saturation limits."),
            spec(:coordination_rules, "Coordination rules", "Margins, grading intervals, and selectivity rules."),
            spec(:manufacturer_curves, "Manufacturer curves", "Device curve libraries and tolerance bands."),
            spec(:standards, "Standards", "IEC/IEEE/company protection requirements."),
        ];
        notes = "Protection can run from device data plus fault duties. A full EMT or power-flow model is optional.",
    ),

    :arc_flash => StudyInputProfile(
        :arc_flash,
        [
            spec(:equipment, "Equipment", "Enclosure/equipment type, voltage, gap, and working distance."),
            spec(:fault_current, "Fault current", "Available arcing or bolted fault current.", unit = "A or kA"),
            spec(:clearing_time, "Clearing time", "Protective-device clearing time.", unit = "s"),
            spec(:standard, "Standard", "Arc-flash method or project standard."),
        ],
        [
            spec(:protective_device, "Protective device", "Device curve used to derive clearing time."),
            spec(:electrode_configuration, "Electrode configuration", "Configuration used by the selected method."),
            spec(:labels, "Labels", "Label formatting and PPE category settings."),
        ];
        notes = "Arc flash can use protection results, but it must also accept direct fault current and clearing time.",
    ),

    :pv_cable_sizing => StudyInputProfile(
        :pv_cable_sizing,
        [
            spec(:circuit_current, "Circuit current", "Design current for the PV DC or AC circuit.", unit = "A"),
            spec(:voltage, "Voltage", "Circuit voltage.", unit = "V"),
            spec(:route_length, "Route length", "One-way cable route length.", unit = "m"),
            spec(:installation_method, "Installation method", "Tray, conduit, buried, duct bank, or free air."),
            spec(:conductor_material, "Conductor material", "Copper or aluminum."),
        ],
        [
            spec(:ambient_temperature, "Ambient temperature", "Ambient or soil temperature.", unit = "degC"),
            spec(:voltage_drop_limit, "Voltage-drop limit", "Maximum permitted voltage drop.", unit = "%"),
            spec(:grouping_derating, "Grouping derating", "Circuit grouping and spacing derating factors."),
            spec(:short_circuit_duty, "Short-circuit duty", "Fault current and clearing time for withstand check."),
            spec(:standard, "Standard", "NEC/IEC/project sizing method."),
        ];
        notes = "PV cable sizing should not require a full network model.",
    ),

    :cable_thermal_capacity => StudyInputProfile(
        :cable_thermal_capacity,
        [
            spec(:cable_construction, "Cable construction", "Conductor, insulation, screen/sheath, and jacket data."),
            spec(:installation_environment, "Installation environment", "Air, soil, duct, tray, or tunnel conditions."),
            spec(:load_profile, "Load profile", "Current versus time.", unit = "A"),
            spec(:thermal_limits, "Thermal limits", "Normal, emergency, and short-time temperature limits.", unit = "degC"),
        ],
        [
            spec(:soil_data, "Soil data", "Thermal resistivity, moisture, and burial details."),
            spec(:ambient_profile, "Ambient profile", "Ambient temperature versus time."),
            spec(:adjacent_heat_sources, "Adjacent heat sources", "Nearby loaded circuits or external heat sources."),
            spec(:standard, "Standard", "IEC/IEEE/CIGRE/project method."),
        ];
        notes = "Thermal cable capacity needs cable and environment data, not protection relay settings.",
    ),

    :grounding_grid => StudyInputProfile(
        :grounding_grid,
        [
            spec(:soil_model, "Soil model", "Single or multilayer soil resistivity model.", unit = "ohm-m"),
            spec(:fault_current, "Ground fault current", "Current injected into the grounding system.", unit = "A"),
            spec(:clearing_time, "Clearing time", "Fault clearing time.", unit = "s"),
            spec(:grid_geometry, "Grid geometry", "Conductor layout, rods, depth, and area."),
        ],
        [
            spec(:surface_layer, "Surface layer", "Crushed rock/asphalt surface layer data."),
            spec(:body_weight, "Body weight", "Assumed body weight for touch/step limit method."),
            spec(:standard, "Standard", "IEEE 80 or project method."),
        ];
        notes = "Grounding studies should not require load-flow or EMT model parameters.",
    ),

    :induction_machine => StudyInputProfile(
        :induction_machine,
        [
            spec(:nameplate, "Nameplate", "Rated power, voltage, frequency, speed, and efficiency."),
            spec(:equivalent_circuit, "Equivalent circuit", "Machine R/X parameters or locked-rotor data."),
            spec(:mechanical_load, "Mechanical load", "Torque-speed or load torque model."),
        ],
        [
            spec(:thermal_model, "Thermal model", "Thermal limits and time constants."),
            spec(:starter, "Starter", "DOL, soft starter, VFD, or reduced-voltage starter."),
            spec(:supply_network, "Supply network", "Network Thevenin or full network connection."),
        ];
        notes = "Machine studies need machine data first; the full system model is optional for standalone checks.",
    ),

    :emt => StudyInputProfile(
        :emt,
        [
            spec(:timestep, "Timestep", "Fixed simulation timestep.", unit = "s"),
            spec(:duration, "Duration", "Simulation duration.", unit = "s"),
            spec(:network_elements, "Network elements", "Branches, sources, switches, and injections."),
            spec(:initial_conditions, "Initial conditions", "Initial voltages, currents, and model states."),
        ],
        [
            spec(:events, "Events", "Switching, faults, references, and disturbances."),
            spec(:controls, "Controls", "Control blocks and converter controllers."),
            spec(:output_channels, "Output channels", "Signals to record."),
        ];
        notes = "EMT requires timestep state; protection and sizing studies should not inherit this requirement.",
    ),
)

available_input_profiles() = copy(PROFILES)

function input_profile(study::Symbol)
    profile = get(PROFILES, study, nothing)
    profile === nothing && error("No input profile is defined for study $(study). Add one before implementing that study.")
    return profile
end

end
