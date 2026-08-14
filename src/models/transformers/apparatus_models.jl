export AbstractTransformerTierModel,
       TransformerTerminalMatrices,
       LowFrequencyTransformerModel,
       BCTRANTransformerModel,
       HybridTransformerModel,
       MagneticEquivalentCircuitModel,
       WidebandTransformerModel,
       TransformerLadderBranch,
       GreyBoxTransformerModel,
       WhiteBoxTransformerModel,
       TransformerRuntimeSettings,
       TransformerInitializationMode,
       DeenergizedTransformerInitialization,
       SpecifiedTransformerInitialization,
       SinusoidalTransformerOperatingPoint,
       TransformerApparatusSpecification,
       TransformerApparatusPreparation,
       TransformerApparatusReadiness,
       TransformerApparatusRefusal,
       prepare_transformer_apparatus,
       transformer_apparatus_readiness,
       transformer_apparatus_signature

abstract type AbstractTransformerTierModel end

@enum TransformerInitializationMode begin
    DeenergizedTransformerInitialization
    SpecifiedTransformerInitialization
    SinusoidalTransformerOperatingPoint
end

const _TRANSFORMER_INITIALIZATION_MODE_IDS = Dict(
    DeenergizedTransformerInitialization => :deenergized,
    SpecifiedTransformerInitialization => :specified_state,
    SinusoidalTransformerOperatingPoint => :sinusoidal_operating_point,
)

struct TransformerApparatusRefusal <: Exception
    code::Symbol
    operation::Symbol
    apparatus::Symbol
    tier::TransformerApparatusTier
    message::String
    diagnostics::NamedTuple
end

Base.showerror(io::IO, refusal::TransformerApparatusRefusal) = print(
    io,
    String(refusal.code),
    " during ",
    String(refusal.operation),
    " for transformer apparatus ",
    String(refusal.apparatus),
    " (",
    String(_TRANSFORMER_TIER_IDS[refusal.tier]),
    "): ",
    refusal.message,
)

function _transformer_refusal(
    code::Symbol,
    operation::Symbol,
    apparatus::Symbol,
    tier::TransformerApparatusTier,
    message::AbstractString;
    diagnostics=NamedTuple(),
)
    throw(TransformerApparatusRefusal(
        code,
        operation,
        apparatus,
        tier,
        String(message),
        diagnostics,
    ))
end

function _transformer_symmetric_matrix(
    values,
    dimension::Int,
    label::AbstractString;
    positive_semidefinite::Bool=true,
    positive_definite::Bool=false,
)
    matrix = Matrix{Float64}(values)
    size(matrix) == (dimension, dimension) || throw(DimensionMismatch(
        "transformer $label matrix must be $dimension by $dimension",
    ))
    all(isfinite, matrix) || throw(ArgumentError(
        "transformer $label matrix entries must be finite",
    ))
    scale = max(maximum(abs, matrix; init=0.0), 1.0)
    tolerance = 256.0 * eps(Float64) * scale
    maximum(abs, matrix - transpose(matrix); init=0.0) <= tolerance ||
        throw(ArgumentError("transformer $label matrix must be symmetric"))
    symmetric_matrix = 0.5 .* (matrix .+ transpose(matrix))
    minimum_eigenvalue = minimum(eigvals(Symmetric(symmetric_matrix)); init=Inf)
    positive_semidefinite && minimum_eigenvalue < -tolerance && throw(ArgumentError(
        "transformer $label matrix must be positive semidefinite",
    ))
    positive_definite && minimum_eigenvalue <= tolerance && throw(ArgumentError(
        "transformer $label matrix must be positive definite",
    ))
    return symmetric_matrix, minimum_eigenvalue
end

"""Complete reciprocal coil-coordinate R-L plus parallel C-G terminal matrices."""
struct TransformerTerminalMatrices
    resistance_ohm::Matrix{Float64}
    inductance_h::Matrix{Float64}
    capacitance_f::Matrix{Float64}
    conductance_s::Matrix{Float64}
    minimum_resistance_eigenvalue_ohm::Float64
    minimum_inductance_eigenvalue_h::Float64
    minimum_capacitance_eigenvalue_f::Float64
    minimum_conductance_eigenvalue_s::Float64
    deterministic_signature_sha256::String
end

function _transformer_matrix_signature(matrices...)
    io = IOBuffer()
    for matrix in matrices
        println(io, size(matrix, 1), 'x', size(matrix, 2))
        for value in matrix
            println(io, bitstring(Float64(value)))
        end
    end
    return bytes2hex(sha256(take!(io)))
end

function TransformerTerminalMatrices(
    resistance_ohm,
    inductance_h;
    capacitance_f=nothing,
    conductance_s=nothing,
)
    dimension = size(resistance_ohm, 1)
    dimension > 0 || throw(ArgumentError(
        "transformer terminal matrices require at least one coil",
    ))
    resistance, minimum_resistance = _transformer_symmetric_matrix(
        resistance_ohm,
        dimension,
        "resistance",
    )
    inductance, minimum_inductance = _transformer_symmetric_matrix(
        inductance_h,
        dimension,
        "inductance",
    )
    maximum(abs, inductance; init=0.0) > 0.0 || throw(ArgumentError(
        "transformer inductance matrix must contain physical storage",
    ))
    capacitance_source = capacitance_f === nothing ? zeros(dimension, dimension) :
        capacitance_f
    conductance_source = conductance_s === nothing ? zeros(dimension, dimension) :
        conductance_s
    capacitance, minimum_capacitance = _transformer_symmetric_matrix(
        capacitance_source,
        dimension,
        "capacitance",
    )
    conductance, minimum_conductance = _transformer_symmetric_matrix(
        conductance_source,
        dimension,
        "conductance",
    )
    signature = _transformer_matrix_signature(
        resistance,
        inductance,
        capacitance,
        conductance,
    )
    return TransformerTerminalMatrices(
        resistance,
        inductance,
        capacitance,
        conductance,
        minimum_resistance,
        minimum_inductance,
        minimum_capacitance,
        minimum_conductance,
        signature,
    )
end

struct LowFrequencyTransformerModel <: AbstractTransformerTierModel
    matrices::TransformerTerminalMatrices
    source_kind::Symbol

    function LowFrequencyTransformerModel(
        matrices::TransformerTerminalMatrices;
        source_kind::Symbol=:explicit_si_matrices,
    )
        source_kind in (:explicit_si_matrices, :short_circuit_tests, :design_geometry, :generic) ||
            throw(ArgumentError("low-frequency transformer source kind is unsupported"))
        return new(matrices, source_kind)
    end
end

function LowFrequencyTransformerModel(result::SaturableTransformerParameterResult)
    result.physical_checks_passed || throw(ArgumentError(
        "saturable transformer parameter result must pass physical checks",
    ))
    return LowFrequencyTransformerModel(
        TransformerTerminalMatrices(
            result.branch_resistance_matrix_ohm,
            result.physical_inductance_matrix_h,
        );
        source_kind=:short_circuit_tests,
    )
end

struct BCTRANTransformerModel <: AbstractTransformerTierModel
    matrices::TransformerTerminalMatrices
    positive_pair_reconstruction_residual::Float64
    zero_pair_reconstruction_residual::Float64
    inverse_reconstruction_residual::Float64

    function BCTRANTransformerModel(
        matrices::TransformerTerminalMatrices;
        positive_pair_reconstruction_residual::Real,
        zero_pair_reconstruction_residual::Real,
        inverse_reconstruction_residual::Real,
    )
        residuals = Float64.(Tuple((
            positive_pair_reconstruction_residual,
            zero_pair_reconstruction_residual,
            inverse_reconstruction_residual,
        )))
        all(value -> isfinite(value) && value >= 0.0, residuals) ||
            throw(ArgumentError("BCTRAN reconstruction residuals must be finite and nonnegative"))
        return new(matrices, residuals...)
    end
end

function BCTRANTransformerModel(result::MultiphaseTransformerParameterResult)
    result.physical_checks_passed || throw(ArgumentError(
        "BCTRAN transformer parameter result must pass physical checks",
    ))
    frequency = result.input.frequency_hz
    frequency > 0.0 || throw(ArgumentError("BCTRAN result frequency must be positive"))
    inductance = result.reactance_matrix_ohm ./ (2.0 * pi * frequency)
    return BCTRANTransformerModel(
        TransformerTerminalMatrices(result.resistance_matrix_ohm, inductance);
        positive_pair_reconstruction_residual=
            result.positive_pair_reconstruction_residual,
        zero_pair_reconstruction_residual=result.zero_pair_reconstruction_residual,
        inverse_reconstruction_residual=result.inverse_reconstruction_residual,
    )
end

struct HybridTransformerModel <: AbstractTransformerTierModel
    leakage::TransformerTerminalMatrices
    magnetic_graph::TransformerMagneticGraph
    winding_loss_state_matrix_per_s::Matrix{Float64}
    winding_loss_input_matrix::Matrix{Float64}
    winding_loss_output_matrix_ohm_per_s::Matrix{Float64}
    winding_loss_direct_ohm::Matrix{Float64}
    winding_loss_storage_matrix_j::Matrix{Float64}
    winding_loss_dissipation_matrix_j_per_s::Matrix{Float64}
    continuous_passivity_margin_ohm::Float64
end

function HybridTransformerModel(
    leakage::TransformerTerminalMatrices,
    magnetic_graph::TransformerMagneticGraph;
    winding_loss_state_matrix_per_s=zeros(0, 0),
    winding_loss_input_matrix=zeros(0, size(leakage.resistance_ohm, 1)),
    winding_loss_output_matrix_ohm_per_s=zeros(
        size(leakage.resistance_ohm, 1),
        0,
    ),
    winding_loss_direct_ohm=zeros(size(leakage.resistance_ohm)),
    winding_loss_storage_matrix_j=zeros(
        size(winding_loss_state_matrix_per_s, 1),
        size(winding_loss_state_matrix_per_s, 1),
    ),
    continuous_passivity_margin_ohm::Real=0.0,
)
    coil_count = size(leakage.resistance_ohm, 1)
    size(magnetic_graph.winding_turns, 2) == coil_count || throw(DimensionMismatch(
        "hybrid transformer magnetic turns must cover every electrical coil",
    ))
    state = Matrix{Float64}(winding_loss_state_matrix_per_s)
    input = Matrix{Float64}(winding_loss_input_matrix)
    output = Matrix{Float64}(winding_loss_output_matrix_ohm_per_s)
    direct, minimum_direct = _transformer_symmetric_matrix(
        winding_loss_direct_ohm,
        coil_count,
        "winding-loss direct resistance",
    )
    state_count = size(state, 1)
    size(state) == (state_count, state_count) &&
        size(input) == (state_count, coil_count) &&
        size(output) == (coil_count, state_count) || throw(DimensionMismatch(
            "hybrid transformer winding-loss state-space dimensions disagree",
        ))
    all(matrix -> all(isfinite, matrix), (state, input, output)) ||
        throw(ArgumentError("hybrid transformer winding-loss matrices must be finite"))
    state_count == 0 || maximum(real, eigvals(state); init=-Inf) < 0.0 ||
        throw(ArgumentError("hybrid transformer winding-loss state must be strictly stable"))
    storage, minimum_storage = _transformer_symmetric_matrix(
        winding_loss_storage_matrix_j,
        state_count,
        "winding-loss storage",
    )
    state_count == 0 || minimum_storage > 256.0 * eps(Float64) ||
        throw(ArgumentError(
            "hybrid transformer winding-loss storage matrix must be positive definite",
        ))
    passivity_scale = max(
        maximum(abs, storage * input; init=0.0),
        maximum(abs, transpose(output); init=0.0),
        1.0,
    )
    maximum(abs, storage * input - transpose(output); init=0.0) <=
        1.0e-10 * passivity_scale || throw(ArgumentError(
            "hybrid transformer winding-loss realization violates the storage/input-output identity",
        ))
    dissipation = -(transpose(state) * storage + storage * state)
    dissipation, minimum_dissipation = _transformer_symmetric_matrix(
        dissipation,
        state_count,
        "winding-loss state dissipation",
    )
    minimum_dissipation >= -256.0 * eps(Float64) || throw(ArgumentError(
        "hybrid transformer winding-loss state dissipation must be passive",
    ))
    margin = Float64(continuous_passivity_margin_ohm)
    isfinite(margin) && margin >= 0.0 || throw(ArgumentError(
        "hybrid transformer winding-loss passivity margin must be finite and nonnegative",
    ))
    minimum_direct >= -256.0 * eps(Float64) || throw(ArgumentError(
        "hybrid transformer winding-loss direct response must be passive",
    ))
    return HybridTransformerModel(
        leakage,
        magnetic_graph,
        state,
        input,
        output,
        direct,
        storage,
        dissipation,
        margin,
    )
end

struct MagneticEquivalentCircuitModel <: AbstractTransformerTierModel
    winding_resistance_ohm::Matrix{Float64}
    leakage_inductance_h::Matrix{Float64}
    terminal_capacitance_f::Matrix{Float64}
    terminal_conductance_s::Matrix{Float64}
    magnetic_graph::TransformerMagneticGraph
end

function MagneticEquivalentCircuitModel(
    winding_resistance_ohm,
    leakage_inductance_h,
    magnetic_graph::TransformerMagneticGraph;
    terminal_capacitance_f=nothing,
    terminal_conductance_s=nothing,
)
    matrices = TransformerTerminalMatrices(
        winding_resistance_ohm,
        leakage_inductance_h;
        capacitance_f=terminal_capacitance_f,
        conductance_s=terminal_conductance_s,
    )
    size(magnetic_graph.winding_turns, 2) == size(matrices.resistance_ohm, 1) ||
        throw(DimensionMismatch(
            "MEC transformer magnetic turns must cover every electrical coil",
        ))
    return MagneticEquivalentCircuitModel(
        matrices.resistance_ohm,
        matrices.inductance_h,
        matrices.capacitance_f,
        matrices.conductance_s,
        magnetic_graph,
    )
end

"""Complete ordered stable positive-real terminal admittance state-space model."""
struct WidebandTransformerModel <: AbstractTransformerTierModel
    state_matrix_per_s::Matrix{Float64}
    input_matrix::Matrix{Float64}
    output_matrix_s_per_s::Matrix{Float64}
    direct_admittance_s::Matrix{Float64}
    storage_matrix_j::Matrix{Float64}
    dissipation_matrix_j_per_s::Matrix{Float64}
    port_order::Vector{Symbol}
    frequency_band_hz::NTuple{2,Float64}
    continuous_passivity_margin_s::Float64
    enforcement_perturbation_relative::Float64
    source_response_sha256::String
    deterministic_signature_sha256::String
end

function WidebandTransformerModel(
    state_matrix_per_s,
    input_matrix,
    output_matrix_s_per_s,
    direct_admittance_s;
    storage_matrix_j,
    port_order,
    frequency_band_hz,
    continuous_passivity_margin_s::Real,
    enforcement_perturbation_relative::Real,
    source_response_sha256::AbstractString,
)
    state = Matrix{Float64}(state_matrix_per_s)
    input = Matrix{Float64}(input_matrix)
    output = Matrix{Float64}(output_matrix_s_per_s)
    direct = Matrix{Float64}(direct_admittance_s)
    state_count = size(state, 1)
    port_count = size(direct, 1)
    2 <= port_count <= 30 && size(direct) == (port_count, port_count) ||
        throw(ArgumentError("wideband transformer requires 2 through 30 square ports"))
    size(state) == (state_count, state_count) && state_count > 0 ||
        throw(DimensionMismatch("wideband transformer state matrix must be nonempty and square"))
    size(input) == (state_count, port_count) &&
        size(output) == (port_count, state_count) || throw(DimensionMismatch(
            "wideband transformer state-space input/output dimensions disagree",
        ))
    all(matrix -> all(isfinite, matrix), (state, input, output, direct)) ||
        throw(ArgumentError("wideband transformer state-space matrices must be finite"))
    maximum(real, eigvals(state); init=-Inf) < 0.0 || throw(ArgumentError(
        "wideband transformer state matrix must be strictly stable",
    ))
    direct, _ = _transformer_symmetric_matrix(
        direct,
        port_count,
        "wideband direct admittance",
    )
    storage, minimum_storage = _transformer_symmetric_matrix(
        storage_matrix_j,
        state_count,
        "wideband storage",
    )
    minimum_storage > 0.0 || throw(ArgumentError(
        "wideband transformer storage matrix must be positive definite",
    ))
    passivity_scale = max(
        maximum(abs, storage * input; init=0.0),
        maximum(abs, transpose(output); init=0.0),
        1.0,
    )
    maximum(abs, storage * input - transpose(output); init=0.0) <=
        1.0e-10 * passivity_scale || throw(ArgumentError(
            "wideband transformer realization violates the storage/input-output identity",
        ))
    dissipation = -(transpose(state) * storage + storage * state)
    dissipation, _ = _transformer_symmetric_matrix(
        dissipation,
        state_count,
        "wideband state dissipation",
    )
    ports = _connection_symbols(port_order, "wideband port order")
    length(ports) == port_count || throw(DimensionMismatch(
        "wideband transformer port identities must cover every port",
    ))
    band_values = Float64.(Tuple(frequency_band_hz))
    length(band_values) == 2 && 0.0 < band_values[1] < band_values[2] <= 2.0e6 ||
        throw(ArgumentError(
            "wideband transformer frequency band must increase inside (0,2 MHz]",
        ))
    margin = Float64(continuous_passivity_margin_s)
    isfinite(margin) && margin >= 0.0 || throw(ArgumentError(
        "wideband transformer continuous passivity margin must be nonnegative",
    ))
    perturbation = Float64(enforcement_perturbation_relative)
    isfinite(perturbation) && 0.0 <= perturbation <= 0.05 || throw(ArgumentError(
        "wideband transformer enforcement perturbation must lie from zero through five percent",
    ))
    source_signature = lowercase(String(source_response_sha256))
    occursin(r"^[0-9a-f]{64}$", source_signature) || throw(ArgumentError(
        "wideband transformer source response identity must be lowercase SHA-256",
    ))
    signature = _transformer_matrix_signature(
        state,
        input,
        output,
        direct,
        storage,
        dissipation,
    )
    io = IOBuffer()
    println(io, signature)
    println(io, join(String.(ports), ','))
    println(io, bitstring(band_values[1]))
    println(io, bitstring(band_values[2]))
    println(io, bitstring(margin))
    println(io, bitstring(perturbation))
    println(io, source_signature)
    deterministic_signature = bytes2hex(sha256(take!(io)))
    return WidebandTransformerModel(
        state,
        input,
        output,
        direct,
        storage,
        dissipation,
        ports,
        (band_values[1], band_values[2]),
        margin,
        perturbation,
        source_signature,
        deterministic_signature,
    )
end

struct TransformerLadderBranch
    id::Symbol
    from_node_index::Int
    to_node_index::Int
    resistance_ohm::Float64
    inductance_h::Float64

    function TransformerLadderBranch(
        id::Symbol,
        from_node_index::Integer,
        to_node_index::Integer;
        resistance_ohm::Real,
        inductance_h::Real,
    )
        id == Symbol("") && throw(ArgumentError(
            "transformer ladder branch identity must not be empty",
        ))
        from_node = Int(from_node_index)
        to_node = Int(to_node_index)
        from_node >= 0 && to_node >= 0 && from_node != to_node ||
            throw(ArgumentError(
                "transformer ladder branch nodes must be distinct and nonnegative",
            ))
        resistance = Float64(resistance_ohm)
        inductance = Float64(inductance_h)
        isfinite(resistance) && resistance >= 0.0 || throw(ArgumentError(
            "transformer ladder resistance must be finite and nonnegative",
        ))
        isfinite(inductance) && inductance >= 0.0 || throw(ArgumentError(
            "transformer ladder inductance must be finite and nonnegative",
        ))
        resistance > 0.0 || inductance > 0.0 || throw(ArgumentError(
            "transformer ladder branch must contain resistance or inductance",
        ))
        return new(id, from_node, to_node, resistance, inductance)
    end
end

struct GreyBoxTransformerModel <: AbstractTransformerTierModel
    node_order::Vector{Symbol}
    terminal_node_indices::Vector{Int}
    branches::Vector{TransformerLadderBranch}
    capacitance_f::Matrix{Float64}
    conductance_s::Matrix{Float64}
    source_response_sha256::String
    identification_residual_relative::Float64
    parameter_nonuniqueness::String
end

function GreyBoxTransformerModel(;
    node_order,
    terminal_node_indices,
    branches,
    capacitance_f,
    conductance_s,
    source_response_sha256::AbstractString,
    identification_residual_relative::Real,
    parameter_nonuniqueness::AbstractString,
)
    nodes = _connection_symbols(node_order, "grey-box node order")
    2 <= length(nodes) <= 128 || throw(ArgumentError(
        "grey-box transformer requires 2 through 128 represented nodes",
    ))
    terminal_indices = Int.(terminal_node_indices)
    length(terminal_indices) >= 2 && length(unique(terminal_indices)) == length(terminal_indices) &&
        all(index -> 1 <= index <= length(nodes), terminal_indices) || throw(ArgumentError(
            "grey-box terminal nodes must be unique represented node indices",
        ))
    branch_rows = TransformerLadderBranch[branches...]
    isempty(branch_rows) && throw(ArgumentError(
        "grey-box transformer requires at least one physical R-L branch",
    ))
    length(unique(getfield.(branch_rows, :id))) == length(branch_rows) ||
        throw(ArgumentError("grey-box branch identities must be unique"))
    all(branch -> max(branch.from_node_index, branch.to_node_index) <= length(nodes), branch_rows) ||
        throw(ArgumentError("grey-box branch node is outside its represented graph"))
    capacitance, _ = _transformer_symmetric_matrix(
        capacitance_f,
        length(nodes),
        "grey-box capacitance",
    )
    conductance, _ = _transformer_symmetric_matrix(
        conductance_s,
        length(nodes),
        "grey-box conductance",
    )
    source_signature = lowercase(String(source_response_sha256))
    occursin(r"^[0-9a-f]{64}$", source_signature) || throw(ArgumentError(
        "grey-box source response identity must be lowercase SHA-256",
    ))
    residual = Float64(identification_residual_relative)
    isfinite(residual) && residual >= 0.0 || throw(ArgumentError(
        "grey-box identification residual must be finite and nonnegative",
    ))
    nonuniqueness = strip(String(parameter_nonuniqueness))
    isempty(nonuniqueness) && throw(ArgumentError(
        "grey-box parameter nonuniqueness disclosure must not be empty",
    ))
    return GreyBoxTransformerModel(
        nodes,
        terminal_indices,
        branch_rows,
        capacitance,
        conductance,
        source_signature,
        residual,
        nonuniqueness,
    )
end

struct WhiteBoxTransformerModel <: AbstractTransformerTierModel
    winding_order::Vector{Symbol}
    section_count_per_winding::Vector{Int}
    section_length_m::Vector{Float64}
    series_resistance_ohm_per_m::Vector{Matrix{Float64}}
    series_inductance_h_per_m::Vector{Matrix{Float64}}
    shunt_conductance_s_per_m::Vector{Matrix{Float64}}
    shunt_capacitance_f_per_m::Vector{Matrix{Float64}}
    geometry_sha256::String
    frequency_band_hz::NTuple{2,Float64}
    section_refinement_residual_relative::Float64
end

function WhiteBoxTransformerModel(;
    winding_order,
    section_count_per_winding,
    section_length_m,
    series_resistance_ohm_per_m,
    series_inductance_h_per_m,
    shunt_conductance_s_per_m,
    shunt_capacitance_f_per_m,
    geometry_sha256::AbstractString,
    frequency_band_hz,
    section_refinement_residual_relative::Real,
)
    windings = _connection_symbols(winding_order, "white-box winding order")
    counts = Int.(section_count_per_winding)
    length(counts) == length(windings) &&
        all(count -> 2 <= count <= 32, counts) ||
        throw(ArgumentError(
            "white-box transformer requires 2 through 32 sections per winding",
        ))
    total_sections = sum(counts)
    total_sections <= 128 || throw(ArgumentError(
        "white-box transformer exceeds 128 coupled conductor sections",
    ))
    lengths = Float64.(section_length_m)
    length(lengths) == maximum(counts) && all(value -> isfinite(value) && value > 0.0, lengths) ||
        throw(ArgumentError(
            "white-box section lengths must cover the maximum section count and be positive",
        ))
    matrix_groups = (
        series_resistance_ohm_per_m,
        series_inductance_h_per_m,
        shunt_conductance_s_per_m,
        shunt_capacitance_f_per_m,
    )
    all(length(group) == length(lengths) for group in matrix_groups) ||
        throw(DimensionMismatch(
            "white-box per-section matrix groups must match section lengths",
        ))
    conductor_count = length(windings)
    checked_groups = Vector{Matrix{Float64}}[]
    for (group, label) in zip(
        matrix_groups,
        ("series resistance", "series inductance", "shunt conductance", "shunt capacitance"),
    )
        checked = Matrix{Float64}[]
        for values in group
            matrix, _ = _transformer_symmetric_matrix(
                values,
                conductor_count,
                "white-box $label",
            )
            push!(checked, matrix)
        end
        push!(checked_groups, checked)
    end
    geometry_signature = lowercase(String(geometry_sha256))
    occursin(r"^[0-9a-f]{64}$", geometry_signature) || throw(ArgumentError(
        "white-box geometry identity must be lowercase SHA-256",
    ))
    band = Float64.(Tuple(frequency_band_hz))
    length(band) == 2 && 0.0 < band[1] < band[2] <= 2.0e6 || throw(ArgumentError(
        "white-box frequency band must increase inside (0,2 MHz]",
    ))
    residual = Float64(section_refinement_residual_relative)
    isfinite(residual) && residual >= 0.0 || throw(ArgumentError(
        "white-box section refinement residual must be finite and nonnegative",
    ))
    return WhiteBoxTransformerModel(
        windings,
        counts,
        lengths,
        checked_groups[1],
        checked_groups[2],
        checked_groups[3],
        checked_groups[4],
        geometry_signature,
        (band[1], band[2]),
        residual,
    )
end

Base.@kwdef struct TransformerRuntimeSettings
    timestep_s::Float64
    initialization_frequency_hz::Float64
    kcl_absolute_tolerance_a::Float64=1.0e-9
    magnetic_continuity_absolute_tolerance_wb::Float64=1.0e-12
    energy_absolute_tolerance_j::Float64=1.0e-10
    nonlinear_residual_relative_tolerance::Float64=1.0e-9
    maximum_local_nonlinear_iterations::Int=40
end

function _validated_transformer_runtime_settings(settings::TransformerRuntimeSettings)
    for (value, label) in (
        (settings.timestep_s, "timestep"),
        (settings.initialization_frequency_hz, "initialization frequency"),
        (settings.kcl_absolute_tolerance_a, "KCL tolerance"),
        (settings.magnetic_continuity_absolute_tolerance_wb, "magnetic continuity tolerance"),
        (settings.energy_absolute_tolerance_j, "energy tolerance"),
        (settings.nonlinear_residual_relative_tolerance, "nonlinear residual tolerance"),
    )
        isfinite(value) && value > 0.0 || throw(ArgumentError(
            "transformer runtime $label must be finite and positive",
        ))
    end
    2.0e-9 <= settings.timestep_s <= 100.0e-6 || throw(ArgumentError(
        "transformer runtime timestep is outside the admitted 2 ns through 100 microsecond domain",
    ))
    16.7 <= settings.initialization_frequency_hz <= 400.0 || throw(ArgumentError(
        "transformer initialization frequency is outside the admitted 16.7 through 400 Hz domain",
    ))
    settings.maximum_local_nonlinear_iterations > 0 || throw(ArgumentError(
        "transformer maximum local nonlinear iterations must be positive",
    ))
    return settings
end

struct TransformerApparatusSpecification{M<:AbstractTransformerTierModel}
    id::Symbol
    tier::TransformerApparatusTier
    phase_count::Int
    rated_power_va::Float64
    rated_voltage_v::Float64
    rated_frequency_hz::Float64
    reactor_definition::Union{Nothing,ReactorApparatusDefinition}
    connection::TransformerConnectionTopology
    model::M
    settings::TransformerRuntimeSettings
    sources::Vector{TransformerSourceRecord}
    uncertainty::String
    validity_domain::String
    deterministic_signature_sha256::String
end

function _expected_transformer_model_type(tier::TransformerApparatusTier)
    tier === LowFrequencyTerminalTier && return LowFrequencyTransformerModel
    tier === BCTRANTerminalTier && return BCTRANTransformerModel
    tier === HybridTransformerTier && return HybridTransformerModel
    tier === MagneticEquivalentCircuitTier && return MagneticEquivalentCircuitModel
    tier === WidebandBlackBoxTier && return WidebandTransformerModel
    tier === GreyBoxLadderTier && return GreyBoxTransformerModel
    tier === WhiteBoxWindingTier && return WhiteBoxTransformerModel
    error("unreachable transformer tier")
end

function _transformer_matrix_has_mutual_coupling(matrix)
    matrix isa AbstractMatrix || return false
    rows, columns = size(matrix)
    rows == columns || return false
    return any(
        row != column && matrix[row, column] != 0.0
        for row in 1:rows, column in 1:columns
    )
end

function _transformer_model_has_mutual_coupling(model::AbstractTransformerTierModel)
    if model isa LowFrequencyTransformerModel || model isa BCTRANTransformerModel
        return _transformer_matrix_has_mutual_coupling(model.matrices.inductance_h) ||
            _transformer_matrix_has_mutual_coupling(model.matrices.capacitance_f)
    elseif model isa HybridTransformerModel
        return _transformer_matrix_has_mutual_coupling(model.leakage.inductance_h) ||
            _transformer_matrix_has_mutual_coupling(
                transpose(model.magnetic_graph.winding_turns) *
                model.magnetic_graph.winding_turns,
            )
    elseif model isa MagneticEquivalentCircuitModel
        return _transformer_matrix_has_mutual_coupling(model.leakage_inductance_h) ||
            _transformer_matrix_has_mutual_coupling(
                transpose(model.magnetic_graph.winding_turns) *
                model.magnetic_graph.winding_turns,
            )
    elseif model isa WidebandTransformerModel
        return _transformer_matrix_has_mutual_coupling(model.direct_admittance_s) ||
            _transformer_matrix_has_mutual_coupling(
                model.output_matrix_s_per_s * model.input_matrix,
            )
    elseif model isa GreyBoxTransformerModel
        return any(
            branch -> _transformer_matrix_has_mutual_coupling(branch.inductance_h),
            model.branches,
        )
    else
        model = model::WhiteBoxTransformerModel
        return any(_transformer_matrix_has_mutual_coupling, model.series_inductance_h_per_m) ||
            any(_transformer_matrix_has_mutual_coupling, model.shunt_capacitance_f_per_m)
    end
end

function _transformer_model_signature(model::AbstractTransformerTierModel)
    if model isa LowFrequencyTransformerModel || model isa BCTRANTransformerModel
        return model.matrices.deterministic_signature_sha256
    elseif model isa HybridTransformerModel
        return bytes2hex(sha256(
            model.leakage.deterministic_signature_sha256 *
            model.magnetic_graph.deterministic_signature_sha256 *
            _transformer_matrix_signature(
                model.winding_loss_state_matrix_per_s,
                model.winding_loss_input_matrix,
                model.winding_loss_output_matrix_ohm_per_s,
                model.winding_loss_direct_ohm,
                model.winding_loss_storage_matrix_j,
                model.winding_loss_dissipation_matrix_j_per_s,
            ),
        ))
    elseif model isa MagneticEquivalentCircuitModel
        return bytes2hex(sha256(
            _transformer_matrix_signature(
                model.winding_resistance_ohm,
                model.leakage_inductance_h,
                model.terminal_capacitance_f,
                model.terminal_conductance_s,
            ) * model.magnetic_graph.deterministic_signature_sha256,
        ))
    elseif model isa WidebandTransformerModel
        return model.deterministic_signature_sha256
    elseif model isa GreyBoxTransformerModel
        io = IOBuffer()
        println(io, join(String.(model.node_order), ','))
        println(io, join(model.terminal_node_indices, ','))
        println(io, model.source_response_sha256)
        println(io, bitstring(model.identification_residual_relative))
        println(io, _transformer_matrix_signature(model.capacitance_f, model.conductance_s))
        for branch in model.branches
            println(io, branch.id, ',', branch.from_node_index, ',', branch.to_node_index)
            println(io, bitstring(branch.resistance_ohm), ',', bitstring(branch.inductance_h))
        end
        return bytes2hex(sha256(take!(io)))
    else
        model = model::WhiteBoxTransformerModel
        io = IOBuffer()
        println(io, join(String.(model.winding_order), ','))
        println(io, join(model.section_count_per_winding, ','))
        println(io, model.geometry_sha256)
        for matrices in (
            model.series_resistance_ohm_per_m,
            model.series_inductance_h_per_m,
            model.shunt_conductance_s_per_m,
            model.shunt_capacitance_f_per_m,
        )
            println(io, _transformer_matrix_signature(matrices...))
        end
        return bytes2hex(sha256(take!(io)))
    end
end

function TransformerApparatusSpecification(
    id::Symbol,
    tier::TransformerApparatusTier,
    connection::TransformerConnectionTopology,
    model::M,
    settings::TransformerRuntimeSettings;
    phase_count::Integer,
    rated_power_va::Real,
    rated_voltage_v::Real,
    rated_frequency_hz::Real,
    reactor_definition::Union{Nothing,ReactorApparatusDefinition}=nothing,
    sources,
    uncertainty::AbstractString,
    validity_domain::AbstractString,
) where {M<:AbstractTransformerTierModel}
    id == Symbol("") && throw(ArgumentError(
        "transformer apparatus identity must not be empty",
    ))
    expected_model = _expected_transformer_model_type(tier)
    model isa expected_model || throw(ArgumentError(
        "transformer tier $(_TRANSFORMER_TIER_IDS[tier]) requires model $expected_model",
    ))
    phases = Int(phase_count)
    phases in (1, 3) && phases == length(connection.phase_order) ||
        throw(ArgumentError(
            "transformer phase count must be one or three and match the connection",
        ))
    rated_power = Float64(rated_power_va)
    rated_voltage = Float64(rated_voltage_v)
    rated_frequency = Float64(rated_frequency_hz)
    100.0 <= rated_power <= 1.5e9 || throw(ArgumentError(
        "transformer rated power is outside the admitted 0.1 kVA through 1.5 GVA domain",
    ))
    10.0 <= rated_voltage <= 765.0e3 || throw(ArgumentError(
        "transformer rated voltage is outside the admitted 10 V through 765 kV domain",
    ))
    16.7 <= rated_frequency <= 400.0 || throw(ArgumentError(
        "transformer rated frequency is outside the admitted 16.7 through 400 Hz domain",
    ))
    checked_settings = _validated_transformer_runtime_settings(settings)
    if reactor_definition !== nothing
        reactor = something(reactor_definition)
        coil_count = length(connection.coil_order)
        if reactor.winding_configuration === SingleCoilReactorWinding
            coil_count == 1 || throw(ArgumentError(
                "single-coil reactor topology requires exactly one declared coil",
            ))
        else
            coil_count >= 2 || throw(ArgumentError(
                "split or mutually coupled reactor topology requires at least two coils",
            ))
        end
        reactor.winding_configuration === MutuallyCoupledReactorWinding &&
            !_transformer_model_has_mutual_coupling(model) && throw(ArgumentError(
                "mutually coupled reactor topology requires explicit cross-coil coupling",
            ))
        if reactor.magnetic_construction === AirCoreReactorConstruction
            tier in (
                HybridTransformerTier,
                MagneticEquivalentCircuitTier,
            ) && throw(ArgumentError(
                "air-core reactor models cannot invent a ferromagnetic core state",
            ))
        else
            tier in (
                HybridTransformerTier,
                MagneticEquivalentCircuitTier,
            ) || throw(ArgumentError(
                "iron-core reactor models require an explicit hybrid or magnetic-equivalent core graph",
            ))
            graph = _transformer_model_magnetic_graph(model)
            graph === nothing && error("unreachable iron-core reactor model")
            graph_gap_length = sum(branch.air_gap_length_m for branch in graph.branches)
            gap_tolerance = 64.0 * eps(Float64) * max(
                graph_gap_length,
                reactor.total_air_gap_length_m,
                1.0,
            )
            abs(graph_gap_length - reactor.total_air_gap_length_m) <= gap_tolerance ||
                throw(ArgumentError(
                    "reactor declared total air-gap length does not match its magnetic graph",
                ))
            for branch in graph.branches
                branch.air_gap_length_m == 0.0 && continue
                abs(
                    branch.air_gap_effective_area_factor -
                    reactor.air_gap_effective_area_factor,
                ) <= 64.0 * eps(Float64) * max(
                    branch.air_gap_effective_area_factor,
                    reactor.air_gap_effective_area_factor,
                    1.0,
                ) || throw(ArgumentError(
                    "reactor gap-fringing factor does not match its magnetic graph",
                ))
            end
        end
    end
    source_rows = TransformerSourceRecord[sources...]
    isempty(source_rows) && throw(ArgumentError(
        "transformer apparatus requires at least one source record",
    ))
    source_ids = getfield.(source_rows, :id)
    length(unique(source_ids)) == length(source_ids) || throw(ArgumentError(
        "transformer source identities must be unique",
    ))
    uncertainty_text = strip(String(uncertainty))
    validity_text = strip(String(validity_domain))
    isempty(uncertainty_text) && throw(ArgumentError(
        "transformer uncertainty must be explicit, including unknown",
    ))
    isempty(validity_text) && throw(ArgumentError(
        "transformer validity domain must not be empty",
    ))
    model_signature = _transformer_model_signature(model)
    io = IOBuffer()
    println(io, id)
    println(io, _TRANSFORMER_TIER_IDS[tier])
    println(io, connection.deterministic_signature_sha256)
    println(io, model_signature)
    if reactor_definition === nothing
        println(io, "transformer")
    else
        reactor = something(reactor_definition)
        println(io, "reactor")
        println(io, Int(reactor.application))
        println(io, Int(reactor.magnetic_construction))
        println(io, Int(reactor.winding_configuration))
        println(io, Int(reactor.control_mode))
        println(io, Int(reactor.gap_model))
        println(io, bitstring(reactor.total_air_gap_length_m))
        println(io, bitstring(reactor.air_gap_effective_area_factor))
    end
    for value in (
        rated_power,
        rated_voltage,
        rated_frequency,
        checked_settings.timestep_s,
        checked_settings.initialization_frequency_hz,
    )
        println(io, bitstring(value))
    end
    for source in source_rows
        println(io, source.id)
        println(io, source.content_sha256)
    end
    println(io, uncertainty_text)
    println(io, validity_text)
    signature = bytes2hex(sha256(take!(io)))
    return TransformerApparatusSpecification{M}(
        id,
        tier,
        phases,
        rated_power,
        rated_voltage,
        rated_frequency,
        reactor_definition,
        connection,
        model,
        checked_settings,
        source_rows,
        uncertainty_text,
        validity_text,
        signature,
    )
end

transformer_apparatus_signature(specification::TransformerApparatusSpecification) =
    specification.deterministic_signature_sha256

struct TransformerApparatusPreparation{S<:TransformerApparatusSpecification}
    specification::S
    terminal_order::Vector{Symbol}
    coil_order::Vector{Symbol}
    effective_linear_inductance_h::Union{Nothing,Matrix{Float64}}
    initialization_mode::TransformerInitializationMode
    initial_time_s::Float64
    initial_node_voltage_v::Vector{Float64}
    initial_node_voltage_derivative_v_per_s::Vector{Float64}
    initial_coil_current_a::Vector{Float64}
    initial_coil_current_derivative_a_per_s::Vector{Float64}
    initial_branch_flux_wb::Union{Nothing,Vector{Float64}}
    initial_branch_flux_derivative_wb_per_s::Union{Nothing,Vector{Float64}}
    residual_flux_projection_correction_wb::Float64
    operating_point_electrical_residual_v::Float64
    preparation_signature_sha256::String
end

function _transformer_model_coil_count(model::AbstractTransformerTierModel)
    model isa LowFrequencyTransformerModel && return size(model.matrices.resistance_ohm, 1)
    model isa BCTRANTransformerModel && return size(model.matrices.resistance_ohm, 1)
    model isa HybridTransformerModel && return size(model.leakage.resistance_ohm, 1)
    model isa MagneticEquivalentCircuitModel && return size(model.winding_resistance_ohm, 1)
    return nothing
end

function _transformer_model_magnetic_graph(model::AbstractTransformerTierModel)
    model isa HybridTransformerModel && return model.magnetic_graph
    model isa MagneticEquivalentCircuitModel && return model.magnetic_graph
    return nothing
end

function _transformer_sinusoidal_frequencies(settings::TransformerRuntimeSettings)
    physical_angular_frequency = 2.0 * pi * settings.initialization_frequency_hz
    half_step_angle = 0.5 * physical_angular_frequency * settings.timestep_s
    abs(half_step_angle) < 0.5 * pi || throw(DomainError(
        half_step_angle,
        "transformer sinusoidal initialization reaches the trapezoidal Nyquist boundary",
    ))
    discrete_angular_frequency =
        (2.0 / settings.timestep_s) * tan(half_step_angle)
    return physical_angular_frequency, discrete_angular_frequency
end

function _transformer_sinusoidal_phasor(
    value,
    physical_derivative,
    physical_angular_frequency::Float64,
)
    return ComplexF64.(value) .-
        im .* ComplexF64.(physical_derivative) ./ physical_angular_frequency
end

function _transformer_initial_discrete_derivative(
    value,
    physical_derivative,
    initialization_mode::TransformerInitializationMode,
    settings::TransformerRuntimeSettings,
)
    initialization_mode === SinusoidalTransformerOperatingPoint ||
        return Float64.(physical_derivative)
    physical_frequency, discrete_frequency =
        _transformer_sinusoidal_frequencies(settings)
    phasor = _transformer_sinusoidal_phasor(
        value,
        physical_derivative,
        physical_frequency,
    )
    return real.(im * discrete_frequency .* phasor)
end

function _transformer_sinusoidal_linear_state(
    state_matrix_per_s,
    input_matrix,
    input_value,
    input_derivative,
    settings::TransformerRuntimeSettings,
)
    state_count = size(state_matrix_per_s, 1)
    state_count == 0 && return (
        state=Float64[],
        derivative=Float64[],
        phasor=ComplexF64[],
        residual=0.0,
    )
    physical_frequency, discrete_frequency =
        _transformer_sinusoidal_frequencies(settings)
    input_phasor = _transformer_sinusoidal_phasor(
        input_value,
        input_derivative,
        physical_frequency,
    )
    state_phasor = (
        im * discrete_frequency .* Matrix{ComplexF64}(I, state_count, state_count) .-
        ComplexF64.(state_matrix_per_s)
    ) \ (ComplexF64.(input_matrix) * input_phasor)
    state = real.(state_phasor)
    derivative = real.(im * discrete_frequency .* state_phasor)
    residual = maximum(
        abs,
        derivative .- state_matrix_per_s * state .- input_matrix * input_value;
        init=0.0,
    )
    return (
        state=state,
        derivative=derivative,
        phasor=state_phasor,
        residual=residual,
    )
end

function _transformer_initial_winding_loss_state(
    model::HybridTransformerModel,
    initialization_mode::TransformerInitializationMode,
    initial_coil_current_a,
    initial_coil_current_derivative_a_per_s,
    settings::TransformerRuntimeSettings,
)
    if initialization_mode === SinusoidalTransformerOperatingPoint
        state = _transformer_sinusoidal_linear_state(
            model.winding_loss_state_matrix_per_s,
            model.winding_loss_input_matrix,
            initial_coil_current_a,
            initial_coil_current_derivative_a_per_s,
            settings,
        )
    else
        state = (
            state=zeros(size(model.winding_loss_state_matrix_per_s, 1)),
            derivative=zeros(size(model.winding_loss_state_matrix_per_s, 1)),
            phasor=zeros(ComplexF64, size(model.winding_loss_state_matrix_per_s, 1)),
            residual=0.0,
        )
    end
    voltage = model.winding_loss_output_matrix_ohm_per_s * state.state .+
        model.winding_loss_direct_ohm * initial_coil_current_a
    return merge(state, (voltage=voltage,))
end

function _project_transformer_magnetic_continuity(
    graph::TransformerMagneticGraph,
    requested_flux_wb,
    tolerance_wb::Float64,
)
    requested = Float64.(requested_flux_wb)
    length(requested) == length(graph.branches) && all(isfinite, requested) ||
        throw(DimensionMismatch(
            "transformer residual-flux vector must be finite and cover every magnetic branch",
        ))
    correction = transpose(graph.incidence) * (
        (graph.incidence * transpose(graph.incidence)) \
        (graph.incidence * requested)
    )
    projected = requested - correction
    maximum_correction = maximum(abs, correction; init=0.0)
    maximum_correction <= tolerance_wb || throw(DomainError(
        maximum_correction,
        "transformer residual-flux request is outside its declared magnetic-continuity projection tolerance",
    ))
    for branch_index in eachindex(graph.branches)
        branch = graph.branches[branch_index]
        material = graph.materials[branch.material_index]
        flux_density = projected[branch_index] / branch.cross_section_m2
        if material isa LinearTransformerMagneticMaterial
            abs(flux_density) <= material.maximum_flux_density_t || throw(DomainError(
                flux_density,
                "transformer residual flux exceeds a linear magnetic-material domain",
            ))
        elseif material isa PiecewiseLinearTransformerMagneticMaterial
            abs(flux_density) <= last(material.flux_density_t) || throw(DomainError(
                flux_density,
                "transformer residual flux exceeds a saturation material domain",
            ))
        else
            material = material::TellinenTransformerMagneticMaterial
            lower, _ = _tellinen_curve_value(material.lower_branch, 0.0)
            upper, _ = _tellinen_curve_value(material.upper_branch, 0.0)
            lower <= flux_density <= upper || throw(DomainError(
                flux_density,
                "transformer residual flux is not remanently admissible at zero field",
            ))
        end
    end
    return projected, maximum_correction
end

function prepare_transformer_apparatus(
    specification::TransformerApparatusSpecification;
    initialization_mode::TransformerInitializationMode=SpecifiedTransformerInitialization,
    initial_time_s::Real=0.0,
    initial_node_voltage_v=zeros(length(specification.connection.node_order)),
    initial_node_voltage_derivative_v_per_s=zeros(
        length(specification.connection.node_order),
    ),
    initial_coil_current_a=zeros(length(specification.connection.coil_order)),
    initial_coil_current_derivative_a_per_s=zeros(
        length(specification.connection.coil_order),
    ),
    initial_branch_flux_wb=nothing,
    initial_branch_flux_derivative_wb_per_s=nothing,
    residual_flux_projection_tolerance_wb::Real=
        specification.settings.magnetic_continuity_absolute_tolerance_wb,
    flux_derivative_projection_tolerance_wb_per_s::Real=
        specification.settings.magnetic_continuity_absolute_tolerance_wb /
        specification.settings.timestep_s,
)
    connection = specification.connection
    model = specification.model
    coil_count = _transformer_model_coil_count(model)
    coil_count === nothing || coil_count == length(connection.coil_order) ||
        _transformer_refusal(
            :coil_count_mismatch,
            :prepare,
            specification.id,
            specification.tier,
            "model coil count does not match the explicit connection";
            diagnostics=(model_coil_count=coil_count, connection_coil_count=length(connection.coil_order)),
        )
    if model isa WidebandTransformerModel
        model.port_order == connection.node_order || _transformer_refusal(
            :terminal_order_mismatch,
            :prepare,
            specification.id,
            specification.tier,
            "wideband port order must equal the explicit physical terminal order",
        )
        last(model.frequency_band_hz) * specification.settings.timestep_s < 0.5 ||
            _transformer_refusal(
                :frequency_timestep_domain,
                :prepare,
                specification.id,
                specification.tier,
                "wideband source band reaches or exceeds fixed-step Nyquist",
            )
        if initialization_mode === SinusoidalTransformerOperatingPoint
            first(model.frequency_band_hz) <=
                specification.settings.initialization_frequency_hz <=
                last(model.frequency_band_hz) || _transformer_refusal(
                    :initialization_frequency_outside_source_band,
                    :initialize,
                    specification.id,
                    specification.tier,
                    "wideband sinusoidal initialization frequency is outside the immutable source band",
                )
        end
    elseif model isa GreyBoxTransformerModel
        length(model.terminal_node_indices) == length(connection.node_order) ||
            _transformer_refusal(
                :terminal_count_mismatch,
                :prepare,
                specification.id,
                specification.tier,
                "grey-box represented terminal count does not match the connection",
            )
    elseif model isa WhiteBoxTransformerModel
        model.winding_order == connection.winding_order || _transformer_refusal(
            :winding_order_mismatch,
            :prepare,
            specification.id,
            specification.tier,
            "white-box winding order does not match the connection",
        )
        last(model.frequency_band_hz) * specification.settings.timestep_s < 0.5 ||
            _transformer_refusal(
                :frequency_timestep_domain,
                :prepare,
                specification.id,
                specification.tier,
                "white-box source band reaches or exceeds fixed-step Nyquist",
            )
        if initialization_mode === SinusoidalTransformerOperatingPoint
            first(model.frequency_band_hz) <=
                specification.settings.initialization_frequency_hz <=
                last(model.frequency_band_hz) || _transformer_refusal(
                    :initialization_frequency_outside_source_band,
                    :initialize,
                    specification.id,
                    specification.tier,
                    "white-box sinusoidal initialization frequency is outside the immutable source band",
                )
        end
    end
    initial_time = Float64(initial_time_s)
    isfinite(initial_time) && initial_time >= 0.0 || _transformer_refusal(
        :initial_time_domain,
        :initialize,
        specification.id,
        specification.tier,
        "transformer initial time must be finite and nonnegative",
    )
    initial_voltage = Float64.(initial_node_voltage_v)
    initial_voltage_derivative = Float64.(initial_node_voltage_derivative_v_per_s)
    initial_current = Float64.(initial_coil_current_a)
    initial_current_derivative = Float64.(initial_coil_current_derivative_a_per_s)
    length(initial_voltage) == length(connection.node_order) && all(isfinite, initial_voltage) ||
        _transformer_refusal(
            :initial_voltage_mismatch,
            :initialize,
            specification.id,
            specification.tier,
            "initial terminal voltages must be finite and match the connection",
        )
    length(initial_current) == length(connection.coil_order) && all(isfinite, initial_current) ||
        _transformer_refusal(
            :initial_current_mismatch,
            :initialize,
            specification.id,
            specification.tier,
            "initial winding currents must be finite and match the connection",
        )
    length(initial_voltage_derivative) == length(connection.node_order) &&
        all(isfinite, initial_voltage_derivative) || _transformer_refusal(
            :initial_voltage_derivative_mismatch,
            :initialize,
            specification.id,
            specification.tier,
            "initial terminal-voltage derivatives must be finite and match the connection",
        )
    length(initial_current_derivative) == length(connection.coil_order) &&
        all(isfinite, initial_current_derivative) || _transformer_refusal(
            :initial_current_derivative_mismatch,
            :initialize,
            specification.id,
            specification.tier,
            "initial winding-current derivatives must be finite and match the connection",
        )
    if initialization_mode === DeenergizedTransformerInitialization
        maximum(abs, initial_voltage; init=0.0) == 0.0 &&
            maximum(abs, initial_voltage_derivative; init=0.0) == 0.0 &&
            maximum(abs, initial_current; init=0.0) == 0.0 &&
            maximum(abs, initial_current_derivative; init=0.0) == 0.0 ||
            _transformer_refusal(
                :energized_deenergized_initialization,
                :initialize,
                specification.id,
                specification.tier,
                "deenergized transformer initialization cannot contain electrical excitation",
            )
    end
    magnetic_graph = _transformer_model_magnetic_graph(model)
    projection_tolerance = Float64(residual_flux_projection_tolerance_wb)
    derivative_projection_tolerance =
        Float64(flux_derivative_projection_tolerance_wb_per_s)
    isfinite(projection_tolerance) && projection_tolerance >= 0.0 ||
        _transformer_refusal(
            :residual_flux_projection_tolerance,
            :initialize,
            specification.id,
            specification.tier,
            "residual-flux projection tolerance must be finite and nonnegative",
        )
    isfinite(derivative_projection_tolerance) &&
        derivative_projection_tolerance >= 0.0 || _transformer_refusal(
            :flux_derivative_projection_tolerance,
            :initialize,
            specification.id,
            specification.tier,
            "flux-derivative projection tolerance must be finite and nonnegative",
        )
    if magnetic_graph === nothing && (
        initial_branch_flux_wb !== nothing ||
        initial_branch_flux_derivative_wb_per_s !== nothing
    )
        _transformer_refusal(
            :unrepresented_initial_magnetic_state,
            :initialize,
            specification.id,
            specification.tier,
            "selected transformer tier does not represent magnetic branch state",
        )
    end
    projected_branch_flux = nothing
    projection_correction = 0.0
    if initial_branch_flux_wb !== nothing
        try
            projected_branch_flux, projection_correction =
                _project_transformer_magnetic_continuity(
                    something(magnetic_graph),
                    initial_branch_flux_wb,
                    projection_tolerance,
                )
        catch error
            _transformer_refusal(
                :inadmissible_residual_flux,
                :initialize,
                specification.id,
                specification.tier,
                sprint(showerror, error),
            )
        end
    end
    projected_branch_flux_derivative = nothing
    if initial_branch_flux_derivative_wb_per_s !== nothing
        graph = something(magnetic_graph)
        requested_derivative = Float64.(initial_branch_flux_derivative_wb_per_s)
        length(requested_derivative) == length(graph.branches) &&
            all(isfinite, requested_derivative) || _transformer_refusal(
                :initial_flux_derivative_mismatch,
                :initialize,
                specification.id,
                specification.tier,
                "initial branch-flux derivatives must be finite and cover every magnetic branch",
            )
        derivative_correction = transpose(graph.incidence) * (
            (graph.incidence * transpose(graph.incidence)) \
            (graph.incidence * requested_derivative)
        )
        maximum(abs, derivative_correction; init=0.0) <=
            derivative_projection_tolerance || _transformer_refusal(
                :inadmissible_flux_derivative,
                :initialize,
                specification.id,
                specification.tier,
                "initial branch-flux derivative violates magnetic continuity beyond tolerance",
            )
        projected_branch_flux_derivative = requested_derivative - derivative_correction
    end
    if initialization_mode === DeenergizedTransformerInitialization &&
       projected_branch_flux_derivative !== nothing &&
       maximum(abs, projected_branch_flux_derivative; init=0.0) > 0.0
        _transformer_refusal(
            :changing_flux_in_deenergized_initialization,
            :initialize,
            specification.id,
            specification.tier,
            "deenergized transformer initialization requires static residual flux",
        )
    end
    coil_voltage = transpose(connection.incidence) * initial_voltage
    coil_voltage_derivative =
        transpose(connection.incidence) * initial_voltage_derivative
    electrical_residual = zeros(length(connection.coil_order))
    complex_electrical_residual = ComplexF64[]
    if model isa LowFrequencyTransformerModel || model isa BCTRANTransformerModel
        if initialization_mode === SinusoidalTransformerOperatingPoint
            physical_frequency, discrete_frequency =
                _transformer_sinusoidal_frequencies(specification.settings)
            voltage_phasor = _transformer_sinusoidal_phasor(
                coil_voltage,
                coil_voltage_derivative,
                physical_frequency,
            )
            current_phasor = _transformer_sinusoidal_phasor(
                initial_current,
                initial_current_derivative,
                physical_frequency,
            )
            complex_electrical_residual = voltage_phasor .-
                (
                    model.matrices.resistance_ohm .+
                    im * discrete_frequency .* model.matrices.inductance_h
                ) * current_phasor
        else
            electrical_residual .= coil_voltage .-
                model.matrices.resistance_ohm * initial_current .-
                model.matrices.inductance_h * initial_current_derivative
        end
    elseif model isa HybridTransformerModel || model isa MagneticEquivalentCircuitModel
        parts_resistance = model isa HybridTransformerModel ?
            model.leakage.resistance_ohm : model.winding_resistance_ohm
        parts_leakage = model isa HybridTransformerModel ?
            model.leakage.inductance_h : model.leakage_inductance_h
        if initialization_mode === SinusoidalTransformerOperatingPoint &&
           projected_branch_flux_derivative === nothing
            _transformer_refusal(
                :missing_operating_point_flux_derivative,
                :initialize,
                specification.id,
                specification.tier,
                "sinusoidal magnetic initialization requires complete branch-flux derivatives",
            )
        end
        if initialization_mode === SinusoidalTransformerOperatingPoint &&
           projected_branch_flux === nothing
            _transformer_refusal(
                :missing_operating_point_flux,
                :initialize,
                specification.id,
                specification.tier,
                "sinusoidal magnetic initialization requires complete branch fluxes",
            )
        end
        flux_derivative = projected_branch_flux_derivative === nothing ?
            zeros(length(something(magnetic_graph).branches)) :
            projected_branch_flux_derivative
        if initialization_mode === SinusoidalTransformerOperatingPoint
            physical_frequency, discrete_frequency =
                _transformer_sinusoidal_frequencies(specification.settings)
            voltage_phasor = _transformer_sinusoidal_phasor(
                coil_voltage,
                coil_voltage_derivative,
                physical_frequency,
            )
            current_phasor = _transformer_sinusoidal_phasor(
                initial_current,
                initial_current_derivative,
                physical_frequency,
            )
            flux_phasor = _transformer_sinusoidal_phasor(
                something(projected_branch_flux),
                flux_derivative,
                physical_frequency,
            )
            winding_loss_voltage_phasor = zeros(
                ComplexF64,
                length(connection.coil_order),
            )
            if model isa HybridTransformerModel
                winding_loss = _transformer_initial_winding_loss_state(
                    model,
                    initialization_mode,
                    initial_current,
                    initial_current_derivative,
                    specification.settings,
                )
                winding_loss_voltage_phasor .=
                    model.winding_loss_output_matrix_ohm_per_s *
                        winding_loss.phasor .+
                    model.winding_loss_direct_ohm * current_phasor
            end
            complex_electrical_residual = voltage_phasor .-
                (
                    parts_resistance .+
                    im * discrete_frequency .* parts_leakage
                ) * current_phasor .-
                im * discrete_frequency .* (
                    transpose(something(magnetic_graph).winding_turns) * flux_phasor
                ) .- winding_loss_voltage_phasor
        else
            winding_loss_voltage = model isa HybridTransformerModel ?
                _transformer_initial_winding_loss_state(
                    model,
                    initialization_mode,
                    initial_current,
                    initial_current_derivative,
                    specification.settings,
                ).voltage : zeros(length(connection.coil_order))
            electrical_residual .= coil_voltage .-
                parts_resistance * initial_current .-
                parts_leakage * initial_current_derivative .-
                transpose(something(magnetic_graph).winding_turns) * flux_derivative .-
                winding_loss_voltage
        end
    end
    maximum_electrical_residual = isempty(complex_electrical_residual) ?
        maximum(abs, electrical_residual; init=0.0) :
        maximum(abs, complex_electrical_residual; init=0.0)
    if initialization_mode === SinusoidalTransformerOperatingPoint
        residual_scale = max(maximum(abs, coil_voltage; init=0.0), 1.0)
        residual_tolerance = max(
            1.0e-7,
            specification.settings.nonlinear_residual_relative_tolerance * residual_scale,
        )
        maximum_electrical_residual <= residual_tolerance || _transformer_refusal(
            :inconsistent_operating_point,
            :initialize,
            specification.id,
            specification.tier,
            "transformer sinusoidal operating point violates its electrical residual";
            diagnostics=(
                maximum_electrical_residual_v=maximum_electrical_residual,
                tolerance_v=residual_tolerance,
            ),
        )
    end
    effective_inductance = if model isa LowFrequencyTransformerModel ||
                              model isa BCTRANTransformerModel
        copy(model.matrices.inductance_h)
    elseif model isa HybridTransformerModel
        all(material -> material isa LinearTransformerMagneticMaterial, model.magnetic_graph.materials) ?
        model.leakage.inductance_h +
            transformer_magnetic_linear_inductance(model.magnetic_graph) : nothing
    elseif model isa MagneticEquivalentCircuitModel
        all(material -> material isa LinearTransformerMagneticMaterial, model.magnetic_graph.materials) ?
        model.leakage_inductance_h +
            transformer_magnetic_linear_inductance(model.magnetic_graph) : nothing
    else
        nothing
    end
    io = IOBuffer()
    println(io, specification.deterministic_signature_sha256)
    println(io, Int(initialization_mode))
    println(io, bitstring(initial_time))
    for values in (
        initial_voltage,
        initial_voltage_derivative,
        initial_current,
        initial_current_derivative,
    )
        for value in values
            println(io, bitstring(value))
        end
    end
    for values in (projected_branch_flux, projected_branch_flux_derivative)
        if values === nothing
            println(io, "nothing")
        else
            println(io, length(values))
            for value in values
                println(io, bitstring(value))
            end
        end
    end
    println(io, bitstring(projection_correction))
    println(io, bitstring(maximum_electrical_residual))
    preparation_signature = bytes2hex(sha256(take!(io)))
    return TransformerApparatusPreparation(
        specification,
        copy(connection.node_order),
        copy(connection.coil_order),
        effective_inductance,
        initialization_mode,
        initial_time,
        initial_voltage,
        initial_voltage_derivative,
        initial_current,
        initial_current_derivative,
        projected_branch_flux,
        projected_branch_flux_derivative,
        projection_correction,
        maximum_electrical_residual,
        preparation_signature,
    )
end

struct TransformerApparatusReadiness
    apparatus::Symbol
    tier::TransformerApparatusTier
    ready::Bool
    production_backend_available::Bool
    unavailable_outputs::Vector{Symbol}
    unsupported_phenomena::Vector{Symbol}
    preparation_signature_sha256::String
end

function transformer_apparatus_readiness(
    preparation::TransformerApparatusPreparation;
    production_backend_available::Bool=false,
)
    tier = preparation.specification.tier
    unavailable = Symbol[]
    tier in (LowFrequencyTerminalTier, BCTRANTerminalTier, WidebandBlackBoxTier) &&
        append!(unavailable, (:core_branch_flux_wb, :winding_internal_voltage_v))
    tier === GreyBoxLadderTier && push!(unavailable, :unrepresented_turn_voltage_v)
    return TransformerApparatusReadiness(
        preparation.specification.id,
        tier,
        true,
        production_backend_available,
        unavailable,
        collect(_TRANSFORMER_TIER_UNSUPPORTED[tier]),
        preparation.preparation_signature_sha256,
    )
end
