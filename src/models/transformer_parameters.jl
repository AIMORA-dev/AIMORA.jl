module TransformerParameters

using LinearAlgebra

export AbstractTransformerParameterCase,
       MultiphaseTransformerParameterCase,
       MultiphaseTransformerParameterResult,
       MultiphaseTransformerWinding,
       TransformerShortCircuitCase,
       TransformerShortCircuitResult,
       TransformerMagnetizingShunt,
       SaturableTransformerParameterCase,
       SaturableTransformerParameterResult,
       TransformerSequenceShortCircuitTest,
       TransformerParameterWinding,
       TransformerGeneratedBranchRow,
       multiphase_transformer_parameters,
       transformer_short_circuit_parameters,
       saturable_transformer_parameters

abstract type AbstractTransformerParameterCase end

struct TransformerParameterWinding
    winding_number::Int
    from_node::Symbol
    to_node::Symbol
    resistance_ohm::Float64
    inductance_or_reactance::Float64
    rated_voltage_kv::Float64

    function TransformerParameterWinding(
        winding_number::Integer,
        from_node,
        to_node,
        resistance_ohm::Real,
        inductance_or_reactance::Real,
        rated_voltage_kv::Real,
    )
        number = Int(winding_number)
        resistance = Float64(resistance_ohm)
        inductance = Float64(inductance_or_reactance)
        voltage = Float64(rated_voltage_kv)
        number > 0 ||
            throw(ArgumentError("transformer winding number must be positive"))
        isfinite(resistance) && resistance >= 0.0 ||
            throw(ArgumentError("transformer winding resistance must be finite and nonnegative"))
        isfinite(inductance) && inductance >= 0.0 ||
            throw(ArgumentError("transformer winding inductance/reactance must be finite and nonnegative"))
        isfinite(voltage) && voltage > 0.0 ||
            throw(ArgumentError("transformer winding voltage must be finite and positive"))
        return new(
            number,
            Symbol(from_node),
            Symbol(to_node),
            resistance,
            inductance,
            voltage,
        )
    end
end

struct TransformerShortCircuitCase <: AbstractTransformerParameterCase
    source_line::Int
    winding_voltages_kv::Vector{Float64}
    pair_losses_kw::Vector{Float64}
    pair_impedance_percent::Vector{Float64}
    pair_ratings_mva::Vector{Float64}
    magnetizing_current_percent::Float64
    magnetizing_rating_mva::Float64
    node_pairs::Vector{NTuple{2,Symbol}}

    function TransformerShortCircuitCase(
        source_line::Integer,
        winding_voltages_kv::AbstractVector{<:Real},
        pair_losses_kw::AbstractVector{<:Real},
        pair_impedance_percent::AbstractVector{<:Real},
        pair_ratings_mva::AbstractVector{<:Real},
        magnetizing_current_percent::Real,
        magnetizing_rating_mva::Real,
        node_pairs::AbstractVector{<:Tuple},
    )
        voltages = Float64.(winding_voltages_kv)
        winding_count = length(voltages)
        winding_count in (2, 3) ||
            throw(ArgumentError("OVER41 short-circuit conversion requires two or three windings"))
        pair_count = winding_count * (winding_count - 1) ÷ 2
        losses = Float64.(pair_losses_kw)
        impedances = Float64.(pair_impedance_percent)
        ratings = Float64.(pair_ratings_mva)
        length(losses) == pair_count ||
            throw(ArgumentError("short-circuit loss count must equal $pair_count"))
        length(impedances) == pair_count ||
            throw(ArgumentError("short-circuit impedance count must equal $pair_count"))
        length(ratings) == pair_count ||
            throw(ArgumentError("short-circuit rating count must equal $pair_count"))
        all(value -> isfinite(value) && value > 0.0, voltages) ||
            throw(ArgumentError("winding voltages must be finite and positive"))
        all(value -> isfinite(value) && value >= 0.0, losses) ||
            throw(ArgumentError("short-circuit losses must be finite and nonnegative"))
        all(value -> isfinite(value) && value != 0.0, impedances) ||
            throw(ArgumentError("short-circuit impedances must be finite and nonzero"))
        all(value -> isfinite(value) && value > 0.0, ratings) ||
            throw(ArgumentError("short-circuit ratings must be finite and positive"))
        excitation = Float64(magnetizing_current_percent)
        excitation_rating = Float64(magnetizing_rating_mva)
        isfinite(excitation) && excitation >= 0.0 ||
            throw(ArgumentError("magnetizing current must be finite and nonnegative"))
        isfinite(excitation_rating) && excitation_rating > 0.0 ||
            throw(ArgumentError("magnetizing rating must be finite and positive"))
        normalized_pairs = NTuple{2,Symbol}[
            (Symbol(pair[1]), Symbol(pair[2])) for pair in node_pairs
        ]
        if isempty(normalized_pairs)
            normalized_pairs = [
                (Symbol("W$(index)A"), Symbol("W$(index)B"))
                for index in 1:winding_count
            ]
        end
        length(normalized_pairs) == winding_count ||
            throw(ArgumentError("node-pair count must equal the winding count"))
        return new(
            Int(source_line),
            voltages,
            losses,
            impedances,
            ratings,
            excitation,
            excitation_rating,
            normalized_pairs,
        )
    end
end

struct SaturableTransformerParameterCase <: AbstractTransformerParameterCase
    source_line::Int
    frequency_hz::Float64
    inductance_frequency_hz::Float64
    magnetizing_current_a::Float64
    magnetizing_flux_wb_turn::Float64
    magnetizing_resistance_ohm::Float64
    windings::Vector{TransformerParameterWinding}

    function SaturableTransformerParameterCase(
        source_line::Integer,
        frequency_hz::Real,
        inductance_frequency_hz::Real,
        magnetizing_current_a::Real,
        magnetizing_flux_wb_turn::Real,
        magnetizing_resistance_ohm::Real,
        windings::AbstractVector{TransformerParameterWinding},
    )
        frequency = Float64(frequency_hz)
        inductance_frequency = Float64(inductance_frequency_hz)
        current = Float64(magnetizing_current_a)
        flux = Float64(magnetizing_flux_wb_turn)
        resistance = Float64(magnetizing_resistance_ohm)
        isfinite(frequency) && frequency > 0.0 ||
            throw(ArgumentError("transformer parameter frequency must be finite and positive"))
        isfinite(inductance_frequency) && inductance_frequency >= 0.0 ||
            throw(ArgumentError("inductance frequency must be finite and nonnegative"))
        isfinite(current) && current > 0.0 ||
            throw(ArgumentError("magnetizing current must be finite and positive"))
        isfinite(flux) && flux > 0.0 ||
            throw(ArgumentError("magnetizing flux must be finite and positive"))
        (isfinite(resistance) && resistance > 0.0) || isinf(resistance) ||
            throw(ArgumentError("magnetizing resistance must be positive or Inf"))
        rows = collect(windings)
        isempty(rows) &&
            throw(ArgumentError("saturable-transformer conversion requires at least one winding"))
        length(rows) <= 20 ||
            throw(ArgumentError("OVER41 accepts at most 20 saturable-transformer windings"))
        getfield.(rows, :winding_number) == collect(1:length(rows)) ||
            throw(ArgumentError("transformer windings must be numbered consecutively from one"))
        return new(
            Int(source_line),
            frequency,
            inductance_frequency,
            current,
            flux,
            resistance,
            rows,
        )
    end
end

struct MultiphaseTransformerWinding
    winding_number::Int
    rated_voltage_kv::Float64
    resistance_ohm::Float64
    node_pairs::NTuple{3,NTuple{2,Symbol}}

    function MultiphaseTransformerWinding(
        winding_number::Integer,
        rated_voltage_kv::Real,
        resistance_ohm::Real,
        node_pairs::NTuple{3,<:Tuple},
    )
        number = Int(winding_number)
        voltage = Float64(rated_voltage_kv)
        resistance = Float64(resistance_ohm)
        number > 0 ||
            throw(ArgumentError("multiphase transformer winding number must be positive"))
        isfinite(voltage) && voltage > 0.0 ||
            throw(ArgumentError("multiphase transformer winding voltage must be finite and positive"))
        isfinite(resistance) && resistance >= 0.0 ||
            throw(ArgumentError("multiphase transformer winding resistance must be finite and nonnegative"))
        pairs = ntuple(Val(3)) do phase
            pair = node_pairs[phase]
            return (Symbol(pair[1]), Symbol(pair[2]))
        end
        return new(number, voltage, resistance, pairs)
    end
end

struct TransformerSequenceShortCircuitTest
    from_winding::Int
    to_winding::Int
    load_loss_kw::Float64
    positive_impedance_percent::Float64
    positive_rating_mva::Float64
    zero_impedance_percent::Float64
    zero_rating_mva::Float64
    closed_delta_winding::Int

    function TransformerSequenceShortCircuitTest(
        from_winding::Integer,
        to_winding::Integer,
        load_loss_kw::Real,
        positive_impedance_percent::Real,
        positive_rating_mva::Real,
        zero_impedance_percent::Real,
        zero_rating_mva::Real,
        closed_delta_winding::Integer = 0,
    )
        from = Int(from_winding)
        to = Int(to_winding)
        from > 0 && to > 0 && from != to ||
            throw(ArgumentError("short-circuit test requires two distinct positive winding numbers"))
        loss = Float64(load_loss_kw)
        positive_impedance = Float64(positive_impedance_percent)
        positive_rating = Float64(positive_rating_mva)
        zero_impedance = Float64(zero_impedance_percent)
        zero_rating = Float64(zero_rating_mva)
        all(isfinite, (
            loss,
            positive_impedance,
            positive_rating,
            zero_impedance,
            zero_rating,
        )) ||
            throw(ArgumentError("short-circuit test values must be finite"))
        loss >= 0.0 ||
            throw(ArgumentError("short-circuit load loss must be nonnegative"))
        positive_impedance > 0.0 && zero_impedance > 0.0 ||
            throw(ArgumentError("sequence impedances must be positive"))
        positive_rating > 0.0 && zero_rating > 0.0 ||
            throw(ArgumentError("sequence ratings must be positive"))
        return new(
            min(from, to),
            max(from, to),
            loss,
            positive_impedance,
            positive_rating,
            zero_impedance,
            zero_rating,
            from <= to ? Int(closed_delta_winding) : -Int(closed_delta_winding),
        )
    end
end

struct MultiphaseTransformerParameterCase <: AbstractTransformerParameterCase
    source_line::Int
    frequency_hz::Float64
    positive_exciting_current_percent::Float64
    positive_excitation_rating_mva::Float64
    positive_excitation_loss_kw::Float64
    zero_exciting_current_percent::Float64
    zero_excitation_rating_mva::Float64
    zero_excitation_loss_kw::Float64
    phase_count::Int
    excitation_test_winding::Int
    magnetizing_output_winding::Int
    output_representation::Symbol
    calculate_winding_resistances_from_losses::Bool
    windings::Vector{MultiphaseTransformerWinding}
    short_circuit_tests::Vector{TransformerSequenceShortCircuitTest}

    function MultiphaseTransformerParameterCase(
        source_line::Integer,
        frequency_hz::Real,
        positive_exciting_current_percent::Real,
        positive_excitation_rating_mva::Real,
        positive_excitation_loss_kw::Real,
        zero_exciting_current_percent::Real,
        zero_excitation_rating_mva::Real,
        zero_excitation_loss_kw::Real,
        phase_count::Integer,
        excitation_test_winding::Integer,
        magnetizing_output_winding::Integer,
        output_representation::Symbol,
        calculate_winding_resistances_from_losses::Bool,
        windings::AbstractVector{MultiphaseTransformerWinding},
        short_circuit_tests::AbstractVector{TransformerSequenceShortCircuitTest},
    )
        frequency = Float64(frequency_hz)
        isfinite(frequency) && frequency > 0.0 ||
            throw(ArgumentError("BCTRAN frequency must be finite and positive"))
        phases = Int(phase_count)
        phases in (1, 3) ||
            throw(ArgumentError("BCTRAN phase count must be one or three"))
        representation = output_representation
        representation in (:inverse_inductance, :reactance) ||
            throw(ArgumentError("BCTRAN output representation must be :inverse_inductance or :reactance"))
        winding_rows = collect(windings)
        winding_count = length(winding_rows)
        2 <= winding_count <= 10 ||
            throw(ArgumentError("BCTRAN accepts between two and ten windings"))
        getfield.(winding_rows, :winding_number) == collect(1:winding_count) ||
            throw(ArgumentError("BCTRAN windings must be numbered consecutively from one"))
        tests = collect(short_circuit_tests)
        expected_test_count = winding_count * (winding_count - 1) ÷ 2
        length(tests) == expected_test_count ||
            throw(ArgumentError("BCTRAN requires $expected_test_count pairwise short-circuit tests"))
        pairs = Set{Tuple{Int,Int}}()
        for test in tests
            test.to_winding <= winding_count ||
                throw(ArgumentError("short-circuit test winding exceeds the winding count"))
            pair = (test.from_winding, test.to_winding)
            pair in pairs &&
                throw(ArgumentError("duplicate short-circuit test for windings $(pair[1])-$(pair[2])"))
            push!(pairs, pair)
            delta = abs(test.closed_delta_winding)
            0 <= delta <= winding_count ||
                throw(ArgumentError("closed-delta winding exceeds the winding count"))
            delta in pair &&
                throw(ArgumentError("closed-delta winding must differ from both tested windings"))
        end
        excitation_test = Int(excitation_test_winding)
        output_winding = Int(magnetizing_output_winding)
        0 <= excitation_test <= winding_count ||
            throw(ArgumentError("excitation-test winding exceeds the winding count"))
        0 <= output_winding <= winding_count ||
            throw(ArgumentError("magnetizing-output winding exceeds the winding count"))
        (excitation_test == 0) == (output_winding == 0) ||
            throw(ArgumentError("excitation-test and magnetizing-output windings must both be zero or both be specified"))
        excitation_values = Float64[
            positive_exciting_current_percent,
            positive_excitation_rating_mva,
            positive_excitation_loss_kw,
            zero_exciting_current_percent,
            zero_excitation_rating_mva,
            zero_excitation_loss_kw,
        ]
        phases == 1 &&
            (excitation_values[4:6] .= excitation_values[1:3])
        all(isfinite, excitation_values) ||
            throw(ArgumentError("BCTRAN excitation values must be finite"))
        all(excitation_values[[1, 3, 4, 6]] .>= 0.0) ||
            throw(ArgumentError("BCTRAN excitation current and losses must be nonnegative"))
        all(excitation_values[[2, 5]] .> 0.0) ||
            throw(ArgumentError("BCTRAN excitation ratings must be positive"))
        return new(
            Int(source_line),
            frequency,
            excitation_values...,
            phases,
            excitation_test,
            output_winding,
            representation,
            calculate_winding_resistances_from_losses,
            winding_rows,
            tests,
        )
    end
end

struct TransformerGeneratedBranchRow
    branch_type::Int
    winding_number::Int
    from_node::Symbol
    to_node::Symbol
    resistance_values_ohm::Vector{Float64}
    inductance_or_reactance_values::Vector{Float64}
end

struct TransformerShortCircuitResult
    input::TransformerShortCircuitCase
    admittance_matrix_s::Matrix{ComplexF64}
    impedance_matrix_ohm::Matrix{ComplexF64}
    generated_branches::Vector{TransformerGeneratedBranchRow}
    inverse_residual::Float64
    symmetry_residual::Float64
    minimum_resistance_eigenvalue_ohm::Float64
    physical_checks_passed::Bool
end

struct SaturableTransformerParameterResult
    input::SaturableTransformerParameterCase
    magnetizing_parallel_inductance_h::Float64
    magnetizing_reactance_ohm::Float64
    equivalent_series_resistance_ohm::Float64
    equivalent_series_inductance_h::Float64
    branch_resistance_matrix_ohm::Matrix{Float64}
    branch_inductance_or_reactance_matrix::Matrix{Float64}
    physical_inductance_matrix_h::Matrix{Float64}
    generated_branches::Vector{TransformerGeneratedBranchRow}
    symmetry_residual::Float64
    minimum_resistance_eigenvalue_ohm::Float64
    minimum_inductance_eigenvalue_h::Float64
    physical_checks_passed::Bool
end

struct TransformerMagnetizingShunt
    winding_number::Int
    self_resistance_ohm::Float64
    mutual_resistance_ohm::Float64
end

struct MultiphaseTransformerParameterResult
    input::MultiphaseTransformerParameterCase
    winding_resistances_ohm::Vector{Float64}
    positive_pair_reactance_normalized::Vector{Float64}
    zero_pair_reactance_normalized::Vector{Float64}
    positive_reduced_reactance_normalized::Matrix{Float64}
    zero_reduced_reactance_normalized::Matrix{Float64}
    positive_winding_admittance_normalized::Matrix{Float64}
    zero_winding_admittance_normalized::Matrix{Float64}
    resistance_matrix_ohm::Matrix{Float64}
    inverse_inductance_matrix_per_h::Matrix{Float64}
    reactance_matrix_ohm::Matrix{Float64}
    magnetizing_shunts::Vector{TransformerMagnetizingShunt}
    generated_branches::Vector{TransformerGeneratedBranchRow}
    positive_pair_reconstruction_residual::Float64
    zero_pair_reconstruction_residual::Float64
    inverse_reconstruction_residual::Float64
    symmetry_residual::Float64
    minimum_resistance_eigenvalue_ohm::Float64
    minimum_inverse_inductance_eigenvalue_per_h::Float64
    physical_checks_passed::Bool
end

function _signed_reactive_component(
    impedance_percent::Float64,
    resistance_per_unit::Float64,
)
    magnitude_squared = (impedance_percent / 100.0)^2
    reactive_squared = magnitude_squared - resistance_per_unit^2
    tolerance = 64.0 * eps(Float64) * max(magnitude_squared, resistance_per_unit^2, 1.0)
    reactive_squared >= -tolerance ||
        throw(ArgumentError("transformer losses exceed the specified impedance magnitude"))
    reactive = sqrt(max(reactive_squared, 0.0))
    return signbit(impedance_percent) ? -reactive : reactive
end

function _matrix_diagnostics(matrix::AbstractMatrix{<:Complex})
    symmetry = maximum(abs.(matrix .- transpose(matrix)); init = 0.0)
    real_part = Symmetric(real.(matrix))
    minimum_resistance = minimum(eigvals(real_part))
    return symmetry, minimum_resistance
end

function _short_circuit_two_winding(case::TransformerShortCircuitCase)
    voltage_1, voltage_2 = case.winding_voltages_kv
    rating = only(case.pair_ratings_mva)
    resistance_per_unit = only(case.pair_losses_kw) / (rating * 1000.0)
    reactance_per_unit = _signed_reactive_component(
        only(case.pair_impedance_percent),
        resistance_per_unit,
    )
    magnetizing = case.magnetizing_current_percent *
        case.magnetizing_rating_mva / 200.0
    denominator = (resistance_per_unit^2 + reactance_per_unit^2) / rating
    denominator > 0.0 ||
        throw(ArgumentError("two-winding transformer test produces singular admittance"))
    real_inverse = resistance_per_unit / denominator
    imaginary_inverse = -reactance_per_unit / denominator
    magnetizing > abs(imaginary_inverse) * eps(Float64) ||
        throw(ArgumentError("two-winding transformer admittance is singular"))
    admittance = ComplexF64[
        complex(-real_inverse / voltage_1^2, (magnetizing - imaginary_inverse) / voltage_1^2) complex(real_inverse / (voltage_1 * voltage_2), imaginary_inverse / (voltage_1 * voltage_2))
        complex(real_inverse / (voltage_1 * voltage_2), imaginary_inverse / (voltage_1 * voltage_2)) complex(-real_inverse / voltage_2^2, (magnetizing - imaginary_inverse) / voltage_2^2)
    ]
    return admittance
end

function _short_circuit_three_winding(case::TransformerShortCircuitCase)
    pair_resistance = Float64[]
    pair_reactance = Float64[]
    for index in eachindex(case.pair_losses_kw)
        rating = case.pair_ratings_mva[index]
        resistance_per_unit = case.pair_losses_kw[index] / (rating * 1000.0)
        push!(pair_resistance, resistance_per_unit / rating)
        push!(
            pair_reactance,
            _signed_reactive_component(
                case.pair_impedance_percent[index],
                resistance_per_unit,
            ) / rating,
        )
    end
    resistance_sum = sum(pair_resistance) / 2.0
    reactance_sum = sum(pair_reactance) / 2.0
    real_star = resistance_sum - pair_resistance[1]
    imaginary_star = reactance_sum - pair_reactance[1]
    denominator_real =
        pair_resistance[2] * pair_resistance[3] -
        pair_reactance[2] * pair_reactance[3] -
        real_star^2 + imaginary_star^2
    denominator_imaginary =
        pair_resistance[2] * pair_reactance[3] +
        pair_reactance[2] * pair_resistance[3] -
        2.0 * real_star * imaginary_star
    denominator_norm = denominator_real^2 + denominator_imaginary^2
    denominator_norm > 0.0 ||
        throw(ArgumentError("three-winding transformer test produces singular admittance"))
    inverse_real = denominator_real / denominator_norm
    inverse_imaginary = -denominator_imaginary / denominator_norm
    magnetizing = case.magnetizing_current_percent *
        case.magnetizing_rating_mva / 300.0
    admittance = zeros(ComplexF64, 3, 3)
    pair_index = 0
    for row in 1:3
        for column in 1:row
            voltage_product =
                case.winding_voltages_kv[row] * case.winding_voltages_kv[column]
            if column != row
                pair_index += 1
                real_component =
                    (pair_resistance[pair_index] - resistance_sum) / voltage_product
                imaginary_component =
                    (pair_reactance[pair_index] - reactance_sum) / voltage_product
                value = complex(
                    imaginary_component * inverse_imaginary -
                        real_component * inverse_real,
                    -imaginary_component * inverse_real -
                        real_component * inverse_imaginary,
                )
            else
                opposite_pair = 4 - row
                real_component = pair_resistance[opposite_pair] / voltage_product
                imaginary_component = pair_reactance[opposite_pair] / voltage_product
                real_value =
                    imaginary_component * inverse_imaginary -
                    real_component * inverse_real
                imaginary_value =
                    -imaginary_component * inverse_real -
                    real_component * inverse_imaginary
                magnetizing_component = magnetizing / voltage_product
                magnetizing_component > abs(imaginary_value) * eps(Float64) ||
                    throw(ArgumentError("three-winding transformer admittance is singular"))
                value = complex(real_value, imaginary_value + magnetizing_component)
            end
            admittance[row, column] = value
            admittance[column, row] = value
        end
    end
    return admittance
end

function _generated_complex_branches(
    matrix::AbstractMatrix{<:Complex},
    node_pairs::AbstractVector{<:Tuple},
)
    return TransformerGeneratedBranchRow[
        TransformerGeneratedBranchRow(
            50 + row,
            row,
            node_pairs[row][1],
            node_pairs[row][2],
            Float64[real(matrix[row, column]) for column in 1:row],
            Float64[imag(matrix[row, column]) for column in 1:row],
        )
        for row in axes(matrix, 1)
    ]
end

function transformer_short_circuit_parameters(case::TransformerShortCircuitCase)
    admittance = length(case.winding_voltages_kv) == 2 ?
        _short_circuit_two_winding(case) :
        _short_circuit_three_winding(case)
    impedance = -inv(admittance)
    identity_residual = maximum(
        abs.(admittance * impedance + I);
        init = 0.0,
    )
    symmetry, minimum_resistance = _matrix_diagnostics(impedance)
    scale = maximum(abs, impedance; init = 1.0)
    tolerance = 5.0e-11 * max(scale, 1.0)
    checks =
        all(value -> isfinite(real(value)) && isfinite(imag(value)), impedance) &&
        identity_residual <= 5.0e-11 &&
        symmetry <= tolerance &&
        minimum_resistance >= -tolerance
    return TransformerShortCircuitResult(
        case,
        admittance,
        impedance,
        _generated_complex_branches(impedance, case.node_pairs),
        identity_residual,
        symmetry,
        minimum_resistance,
        checks,
    )
end

function _saturable_branch_inductance_value(
    equivalent_inductance_h::Float64,
    inductance_frequency_hz::Float64,
)
    if inductance_frequency_hz == 0.0
        return equivalent_inductance_h * 1000.0
    end
    return equivalent_inductance_h * 2.0 * pi * inductance_frequency_hz
end

function saturable_transformer_parameters(case::SaturableTransformerParameterCase)
    omega = 2.0 * pi * case.frequency_hz
    parallel_inductance =
        case.magnetizing_flux_wb_turn / case.magnetizing_current_a
    reactance = omega * parallel_inductance
    if isinf(case.magnetizing_resistance_ohm)
        series_resistance = 0.0
        series_inductance = parallel_inductance
    else
        resistance_squared = case.magnetizing_resistance_ohm^2
        reactance_squared = reactance^2
        denominator = resistance_squared + reactance_squared
        series_resistance =
            case.magnetizing_resistance_ohm * reactance_squared / denominator
        series_inductance =
            resistance_squared * parallel_inductance / denominator
    end
    branch_inductance = _saturable_branch_inductance_value(
        series_inductance,
        case.inductance_frequency_hz,
    )
    turns = Float64[
        winding.rated_voltage_kv / first(case.windings).rated_voltage_kv
        for winding in case.windings
    ]
    resistance_matrix = series_resistance .* (turns * transpose(turns))
    branch_inductance_matrix =
        branch_inductance .* (turns * transpose(turns))
    for (index, winding) in enumerate(case.windings)
        resistance_matrix[index, index] += winding.resistance_ohm
        branch_inductance_matrix[index, index] +=
            winding.inductance_or_reactance
    end
    physical_inductance_matrix =
        case.inductance_frequency_hz == 0.0 ?
        branch_inductance_matrix ./ 1000.0 :
        branch_inductance_matrix ./
        (2.0 * pi * case.inductance_frequency_hz)
    symmetry = max(
        maximum(abs.(resistance_matrix .- transpose(resistance_matrix)); init = 0.0),
        maximum(
            abs.(physical_inductance_matrix .- transpose(physical_inductance_matrix));
            init = 0.0,
        ),
    )
    minimum_resistance = minimum(eigvals(Symmetric(resistance_matrix)))
    minimum_inductance = minimum(eigvals(Symmetric(physical_inductance_matrix)))
    scale = max(
        maximum(abs, resistance_matrix; init = 1.0),
        maximum(abs, physical_inductance_matrix; init = 1.0),
    )
    tolerance = 5.0e-12 * max(scale, 1.0)
    checks =
        all(isfinite, resistance_matrix) &&
        all(isfinite, branch_inductance_matrix) &&
        all(isfinite, physical_inductance_matrix) &&
        symmetry <= tolerance &&
        minimum_resistance >= -tolerance &&
        minimum_inductance >= -tolerance
    branches = TransformerGeneratedBranchRow[
        TransformerGeneratedBranchRow(
            50 + row,
            row,
            case.windings[row].from_node,
            case.windings[row].to_node,
            copy(resistance_matrix[row, 1:row]),
            copy(branch_inductance_matrix[row, 1:row]),
        )
        for row in eachindex(case.windings)
    ]
    return SaturableTransformerParameterResult(
        case,
        parallel_inductance,
        reactance,
        series_resistance,
        series_inductance,
        resistance_matrix,
        branch_inductance_matrix,
        physical_inductance_matrix,
        branches,
        symmetry,
        minimum_resistance,
        minimum_inductance,
        checks,
    )
end

function _bctran_pair_index(
    tests::AbstractVector{TransformerSequenceShortCircuitTest},
    first_winding::Int,
    second_winding::Int,
)
    first, second = minmax(first_winding, second_winding)
    index = findfirst(tests) do test
        test.from_winding == first && test.to_winding == second
    end
    index === nothing &&
        throw(ArgumentError("missing short-circuit test for windings $first-$second"))
    return index
end

function _bctran_winding_resistances(
    case::MultiphaseTransformerParameterCase,
)
    winding_count = length(case.windings)
    voltages = getfield.(case.windings, :rated_voltage_kv)
    normalized = Float64[
        winding.resistance_ohm / (3.0 * winding.rated_voltage_kv^2)
        for winding in case.windings
    ]
    derive_from_losses =
        case.calculate_winding_resistances_from_losses &&
        winding_count <= 3 &&
        all(test -> test.load_loss_kw > 0.0, case.short_circuit_tests)
    if derive_from_losses
        loss_resistance = zeros(Float64, max(4, winding_count))
        for test in case.short_circuit_tests
            value =
                test.load_loss_kw * 1.0e-3 /
                test.positive_rating_mva^2 / 2.0
            loss_resistance[test.from_winding] += value
            loss_resistance[test.to_winding] += value
            opposite = 6 - test.from_winding - test.to_winding
            opposite <= 0 && (opposite = 4)
            loss_resistance[opposite] -= value
        end
        normalized .= max.(loss_resistance[1:winding_count], 0.0)
    end
    physical = normalized .* (3.0 .* voltages .^ 2)
    return normalized, physical
end

function _bctran_pair_reactances(
    case::MultiphaseTransformerParameterCase,
    winding_resistance_normalized::AbstractVector{<:Real},
)
    positive = Float64[]
    zero = Float64[]
    for test in case.short_circuit_tests
        resistance =
            winding_resistance_normalized[test.from_winding] +
            winding_resistance_normalized[test.to_winding]
        resistance_squared = resistance^2
        positive_magnitude =
            test.positive_impedance_percent * 0.01 /
            test.positive_rating_mva
        positive_squared = positive_magnitude^2 - resistance_squared
        positive_tolerance =
            64.0 * eps(Float64) *
            max(positive_magnitude^2, resistance_squared, 1.0)
        positive_squared >= -positive_tolerance ||
            throw(ArgumentError(
                "positive-sequence load losses or winding resistances exceed the impedance for windings $(test.from_winding)-$(test.to_winding)",
            ))
        push!(positive, sqrt(max(positive_squared, 0.0)))
        zero_magnitude =
            test.zero_impedance_percent * 0.01 /
            test.zero_rating_mva
        if test.closed_delta_winding == 0
            zero_squared = zero_magnitude^2 - resistance_squared
            zero_tolerance =
                64.0 * eps(Float64) *
                max(zero_magnitude^2, resistance_squared, 1.0)
            zero_squared >= -zero_tolerance ||
                throw(ArgumentError(
                    "zero-sequence winding resistances exceed the impedance for windings $(test.from_winding)-$(test.to_winding)",
                ))
            push!(zero, sqrt(max(zero_squared, 0.0)))
        else
            push!(zero, zero_magnitude)
        end
    end
    return positive, zero
end

function _bctran_closed_delta_newton(
    initial_reactance::Float64,
    winding_resistance_normalized::AbstractVector{<:Real},
    test_from::Int,
    test_to::Int,
    delta_winding::Int,
    measured_test_reactance::Float64,
    first_adjacent_reactance::Float64,
    second_adjacent_reactance::Float64,
)
    reactance = initial_reactance
    tolerance = 32.0 * eps(Float64)
    for _ in 1:200
        resistance_test_from = winding_resistance_normalized[test_from]
        resistance_test_to = winding_resistance_normalized[test_to]
        resistance_delta = winding_resistance_normalized[delta_winding]
        resistance_total = resistance_test_to + resistance_delta
        reactance_difference =
            second_adjacent_reactance - first_adjacent_reactance
        real_numerator =
            reactance^2 -
            2.0 * first_adjacent_reactance * reactance +
            resistance_test_from * resistance_total +
            resistance_test_to * resistance_delta -
            first_adjacent_reactance * reactance_difference
        imaginary_numerator =
            2.0 * resistance_delta * reactance +
            resistance_test_from * second_adjacent_reactance +
            resistance_delta * reactance_difference +
            resistance_test_to * first_adjacent_reactance
        target =
            measured_test_reactance^2 *
            (resistance_total^2 + second_adjacent_reactance^2)
        residual =
            real_numerator^2 + imaginary_numerator^2 - target
        derivative =
            2.0 * real_numerator *
            (2.0 * reactance - 2.0 * first_adjacent_reactance) +
            4.0 * resistance_delta * imaginary_numerator
        abs(derivative) > eps(Float64) ||
            throw(ArgumentError("closed-delta correction has a zero Newton derivative"))
        increment = -residual / derivative
        reactance += increment
        abs(increment) <= tolerance && return reactance
    end
    throw(ArgumentError("closed-delta correction did not converge in 200 iterations"))
end

function _bctran_correct_closed_delta!(
    zero_pair_reactance::Vector{Float64},
    case::MultiphaseTransformerParameterCase,
    winding_resistance_normalized::AbstractVector{<:Real},
)
    delta_designations =
        getfield.(case.short_circuit_tests, :closed_delta_winding)
    for test_index in eachindex(case.short_circuit_tests)
        signed_delta = delta_designations[test_index]
        signed_delta == 0 && continue
        test = case.short_circuit_tests[test_index]
        delta_winding = abs(signed_delta)
        test_from = test.from_winding
        test_to = test.to_winding
        first_adjacent_index = _bctran_pair_index(
            case.short_circuit_tests,
            test_from,
            delta_winding,
        )
        second_adjacent_index = _bctran_pair_index(
            case.short_circuit_tests,
            test_to,
            delta_winding,
        )
        if signed_delta < 0
            test_from, test_to = test_to, test_from
            first_adjacent_index, second_adjacent_index =
                second_adjacent_index, first_adjacent_index
        end
        measured_test_reactance = zero_pair_reactance[test_index]
        first_adjacent_reactance =
            zero_pair_reactance[first_adjacent_index]
        second_adjacent_reactance =
            zero_pair_reactance[second_adjacent_index]
        first_component =
            (
                measured_test_reactance +
                first_adjacent_reactance -
                second_adjacent_reactance
            ) / 2.0
        first_adjacent_delta = delta_designations[first_adjacent_index]
        second_adjacent_delta = delta_designations[second_adjacent_index]
        if first_adjacent_delta == 0 && second_adjacent_delta == 0
            first_component = _bctran_closed_delta_newton(
                first_component,
                winding_resistance_normalized,
                test_from,
                test_to,
                delta_winding,
                measured_test_reactance,
                first_adjacent_reactance,
                second_adjacent_reactance,
            )
            second_component =
                second_adjacent_reactance -
                first_adjacent_reactance +
                first_component
            third_component = first_adjacent_reactance - first_component
        else
            first_adjacent_delta != 0 && second_adjacent_delta == 0 ||
                throw(ArgumentError(
                    "inconsistent closed-delta designations around winding $delta_winding",
                ))
            adjacent_test = case.short_circuit_tests[first_adjacent_index]
            adjacent_open_winding =
                first_adjacent_delta < 0 ?
                adjacent_test.to_winding :
                adjacent_test.from_winding
            adjacent_open_winding == test_from ||
                throw(ArgumentError(
                    "closed-delta designation is inconsistent with the adjacent test",
                ))
            second_component = second_adjacent_reactance / 2.0
            third_component = second_component
            resistance_delta =
                winding_resistance_normalized[delta_winding]
            resistance_test_to =
                winding_resistance_normalized[test_to]
            resistance_total = resistance_delta + resistance_test_to
            denominator =
                resistance_total^2 + second_adjacent_reactance^2
            denominator > 0.0 ||
                throw(ArgumentError("closed-delta correction has a zero denominator"))
            real_product =
                (
                    resistance_delta * resistance_test_to -
                    second_component^2
                )
            imaginary_product =
                second_component * resistance_total
            equivalent_real =
                (
                    real_product * resistance_total +
                    imaginary_product * second_adjacent_reactance
                ) / denominator
            equivalent_imaginary =
                (
                    resistance_total * imaginary_product -
                    real_product * second_adjacent_reactance
                ) / denominator
            remaining_squared =
                measured_test_reactance^2 -
                (
                    winding_resistance_normalized[test_from] +
                    equivalent_real
                )^2
            remaining_squared >= -64.0 * eps(Float64) ||
                throw(ArgumentError("closed-delta correction has no physical solution"))
            first_component =
                sqrt(max(remaining_squared, 0.0)) -
                equivalent_imaginary
        end
        zero_pair_reactance[test_index] =
            first_component + second_component
        delta_designations[test_index] = 0
        if first_adjacent_delta != 0
            zero_pair_reactance[first_adjacent_index] =
                first_component + third_component
            delta_designations[first_adjacent_index] = 0
        end
    end
    return zero_pair_reactance
end

function _bctran_reduced_reactance(
    pair_reactance::AbstractVector{<:Real},
    tests::AbstractVector{TransformerSequenceShortCircuitTest},
    winding_count::Int,
)
    reduced = zeros(Float64, winding_count - 1, winding_count - 1)
    for winding in 1:(winding_count - 1)
        reduced[winding, winding] =
            pair_reactance[_bctran_pair_index(tests, winding, winding_count)]
    end
    for row in 2:(winding_count - 1)
        for column in 1:(row - 1)
            pair =
                pair_reactance[_bctran_pair_index(tests, column, row)]
            value =
                (
                    reduced[row, row] +
                    reduced[column, column] -
                    pair
                ) / 2.0
            reduced[row, column] = value
            reduced[column, row] = value
        end
    end
    return reduced
end

function _bctran_expand_admittance(
    reduced_reactance::AbstractMatrix{<:Real},
)
    reduced_admittance = inv(Symmetric(Matrix(reduced_reactance)))
    winding_count = size(reduced_reactance, 1) + 1
    expanded = zeros(Float64, winding_count, winding_count)
    expanded[1:(end - 1), 1:(end - 1)] .= reduced_admittance
    expanded[1:(end - 1), end] .= -sum(reduced_admittance; dims = 2)[:]
    expanded[end, 1:(end - 1)] .= -sum(reduced_admittance; dims = 1)[:]
    expanded[end, end] = sum(reduced_admittance)
    return expanded
end

function _bctran_pair_reconstruction_residual(
    reduced::AbstractMatrix{<:Real},
    pair_reactance::AbstractVector{<:Real},
    tests::AbstractVector{TransformerSequenceShortCircuitTest},
)
    winding_count = size(reduced, 1) + 1
    maximum_error = 0.0
    for (index, test) in enumerate(tests)
        from = test.from_winding
        to = test.to_winding
        reconstructed = if to == winding_count
            reduced[from, from]
        else
            reduced[from, from] +
            reduced[to, to] -
            2.0 * reduced[from, to]
        end
        maximum_error =
            max(maximum_error, abs(reconstructed - pair_reactance[index]))
    end
    return maximum_error
end

function _bctran_excitation(
    case::MultiphaseTransformerParameterCase,
    winding_resistance_normalized::AbstractVector{<:Real},
    positive_pair_reactance::AbstractVector{<:Real},
    zero_pair_reactance::AbstractVector{<:Real},
)
    positive_loss = case.positive_excitation_loss_kw * 1.0e-3
    zero_loss = case.zero_excitation_loss_kw * 1.0e-3
    positive_apparent =
        case.positive_exciting_current_percent *
        case.positive_excitation_rating_mva * 0.01
    zero_apparent =
        case.zero_exciting_current_percent *
        case.zero_excitation_rating_mva * 0.01
    positive_reactive_squared = positive_apparent^2 - positive_loss^2
    zero_reactive_squared = zero_apparent^2 - zero_loss^2
    positive_reactive_squared >= -64.0 * eps(Float64) ||
        throw(ArgumentError("positive-sequence excitation loss exceeds exciting apparent power"))
    zero_reactive_squared >= -64.0 * eps(Float64) ||
        throw(ArgumentError("zero-sequence excitation loss exceeds exciting apparent power"))
    positive_reactive = sqrt(max(positive_reactive_squared, 0.0))
    zero_reactive = sqrt(max(zero_reactive_squared, 0.0))
    if case.phase_count == 1
        zero_loss = positive_loss
        zero_reactive = positive_reactive
    end
    if case.excitation_test_winding > 0
        same_winding =
            case.excitation_test_winding ==
            case.magnetizing_output_winding
        pair_index =
            same_winding ?
            nothing :
            _bctran_pair_index(
                case.short_circuit_tests,
                case.excitation_test_winding,
                case.magnetizing_output_winding,
            )
        series_resistance =
            winding_resistance_normalized[case.excitation_test_winding]
        function adjust(loss, reactive, series_reactance)
            magnitude_squared = loss^2 + reactive^2
            minimum_loss = magnitude_squared * series_resistance
            adjusted_loss = max(loss, minimum_loss)
            magnitude_squared = adjusted_loss^2 + reactive^2
            magnitude_squared <= eps(Float64) &&
                return adjusted_loss, reactive
            real_component =
                adjusted_loss / magnitude_squared - series_resistance
            imaginary_component =
                reactive / magnitude_squared - series_reactance
            denominator = real_component^2 + imaginary_component^2
            denominator > 0.0 ||
                throw(ArgumentError("excitation correction has a zero denominator"))
            return (
                real_component / denominator,
                imaginary_component / denominator,
            )
        end
        positive_loss, positive_reactive = adjust(
            positive_loss,
            positive_reactive,
            pair_index === nothing ? 0.0 : positive_pair_reactance[pair_index],
        )
        zero_loss, zero_reactive = adjust(
            zero_loss,
            zero_reactive,
            pair_index === nothing ? 0.0 : zero_pair_reactance[pair_index],
        )
    end
    zero_loss = max(zero_loss, positive_loss)
    shunt_self_normalized = 0.0
    shunt_mutual_normalized = 0.0
    if abs(positive_loss) >= eps(Float64) &&
       abs(zero_loss) >= eps(Float64)
        positive_resistance = 1.0 / positive_loss
        shunt_mutual_normalized =
            (positive_loss / zero_loss - 1.0) *
            positive_resistance / 3.0
        shunt_self_normalized =
            shunt_mutual_normalized + positive_resistance
    end
    reactive_mutual = (zero_reactive - positive_reactive) / 3.0
    reactive_self = reactive_mutual + positive_reactive
    return (
        reactive_self = reactive_self,
        reactive_mutual = reactive_mutual,
        shunt_self_normalized = shunt_self_normalized,
        shunt_mutual_normalized = shunt_mutual_normalized,
    )
end

function multiphase_transformer_parameters(
    case::MultiphaseTransformerParameterCase,
)
    winding_count = length(case.windings)
    normalized_resistance, physical_resistance =
        _bctran_winding_resistances(case)
    positive_pair, zero_pair =
        _bctran_pair_reactances(case, normalized_resistance)
    _bctran_correct_closed_delta!(
        zero_pair,
        case,
        normalized_resistance,
    )
    positive_reduced = _bctran_reduced_reactance(
        positive_pair,
        case.short_circuit_tests,
        winding_count,
    )
    zero_reduced = _bctran_reduced_reactance(
        zero_pair,
        case.short_circuit_tests,
        winding_count,
    )
    positive_admittance = _bctran_expand_admittance(positive_reduced)
    zero_admittance = _bctran_expand_admittance(zero_reduced)
    phase_self =
        (2.0 .* positive_admittance .+ zero_admittance) ./ 3.0
    phase_mutual =
        (zero_admittance .- positive_admittance) ./ 3.0
    excitation = _bctran_excitation(
        case,
        normalized_resistance,
        positive_pair,
        zero_pair,
    )
    if case.magnetizing_output_winding > 0
        winding = case.magnetizing_output_winding
        phase_self[winding, winding] += excitation.reactive_self
        phase_mutual[winding, winding] += excitation.reactive_mutual
    else
        for winding in 1:winding_count
            phase_self[winding, winding] +=
                excitation.reactive_self / winding_count
            phase_mutual[winding, winding] +=
                excitation.reactive_mutual / winding_count
        end
    end
    voltages = getfield.(case.windings, :rated_voltage_kv)
    omega = 2.0 * pi * case.frequency_hz
    scaling = omega ./ (3.0 .* (voltages * transpose(voltages)))
    phase_self .*= scaling
    phase_mutual .*= scaling
    matrix_size = case.phase_count * winding_count
    resistance_matrix = zeros(Float64, matrix_size, matrix_size)
    inverse_inductance_matrix = zeros(Float64, matrix_size, matrix_size)
    for phase_row in 1:case.phase_count
        for phase_column in 1:case.phase_count
            row_range =
                ((phase_row - 1) * winding_count + 1):(phase_row * winding_count)
            column_range =
                ((phase_column - 1) * winding_count + 1):(phase_column * winding_count)
            inverse_inductance_matrix[row_range, column_range] .=
                phase_row == phase_column ? phase_self : phase_mutual
        end
    end
    for phase in 1:case.phase_count
        offset = (phase - 1) * winding_count
        for winding in 1:winding_count
            resistance_matrix[offset + winding, offset + winding] =
                physical_resistance[winding]
        end
    end
    inverse_scale =
        maximum(abs, inverse_inductance_matrix; init = 1.0)
    inverse_tolerance =
        128.0 * eps(Float64) * max(inverse_scale, 1.0)
    minimum_inverse_inductance =
        minimum(eigvals(Symmetric(inverse_inductance_matrix)))
    inductance_matrix =
        minimum_inverse_inductance > inverse_tolerance ?
        inv(Symmetric(inverse_inductance_matrix)) :
        pinv(Symmetric(inverse_inductance_matrix))
    reactance_matrix = omega .* Matrix(inductance_matrix)
    inverse_reconstruction = maximum(
        abs.(
            inverse_inductance_matrix *
            Matrix(inductance_matrix) *
            inverse_inductance_matrix -
            inverse_inductance_matrix
        );
        init = 0.0,
    )
    shunts = TransformerMagnetizingShunt[]
    if excitation.shunt_self_normalized != 0.0 ||
       excitation.shunt_mutual_normalized != 0.0
        target_windings =
            case.magnetizing_output_winding > 0 ?
            (case.magnetizing_output_winding:case.magnetizing_output_winding) :
            (1:winding_count)
        multiplier =
            case.magnetizing_output_winding > 0 ? 1.0 : winding_count
        for winding in target_windings
            normalized_voltage_squared =
                1.0 / (3.0 * case.windings[winding].rated_voltage_kv^2)
            push!(
                shunts,
                TransformerMagnetizingShunt(
                    winding,
                    excitation.shunt_self_normalized *
                    multiplier / normalized_voltage_squared,
                    excitation.shunt_mutual_normalized *
                    multiplier / normalized_voltage_squared,
                ),
            )
        end
    end
    output_matrix =
        case.output_representation == :inverse_inductance ?
        inverse_inductance_matrix :
        reactance_matrix
    branches = TransformerGeneratedBranchRow[]
    for row in 1:matrix_size
        phase = fld(row - 1, winding_count) + 1
        winding = mod(row - 1, winding_count) + 1
        nodes = case.windings[winding].node_pairs[phase]
        push!(
            branches,
            TransformerGeneratedBranchRow(
                row,
                winding,
                nodes[1],
                nodes[2],
                copy(resistance_matrix[row, 1:row]),
                copy(output_matrix[row, 1:row]),
            ),
        )
    end
    positive_residual = _bctran_pair_reconstruction_residual(
        positive_reduced,
        positive_pair,
        case.short_circuit_tests,
    )
    zero_residual = _bctran_pair_reconstruction_residual(
        zero_reduced,
        zero_pair,
        case.short_circuit_tests,
    )
    symmetry = max(
        maximum(
            abs.(inverse_inductance_matrix .- transpose(inverse_inductance_matrix));
            init = 0.0,
        ),
        maximum(
            abs.(reactance_matrix .- transpose(reactance_matrix));
            init = 0.0,
        ),
    )
    minimum_resistance =
        minimum(eigvals(Symmetric(resistance_matrix)))
    scale = max(
        maximum(abs, reactance_matrix; init = 1.0),
        maximum(abs, inverse_inductance_matrix; init = 1.0),
    )
    tolerance = 5.0e-10 * max(scale, 1.0)
    checks =
        all(isfinite, resistance_matrix) &&
        all(isfinite, inverse_inductance_matrix) &&
        all(isfinite, reactance_matrix) &&
        positive_residual <= 5.0e-12 &&
        zero_residual <= 5.0e-12 &&
        inverse_reconstruction <= tolerance &&
        symmetry <= tolerance &&
        minimum_resistance >= -tolerance &&
        minimum_inverse_inductance >= -tolerance
    return MultiphaseTransformerParameterResult(
        case,
        physical_resistance,
        positive_pair,
        zero_pair,
        positive_reduced,
        zero_reduced,
        positive_admittance,
        zero_admittance,
        resistance_matrix,
        inverse_inductance_matrix,
        reactance_matrix,
        shunts,
        branches,
        positive_residual,
        zero_residual,
        inverse_reconstruction,
        symmetry,
        minimum_resistance,
        minimum_inverse_inductance,
        checks,
    )
end

end
