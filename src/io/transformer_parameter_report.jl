module TransformerParameterReport

using Printf

using ..TransformerParameterStudy:
    TransformerParameterStudyResult
using ..TransformerParameters:
    MultiphaseTransformerParameterResult,
    SaturableTransformerParameterResult,
    TransformerShortCircuitResult

export transformer_parameter_report_text,
       write_transformer_parameter_report

function _write_matrix(
    io::IO,
    label::AbstractString,
    matrix::AbstractMatrix{<:Real},
)
    println(io, label)
    for row in axes(matrix, 1)
        @printf(io, "%3d", row)
        for column in axes(matrix, 2)
            @printf(io, " %22.13e", matrix[row, column])
        end
        println(io)
    end
end

function _write_matrix(
    io::IO,
    label::AbstractString,
    matrix::AbstractMatrix{<:Complex},
)
    println(io, label)
    for row in axes(matrix, 1)
        @printf(io, "%3d", row)
        for column in axes(matrix, 2)
            @printf(
                io,
                " %22.13e%+22.13ej",
                real(matrix[row, column]),
                imag(matrix[row, column]),
            )
        end
        println(io)
    end
end

function _write_branches(io::IO, result)
    println(io, "GENERATED BRANCH ROWS")
    for row in result.generated_branches
        @printf(
            io,
            "%2d %-6s %-6s",
            row.branch_type,
            String(row.from_node),
            String(row.to_node),
        )
        for index in eachindex(row.resistance_values_ohm)
            @printf(
                io,
                " %22.13e %22.13e",
                row.resistance_values_ohm[index],
                row.inductance_or_reactance_values[index],
            )
        end
        println(io)
    end
end

function _write_case(
    io::IO,
    case_index::Int,
    result::MultiphaseTransformerParameterResult,
)
    println(io)
    println(io, "CASE $case_index MULTIPHASE_TRANSFORMER")
    println(io, "WINDING_COUNT ", length(result.input.windings))
    println(io, "PHASE_COUNT ", result.input.phase_count)
    println(io, "OUTPUT_REPRESENTATION ", result.input.output_representation)
    _write_matrix(
        io,
        "POSITIVE_REDUCED_REACTANCE_NORMALIZED",
        result.positive_reduced_reactance_normalized,
    )
    _write_matrix(
        io,
        "ZERO_REDUCED_REACTANCE_NORMALIZED",
        result.zero_reduced_reactance_normalized,
    )
    _write_matrix(
        io,
        "RESISTANCE_MATRIX_OHM",
        result.resistance_matrix_ohm,
    )
    _write_matrix(
        io,
        "INVERSE_INDUCTANCE_MATRIX_PER_H",
        result.inverse_inductance_matrix_per_h,
    )
    _write_matrix(
        io,
        "REACTANCE_MATRIX_OHM",
        result.reactance_matrix_ohm,
    )
    println(io, "MAGNETIZING_SHUNT_COUNT ", length(result.magnetizing_shunts))
    for shunt in result.magnetizing_shunts
        @printf(
            io,
            "MAGNETIZING_SHUNT %d %.17g %.17g\n",
            shunt.winding_number,
            shunt.self_resistance_ohm,
            shunt.mutual_resistance_ohm,
        )
    end
    _write_branches(io, result)
    @printf(
        io,
        "POSITIVE_PAIR_RECONSTRUCTION_RESIDUAL %.17g\n",
        result.positive_pair_reconstruction_residual,
    )
    @printf(
        io,
        "ZERO_PAIR_RECONSTRUCTION_RESIDUAL %.17g\n",
        result.zero_pair_reconstruction_residual,
    )
    @printf(
        io,
        "INVERSE_RECONSTRUCTION_RESIDUAL %.17g\n",
        result.inverse_reconstruction_residual,
    )
    @printf(io, "SYMMETRY_RESIDUAL %.17g\n", result.symmetry_residual)
    @printf(
        io,
        "MINIMUM_RESISTANCE_EIGENVALUE_OHM %.17g\n",
        result.minimum_resistance_eigenvalue_ohm,
    )
    @printf(
        io,
        "MINIMUM_INVERSE_INDUCTANCE_EIGENVALUE_PER_H %.17g\n",
        result.minimum_inverse_inductance_eigenvalue_per_h,
    )
    println(io, "PHYSICAL_CHECKS_PASSED ", result.physical_checks_passed)
end

function _write_case(
    io::IO,
    case_index::Int,
    result::TransformerShortCircuitResult,
)
    println(io)
    println(io, "CASE $case_index SHORT_CIRCUIT")
    println(io, "WINDING_COUNT ", length(result.input.winding_voltages_kv))
    _write_matrix(io, "ADMITTANCE_MATRIX_S", result.admittance_matrix_s)
    _write_matrix(io, "IMPEDANCE_MATRIX_OHM", result.impedance_matrix_ohm)
    _write_branches(io, result)
    @printf(io, "INVERSE_RESIDUAL %.17g\n", result.inverse_residual)
    @printf(io, "SYMMETRY_RESIDUAL %.17g\n", result.symmetry_residual)
    @printf(
        io,
        "MINIMUM_RESISTANCE_EIGENVALUE_OHM %.17g\n",
        result.minimum_resistance_eigenvalue_ohm,
    )
    println(io, "PHYSICAL_CHECKS_PASSED ", result.physical_checks_passed)
end

function _write_case(
    io::IO,
    case_index::Int,
    result::SaturableTransformerParameterResult,
)
    println(io)
    println(io, "CASE $case_index SATURABLE_CONVERSION")
    println(io, "WINDING_COUNT ", length(result.input.windings))
    @printf(
        io,
        "MAGNETIZING_PARALLEL_INDUCTANCE_H %.17g\n",
        result.magnetizing_parallel_inductance_h,
    )
    @printf(
        io,
        "MAGNETIZING_REACTANCE_OHM %.17g\n",
        result.magnetizing_reactance_ohm,
    )
    @printf(
        io,
        "EQUIVALENT_SERIES_RESISTANCE_OHM %.17g\n",
        result.equivalent_series_resistance_ohm,
    )
    @printf(
        io,
        "EQUIVALENT_SERIES_INDUCTANCE_H %.17g\n",
        result.equivalent_series_inductance_h,
    )
    _write_matrix(
        io,
        "BRANCH_RESISTANCE_MATRIX_OHM",
        result.branch_resistance_matrix_ohm,
    )
    _write_matrix(
        io,
        "BRANCH_INDUCTANCE_OR_REACTANCE_MATRIX",
        result.branch_inductance_or_reactance_matrix,
    )
    _write_matrix(
        io,
        "PHYSICAL_INDUCTANCE_MATRIX_H",
        result.physical_inductance_matrix_h,
    )
    _write_branches(io, result)
    @printf(io, "SYMMETRY_RESIDUAL %.17g\n", result.symmetry_residual)
    @printf(
        io,
        "MINIMUM_RESISTANCE_EIGENVALUE_OHM %.17g\n",
        result.minimum_resistance_eigenvalue_ohm,
    )
    @printf(
        io,
        "MINIMUM_INDUCTANCE_EIGENVALUE_H %.17g\n",
        result.minimum_inductance_eigenvalue_h,
    )
    println(io, "PHYSICAL_CHECKS_PASSED ", result.physical_checks_passed)
end

function transformer_parameter_report_text(
    study::TransformerParameterStudyResult,
)
    return sprint() do io
        println(io, "AIMORA TRANSFORMER PARAMETER STUDY")
        println(io, "SOURCE ", study.source)
        println(io, "CASE_COUNT ", length(study.case_results))
        for (index, result) in enumerate(study.case_results)
            _write_case(io, index, result)
        end
        println(io)
        println(io, "GENERATED_BRANCH_COUNT ", study.generated_branch_count)
        println(io, "PHYSICAL_CHECKS_PASSED ", study.physical_checks_passed)
    end
end

function write_transformer_parameter_report(
    path::AbstractString,
    study::TransformerParameterStudyResult,
)
    study.physical_checks_passed ||
        throw(ArgumentError("transformer parameter study physical checks did not pass"))
    mkpath(dirname(path))
    write(path, transformer_parameter_report_text(study))
    return abspath(path)
end

end
