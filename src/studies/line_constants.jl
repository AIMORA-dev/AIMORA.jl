module LineConstantsStudy

using Printf

using ..DeckParser:
    DeckLineConstantsPhysicalConductor,
    DeckParseResult,
    assert_deck_valid!,
    deck_line_constants_frequency_cards,
    deck_line_constants_physical_conductors
using ..Lines:
    OverheadLineConductor,
    OverheadLineConstants,
    overhead_line_constants

export LineConstantsStudyResult,
       run_line_constants_study,
       write_line_constants_report

struct LineConstantsStudyResult
    source::String
    physical_conductors::Vector{DeckLineConstantsPhysicalConductor}
    frequency_results::Vector{OverheadLineConstants}
    physical_checks_passed::Bool
end

function _overhead_line_conductor(row::DeckLineConstantsPhysicalConductor)
    return OverheadLineConductor(
        row.phase_number,
        row.skin_effect_type,
        row.resistance_ohm_per_mile,
        row.reactance_type,
        row.reactance_or_gmr,
        row.diameter_inches * 0.0254,
        row.horizontal_ft * 0.3048,
        row.average_height_ft * 0.3048;
        name = row.conductor_name,
    )
end

function run_line_constants_study(parsed::DeckParseResult)
    assert_deck_valid!(parsed)
    physical_rows = deck_line_constants_physical_conductors(parsed)
    frequency_cards = deck_line_constants_frequency_cards(parsed)
    isempty(physical_rows) &&
        throw(ArgumentError("deck contains no accepted LINE CONSTANTS conductor rows"))
    isempty(frequency_cards) &&
        throw(ArgumentError("deck contains no accepted LINE CONSTANTS frequency rows"))
    conductors = OverheadLineConductor[
        _overhead_line_conductor(row) for row in physical_rows
    ]
    results = OverheadLineConstants[
        overhead_line_constants(
            conductors,
            card.earth_resistivity_ohm_m,
            card.frequency_hz;
            conductance_mho_per_mile = card.conductance_mho_per_mile,
            earth_return_correction_tolerance = card.carson_correction_factor,
        )
        for card in frequency_cards
    ]
    checks =
        all(result -> result.physical_checks_passed, results) &&
        all(
            index -> results[index].frequency_hz == frequency_cards[index].frequency_hz,
            eachindex(results),
        )
    return LineConstantsStudyResult(
        parsed.source,
        physical_rows,
        results,
        checks,
    )
end

function _write_real_matrix(
    io::IO,
    title::AbstractString,
    matrix::AbstractMatrix{<:Real},
    unit::AbstractString,
)
    println(io)
    println(io, title, " (", unit, ")")
    for row in axes(matrix, 1)
        @printf(io, "%3d", row)
        for column in axes(matrix, 2)
            @printf(io, " %16.8e", matrix[row, column])
        end
        println(io)
    end
    return nothing
end

function _write_complex_matrix(
    io::IO,
    title::AbstractString,
    matrix::AbstractMatrix{<:Complex},
    unit::AbstractString,
)
    println(io)
    println(io, title, " (", unit, ")")
    for row in axes(matrix, 1)
        @printf(io, "%3d", row)
        for column in axes(matrix, 2)
            @printf(
                io,
                " %16.8e%+16.8ej",
                real(matrix[row, column]),
                imag(matrix[row, column]),
            )
        end
        println(io)
    end
    return nothing
end

function _write_conductor_table(io::IO, result::LineConstantsStudyResult)
    println(io)
    println(io, "PHYSICAL CONDUCTORS")
    println(
        io,
        "ROW PHASE BUNDLE R(OHM/MILE) DIAMETER(IN) X(FT) HEIGHT(FT) NAME",
    )
    for (index, row) in enumerate(result.physical_conductors)
        @printf(
            io,
            "%3d %5d %6d %11.6f %12.6f %10.4f %10.4f",
            index,
            row.phase_number,
            row.bundle_ordinal,
            row.resistance_ohm_per_mile,
            row.diameter_inches,
            row.horizontal_ft,
            row.average_height_ft,
        )
        isempty(row.conductor_name) ? println(io) : println(io, ' ', row.conductor_name)
    end
    return nothing
end

function _write_frequency_result(io::IO, result::OverheadLineConstants)
    println(io)
    println(io, "FREQUENCY RESULT")
    @printf(io, "FREQUENCY_HZ %.12g\n", result.frequency_hz)
    @printf(io, "EARTH_RESISTIVITY_OHM_M %.12g\n", result.earth_resistivity_ohm_m)
    @printf(
        io,
        "EARTH_RETURN_CORRECTION_TOLERANCE %.12g\n",
        result.earth_return_correction_tolerance,
    )
    println(
        io,
        "EARTH_RETURN_CORRECTION_APPLIED ",
        result.earth_return_correction_applied,
    )
    @printf(io, "CONDUCTANCE_MHO_PER_MILE %.12g\n", result.conductance_mho_per_mile)
    println(io, "PHASE_COUNT ", result.phase_count)
    println(io, "GROUNDED_CONDUCTOR_COUNT ", result.grounded_conductor_count)

    _write_real_matrix(
        io,
        "CAPACITANCE MATRIX FOR PHYSICAL CONDUCTORS",
        result.physical_capacitance_matrix_f_per_mile,
        "F/MILE",
    )
    _write_real_matrix(
        io,
        "CAPACITANCE MATRIX FOR EQUIVALENT PHASE CONDUCTORS",
        result.equivalent_phase_capacitance_matrix_f_per_mile,
        "F/MILE",
    )
    _write_complex_matrix(
        io,
        "CAPACITANCE MATRIX FOR SEQUENCE COMPONENTS",
        result.sequence_capacitance_matrix_f_per_mile,
        "F/MILE",
    )
    _write_complex_matrix(
        io,
        "IMPEDANCE MATRIX FOR PHYSICAL CONDUCTORS",
        result.physical_impedance_matrix_ohm_per_mile,
        "OHM/MILE",
    )
    _write_complex_matrix(
        io,
        "IMPEDANCE MATRIX FOR EQUIVALENT PHASE CONDUCTORS",
        result.equivalent_phase_impedance_matrix_ohm_per_mile,
        "OHM/MILE",
    )
    _write_complex_matrix(
        io,
        "IMPEDANCE MATRIX FOR SEQUENCE COMPONENTS",
        result.sequence_impedance_matrix_ohm_per_mile,
        "OHM/MILE",
    )
    _write_complex_matrix(
        io,
        "INVERTED IMPEDANCE MATRIX FOR PHYSICAL CONDUCTORS",
        result.physical_inverse_impedance_matrix_mho_mile,
        "MHO-MILE",
    )

    println(io)
    println(
        io,
        "SEQUENCE |ZC|(OHM) ANGLE(DEG) ATTEN(DB/MILE) VELOCITY(MILE/S) " *
        "WAVELENGTH(MILE) R(OHM/MILE) X(OHM/MILE) B(MHO/MILE)",
    )
    for row in result.sequence_constants
        @printf(
            io,
            "%-8s %12.5f %11.5f %15.7g %16.7g %16.7g %13.7g %13.7g %13.7g\n",
            String(row.sequence),
            row.surge_impedance_magnitude_ohm,
            row.surge_impedance_angle_deg,
            row.attenuation_db_per_mile,
            row.velocity_miles_per_s,
            row.wavelength_miles,
            row.resistance_ohm_per_mile,
            row.reactance_ohm_per_mile,
            row.susceptance_mho_per_mile,
        )
    end
    println(io, "PHYSICAL_CHECKS_PASSED ", result.physical_checks_passed)
    return nothing
end

function write_line_constants_report(
    path::AbstractString,
    result::LineConstantsStudyResult,
)
    result.physical_checks_passed ||
        throw(ArgumentError("line-constants study physical checks did not pass"))
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "AIMORA OVERHEAD LINE CONSTANTS")
        println(io, "SOURCE ", result.source)
        _write_conductor_table(io, result)
        for frequency_result in result.frequency_results
            _write_frequency_result(io, frequency_result)
        end
    end
    return abspath(path)
end

end
