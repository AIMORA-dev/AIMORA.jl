export OverheadLineConductor,
       OverheadLineSequenceConstants,
       OverheadLineConstants,
       overhead_line_constants

const LINE_MILE_LENGTH_M = 1609.344
const LINE_FOOT_LENGTH_M = 0.3048
const LINE_INCH_LENGTH_M = 0.0254

"""
One physical overhead conductor after bundle expansion.

`phase_number == 0` denotes a grounded shield conductor. Positive phase
numbers identify conductors constrained to the same phase potential.
"""
struct OverheadLineConductor
    phase_number::Int
    skin_depth_parameter::Float64
    resistance_ohm_per_mile::Float64
    reactance_specification::Int
    reactance_or_gmr::Float64
    diameter_m::Float64
    horizontal_position_m::Float64
    average_height_m::Float64
    name::String

    function OverheadLineConductor(
        phase_number::Integer,
        skin_depth_parameter::Real,
        resistance_ohm_per_mile::Real,
        reactance_specification::Integer,
        reactance_or_gmr::Real,
        diameter_m::Real,
        horizontal_position_m::Real,
        average_height_m::Real;
        name::AbstractString = "",
    )
        phase = Int(phase_number)
        phase >= 0 ||
            throw(ArgumentError("overhead conductor phase_number must be nonnegative"))
        skin = Float64(skin_depth_parameter)
        isfinite(skin) && 0.0 <= skin <= 0.5 ||
            throw(ArgumentError("skin_depth_parameter must be between zero and 0.5"))
        resistance = Float64(resistance_ohm_per_mile)
        isfinite(resistance) && resistance >= 0.0 ||
            throw(ArgumentError("conductor resistance must be finite and nonnegative"))
        specification = Int(reactance_specification)
        0 <= specification <= 4 ||
            throw(ArgumentError("reactance_specification must be between 0 and 4"))
        reactance_value = Float64(reactance_or_gmr)
        isfinite(reactance_value) ||
            throw(ArgumentError("reactance_or_gmr must be finite"))
        diameter = Float64(diameter_m)
        isfinite(diameter) && diameter > 0.0 ||
            throw(ArgumentError("conductor diameter_m must be finite and positive"))
        horizontal = Float64(horizontal_position_m)
        height = Float64(average_height_m)
        isfinite(horizontal) ||
            throw(ArgumentError("conductor horizontal_position_m must be finite"))
        isfinite(height) && height > diameter / 2.0 ||
            throw(ArgumentError("conductor average_height_m must exceed its radius"))
        if specification in (2, 3)
            reactance_value > 0.0 ||
                throw(ArgumentError("reactance specification $specification requires a positive GMR value"))
        end
        return new(
            phase,
            skin,
            resistance,
            specification,
            reactance_value,
            diameter,
            horizontal,
            height,
            String(name),
        )
    end
end

struct OverheadLineSequenceConstants
    sequence::Symbol
    resistance_ohm_per_mile::Float64
    reactance_ohm_per_mile::Float64
    conductance_mho_per_mile::Float64
    susceptance_mho_per_mile::Float64
    surge_impedance_magnitude_ohm::Float64
    surge_impedance_angle_deg::Float64
    attenuation_db_per_mile::Float64
    velocity_miles_per_s::Float64
    wavelength_miles::Float64
end

struct OverheadLineConstants
    frequency_hz::Float64
    angular_frequency_rad_s::Float64
    earth_resistivity_ohm_m::Float64
    earth_return_correction_tolerance::Float64
    earth_return_correction_applied::Bool
    conductance_mho_per_mile::Float64
    conductors::Vector{OverheadLineConductor}
    phase_count::Int
    grounded_conductor_count::Int
    potential_coefficient_matrix_v_m_per_c::Matrix{Float64}
    physical_capacitance_matrix_f_per_mile::Matrix{Float64}
    physical_impedance_matrix_ohm_per_mile::Matrix{ComplexF64}
    physical_inverse_impedance_matrix_mho_mile::Matrix{ComplexF64}
    equivalent_phase_capacitance_matrix_f_per_mile::Matrix{Float64}
    equivalent_phase_impedance_matrix_ohm_per_mile::Matrix{ComplexF64}
    sequence_capacitance_matrix_f_per_mile::Matrix{ComplexF64}
    sequence_impedance_matrix_ohm_per_mile::Matrix{ComplexF64}
    sequence_constants::Vector{OverheadLineSequenceConstants}
    capacitance_symmetry_max_abs_error::Float64
    impedance_symmetry_max_abs_error::Float64
    minimum_capacitance_eigenvalue_f_per_mile::Float64
    minimum_resistance_eigenvalue_ohm_per_mile::Float64
    physical_checks_passed::Bool
end

function _overhead_line_phase_count(conductors)
    phase_count = maximum(getfield.(conductors, :phase_number); init = 0)
    phase_count > 0 ||
        throw(ArgumentError("overhead line requires at least one phase conductor"))
    phases = Set(
        conductor.phase_number for conductor in conductors
        if conductor.phase_number > 0
    )
    phases == Set(1:phase_count) ||
        throw(ArgumentError("positive phase numbers must be contiguous from one"))
    return phase_count
end

function _overhead_line_geometry_checks(conductors)
    for row in eachindex(conductors), column in 1:(row - 1)
        first = conductors[row]
        second = conductors[column]
        hypot(
            first.horizontal_position_m - second.horizontal_position_m,
            first.average_height_m - second.average_height_m,
        ) > 0.0 || throw(ArgumentError("physical overhead conductors cannot share coordinates"))
    end
    return nothing
end

function _overhead_line_potential_coefficients(conductors)
    count = length(conductors)
    scale = inv(2.0 * pi * LINE_VACUUM_PERMITTIVITY_F_PER_M)
    potential = Matrix{Float64}(undef, count, count)
    for column in 1:count, row in 1:column
        first = conductors[row]
        second = conductors[column]
        dx = first.horizontal_position_m - second.horizontal_position_m
        image_distance = hypot(dx, first.average_height_m + second.average_height_m)
        value = if row == column
            radius = first.diameter_m / 2.0
            scale * log(2.0 * first.average_height_m / radius)
        else
            direct_distance = hypot(dx, first.average_height_m - second.average_height_m)
            direct_distance > 0.0 ||
                throw(ArgumentError("physical overhead conductor spacing must be positive"))
            scale * log(image_distance / direct_distance)
        end
        potential[row, column] = value
        potential[column, row] = value
    end
    return potential
end

function _overhead_line_internal_impedance(
    conductor::OverheadLineConductor,
    frequency_hz::Float64,
)
    resistance_per_m = conductor.resistance_ohm_per_mile / LINE_MILE_LENGTH_M
    conductor.skin_depth_parameter > 0.0 || return ComplexF64(resistance_per_m)
    resistance_per_m > 0.0 ||
        throw(ArgumentError("skin-effect calculation requires positive conductor resistance"))
    radius = conductor.diameter_m / 2.0
    inner_radius_ratio = 1.0 - 2.0 * conductor.skin_depth_parameter
    conducting_area = pi * radius^2 * (1.0 - inner_radius_ratio^2)
    resistivity = resistance_per_m * conducting_area
    outer_diffusion_factor =
        radius * sqrt(LINE_VACUUM_PERMEABILITY_H_PER_M / resistivity)
    return cable_skin_effect_internal_impedance(
        inner_radius_ratio * outer_diffusion_factor,
        outer_diffusion_factor,
        1.0,
        frequency_hz,
    )
end

function _overhead_line_effective_radius(
    conductor::OverheadLineConductor,
    frequency_hz::Float64,
)
    specification = conductor.reactance_specification
    if specification == 4
        return conductor.diameter_m / 2.0
    elseif specification == 3
        return conductor.reactance_or_gmr * conductor.diameter_m / 2.0
    elseif specification == 2
        return conductor.reactance_or_gmr * LINE_INCH_LENGTH_M
    end

    omega = 2.0 * pi * frequency_hz
    reactance_scale = omega * 0.00064373888
    specified_reactance = conductor.reactance_or_gmr
    specification == 1 && (specified_reactance *= frequency_hz / 60.0)
    logarithmic_radius =
        specified_reactance / reactance_scale +
        0.5 * log(2.0 * conductor.average_height_m / LINE_FOOT_LENGTH_M)
    radius_inches =
        24.0 * (conductor.average_height_m / LINE_FOOT_LENGTH_M) /
        exp(2.0 * logarithmic_radius)
    radius_inches > 0.0 ||
        throw(ArgumentError("reactance specification produces a nonpositive effective radius"))
    return radius_inches * LINE_INCH_LENGTH_M
end

function _overhead_line_physical_impedance(
    conductors,
    frequency_hz::Float64,
    earth_resistivity_ohm_m::Float64,
    earth_return_correction_applied::Bool,
)
    count = length(conductors)
    omega = 2.0 * pi * frequency_hz
    air_scale =
        omega * LINE_VACUUM_PERMEABILITY_H_PER_M / (2.0 * pi)
    impedance = Matrix{ComplexF64}(undef, count, count)
    for column in 1:count, row in 1:column
        first = conductors[row]
        second = conductors[column]
        dx = abs(first.horizontal_position_m - second.horizontal_position_m)
        image_distance = hypot(dx, first.average_height_m + second.average_height_m)
        value = if row == column
            effective_radius = _overhead_line_effective_radius(first, frequency_hz)
            internal = first.reactance_specification == 4 ?
                _overhead_line_internal_impedance(first, frequency_hz) :
                ComplexF64(first.resistance_ohm_per_mile / LINE_MILE_LENGTH_M)
            internal +
            ComplexF64(0.0, air_scale * log(2.0 * first.average_height_m / effective_radius)) +
            (earth_return_correction_applied ? cable_homogeneous_earth_return_impedance(
                image_distance,
                0.0,
                earth_resistivity_ohm_m,
                frequency_hz,
            ) : 0.0im)
        else
            direct_distance = hypot(
                dx,
                first.average_height_m - second.average_height_m,
            )
            image_angle = asin(clamp(dx / image_distance, 0.0, 1.0))
            ComplexF64(0.0, air_scale * log(image_distance / direct_distance)) +
            (earth_return_correction_applied ? cable_homogeneous_earth_return_impedance(
                image_distance,
                image_angle,
                earth_resistivity_ohm_m,
                frequency_hz,
            ) : 0.0im)
        end
        impedance[row, column] = value
        impedance[column, row] = value
    end
    return impedance .* LINE_MILE_LENGTH_M
end

function _overhead_line_phase_matrices(
    conductors,
    physical_capacitance,
    physical_impedance,
    phase_count,
)
    active_indices = findall(conductor -> conductor.phase_number > 0, conductors)
    grounded_indices = findall(conductor -> conductor.phase_number == 0, conductors)
    active_impedance = if isempty(grounded_indices)
        physical_impedance[active_indices, active_indices]
    else
        zaa = physical_impedance[active_indices, active_indices]
        zag = physical_impedance[active_indices, grounded_indices]
        zga = physical_impedance[grounded_indices, active_indices]
        zgg = physical_impedance[grounded_indices, grounded_indices]
        zaa - zag * (zgg \ zga)
    end
    conductor_to_phase = zeros(Float64, length(active_indices), phase_count)
    for (row, conductor_index) in enumerate(active_indices)
        conductor_to_phase[row, conductors[conductor_index].phase_number] = 1.0
    end
    transpose_map = transpose(conductor_to_phase)
    phase_impedance = inv(transpose_map * (active_impedance \ conductor_to_phase))
    phase_capacitance =
        transpose_map *
        physical_capacitance[active_indices, active_indices] *
        conductor_to_phase
    return Matrix{Float64}(phase_capacitance), Matrix{ComplexF64}(phase_impedance)
end

function _overhead_line_sequence_transform(phase_count::Int)
    if phase_count == 1
        return ones(ComplexF64, 1, 1)
    elseif phase_count == 2
        return ComplexF64[1.0 1.0; 1.0 -1.0]
    end
    width = 3 * div(phase_count, 3)
    width > 0 ||
        throw(ArgumentError("sequence transformation requires one, two, or at least three phases"))
    transform = zeros(ComplexF64, width, width)
    rotation = cis(2.0 * pi / 3.0)
    block = ComplexF64[
        1.0 1.0 1.0
        1.0 rotation^2 rotation
        1.0 rotation rotation^2
    ]
    for first in 1:3:width
        transform[first:(first + 2), first:(first + 2)] .= block
    end
    return transform
end

function _overhead_line_sequence_matrices(phase_capacitance, phase_impedance)
    transform = _overhead_line_sequence_transform(size(phase_impedance, 1))
    width = size(transform, 1)
    inverse_transform = inv(transform)
    sequence_impedance =
        inverse_transform * phase_impedance[1:width, 1:width] * transform
    sequence_capacitance =
        inverse_transform * phase_capacitance[1:width, 1:width] * transform
    return Matrix{ComplexF64}(sequence_capacitance),
        Matrix{ComplexF64}(sequence_impedance)
end

function _overhead_line_sequence_row(
    sequence::Symbol,
    impedance::ComplexF64,
    capacitance_f_per_mile::Float64,
    conductance_mho_per_mile::Float64,
    frequency_hz::Float64,
)
    omega = 2.0 * pi * frequency_hz
    capacitance_f_per_mile > 0.0 ||
        throw(ArgumentError("$sequence-sequence capacitance must be positive"))
    point = frequency_dependent_line_point(
        real(impedance),
        imag(impedance) / omega,
        conductance_mho_per_mile,
        capacitance_f_per_mile,
        1.0,
        frequency_hz,
    )
    propagation = point.propagation_constant
    phase_constant = imag(propagation)
    phase_constant > 0.0 ||
        throw(ArgumentError("$sequence-sequence phase constant must be positive"))
    wavelength = 2.0 * pi / phase_constant
    return OverheadLineSequenceConstants(
        sequence,
        real(impedance),
        imag(impedance),
        conductance_mho_per_mile,
        omega * capacitance_f_per_mile,
        abs(point.characteristic_impedance),
        rad2deg(angle(point.characteristic_impedance)),
        20.0 / log(10.0) * real(propagation),
        frequency_hz * wavelength,
        wavelength,
    )
end

function _overhead_line_sequence_constants(
    sequence_capacitance,
    sequence_impedance,
    conductance_mho_per_mile,
    frequency_hz,
)
    zero_impedance = sequence_impedance[1, 1]
    zero_capacitance = sequence_capacitance[1, 1]
    rows = OverheadLineSequenceConstants[
        _overhead_line_sequence_row(
            :zero,
            ComplexF64(zero_impedance),
            Float64(real(zero_capacitance)),
            conductance_mho_per_mile,
            frequency_hz,
        ),
    ]
    if size(sequence_impedance, 1) > 1
        positive_impedance = sequence_impedance[2, 2]
        positive_capacitance = sequence_capacitance[2, 2]
        push!(
            rows,
            _overhead_line_sequence_row(
                :positive,
                ComplexF64(positive_impedance),
                Float64(real(positive_capacitance)),
                conductance_mho_per_mile,
                frequency_hz,
            ),
        )
    end
    return rows
end

function overhead_line_constants(
    conductors::AbstractVector{OverheadLineConductor},
    earth_resistivity_ohm_m::Real,
    frequency_hz::Real;
    conductance_mho_per_mile::Real = 3.22e-9,
    earth_return_correction_tolerance::Real = 1.0e-6,
)
    isempty(conductors) &&
        throw(ArgumentError("overhead line requires physical conductors"))
    conductor_rows = collect(conductors)
    _overhead_line_geometry_checks(conductor_rows)
    phase_count = _overhead_line_phase_count(conductor_rows)
    earth_resistivity = Float64(earth_resistivity_ohm_m)
    isfinite(earth_resistivity) && earth_resistivity > 0.0 ||
        throw(ArgumentError("earth_resistivity_ohm_m must be finite and positive"))
    correction_tolerance = Float64(earth_return_correction_tolerance)
    isfinite(correction_tolerance) && correction_tolerance >= 0.0 ||
        throw(ArgumentError(
            "earth_return_correction_tolerance must be finite and nonnegative",
        ))
    correction_applied = correction_tolerance > 0.0
    frequency = Float64(frequency_hz)
    isfinite(frequency) && frequency > 0.0 ||
        throw(ArgumentError("frequency_hz must be finite and positive"))
    conductance = Float64(conductance_mho_per_mile)
    isfinite(conductance) && conductance >= 0.0 ||
        throw(ArgumentError("conductance_mho_per_mile must be finite and nonnegative"))

    potential = _overhead_line_potential_coefficients(conductor_rows)
    physical_capacitance = inv(potential) .* LINE_MILE_LENGTH_M
    physical_impedance = _overhead_line_physical_impedance(
        conductor_rows,
        frequency,
        earth_resistivity,
        correction_applied,
    )
    phase_capacitance, phase_impedance = _overhead_line_phase_matrices(
        conductor_rows,
        physical_capacitance,
        physical_impedance,
        phase_count,
    )
    sequence_capacitance, sequence_impedance =
        _overhead_line_sequence_matrices(phase_capacitance, phase_impedance)
    sequence_rows = _overhead_line_sequence_constants(
        sequence_capacitance,
        sequence_impedance,
        0.0,
        frequency,
    )

    capacitance_symmetry_error =
        maximum(abs, phase_capacitance - transpose(phase_capacitance))
    impedance_symmetry_error =
        maximum(abs, phase_impedance - transpose(phase_impedance))
    minimum_capacitance =
        eigmin(Symmetric((phase_capacitance + transpose(phase_capacitance)) / 2.0))
    resistance_matrix = real.(phase_impedance)
    minimum_resistance =
        eigmin(Symmetric((resistance_matrix + transpose(resistance_matrix)) / 2.0))
    finite_matrices =
        all(isfinite, potential) &&
        all(isfinite, physical_capacitance) &&
        all(value -> isfinite(real(value)) && isfinite(imag(value)), physical_impedance) &&
        all(isfinite, phase_capacitance) &&
        all(value -> isfinite(real(value)) && isfinite(imag(value)), phase_impedance)
    scale_capacitance = maximum(abs, phase_capacitance)
    scale_impedance = maximum(abs, phase_impedance)
    physical_checks =
        finite_matrices &&
        capacitance_symmetry_error <= max(1.0e-15, 1.0e-12 * scale_capacitance) &&
        impedance_symmetry_error <= max(1.0e-12, 1.0e-12 * scale_impedance) &&
        minimum_capacitance >= -max(1.0e-15, 1.0e-12 * scale_capacitance) &&
        minimum_resistance >= -max(1.0e-12, 1.0e-12 * scale_impedance) &&
        all(row -> row.susceptance_mho_per_mile > 0.0, sequence_rows) &&
        all(row -> row.velocity_miles_per_s > 0.0, sequence_rows)

    return OverheadLineConstants(
        frequency,
        2.0 * pi * frequency,
        earth_resistivity,
        correction_tolerance,
        correction_applied,
        conductance,
        conductor_rows,
        phase_count,
        count(conductor -> conductor.phase_number == 0, conductor_rows),
        Matrix{Float64}(potential),
        Matrix{Float64}(physical_capacitance),
        Matrix{ComplexF64}(physical_impedance),
        Matrix{ComplexF64}(inv(physical_impedance)),
        phase_capacitance,
        phase_impedance,
        sequence_capacitance,
        sequence_impedance,
        sequence_rows,
        capacitance_symmetry_error,
        impedance_symmetry_error,
        minimum_capacitance,
        minimum_resistance,
        physical_checks,
    )
end
