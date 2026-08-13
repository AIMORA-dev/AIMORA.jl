using SHA
using TOML

using ..StudyCore: ParameterProvenance, PhysicalModelParameter

export LineParameterSource,
       LineSoilLayer,
       LineSoilProfile,
       LineLayeredEarthReturnPoint,
       LineParameterSegment,
       LineParameterDiagnostics,
       WidebandLineParameterSet,
       LineParameterUncertaintyEnvelope,
       generic_line_parameter_source,
       line_layered_earth_return_impedance,
       line_layered_earth_return_impedance_matrix,
       overhead_line_parameter_segment,
       cable_line_parameter_segment,
       line_parameter_set,
       line_parameter_uncertainty_envelope,
       write_line_parameter_set,
       read_line_parameter_set,
       line_parameter_report_text,
       write_line_parameter_report

const LINE_PARAMETER_SCHEMA_VERSION = 1
const LINE_PARAMETER_MINIMUM_FREQUENCY_HZ = 0.1
const LINE_PARAMETER_MAXIMUM_FREQUENCY_HZ = 1.0e6
const LINE_PARAMETER_MAXIMUM_FREQUENCY_COUNT = 401
const LINE_PARAMETER_MAXIMUM_ROUTE_SEGMENTS = 32
const LINE_PARAMETER_MAXIMUM_ROUTE_LENGTH_M = 2.0e6
const LINE_PARAMETER_LOSS_EIGENVALUE_FLOOR = -1.0e-12

struct LineParameterSource
    provenance::ParameterProvenance
    rights::String
    content_sha256::String
    data_class::Symbol

    function LineParameterSource(
        provenance::ParameterProvenance,
        rights::AbstractString,
        content_sha256::AbstractString,
        data_class::Symbol,
    )
        provenance.nature === PhysicalModelParameter || throw(ArgumentError(
            "line parameter source provenance must describe a physical model parameter",
        ))
        rights_value = String(rights)
        isempty(strip(rights_value)) && throw(ArgumentError(
            "line parameter source rights must not be empty",
        ))
        content_hash = lowercase(String(content_sha256))
        occursin(r"^[0-9a-f]{64}$", content_hash) || throw(ArgumentError(
            "line parameter source content_sha256 must be lowercase 64-hex",
        ))
        data_class in (:generic, :measured, :fitted, :certified, :imported) ||
            throw(ArgumentError("unsupported line parameter source data class"))
        return new(provenance, rights_value, content_hash, data_class)
    end
end

function generic_line_parameter_source(owner::AbstractString; units::AbstractString)
    owner_value = strip(String(owner))
    isempty(owner_value) && throw(ArgumentError("line parameter source owner must not be empty"))
    identity = bytes2hex(sha256(codeunits("AIMORA generic line parameter|" * owner_value * "|" * String(units))))
    return LineParameterSource(
        ParameterProvenance(
            "AIMORA-authored synthetic generic " * owner_value,
            units,
            "direct SI values with explicit phase and terminal orientation",
            "illustrative bounded input or explicit unknown physical uncertainty",
            "REQ-LINE-PARAMETERS-001 frozen generic validity domain",
            PhysicalModelParameter,
        ),
        "Redistributable under the controlling AIMORA repository licence",
        identity,
        :generic,
    )
end

struct LineSoilLayer
    resistivity_ohm_m::Float64
    relative_permittivity::Float64
    relative_permeability::Float64
    thickness_m::Union{Nothing,Float64}
    source::LineParameterSource

    function LineSoilLayer(
        resistivity_ohm_m::Real,
        relative_permittivity::Real,
        relative_permeability::Real,
        thickness_m::Union{Nothing,Real},
        source::LineParameterSource,
    )
        resistivity = _checked_line_positive(resistivity_ohm_m, "soil resistivity_ohm_m")
        permittivity = _checked_line_positive(relative_permittivity, "soil relative_permittivity")
        permeability = _checked_line_positive(relative_permeability, "soil relative_permeability")
        thickness = thickness_m === nothing ? nothing :
            _checked_line_positive(thickness_m, "soil thickness_m")
        return new(resistivity, permittivity, permeability, thickness, source)
    end
end

struct LineSoilProfile
    layers::Vector{LineSoilLayer}
    profile_id::Symbol
    deterministic_signature_sha256::String

    function LineSoilProfile(layers::AbstractVector{LineSoilLayer}; profile_id::Symbol=:soil)
        rows = collect(layers)
        1 <= length(rows) <= 4 || throw(ArgumentError(
            "line soil profile requires one through four layers",
        ))
        for index in eachindex(rows)
            expected_half_space = index == lastindex(rows)
            (rows[index].thickness_m === nothing) == expected_half_space || throw(ArgumentError(
                expected_half_space ?
                "the terminal soil layer must be a half-space" :
                "every nonterminal soil layer requires positive finite thickness",
            ))
        end
        isempty(String(profile_id)) && throw(ArgumentError("soil profile_id must not be empty"))
        io = IOBuffer()
        println(io, profile_id)
        for layer in rows
            println(
                io,
                layer.resistivity_ohm_m,
                '|',
                layer.relative_permittivity,
                '|',
                layer.relative_permeability,
                '|',
                layer.thickness_m,
                '|',
                layer.source.content_sha256,
            )
        end
        return new(rows, profile_id, bytes2hex(sha256(take!(io))))
    end
end

function LineSoilProfile(
    earth_resistivity_ohm_m::Real;
    profile_id::Symbol=:homogeneous_soil,
    source::LineParameterSource=generic_line_parameter_source(
        "homogeneous soil";
        units="ohm metre, dimensionless relative permittivity and permeability",
    ),
)
    return LineSoilProfile(
        [LineSoilLayer(earth_resistivity_ohm_m, 10.0, 1.0, nothing, source)];
        profile_id,
    )
end

struct LineLayeredEarthReturnPoint
    impedance_ohm_per_m::ComplexF64
    quadrature_error_ohm_per_m::Float64
    frequency_hz::Float64
    image_distance_m::Float64
    image_angle_rad::Float64
    soil_signature_sha256::String
    coarse_order::Int
    refined_order::Int
end

function _line_layer_vertical_wavenumber(
    layer::LineSoilLayer,
    spectral_wavenumber_per_m::Float64,
    angular_frequency_rad_s::Float64,
)
    conductivity = inv(layer.resistivity_ohm_m)
    permeability = LINE_VACUUM_PERMEABILITY_H_PER_M * layer.relative_permeability
    permittivity = LINE_VACUUM_PERMITTIVITY_F_PER_M * layer.relative_permittivity
    value = sqrt(ComplexF64(
        spectral_wavenumber_per_m^2 -
        angular_frequency_rad_s^2 * permeability * permittivity,
        angular_frequency_rad_s * permeability * conductivity,
    ))
    real(value) < 0.0 && (value = -value)
    return value, permeability
end

function _line_layered_surface_admittance(
    profile::LineSoilProfile,
    spectral_wavenumber_per_m::Float64,
    angular_frequency_rad_s::Float64,
)
    bottom = last(profile.layers)
    vertical, permeability = _line_layer_vertical_wavenumber(
        bottom,
        spectral_wavenumber_per_m,
        angular_frequency_rad_s,
    )
    input_admittance = vertical / permeability
    for index in (length(profile.layers) - 1):-1:1
        layer = profile.layers[index]
        vertical, permeability = _line_layer_vertical_wavenumber(
            layer,
            spectral_wavenumber_per_m,
            angular_frequency_rad_s,
        )
        characteristic_admittance = vertical / permeability
        propagation = tanh(vertical * something(layer.thickness_m))
        denominator =
            characteristic_admittance + input_admittance * propagation
        abs(denominator) > eps(Float64) || throw(ArgumentError(
            "layered-soil surface-admittance recursion is singular",
        ))
        input_admittance = characteristic_admittance * (
            input_admittance + characteristic_admittance * propagation
        ) / denominator
    end
    return input_admittance
end

function _line_adaptive_simpson_complex(
    integrand,
    left::Float64,
    right::Float64,
    relative_tolerance::Float64;
    maximum_depth::Int=24,
)
    value_scale = Ref(1.0)
    midpoint = (left + right) / 2.0
    left_value = integrand(left)
    midpoint_value = integrand(midpoint)
    right_value = integrand(right)
    value_scale[] = max(
        abs(left_value),
        abs(midpoint_value),
        abs(right_value),
        1.0,
    )
    whole = (right - left) * (left_value + 4.0 * midpoint_value + right_value) / 6.0
    tolerance = relative_tolerance * max(abs(whole), 1.0)
    function refine(a, b, fa, fm, fb, estimate, allowed_error, depth)
        center = (a + b) / 2.0
        left_center = (a + center) / 2.0
        right_center = (center + b) / 2.0
        left_center_value = integrand(left_center)
        right_center_value = integrand(right_center)
        value_scale[] = max(
            value_scale[],
            abs(left_center_value),
            abs(right_center_value),
        )
        left_estimate = (center - a) * (
            fa + 4.0 * left_center_value + fm
        ) / 6.0
        right_estimate = (b - center) * (
            fm + 4.0 * right_center_value + fb
        ) / 6.0
        combined = left_estimate + right_estimate
        correction = combined - estimate
        if depth == 0 || abs(correction) <= 15.0 * allowed_error
            return combined + correction / 15.0, abs(correction) / 15.0
        end
        left_result, left_error = refine(
            a,
            center,
            fa,
            left_center_value,
            fm,
            left_estimate,
            allowed_error / 2.0,
            depth - 1,
        )
        right_result, right_error = refine(
            center,
            b,
            fm,
            right_center_value,
            fb,
            right_estimate,
            allowed_error / 2.0,
            depth - 1,
        )
        return left_result + right_result, left_error + right_error
    end
    result, error = refine(
        left,
        right,
        left_value,
        midpoint_value,
        right_value,
        whole,
        tolerance,
        maximum_depth,
    )
    return result, error, value_scale[]
end

function _line_layered_earth_integral(
    profile::LineSoilProfile,
    angular_frequency_rad_s::Float64,
    image_vertical_distance_m::Float64,
    lateral_distance_m::Float64,
    relative_tolerance::Float64,
)
    function integrand(scaled_spectral_wavenumber)
        spectral = scaled_spectral_wavenumber / image_vertical_distance_m
        surface_admittance = _line_layered_surface_admittance(
            profile,
            spectral,
            angular_frequency_rad_s,
        )
        effective_vertical = LINE_VACUUM_PERMEABILITY_H_PER_M * surface_admittance
        denominator = spectral + effective_vertical
        abs(denominator) > eps(Float64) || throw(ArgumentError(
            "layered-soil earth-return integrand is singular",
        ))
        return exp(-scaled_spectral_wavenumber) *
            cos(lateral_distance_m * spectral) /
            denominator /
            image_vertical_distance_m
    end
    integral, integration_error, integrand_scale = _line_adaptive_simpson_complex(
        integrand,
        0.0,
        48.0,
        relative_tolerance,
    )
    scale = ComplexF64(
        0.0,
        angular_frequency_rad_s * LINE_VACUUM_PERMEABILITY_H_PER_M / pi,
    )
    return scale * integral,
        abs(scale) * integration_error,
        abs(scale) * integrand_scale
end

function line_layered_earth_return_impedance(
    image_distance_m::Real,
    image_angle_rad::Real,
    profile::LineSoilProfile,
    frequency_hz::Real;
    quadrature_relative_tolerance::Real=1.0e-5,
)
    image_distance = _checked_line_positive(
        image_distance_m,
        "layered-earth image_distance_m",
    )
    angle = Float64(image_angle_rad)
    isfinite(angle) && 0.0 <= angle < pi / 2.0 || throw(ArgumentError(
        "layered-earth image_angle_rad must be finite in [0, pi/2)",
    ))
    frequency = _checked_line_positive(frequency_hz, "layered-earth frequency_hz")
    LINE_PARAMETER_MINIMUM_FREQUENCY_HZ <= frequency <= LINE_PARAMETER_MAXIMUM_FREQUENCY_HZ ||
        throw(ArgumentError("layered-earth frequency is outside the released domain"))
    tolerance = Float64(quadrature_relative_tolerance)
    isfinite(tolerance) && tolerance > 0.0 || throw(ArgumentError(
        "layered-earth quadrature tolerance must be finite and positive",
    ))
    vertical = image_distance * cos(angle)
    lateral = image_distance * sin(angle)
    vertical > 0.0 || throw(ArgumentError(
        "layered-earth image vertical distance must be positive",
    ))
    angular_frequency = 2.0 * pi * frequency
    if length(profile.layers) == 1
        homogeneous = cable_homogeneous_earth_return_impedance(
            image_distance,
            angle,
            only(profile.layers).resistivity_ohm_m,
            frequency,
        )
        return LineLayeredEarthReturnPoint(
            homogeneous,
            0.0,
            frequency,
            image_distance,
            angle,
            profile.deterministic_signature_sha256,
            0,
            0,
        )
    end
    refined, error, scale = _line_layered_earth_integral(
        profile,
        angular_frequency,
        vertical,
        lateral,
        tolerance,
    )
    error <= tolerance * max(abs(refined), scale, 1.0e-12) || throw(ArgumentError(
        "layered-soil earth-return quadrature did not converge: error=$(error)",
    ))
    isfinite(real(refined)) && isfinite(imag(refined)) || throw(ArgumentError(
        "layered-soil earth-return impedance must be finite",
    ))
    real(refined) >= LINE_PARAMETER_LOSS_EIGENVALUE_FLOOR || throw(ArgumentError(
        "layered-soil earth-return impedance has negative passive loss",
    ))
    return LineLayeredEarthReturnPoint(
        refined,
        error,
        frequency,
        image_distance,
        angle,
        profile.deterministic_signature_sha256,
        0,
        24,
    )
end

function line_layered_earth_return_impedance_matrix(
    constants::CableGeometryConstants,
    profile::LineSoilProfile,
    frequency_hz::Real;
    quadrature_relative_tolerance::Real=1.0e-5,
)
    matrix = Matrix{ComplexF64}(undef, constants.conductor_count, constants.conductor_count)
    maximum_error = 0.0
    for column in 1:constants.conductor_count, row in 1:column
        point = line_layered_earth_return_impedance(
            constants.image_distance_m[row, column],
            constants.angle_rad[row, column],
            profile,
            frequency_hz;
            quadrature_relative_tolerance,
        )
        matrix[row, column] = point.impedance_ohm_per_m
        matrix[column, row] = point.impedance_ohm_per_m
        maximum_error = max(maximum_error, point.quadrature_error_ohm_per_m)
    end
    return (matrix=matrix, quadrature_error_ohm_per_m=maximum_error)
end

struct LineParameterSegment
    id::Symbol
    kind::Symbol
    length_m::Float64
    phase_order::Vector{Symbol}
    frequencies_hz::Vector{Float64}
    series_impedance_matrices_ohm_per_m::Vector{Matrix{ComplexF64}}
    shunt_admittance_matrices_s_per_m::Vector{Matrix{ComplexF64}}
    source::LineParameterSource
    soil_signature_sha256::String
    quadrature_error_ohm_per_m::Float64
    input_signature_sha256::String
end

function _line_parameter_frequencies(values)
    frequencies = Float64.(values)
    1 <= length(frequencies) <= LINE_PARAMETER_MAXIMUM_FREQUENCY_COUNT ||
        throw(ArgumentError("line parameter frequency count is outside the released domain"))
    all(value -> isfinite(value) &&
        LINE_PARAMETER_MINIMUM_FREQUENCY_HZ <= value <= LINE_PARAMETER_MAXIMUM_FREQUENCY_HZ,
        frequencies) || throw(ArgumentError(
        "line parameter frequencies must be finite inside the released domain",
    ))
    issorted(frequencies) && all(diff(frequencies) .> 0.0) || throw(ArgumentError(
        "line parameter frequencies must be strictly increasing and unique",
    ))
    return frequencies
end

function _line_parameter_matrix(
    matrix,
    phase_count::Int,
    label::AbstractString,
)
    size(matrix) == (phase_count, phase_count) || throw(ArgumentError(
        "$label dimensions must match phase order",
    ))
    result = Matrix{ComplexF64}(matrix)
    all(value -> isfinite(real(value)) && isfinite(imag(value)), result) ||
        throw(ArgumentError("$label entries must be finite"))
    symmetry_error = maximum(abs, result - transpose(result))
    symmetry_error <= 1.0e-10 * max(maximum(abs, result), 1.0) || throw(ArgumentError(
        "$label must be complex symmetric in the declared reciprocal domain",
    ))
    return result
end

function _line_parameter_segment_signature(
    id,
    kind,
    length_m,
    phase_order,
    frequencies,
    series,
    shunt,
    source,
    soil_signature,
)
    io = IOBuffer()
    println(io, id, '|', kind, '|', length_m, '|', join(phase_order, ','), '|', source.content_sha256)
    println(io, soil_signature)
    for index in eachindex(frequencies)
        println(io, frequencies[index])
        for matrix in (series[index], shunt[index]), value in matrix
            println(io, bitstring(real(value)), '|', bitstring(imag(value)))
        end
    end
    return bytes2hex(sha256(take!(io)))
end

function LineParameterSegment(
    id::Symbol,
    kind::Symbol,
    length_m::Real,
    phase_order::AbstractVector{Symbol},
    frequencies_hz,
    series_impedance_matrices,
    shunt_admittance_matrices,
    source::LineParameterSource;
    soil_signature_sha256::AbstractString=repeat("0", 64),
    quadrature_error_ohm_per_m::Real=0.0,
)
    isempty(String(id)) && throw(ArgumentError("line parameter segment id must not be empty"))
    kind in (:overhead, :cable, :mixed, :imported) ||
        throw(ArgumentError("unsupported line parameter segment kind"))
    length_value = _checked_line_positive(length_m, "line parameter segment length_m")
    phases = collect(phase_order)
    isempty(phases) && throw(ArgumentError("line parameter segment requires phase order"))
    length(unique(phases)) == length(phases) || throw(ArgumentError(
        "line parameter segment phase identities must be unique",
    ))
    frequencies = _line_parameter_frequencies(frequencies_hz)
    length(series_impedance_matrices) == length(frequencies) ==
        length(shunt_admittance_matrices) || throw(ArgumentError(
        "line parameter segment matrix rows must match frequency count",
    ))
    series = [
        _line_parameter_matrix(matrix, length(phases), "series impedance matrix")
        for matrix in series_impedance_matrices
    ]
    shunt = [
        _line_parameter_matrix(matrix, length(phases), "shunt admittance matrix")
        for matrix in shunt_admittance_matrices
    ]
    soil_signature = lowercase(String(soil_signature_sha256))
    occursin(r"^[0-9a-f]{64}$", soil_signature) || throw(ArgumentError(
        "line parameter soil signature must be lowercase 64-hex",
    ))
    quadrature_error = Float64(quadrature_error_ohm_per_m)
    isfinite(quadrature_error) && quadrature_error >= 0.0 || throw(ArgumentError(
        "line parameter quadrature error must be finite and nonnegative",
    ))
    signature = _line_parameter_segment_signature(
        id,
        kind,
        length_value,
        phases,
        frequencies,
        series,
        shunt,
        source,
        soil_signature,
    )
    return LineParameterSegment(
        id,
        kind,
        length_value,
        phases,
        frequencies,
        series,
        shunt,
        source,
        soil_signature,
        quadrature_error,
        signature,
    )
end

function _line_default_phase_order(count::Int)
    return [Symbol("phase_", index) for index in 1:count]
end

function _overhead_line_layered_phase_matrices(
    conductors,
    profile::LineSoilProfile,
    frequency_hz::Float64,
    conductance_s_per_m::Float64,
    quadrature_relative_tolerance::Float64,
)
    _overhead_line_geometry_checks(conductors)
    phase_count = _overhead_line_phase_count(conductors)
    potential = _overhead_line_potential_coefficients(conductors)
    physical_capacitance = inv(potential)
    omega = 2.0 * pi * frequency_hz
    air_scale = omega * LINE_VACUUM_PERMEABILITY_H_PER_M / (2.0 * pi)
    physical_impedance = Matrix{ComplexF64}(undef, length(conductors), length(conductors))
    maximum_error = 0.0
    for column in eachindex(conductors), row in 1:column
        first = conductors[row]
        second = conductors[column]
        dx = abs(first.horizontal_position_m - second.horizontal_position_m)
        image_distance = hypot(dx, first.average_height_m + second.average_height_m)
        earth = line_layered_earth_return_impedance(
            image_distance,
            asin(clamp(dx / image_distance, 0.0, 1.0)),
            profile,
            frequency_hz;
            quadrature_relative_tolerance,
        )
        maximum_error = max(maximum_error, earth.quadrature_error_ohm_per_m)
        value = if row == column
            effective_radius = _overhead_line_effective_radius(first, frequency_hz)
            internal = first.reactance_specification == 4 ?
                _overhead_line_internal_impedance(first, frequency_hz) :
                ComplexF64(first.resistance_ohm_per_mile / LINE_MILE_LENGTH_M)
            internal +
            ComplexF64(0.0, air_scale * log(2.0 * first.average_height_m / effective_radius)) +
            earth.impedance_ohm_per_m
        else
            direct_distance = hypot(
                dx,
                first.average_height_m - second.average_height_m,
            )
            ComplexF64(0.0, air_scale * log(image_distance / direct_distance)) +
            earth.impedance_ohm_per_m
        end
        physical_impedance[row, column] = value
        physical_impedance[column, row] = value
    end
    phase_capacitance, phase_impedance = _overhead_line_phase_matrices(
        conductors,
        physical_capacitance,
        physical_impedance,
        phase_count,
    )
    phase_admittance = ComplexF64.(
        conductance_s_per_m .* Matrix{Float64}(I, phase_count, phase_count),
        omega .* phase_capacitance,
    )
    return phase_impedance, phase_admittance, maximum_error
end

function overhead_line_parameter_segment(
    id::Symbol,
    conductors::AbstractVector{OverheadLineConductor},
    profile::LineSoilProfile,
    frequencies_hz;
    length_m::Real,
    phase_order::AbstractVector{Symbol}=_line_default_phase_order(
        _overhead_line_phase_count(conductors),
    ),
    conductance_s_per_m::Real=0.0,
    source::LineParameterSource=generic_line_parameter_source(
        "overhead line geometry and material";
        units="metre, ohm per metre, siemens per metre, hertz",
    ),
    quadrature_relative_tolerance::Real=1.0e-5,
)
    frequencies = _line_parameter_frequencies(frequencies_hz)
    conductance = Float64(conductance_s_per_m)
    isfinite(conductance) && conductance >= 0.0 || throw(ArgumentError(
        "overhead line conductance_s_per_m must be finite and nonnegative",
    ))
    tolerance = Float64(quadrature_relative_tolerance)
    series = Matrix{ComplexF64}[]
    shunt = Matrix{ComplexF64}[]
    maximum_error = 0.0
    for frequency in frequencies
        impedance, admittance, error = _overhead_line_layered_phase_matrices(
            collect(conductors),
            profile,
            frequency,
            conductance,
            tolerance,
        )
        push!(series, impedance)
        push!(shunt, admittance)
        maximum_error = max(maximum_error, error)
    end
    return LineParameterSegment(
        id,
        :overhead,
        length_m,
        phase_order,
        frequencies,
        series,
        shunt,
        source;
        soil_signature_sha256=profile.deterministic_signature_sha256,
        quadrature_error_ohm_per_m=maximum_error,
    )
end

function _cable_layered_phase_series(
    constants::CableGeometryConstants,
    profile::LineSoilProfile,
    frequency_hz::Float64;
    include_external_inductance::Bool,
    include_internal_inductance::Bool,
    include_bounded_skin_effect::Bool,
    skin_effect_diffusion_factors,
    reduce_grounded_conductors::Bool,
    quadrature_relative_tolerance::Float64,
)
    base = cable_phase_series_impedance(
        constants,
        frequency_hz;
        include_external_inductance,
        include_internal_inductance,
        include_bounded_skin_effect,
        skin_effect_diffusion_factors,
        earth_resistivity_ohm_m=nothing,
        reduce_grounded_conductors,
    )
    earth = line_layered_earth_return_impedance_matrix(
        constants,
        profile,
        frequency_hz;
        quadrature_relative_tolerance,
    )
    conductor_series =
        base.conductor_series_impedance_matrix_ohm_per_m + earth.matrix
    active_series, _ = _cable_grounded_reduced_series_impedance(
        conductor_series,
        constants,
        reduce_grounded_conductors,
    )
    active_count = _cable_active_conductor_count(constants)
    active_average = _cable_phase_average_matrix(constants)[:, 1:active_count]
    phase_series = active_average * active_series * transpose(active_average)
    return Matrix{ComplexF64}(phase_series), earth.quadrature_error_ohm_per_m
end

function cable_line_parameter_segment(
    id::Symbol,
    constants::CableGeometryConstants,
    profile::LineSoilProfile,
    frequencies_hz;
    length_m::Real,
    phase_order::AbstractVector{Symbol}=_line_default_phase_order(constants.phase_count),
    relative_permittivity=nothing,
    include_external_inductance::Bool=true,
    include_internal_inductance::Bool=true,
    include_bounded_skin_effect::Bool=true,
    skin_effect_diffusion_factors=nothing,
    reduce_grounded_conductors::Bool=true,
    source::LineParameterSource=generic_line_parameter_source(
        "cable geometry and material";
        units="metre, ohm metre, henry per metre, farad per metre, hertz",
    ),
    quadrature_relative_tolerance::Real=1.0e-5,
)
    frequencies = _line_parameter_frequencies(frequencies_hz)
    tolerance = Float64(quadrature_relative_tolerance)
    series = Matrix{ComplexF64}[]
    shunt = Matrix{ComplexF64}[]
    maximum_error = 0.0
    for frequency in frequencies
        impedance, error = _cable_layered_phase_series(
            constants,
            profile,
            frequency;
            include_external_inductance,
            include_internal_inductance,
            include_bounded_skin_effect,
            skin_effect_diffusion_factors,
            reduce_grounded_conductors,
            quadrature_relative_tolerance=tolerance,
        )
        admittance = cable_phase_electrostatic_admittance(
            constants,
            frequency;
            relative_permittivity,
        ).phase_shunt_admittance_matrix_s_per_m
        push!(series, impedance)
        push!(shunt, admittance)
        maximum_error = max(maximum_error, error)
    end
    return LineParameterSegment(
        id,
        :cable,
        length_m,
        phase_order,
        frequencies,
        series,
        shunt,
        source;
        soil_signature_sha256=profile.deterministic_signature_sha256,
        quadrature_error_ohm_per_m=maximum_error,
    )
end

struct LineParameterDiagnostics
    series_symmetry_max_abs_error::Vector{Float64}
    shunt_symmetry_max_abs_error::Vector{Float64}
    minimum_series_loss_eigenvalue_ohm::Vector{Float64}
    minimum_shunt_loss_eigenvalue_s::Vector{Float64}
    series_condition_number::Vector{Float64}
    shunt_condition_number::Vector{Float64}
    quadrature_error_ohm_per_m::Float64
    physical_checks_passed::Bool
end

struct WidebandLineParameterSet
    schema_version::Int
    route_id::Symbol
    origin::Symbol
    phase_order::Vector{Symbol}
    frequencies_hz::Vector{Float64}
    segments::Vector{LineParameterSegment}
    total_length_m::Float64
    route_series_impedance_matrices_ohm::Vector{Matrix{ComplexF64}}
    route_shunt_admittance_matrices_s::Vector{Matrix{ComplexF64}}
    average_series_impedance_matrices_ohm_per_m::Vector{Matrix{ComplexF64}}
    average_shunt_admittance_matrices_s_per_m::Vector{Matrix{ComplexF64}}
    diagnostics::LineParameterDiagnostics
    deterministic_signature_sha256::String
end

function _line_parameter_permutation(segment_order, route_order)
    length(segment_order) == length(route_order) ||
        throw(ArgumentError("line parameter segment phase count differs from route"))
    indices = Int[]
    for phase in route_order
        index = findfirst(==(phase), segment_order)
        index === nothing && throw(ArgumentError(
            "line parameter segment phase basis is incompatible with route phase $(phase)",
        ))
        push!(indices, index)
    end
    return indices
end

function _line_parameter_diagnostics(series, shunt, quadrature_error)
    series_symmetry = Float64[]
    shunt_symmetry = Float64[]
    series_loss = Float64[]
    shunt_loss = Float64[]
    series_condition = Float64[]
    shunt_condition = Float64[]
    for index in eachindex(series)
        z = series[index]
        y = shunt[index]
        push!(series_symmetry, maximum(abs, z - transpose(z)))
        push!(shunt_symmetry, maximum(abs, y - transpose(y)))
        push!(series_loss, eigmin(Hermitian((z + adjoint(z)) / 2.0)))
        push!(shunt_loss, eigmin(Hermitian((y + adjoint(y)) / 2.0)))
        push!(series_condition, cond(z))
        push!(shunt_condition, iszero(y) ? Inf : cond(y))
    end
    finite_values = all(isfinite, series_symmetry) &&
        all(isfinite, shunt_symmetry) &&
        all(isfinite, series_loss) &&
        all(isfinite, shunt_loss) &&
        all(isfinite, series_condition)
    physical = finite_values &&
        maximum(series_symmetry) <= 1.0e-9 &&
        maximum(shunt_symmetry) <= 1.0e-12 &&
        minimum(series_loss) >= LINE_PARAMETER_LOSS_EIGENVALUE_FLOOR &&
        minimum(shunt_loss) >= LINE_PARAMETER_LOSS_EIGENVALUE_FLOOR
    return LineParameterDiagnostics(
        series_symmetry,
        shunt_symmetry,
        series_loss,
        shunt_loss,
        series_condition,
        shunt_condition,
        quadrature_error,
        physical,
    )
end

function _line_parameter_set_signature(route_id, phase_order, frequencies, segments, series, shunt)
    io = IOBuffer()
    println(io, LINE_PARAMETER_SCHEMA_VERSION, '|', route_id, '|', join(phase_order, ','))
    for segment in segments
        println(io, segment.input_signature_sha256)
    end
    for index in eachindex(frequencies)
        println(io, frequencies[index])
        for matrix in (series[index], shunt[index]), value in matrix
            println(io, bitstring(real(value)), '|', bitstring(imag(value)))
        end
    end
    return bytes2hex(sha256(take!(io)))
end

function line_parameter_set(
    route_id::Symbol,
    segments::AbstractVector{LineParameterSegment};
    phase_order::AbstractVector{Symbol}=isempty(segments) ? Symbol[] :
        first(segments).phase_order,
    origin::Symbol=:native,
)
    rows = collect(segments)
    1 <= length(rows) <= LINE_PARAMETER_MAXIMUM_ROUTE_SEGMENTS || throw(ArgumentError(
        "line parameter route requires one through 32 segments",
    ))
    origin in (:native, :imported) || throw(ArgumentError(
        "line parameter set origin must be :native or :imported",
    ))
    phases = collect(phase_order)
    isempty(phases) && throw(ArgumentError("line parameter route requires phase order"))
    frequencies = copy(first(rows).frequencies_hz)
    total_length = sum(segment.length_m for segment in rows)
    total_length <= LINE_PARAMETER_MAXIMUM_ROUTE_LENGTH_M || throw(ArgumentError(
        "line parameter route length exceeds the released domain",
    ))
    series = [zeros(ComplexF64, length(phases), length(phases)) for _ in frequencies]
    shunt = [zeros(ComplexF64, length(phases), length(phases)) for _ in frequencies]
    maximum_quadrature_error = 0.0
    for segment in rows
        segment.frequencies_hz == frequencies || throw(ArgumentError(
            "line parameter route segments require identical frequency grids",
        ))
        indices = _line_parameter_permutation(segment.phase_order, phases)
        for frequency_index in eachindex(frequencies)
            series[frequency_index] .+= segment.length_m .*
                segment.series_impedance_matrices_ohm_per_m[frequency_index][indices, indices]
            shunt[frequency_index] .+= segment.length_m .*
                segment.shunt_admittance_matrices_s_per_m[frequency_index][indices, indices]
        end
        maximum_quadrature_error = max(
            maximum_quadrature_error,
            segment.quadrature_error_ohm_per_m,
        )
    end
    average_series = [matrix ./ total_length for matrix in series]
    average_shunt = [matrix ./ total_length for matrix in shunt]
    diagnostics = _line_parameter_diagnostics(
        series,
        shunt,
        maximum_quadrature_error,
    )
    diagnostics.physical_checks_passed || throw(ArgumentError(
        "line parameter route failed symmetry or passive-loss checks",
    ))
    signature = _line_parameter_set_signature(
        route_id,
        phases,
        frequencies,
        rows,
        series,
        shunt,
    )
    return WidebandLineParameterSet(
        LINE_PARAMETER_SCHEMA_VERSION,
        route_id,
        origin,
        phases,
        frequencies,
        rows,
        total_length,
        series,
        shunt,
        average_series,
        average_shunt,
        diagnostics,
        signature,
    )
end

struct LineParameterUncertaintyEnvelope
    nominal_signature_sha256::String
    alternative_signatures_sha256::Vector{String}
    series_absolute_radius_ohm::Vector{Matrix{Float64}}
    shunt_absolute_radius_s::Vector{Matrix{Float64}}
    series_relative_max::Float64
    shunt_relative_max::Float64
end

function line_parameter_uncertainty_envelope(
    nominal::WidebandLineParameterSet,
    alternatives::AbstractVector{WidebandLineParameterSet},
)
    rows = collect(alternatives)
    isempty(rows) && throw(ArgumentError(
        "line parameter uncertainty envelope requires at least one alternative",
    ))
    series_radius = [zeros(Float64, size(matrix)) for matrix in nominal.route_series_impedance_matrices_ohm]
    shunt_radius = [zeros(Float64, size(matrix)) for matrix in nominal.route_shunt_admittance_matrices_s]
    series_relative = 0.0
    shunt_relative = 0.0
    for alternative in rows
        alternative.phase_order == nominal.phase_order &&
            alternative.frequencies_hz == nominal.frequencies_hz ||
            throw(ArgumentError("line parameter uncertainty alternatives must share basis and grid"))
        for index in eachindex(nominal.frequencies_hz)
            series_delta = abs.(
                alternative.route_series_impedance_matrices_ohm[index] -
                nominal.route_series_impedance_matrices_ohm[index],
            )
            shunt_delta = abs.(
                alternative.route_shunt_admittance_matrices_s[index] -
                nominal.route_shunt_admittance_matrices_s[index],
            )
            series_radius[index] .= max.(series_radius[index], series_delta)
            shunt_radius[index] .= max.(shunt_radius[index], shunt_delta)
            series_relative = max(
                series_relative,
                maximum(series_delta) /
                max(maximum(abs, nominal.route_series_impedance_matrices_ohm[index]), 1.0e-15),
            )
            shunt_relative = max(
                shunt_relative,
                maximum(shunt_delta) /
                max(maximum(abs, nominal.route_shunt_admittance_matrices_s[index]), 1.0e-18),
            )
        end
    end
    return LineParameterUncertaintyEnvelope(
        nominal.deterministic_signature_sha256,
        getfield.(rows, :deterministic_signature_sha256),
        series_radius,
        shunt_radius,
        series_relative,
        shunt_relative,
    )
end

function _line_source_dictionary(source::LineParameterSource)
    provenance = source.provenance
    return Dict{String,Any}(
        "source" => provenance.source,
        "units" => provenance.units,
        "transformation" => provenance.transformation,
        "uncertainty" => provenance.uncertainty,
        "validity_domain" => provenance.validity_domain,
        "nature" => "PhysicalModelParameter",
        "rights" => source.rights,
        "content_sha256" => source.content_sha256,
        "data_class" => String(source.data_class),
    )
end

function _line_source_from_dictionary(data)
    get(data, "nature", "") == "PhysicalModelParameter" || throw(ArgumentError(
        "imported line parameter source nature is unsupported",
    ))
    return LineParameterSource(
        ParameterProvenance(
            data["source"],
            data["units"],
            data["transformation"],
            data["uncertainty"],
            data["validity_domain"],
            PhysicalModelParameter,
        ),
        data["rights"],
        data["content_sha256"],
        Symbol(data["data_class"]),
    )
end

function _line_complex_matrix_data(matrix)
    return [
        [[real(matrix[row, column]), imag(matrix[row, column])] for column in axes(matrix, 2)]
        for row in axes(matrix, 1)
    ]
end

function _line_complex_matrix_from_data(data)
    rows = length(data)
    rows > 0 || throw(ArgumentError("imported complex matrix must not be empty"))
    columns = length(first(data))
    rows == columns || throw(ArgumentError("imported complex matrix must be square"))
    matrix = Matrix{ComplexF64}(undef, rows, columns)
    for row in 1:rows
        length(data[row]) == columns || throw(ArgumentError(
            "imported complex matrix rows must have equal length",
        ))
        for column in 1:columns
            pair = data[row][column]
            length(pair) == 2 || throw(ArgumentError(
                "imported complex matrix entries require real and imaginary values",
            ))
            matrix[row, column] = ComplexF64(pair[1], pair[2])
        end
    end
    return matrix
end

function write_line_parameter_set(path::AbstractString, parameters::WidebandLineParameterSet)
    data = Dict{String,Any}(
        "schema" => "aimora-wideband-line-parameters-v1",
        "schema_version" => parameters.schema_version,
        "route_id" => String(parameters.route_id),
        "origin" => String(parameters.origin),
        "phase_order" => String.(parameters.phase_order),
        "frequencies_hz" => parameters.frequencies_hz,
        "total_length_m" => parameters.total_length_m,
        "deterministic_signature_sha256" => parameters.deterministic_signature_sha256,
        "segment" => [
            Dict{String,Any}(
                "id" => String(segment.id),
                "kind" => String(segment.kind),
                "length_m" => segment.length_m,
                "phase_order" => String.(segment.phase_order),
                "frequencies_hz" => segment.frequencies_hz,
                "soil_signature_sha256" => segment.soil_signature_sha256,
                "quadrature_error_ohm_per_m" => segment.quadrature_error_ohm_per_m,
                "input_signature_sha256" => segment.input_signature_sha256,
                "source_record" => _line_source_dictionary(segment.source),
                "series_impedance_matrices_ohm_per_m" => [
                    _line_complex_matrix_data(matrix)
                    for matrix in segment.series_impedance_matrices_ohm_per_m
                ],
                "shunt_admittance_matrices_s_per_m" => [
                    _line_complex_matrix_data(matrix)
                    for matrix in segment.shunt_admittance_matrices_s_per_m
                ],
            )
            for segment in parameters.segments
        ],
    )
    open(path, "w") do io
        TOML.print(io, data; sorted=true)
    end
    return String(path)
end

function read_line_parameter_set(path::AbstractString)
    data = TOML.parsefile(path)
    get(data, "schema", "") == "aimora-wideband-line-parameters-v1" ||
        throw(ArgumentError("unsupported line parameter interchange schema"))
    get(data, "schema_version", 0) == LINE_PARAMETER_SCHEMA_VERSION ||
        throw(ArgumentError("unsupported line parameter interchange version"))
    segments = LineParameterSegment[]
    for record in data["segment"]
        segment = LineParameterSegment(
            Symbol(record["id"]),
            Symbol(record["kind"]),
            record["length_m"],
            Symbol.(record["phase_order"]),
            record["frequencies_hz"],
            [_line_complex_matrix_from_data(matrix) for matrix in
                record["series_impedance_matrices_ohm_per_m"]],
            [_line_complex_matrix_from_data(matrix) for matrix in
                record["shunt_admittance_matrices_s_per_m"]],
            _line_source_from_dictionary(record["source_record"]);
            soil_signature_sha256=record["soil_signature_sha256"],
            quadrature_error_ohm_per_m=record["quadrature_error_ohm_per_m"],
        )
        push!(segments, segment)
    end
    parameters = line_parameter_set(
        Symbol(data["route_id"]),
        segments;
        phase_order=Symbol.(data["phase_order"]),
        origin=:imported,
    )
    parameters.deterministic_signature_sha256 ==
        data["deterministic_signature_sha256"] || throw(ArgumentError(
        "imported line parameter deterministic signature mismatch",
    ))
    return parameters
end

function line_parameter_report_text(parameters::WidebandLineParameterSet)
    io = IOBuffer()
    println(io, "AIMORA WIDEBAND LINE PARAMETERS")
    println(io, "schema_version ", parameters.schema_version)
    println(io, "route_id ", parameters.route_id)
    println(io, "origin ", parameters.origin)
    println(io, "phase_order ", join(parameters.phase_order, ","))
    println(io, "segment_count ", length(parameters.segments))
    println(io, "frequency_count ", length(parameters.frequencies_hz))
    println(io, "total_length_m ", parameters.total_length_m)
    println(io, "physical_checks_passed ", parameters.diagnostics.physical_checks_passed)
    println(io, "quadrature_error_ohm_per_m ", parameters.diagnostics.quadrature_error_ohm_per_m)
    for segment in parameters.segments
        println(
            io,
            "segment ",
            segment.id,
            " kind=",
            segment.kind,
            " length_m=",
            segment.length_m,
            " source=",
            segment.source.content_sha256,
        )
    end
    for frequency_index in eachindex(parameters.frequencies_hz)
        println(io, "frequency_hz ", parameters.frequencies_hz[frequency_index])
        for (label, matrix) in (
            ("route_series_impedance_ohm", parameters.route_series_impedance_matrices_ohm[frequency_index]),
            ("route_shunt_admittance_s", parameters.route_shunt_admittance_matrices_s[frequency_index]),
        )
            println(io, label)
            for row in axes(matrix, 1)
                println(
                    io,
                    join(
                        (
                            string(real(matrix[row, column]), ",", imag(matrix[row, column]))
                            for column in axes(matrix, 2)
                        ),
                        " ",
                    ),
                )
            end
        end
    end
    println(io, "deterministic_signature_sha256 ", parameters.deterministic_signature_sha256)
    println(io, "END WIDEBAND LINE PARAMETERS")
    return String(take!(io))
end

function write_line_parameter_report(path::AbstractString, parameters::WidebandLineParameterSet)
    write(path, line_parameter_report_text(parameters))
    return String(path)
end
