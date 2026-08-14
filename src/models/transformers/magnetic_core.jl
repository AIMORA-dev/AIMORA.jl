export AbstractTransformerMagneticMaterial,
       LinearTransformerMagneticMaterial,
       PiecewiseLinearTransformerMagneticMaterial,
       TellinenLimitingCurve,
       TellinenTransformerMagneticMaterial,
       TellinenMagneticState,
       TransformerDynamicCoreLossModel,
       MagneticBranchGeometry,
       TransformerMagneticGraph,
       MagneticGraphLinearResponse,
       magnetic_material_field,
       magnetic_material_differential_reluctivity,
       tellinen_state,
       tellinen_trial_from_flux_density,
       transformer_magnetic_linear_response,
       transformer_magnetic_linear_inductance,
       transformer_magnetic_graph_signature

const TRANSFORMER_VACUUM_PERMEABILITY_H_PER_M = 4.0e-7 * pi

abstract type AbstractTransformerMagneticMaterial end

"""Rate-dependent Bertotti-style core-loss field kept separate from static hysteresis.

The classical term is proportional to `dB/dt`. The excess term approaches
`sign(dB/dt)*sqrt(abs(dB/dt))` outside a small declared regularization rate and
remains differentiable at zero for the coupled Newton solve. Both coefficients
therefore produce nonnegative `H_dynamic*dB/dt` and neither changes the static
magnetization or Tellinen loop law.
"""
struct TransformerDynamicCoreLossModel
    classical_eddy_coefficient_a_s_per_m_t::Float64
    excess_loss_coefficient_a_sqrt_s_per_m_sqrt_t::Float64
    excess_rate_regularization_t_per_s::Float64
    source::TransformerSourceRecord

    function TransformerDynamicCoreLossModel(
        classical_eddy_coefficient_a_s_per_m_t::Real,
        excess_loss_coefficient_a_sqrt_s_per_m_sqrt_t::Real,
        source::TransformerSourceRecord;
        excess_rate_regularization_t_per_s::Real=1.0e-9,
    )
        classical = Float64(classical_eddy_coefficient_a_s_per_m_t)
        excess = Float64(excess_loss_coefficient_a_sqrt_s_per_m_sqrt_t)
        regularization = Float64(excess_rate_regularization_t_per_s)
        isfinite(classical) && classical >= 0.0 || throw(ArgumentError(
            "transformer classical eddy-loss coefficient must be finite and nonnegative",
        ))
        isfinite(excess) && excess >= 0.0 || throw(ArgumentError(
            "transformer excess core-loss coefficient must be finite and nonnegative",
        ))
        classical > 0.0 || excess > 0.0 || throw(ArgumentError(
            "transformer dynamic core-loss model must identify a classical or excess term",
        ))
        isfinite(regularization) && regularization > 0.0 || throw(ArgumentError(
            "transformer excess-loss rate regularization must be finite and positive",
        ))
        return new(classical, excess, regularization, source)
    end
end

function _transformer_dynamic_core_loss_field(
    loss::TransformerDynamicCoreLossModel,
    flux_density_rate_t_per_s::Float64,
)
    rate = flux_density_rate_t_per_s
    regularized_rate_squared = muladd(
        rate,
        rate,
        loss.excess_rate_regularization_t_per_s^2,
    )
    fourth_root = sqrt(sqrt(regularized_rate_squared))
    classical_field = loss.classical_eddy_coefficient_a_s_per_m_t * rate
    excess_field = loss.excess_loss_coefficient_a_sqrt_s_per_m_sqrt_t *
        rate / fourth_root
    classical_derivative = loss.classical_eddy_coefficient_a_s_per_m_t
    excess_derivative = loss.excess_loss_coefficient_a_sqrt_s_per_m_sqrt_t * (
        inv(fourth_root) - 0.5 * rate^2 / (regularized_rate_squared * fourth_root)
    )
    return (
        classical_field_a_per_m=classical_field,
        excess_field_a_per_m=excess_field,
        differential_field_a_s_per_m_t=classical_derivative + excess_derivative,
        classical_loss_density_w_per_m3=classical_field * rate,
        excess_loss_density_w_per_m3=excess_field * rate,
    )
end

"""Linear isotropic magnetic material for air-core, linear-core, and analytic limits."""
struct LinearTransformerMagneticMaterial <: AbstractTransformerMagneticMaterial
    relative_permeability::Float64
    maximum_flux_density_t::Float64
    source::TransformerSourceRecord

    function LinearTransformerMagneticMaterial(
        relative_permeability::Real,
        maximum_flux_density_t::Real,
        source::TransformerSourceRecord,
    )
        permeability = Float64(relative_permeability)
        maximum_flux_density = Float64(maximum_flux_density_t)
        isfinite(permeability) && permeability >= 1.0 || throw(ArgumentError(
            "transformer relative permeability must be finite and at least one",
        ))
        isfinite(maximum_flux_density) && maximum_flux_density > 0.0 ||
            throw(ArgumentError(
                "transformer maximum flux density must be finite and positive",
            ))
        return new(permeability, maximum_flux_density, source)
    end
end

"""Odd-symmetric monotone B-H law used for anhysteretic saturation without hidden extrapolation."""
struct PiecewiseLinearTransformerMagneticMaterial <: AbstractTransformerMagneticMaterial
    flux_density_t::Vector{Float64}
    field_strength_a_per_m::Vector{Float64}
    source::TransformerSourceRecord

    function PiecewiseLinearTransformerMagneticMaterial(
        flux_density_t,
        field_strength_a_per_m,
        source::TransformerSourceRecord,
    )
        flux_density = Float64.(flux_density_t)
        field_strength = Float64.(field_strength_a_per_m)
        length(flux_density) == length(field_strength) >= 3 || throw(ArgumentError(
            "transformer saturation material requires at least three aligned B-H points",
        ))
        first(flux_density) == 0.0 && first(field_strength) == 0.0 ||
            throw(ArgumentError("transformer saturation material must begin at B=H=0"))
        all(isfinite, flux_density) && all(isfinite, field_strength) ||
            throw(ArgumentError("transformer saturation material points must be finite"))
        all(diff(flux_density) .> 0.0) && all(diff(field_strength) .> 0.0) ||
            throw(ArgumentError(
                "transformer saturation B-H magnitudes must increase strictly",
            ))
        last(flux_density) <= 2.5 || throw(ArgumentError(
            "transformer public magnetic material exceeds the admitted 2.5 T domain",
        ))
        return new(flux_density, field_strength, source)
    end
end

function _linear_interpolation_and_slope(x::Float64, grid_x, grid_y)
    first(grid_x) <= x <= last(grid_x) || throw(DomainError(
        x,
        "magnetic material evaluation is outside its declared source domain",
    ))
    index = x == last(grid_x) ? length(grid_x) - 1 : searchsortedlast(grid_x, x)
    index = clamp(index, 1, length(grid_x) - 1)
    fraction = (x - grid_x[index]) / (grid_x[index + 1] - grid_x[index])
    slope = (grid_y[index + 1] - grid_y[index]) /
        (grid_x[index + 1] - grid_x[index])
    return muladd(fraction, grid_y[index + 1] - grid_y[index], grid_y[index]), slope
end

function magnetic_material_field(
    material::LinearTransformerMagneticMaterial,
    flux_density_t::Real,
)
    flux_density = Float64(flux_density_t)
    isfinite(flux_density) && abs(flux_density) <= material.maximum_flux_density_t ||
        throw(DomainError(
            flux_density,
            "linear transformer material flux density exceeds its declared domain",
        ))
    return flux_density /
        (TRANSFORMER_VACUUM_PERMEABILITY_H_PER_M * material.relative_permeability)
end

magnetic_material_differential_reluctivity(
    material::LinearTransformerMagneticMaterial,
    _flux_density_t::Real,
) = inv(TRANSFORMER_VACUUM_PERMEABILITY_H_PER_M * material.relative_permeability)

function magnetic_material_field(
    material::PiecewiseLinearTransformerMagneticMaterial,
    flux_density_t::Real,
)
    flux_density = Float64(flux_density_t)
    isfinite(flux_density) || throw(ArgumentError(
        "transformer magnetic flux density must be finite",
    ))
    field, _ = _linear_interpolation_and_slope(
        abs(flux_density),
        material.flux_density_t,
        material.field_strength_a_per_m,
    )
    return copysign(field, flux_density)
end

function magnetic_material_differential_reluctivity(
    material::PiecewiseLinearTransformerMagneticMaterial,
    flux_density_t::Real,
)
    flux_density = Float64(flux_density_t)
    isfinite(flux_density) || throw(ArgumentError(
        "transformer magnetic flux density must be finite",
    ))
    _, slope = _linear_interpolation_and_slope(
        abs(flux_density),
        material.flux_density_t,
        material.field_strength_a_per_m,
    )
    return slope
end

"""One monotone limiting B(H) curve of the scalar Tellinen model."""
struct TellinenLimitingCurve
    field_strength_a_per_m::Vector{Float64}
    flux_density_t::Vector{Float64}

    function TellinenLimitingCurve(field_strength_a_per_m, flux_density_t)
        field = Float64.(field_strength_a_per_m)
        flux_density = Float64.(flux_density_t)
        length(field) == length(flux_density) >= 5 || throw(ArgumentError(
            "Tellinen limiting curve requires at least five aligned points",
        ))
        all(isfinite, field) && all(isfinite, flux_density) ||
            throw(ArgumentError("Tellinen limiting-curve points must be finite"))
        all(diff(field) .> 0.0) && all(diff(flux_density) .> 0.0) ||
            throw(ArgumentError(
                "Tellinen limiting-curve field and flux density must increase strictly",
            ))
        first(field) < 0.0 < last(field) || throw(ArgumentError(
            "Tellinen limiting curve must span negative and positive field",
        ))
        maximum(abs, flux_density) <= 2.5 || throw(ArgumentError(
            "Tellinen limiting curve exceeds the admitted 2.5 T public domain",
        ))
        return new(field, flux_density)
    end
end

function _tellinen_curve_value(curve::TellinenLimitingCurve, field_strength::Float64)
    return _linear_interpolation_and_slope(
        field_strength,
        curve.field_strength_a_per_m,
        curve.flux_density_t,
    )
end

function _tellinen_curve_value(
    curve::TellinenLimitingCurve,
    field_strength::Float64,
    direction::Int,
)
    direction >= 0 && return _tellinen_curve_value(curve, field_strength)
    field = curve.field_strength_a_per_m
    flux_density = curve.flux_density_t
    first(field) <= field_strength <= last(field) || throw(DomainError(
        field_strength,
        "magnetic material evaluation is outside its declared source domain",
    ))
    index = field_strength == first(field) ? 1 :
        searchsortedfirst(field, field_strength) - 1
    index = clamp(index, 1, length(field) - 1)
    fraction = (field_strength - field[index]) / (field[index + 1] - field[index])
    slope = (flux_density[index + 1] - flux_density[index]) /
        (field[index + 1] - field[index])
    value = muladd(
        fraction,
        flux_density[index + 1] - flux_density[index],
        flux_density[index],
    )
    return value, slope
end

"""Scalar rate-independent Tellinen material with ordered lower/rising and upper/falling envelopes."""
struct TellinenTransformerMagneticMaterial <: AbstractTransformerMagneticMaterial
    lower_branch::TellinenLimitingCurve
    upper_branch::TellinenLimitingCurve
    minimum_branch_fraction::Float64
    integration_field_increment_a_per_m::Float64
    source::TransformerSourceRecord

    function TellinenTransformerMagneticMaterial(
        lower_branch::TellinenLimitingCurve,
        upper_branch::TellinenLimitingCurve,
        source::TransformerSourceRecord;
        minimum_branch_fraction::Real=1.0e-9,
        integration_field_increment_a_per_m::Real=1.0,
    )
        lower_branch.field_strength_a_per_m == upper_branch.field_strength_a_per_m ||
            throw(ArgumentError("Tellinen limiting curves must share one H grid"))
        separation = upper_branch.flux_density_t .- lower_branch.flux_density_t
        all(separation .> 0.0) || throw(ArgumentError(
            "Tellinen upper limiting branch must remain above its lower branch",
        ))
        anhysteretic_flux_density = 0.5 .* (
            upper_branch.flux_density_t .+ lower_branch.flux_density_t
        )
        all(diff(anhysteretic_flux_density) .> 0.0) || throw(ArgumentError(
            "Tellinen limiting curves must define a monotone anhysteretic midpoint",
        ))
        first(anhysteretic_flux_density) < 0.0 < last(anhysteretic_flux_density) ||
            throw(ArgumentError(
                "Tellinen anhysteretic midpoint must span zero flux density",
            ))
        minimum_fraction = Float64(minimum_branch_fraction)
        isfinite(minimum_fraction) && 0.0 < minimum_fraction <= 1.0 ||
            throw(ArgumentError(
                "Tellinen minimum branch fraction must lie in (0,1]",
            ))
        field_increment = Float64(integration_field_increment_a_per_m)
        isfinite(field_increment) && field_increment > 0.0 || throw(ArgumentError(
            "Tellinen integration field increment must be finite and positive",
        ))
        return new(
            lower_branch,
            upper_branch,
            minimum_fraction,
            field_increment,
            source,
        )
    end
end

"""Accepted scalar Tellinen state; reversal fields mutate only after a converged apparatus step."""
struct TellinenMagneticState
    flux_density_t::Float64
    field_strength_a_per_m::Float64
    direction::Int
    reversal_flux_density_t::Float64
    reversal_field_strength_a_per_m::Float64
    loop_integral_j_per_m3::Float64
    reversal_count::Int
end

function tellinen_state(
    material::TellinenTransformerMagneticMaterial;
    field_strength_a_per_m::Real=0.0,
    flux_density_t::Union{Nothing,Real}=nothing,
    direction::Integer=1,
)
    field = Float64(field_strength_a_per_m)
    lower, _ = _tellinen_curve_value(material.lower_branch, field)
    upper, _ = _tellinen_curve_value(material.upper_branch, field)
    flux_density = flux_density_t === nothing ? 0.5 * (lower + upper) :
        Float64(flux_density_t)
    lower <= flux_density <= upper || throw(ArgumentError(
        "initial Tellinen flux density must lie between the limiting branches",
    ))
    direction_value = Int(direction)
    direction_value in (-1, 1) || throw(ArgumentError(
        "Tellinen direction must be -1 or +1",
    ))
    return TellinenMagneticState(
        flux_density,
        field,
        direction_value,
        flux_density,
        field,
        0.0,
        0,
    )
end

function _tellinen_flux_derivative(
    material::TellinenTransformerMagneticMaterial,
    field_strength::Float64,
    flux_density::Float64,
    direction::Int,
    slope_direction::Int=direction,
)
    lower, lower_slope =
        _tellinen_curve_value(material.lower_branch, field_strength, slope_direction)
    upper, upper_slope =
        _tellinen_curve_value(material.upper_branch, field_strength, slope_direction)
    separation = upper - lower
    separation > 0.0 || throw(ArgumentError(
        "Tellinen limiting branches have nonpositive separation",
    ))
    vacuum = TRANSFORMER_VACUUM_PERMEABILITY_H_PER_M
    if direction > 0
        fraction = max(
            material.minimum_branch_fraction,
            (upper - flux_density) / separation,
        )
        return fraction * (lower_slope - vacuum) + vacuum
    end
    fraction = max(
        material.minimum_branch_fraction,
        (flux_density - lower) / separation,
    )
    return fraction * (upper_slope - vacuum) + vacuum
end

function _tellinen_integrate_to_field(
    material::TellinenTransformerMagneticMaterial,
    state::TellinenMagneticState,
    target_field::Float64,
)
    field_delta = target_field - state.field_strength_a_per_m
    direction = field_delta == 0.0 ? state.direction : (field_delta > 0.0 ? 1 : -1)
    field_delta == 0.0 && return state.flux_density_t, direction
    maximum_field_step = material.integration_field_increment_a_per_m
    field_grid = material.lower_branch.field_strength_a_per_m
    field = state.field_strength_a_per_m
    flux_density = state.flux_density_t
    while direction > 0 ? field < target_field : field > target_field
        segment_target = target_field
        if direction > 0
            next_grid_index = searchsortedlast(field_grid, field) + 1
            if next_grid_index <= length(field_grid) &&
               field_grid[next_grid_index] < target_field
                segment_target = field_grid[next_grid_index]
            end
        else
            previous_grid_index = searchsortedfirst(field_grid, field) - 1
            if previous_grid_index >= 1 && field_grid[previous_grid_index] > target_field
                segment_target = field_grid[previous_grid_index]
            end
        end
        segment_start = field
        segment_delta = segment_target - segment_start
        substep_count = ceil(Int, abs(segment_delta) / maximum_field_step)
        uniform_field_step = segment_delta / substep_count
        for substep in 1:substep_count
            next_field = substep == substep_count ? segment_target :
                segment_start + substep * uniform_field_step
            field_step = next_field - field
            derivative_1 = _tellinen_flux_derivative(
                material,
                field,
                flux_density,
                direction,
            )
            midpoint_field = field + 0.5 * field_step
            midpoint_flux_density = flux_density + 0.5 * field_step * derivative_1
            derivative_2 = _tellinen_flux_derivative(
                material,
                midpoint_field,
                midpoint_flux_density,
                direction,
            )
            derivative_3 = _tellinen_flux_derivative(
                material,
                midpoint_field,
                flux_density + 0.5 * field_step * derivative_2,
                direction,
            )
            derivative_4 = _tellinen_flux_derivative(
                material,
                next_field,
                flux_density + field_step * derivative_3,
                direction,
                -direction,
            )
            flux_density += field_step * (
                derivative_1 + 2.0 * derivative_2 + 2.0 * derivative_3 + derivative_4
            ) / 6.0
            field = next_field
            lower, _ = _tellinen_curve_value(material.lower_branch, field)
            upper, _ = _tellinen_curve_value(material.upper_branch, field)
            flux_density = clamp(flux_density, lower, upper)
        end
    end
    return flux_density, direction
end

"""Pure inverse Tellinen trial: find H for a voltage-integrated target B without mutating accepted state."""
function tellinen_trial_from_flux_density(
    material::TellinenTransformerMagneticMaterial,
    state::TellinenMagneticState,
    target_flux_density_t::Real;
    flux_density_tolerance_t::Real=1.0e-11,
    maximum_iterations::Integer=80,
)
    target_flux_density = Float64(target_flux_density_t)
    tolerance = Float64(flux_density_tolerance_t)
    isfinite(target_flux_density) && isfinite(tolerance) && tolerance > 0.0 ||
        throw(ArgumentError("Tellinen trial target and tolerance must be finite"))
    iterations = Int(maximum_iterations)
    iterations > 0 || throw(ArgumentError(
        "Tellinen maximum iterations must be positive",
    ))
    minimum_field = first(material.lower_branch.field_strength_a_per_m)
    maximum_field = last(material.lower_branch.field_strength_a_per_m)
    minimum_flux, _ = _tellinen_integrate_to_field(material, state, minimum_field)
    maximum_flux, _ = _tellinen_integrate_to_field(material, state, maximum_field)
    minimum_flux - tolerance <= target_flux_density <= maximum_flux + tolerance ||
        throw(DomainError(
            target_flux_density,
            "Tellinen target flux density lies outside the reachable limiting-curve domain",
        ))
    lower_field = minimum_field
    upper_field = maximum_field
    field = state.field_strength_a_per_m
    flux_density = state.flux_density_t
    direction = state.direction
    for iteration in 1:iterations
        field = 0.5 * (lower_field + upper_field)
        flux_density, direction = _tellinen_integrate_to_field(material, state, field)
        error = flux_density - target_flux_density
        if abs(error) <= tolerance
            reversal = direction != state.direction
            loop_increment = 0.5 * (state.field_strength_a_per_m + field) *
                (target_flux_density - state.flux_density_t)
            trial_state = TellinenMagneticState(
                target_flux_density,
                field,
                direction,
                reversal ? state.flux_density_t : state.reversal_flux_density_t,
                reversal ? state.field_strength_a_per_m : state.reversal_field_strength_a_per_m,
                state.loop_integral_j_per_m3 + loop_increment,
                state.reversal_count + Int(reversal),
            )
            derivative = _tellinen_flux_derivative(
                material,
                field,
                target_flux_density,
                direction,
            )
            derivative > 0.0 || throw(ArgumentError(
                "Tellinen differential permeability became nonpositive",
            ))
            return (
                state=trial_state,
                differential_reluctivity_m_per_h=inv(derivative),
                iterations=iteration,
                residual_t=error,
            )
        elseif error > 0.0
            upper_field = field
        else
            lower_field = field
        end
    end
    throw(ArgumentError("Tellinen inverse flux-density trial did not converge"))
end

"""One oriented core/yoke/return/gap branch of an explicit magnetic graph."""
struct MagneticBranchGeometry
    id::Symbol
    length_m::Float64
    cross_section_m2::Float64
    air_gap_length_m::Float64
    air_gap_effective_area_factor::Float64
    material_index::Int

    function MagneticBranchGeometry(
        id::Symbol;
        length_m::Real,
        cross_section_m2::Real,
        air_gap_length_m::Real=0.0,
        air_gap_effective_area_factor::Real=1.0,
        material_index::Integer=1,
    )
        id == Symbol("") && throw(ArgumentError(
            "transformer magnetic branch identity must not be empty",
        ))
        length = Float64(length_m)
        area = Float64(cross_section_m2)
        gap = Float64(air_gap_length_m)
        gap_area_factor = Float64(air_gap_effective_area_factor)
        isfinite(length) && length > 0.0 || throw(ArgumentError(
            "transformer magnetic branch length must be finite and positive",
        ))
        isfinite(area) && area > 0.0 || throw(ArgumentError(
            "transformer magnetic branch cross-section must be finite and positive",
        ))
        isfinite(gap) && 0.0 <= gap < length || throw(ArgumentError(
            "transformer magnetic air gap must be finite, nonnegative, and shorter than the branch",
        ))
        isfinite(gap_area_factor) && gap_area_factor >= 1.0 || throw(ArgumentError(
            "transformer magnetic air-gap effective-area factor must be finite and at least one",
        ))
        gap == 0.0 && gap_area_factor != 1.0 && throw(ArgumentError(
            "transformer magnetic branch without an air gap requires unit effective-area factor",
        ))
        material = Int(material_index)
        material > 0 || throw(ArgumentError(
            "transformer magnetic material index must be positive",
        ))
        return new(id, length, area, gap, gap_area_factor, material)
    end
end

struct TransformerMagneticGraph
    node_order::Vector{Symbol}
    branch_order::Vector{Symbol}
    incidence::Matrix{Float64}
    branches::Vector{MagneticBranchGeometry}
    winding_turns::Matrix{Float64}
    materials::Vector{AbstractTransformerMagneticMaterial}
    dynamic_core_loss::Vector{Union{Nothing,TransformerDynamicCoreLossModel}}
    rank::Int
    deterministic_signature_sha256::String
end

function _magnetic_graph_signature(
    nodes,
    branches,
    incidence,
    winding_turns,
    materials,
    dynamic_core_loss,
)
    io = IOBuffer()
    println(io, join(String.(nodes), ','))
    println(io, join(String.(getfield.(branches, :id)), ','))
    for branch in branches
        for value in (
            branch.length_m,
            branch.cross_section_m2,
            branch.air_gap_length_m,
            branch.air_gap_effective_area_factor,
        )
            println(io, bitstring(value))
        end
        println(io, branch.material_index)
    end
    for matrix in (incidence, winding_turns)
        for value in matrix
            println(io, bitstring(value))
        end
    end
    for material in materials
        println(io, string(typeof(material)))
        println(io, material.source.id)
        println(io, material.source.content_sha256)
    end
    for loss in dynamic_core_loss
        if loss === nothing
            println(io, "no_dynamic_core_loss")
        else
            println(io, bitstring(loss.classical_eddy_coefficient_a_s_per_m_t))
            println(io, bitstring(loss.excess_loss_coefficient_a_sqrt_s_per_m_sqrt_t))
            println(io, bitstring(loss.excess_rate_regularization_t_per_s))
            println(io, loss.source.id)
            println(io, loss.source.content_sha256)
        end
    end
    return bytes2hex(sha256(take!(io)))
end

function TransformerMagneticGraph(;
    node_order,
    branches,
    incidence,
    winding_turns,
    materials,
    dynamic_core_loss=nothing,
)
    nodes = _connection_symbols(node_order, "magnetic node order")
    branch_rows = MagneticBranchGeometry[branches...]
    isempty(branch_rows) && throw(ArgumentError(
        "transformer magnetic graph requires at least one branch",
    ))
    branch_ids = getfield.(branch_rows, :id)
    length(unique(branch_ids)) == length(branch_ids) || throw(ArgumentError(
        "transformer magnetic branch identities must be unique",
    ))
    material_rows = AbstractTransformerMagneticMaterial[materials...]
    isempty(material_rows) && throw(ArgumentError(
        "transformer magnetic graph requires at least one material",
    ))
    all(branch -> branch.material_index <= length(material_rows), branch_rows) ||
        throw(ArgumentError("transformer magnetic branch material index is unavailable"))
    dynamic_loss_rows = if dynamic_core_loss === nothing
        Union{Nothing,TransformerDynamicCoreLossModel}[nothing for _ in branch_rows]
    else
        Union{Nothing,TransformerDynamicCoreLossModel}[dynamic_core_loss...]
    end
    length(dynamic_loss_rows) == length(branch_rows) || throw(DimensionMismatch(
        "transformer dynamic core-loss ownership must cover every magnetic branch",
    ))
    graph_incidence = Matrix{Float64}(incidence)
    size(graph_incidence) == (length(nodes), length(branch_rows)) ||
        throw(DimensionMismatch(
            "transformer magnetic incidence size must be node_count by branch_count",
        ))
    all(value -> value in (-1.0, 0.0, 1.0), graph_incidence) ||
        throw(ArgumentError("transformer magnetic incidence entries must be -1, 0, or 1"))
    for branch in axes(graph_incidence, 2)
        nonzero = findall(!iszero, @view graph_incidence[:, branch])
        length(nonzero) in (1, 2) || throw(ArgumentError(
            "transformer magnetic branch must connect a node to reference or two nodes",
        ))
        length(nonzero) == 1 || sort(graph_incidence[nonzero, branch]) == [-1.0, 1.0] ||
            throw(ArgumentError(
                "floating transformer magnetic branch must contain one +1 and one -1",
            ))
    end
    graph_rank = rank(graph_incidence)
    graph_rank == length(nodes) || throw(ArgumentError(
        "transformer magnetic graph must be connected to one declared reference",
    ))
    turns = Matrix{Float64}(winding_turns)
    size(turns, 1) == length(branch_rows) && size(turns, 2) > 0 ||
        throw(DimensionMismatch(
            "transformer magnetic turns matrix must be branch_count by nonzero coil_count",
        ))
    all(isfinite, turns) || throw(ArgumentError(
        "transformer magnetic turns entries must be finite",
    ))
    all(column -> any(!iszero, @view turns[:, column]), axes(turns, 2)) ||
        throw(ArgumentError("every transformer coil must link one magnetic branch"))
    signature = _magnetic_graph_signature(
        nodes,
        branch_rows,
        graph_incidence,
        turns,
        material_rows,
        dynamic_loss_rows,
    )
    return TransformerMagneticGraph(
        nodes,
        branch_ids,
        graph_incidence,
        branch_rows,
        turns,
        material_rows,
        dynamic_loss_rows,
        graph_rank,
        signature,
    )
end

transformer_magnetic_graph_signature(graph::TransformerMagneticGraph) =
    graph.deterministic_signature_sha256

function _linear_branch_reluctances(graph::TransformerMagneticGraph)
    reluctance = zeros(Float64, length(graph.branches))
    for (index, branch) in pairs(graph.branches)
        material = graph.materials[branch.material_index]
        material isa LinearTransformerMagneticMaterial || throw(ArgumentError(
            "linear magnetic response requires a linear material on every branch",
        ))
        core_length = branch.length_m - branch.air_gap_length_m
        reluctance[index] = core_length /
            (
                TRANSFORMER_VACUUM_PERMEABILITY_H_PER_M *
                material.relative_permeability *
                branch.cross_section_m2
            ) + branch.air_gap_length_m /
            (
                TRANSFORMER_VACUUM_PERMEABILITY_H_PER_M *
                branch.cross_section_m2 * branch.air_gap_effective_area_factor
            )
    end
    return reluctance
end

struct MagneticGraphLinearResponse
    winding_current_a::Vector{Float64}
    branch_flux_wb::Vector{Float64}
    branch_flux_density_t::Vector{Float64}
    branch_mmf_drop_at::Vector{Float64}
    magnetic_node_potential_at::Vector{Float64}
    winding_flux_linkage_wb_turn::Vector{Float64}
    stored_energy_j::Float64
    continuity_residual_wb::Float64
    constitutive_residual_at::Float64
end

function transformer_magnetic_linear_response(
    graph::TransformerMagneticGraph,
    winding_current_a,
)
    current = Float64.(winding_current_a)
    length(current) == size(graph.winding_turns, 2) || throw(DimensionMismatch(
        "transformer winding current count must match magnetic turns columns",
    ))
    all(isfinite, current) || throw(ArgumentError(
        "transformer winding currents must be finite",
    ))
    reluctance = _linear_branch_reluctances(graph)
    branch_count = length(graph.branches)
    node_count = length(graph.node_order)
    system = [
        Diagonal(reluctance) transpose(graph.incidence)
        graph.incidence zeros(Float64, node_count, node_count)
    ]
    right_hand_side = vcat(graph.winding_turns * current, zeros(Float64, node_count))
    solution = system \ right_hand_side
    flux = solution[1:branch_count]
    potential = solution[(branch_count + 1):end]
    mmf_drop = reluctance .* flux
    linkage = transpose(graph.winding_turns) * flux
    flux_density = [
        flux[index] / graph.branches[index].cross_section_m2
        for index in eachindex(flux)
    ]
    for (index, density) in pairs(flux_density)
        material = graph.materials[graph.branches[index].material_index]
        abs(density) <= material.maximum_flux_density_t || throw(DomainError(
            density,
            "linear transformer magnetic branch exceeds its declared material domain",
        ))
    end
    continuity = maximum(abs, graph.incidence * flux; init=0.0)
    constitutive = maximum(
        abs,
        mmf_drop + transpose(graph.incidence) * potential -
            graph.winding_turns * current;
        init=0.0,
    )
    energy = 0.5 * dot(flux, mmf_drop)
    return MagneticGraphLinearResponse(
        current,
        flux,
        flux_density,
        mmf_drop,
        potential,
        linkage,
        energy,
        continuity,
        constitutive,
    )
end

function transformer_magnetic_linear_inductance(graph::TransformerMagneticGraph)
    coil_count = size(graph.winding_turns, 2)
    inductance = zeros(Float64, coil_count, coil_count)
    basis = zeros(Float64, coil_count)
    for coil in 1:coil_count
        fill!(basis, 0.0)
        basis[coil] = 1.0
        inductance[:, coil] .=
            transformer_magnetic_linear_response(graph, basis).
            winding_flux_linkage_wb_turn
    end
    symmetry_scale = max(maximum(abs, inductance; init=0.0), 1.0)
    maximum(abs, inductance - transpose(inductance); init=0.0) <=
        256.0 * eps(Float64) * symmetry_scale || throw(ArgumentError(
            "transformer magnetic inductance is not reciprocal",
        ))
    return 0.5 .* (inductance .+ transpose(inductance))
end
