using ..NonlinearNetwork

export FittedZincOxideCurrentBranch

"""Two-terminal odd-symmetric ZnO current branch backed by an accepted fitted characteristic."""
struct FittedZincOxideCurrentBranch <: NonlinearNetwork.AbstractNonlinearCurrentDevice
    positive_node::Int
    negative_node::Int
    characteristic::ZincOxideFitResult
    provenance::NonlinearNetwork.NonlinearParameterProvenance

    function FittedZincOxideCurrentBranch(
        positive_node::Integer,
        negative_node::Integer,
        characteristic::ZincOxideFitResult,
        ;
        provenance::NonlinearNetwork.NonlinearParameterProvenance=
            NonlinearNetwork.NonlinearParameterProvenance(
                "caller-supplied positive current-voltage samples captured by ZincOxideFitResult",
                "ampere and volt",
                "positive continuous piecewise power-law fit with odd-symmetric voltage extension",
                "not supplied; caller owns sample uncertainty and the fit records interpolation error",
                "accepted positive sample range and its declared extrapolation behavior",
                NonlinearNetwork.PhysicalModelParameter,
            ),
    )
        positive = Int(positive_node)
        negative = Int(negative_node)
        positive >= 0 || throw(ArgumentError("ZnO positive terminal must be nonnegative"))
        negative >= 0 || throw(ArgumentError("ZnO negative terminal must be nonnegative"))
        positive != negative || throw(ArgumentError("ZnO terminals must be distinct"))
        characteristic.fit_checks_passed || throw(ArgumentError(
            "ZnO network branch requires an accepted positive continuous fit",
        ))
        provenance.nature === NonlinearNetwork.PhysicalModelParameter || throw(ArgumentError(
            "ZnO branch provenance must describe physical model parameters",
        ))
        return new(positive, negative, characteristic, provenance)
    end
end

NonlinearNetwork.nonlinear_terminal_nodes(device::FittedZincOxideCurrentBranch) =
    (device.positive_node, device.negative_node)

NonlinearNetwork.nonlinear_device_formulation(::FittedZincOxideCurrentBranch) =
    NonlinearNetwork.PhysicalConstitutiveCurrent

NonlinearNetwork.nonlinear_device_provenance(device::FittedZincOxideCurrentBranch) =
    device.provenance

function NonlinearNetwork.nonlinear_current_jacobian!(
    terminal_current_a::AbstractVector{Float64},
    terminal_jacobian_s::AbstractMatrix{Float64},
    device::FittedZincOxideCurrentBranch,
    terminal_voltage_v::AbstractVector{Float64},
    time_s::Float64,
)
    length(terminal_current_a) >= 2 || throw(DimensionMismatch(
        "ZnO current workspace must contain two terminals",
    ))
    size(terminal_jacobian_s, 1) >= 2 && size(terminal_jacobian_s, 2) >= 2 ||
        throw(DimensionMismatch("ZnO Jacobian workspace must be at least 2x2"))
    length(terminal_voltage_v) >= 2 || throw(DimensionMismatch(
        "ZnO voltage workspace must contain two terminals",
    ))
    isfinite(time_s) || throw(ArgumentError("ZnO evaluation time must be finite"))
    branch_voltage_v = terminal_voltage_v[1] - terminal_voltage_v[2]
    evaluation = zinc_oxide_fitted_current_and_derivative(
        device.characteristic,
        branch_voltage_v,
    )
    branch_current_a = evaluation.current_a
    differential_conductance_s = evaluation.derivative_s
    terminal_current_a[1] = branch_current_a
    terminal_current_a[2] = -branch_current_a
    terminal_jacobian_s[1, 1] = differential_conductance_s
    terminal_jacobian_s[1, 2] = -differential_conductance_s
    terminal_jacobian_s[2, 1] = -differential_conductance_s
    terminal_jacobian_s[2, 2] = differential_conductance_s
    return nothing
end
