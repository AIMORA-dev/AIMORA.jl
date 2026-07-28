module AIMORACUDAExt

using AIMORA
using CUDA

if AIMORA.solver_available()
    function AIMORA.Nodal.cuda_fixed_admittance_batch(
        admittance::AbstractMatrix{<:Real},
        batch_count::Int,
    )
        CUDA.functional() || throw(ArgumentError(
            "CUDA.jl loaded, but no functional CUDA device is available",
        ))
        device_admittance = CUDA.CuArray(Float64.(admittance))
        return AIMORA.Nodal.FixedAdmittanceBatchWorkspace(
            device_admittance,
            batch_count,
        )
    end

    function AIMORA.Nodal.accelerator_synchronize(
        workspace::AIMORA.Nodal.FixedAdmittanceBatchWorkspace{A},
    ) where {A<:CUDA.CuArray{Float64,2}}
        CUDA.synchronize()
        return workspace
    end

    function AIMORA.Nodal.fixed_admittance_backend_name(
        ::AIMORA.Nodal.FixedAdmittanceBatchWorkspace{A},
    ) where {A<:CUDA.CuArray{Float64,2}}
        return :cuda
    end

    function AIMORA.Nodal.fixed_admittance_device_resident(
        ::AIMORA.Nodal.FixedAdmittanceBatchWorkspace{A},
    ) where {A<:CUDA.CuArray{Float64,2}}
        return true
    end
end

end
