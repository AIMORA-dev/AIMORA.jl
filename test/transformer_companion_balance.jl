using Test, AIMORA, LinearAlgebra
using AIMORA.NonlinearNetwork

@testset "affine transformer companion preserves an exact resistive operating point" begin
    transformer = AIMORA.TransformerApparatus
    for bctran in (false, true), method in (:backward_euler, :trapezoidal), duration in (1e-5, 4.126e-10)
        runtime = compensated_history_transformer_runtime(; bctran)
        runtime.accepted_state.coil_current_a[1] = 8.0
        runtime.accepted_state.coil_voltage_v[1] = 2.0
        # Exact SI constitutive limit: R=1/4 ohm, I=8 A, V=2 V.
        # No stored-energy change or capacitor current is possible at this point.
        prepare_nonlinear_device_event_substep!(runtime, duration, duration, method)
        current, jacobian, low = zeros(2), zeros(2, 2), zeros(2)
        nonlinear_current_jacobian_compensated!(current, jacobian, low, runtime,
            [2.0, 0.0], [0.0, 0.0], duration)
        @test only(runtime.candidate.coil_current_a) == BigFloat(8)
        @test only(runtime.candidate.capacitor_current_a) == BigFloat(0)
        @test current == [8.0, -8.0]
        @test all(iszero, low)
        @test runtime.accepted_state.accepted_step_count == 0
    end
end

@testset "affine transformer current increments obey the independent R-L balance" begin
    for bctran in (false, true), method in (:backward_euler, :trapezoidal)
        duration = 4.126e-10
        runtime = compensated_history_transformer_runtime(; bctran)
        runtime.accepted_state.coil_current_a[1] = 8.0
        runtime.accepted_state.coil_voltage_v[1] = 2.0
        voltage, voltage_low = [2.0, 0.0], [2.0^-70, 0.0]
        prepare_nonlinear_device_event_substep!(runtime, duration, duration, method)
        current, jacobian, low = zeros(2), zeros(2, 2), zeros(2)
        nonlinear_current_jacobian_compensated!(current, jacobian, low, runtime,
            voltage, voltage_low, duration)
        expected = setprecision(BigFloat, 256) do
            alpha = BigFloat((method === :backward_euler ? 1.0 : 2.0) / duration)
            BigFloat(8) + BigFloat(voltage_low[1]) / (BigFloat(0.25) + alpha * BigFloat(8e-3))
        end
        @test abs(BigFloat(only(runtime.candidate.coil_current_a)) - expected) <= BigFloat(2)^(-101)
        @test only(runtime.candidate.coil_current_a) > BigFloat(8)
        @test runtime.accepted_state.coil_current_a == [8.0]
    end
end
