using Test, AIMORA, LinearAlgebra
using AIMORA.NonlinearNetwork

@testset "affine transformer companion preserves small current beside large history" begin
    transformer=AIMORA.TransformerApparatus
    provenance=AIMORA.StudyCore.ParameterProvenance(
        "Analytic scalar RL companion cancellation fixture","SI",
        "exact supplied scalar parameters","synthetic algebraic precision limit",
        "public transformer numerical contract",AIMORA.StudyCore.PhysicalModelParameter)
    source=transformer.TransformerSourceRecord(:affine_history,repeat("a",64),provenance)
    connection=transformer.TransformerConnectionTopology(
        node_order=[:line,:neutral],coil_order=[:coil],winding_order=[:winding],
        phase_order=[:phase_a],coil_winding=[:winding],coil_phase=[:phase_a],
        incidence=reshape([1.0,-1.0],2,1),vector_group="I0",clock_number=0,phase_shift_rad=0.0)
    resistance,inductance,duration=0.25,8e-3,1e-5
    matrices=transformer.TransformerTerminalMatrices(
        reshape([resistance],1,1),reshape([inductance],1,1))
    for bctran in (false,true),method in (:trapezoidal,:backward_euler)
        model=bctran ? transformer.BCTRANTransformerModel(matrices;
            positive_pair_reconstruction_residual=0.0,zero_pair_reconstruction_residual=0.0,
            inverse_reconstruction_residual=0.0) : transformer.LowFrequencyTransformerModel(matrices)
        tier=bctran ? transformer.BCTRANTerminalTier : transformer.LowFrequencyTerminalTier
        specification=transformer.TransformerApparatusSpecification(
            :affine_companion,tier,connection,model,
            transformer.TransformerRuntimeSettings(timestep_s=duration,initialization_frequency_hz=60.0);
            phase_count=1,rated_power_va=1e6,rated_voltage_v=13.8e3,rated_frequency_hz=60.0,
            sources=[source],uncertainty="exact synthetic parameters",
            validity_domain="scalar affine RL numerical cancellation, not equipment ratings")
        runtime=transformer.transformer_apparatus_runtime(
            transformer.prepare_transformer_apparatus(specification),[1,2])
        incoming_current,incoming_voltage=3.23456789e7,777.12
        runtime.accepted_state.coil_current_a[1]=incoming_current
        runtime.accepted_state.coil_voltage_v[1]=incoming_voltage
        prepare_nonlinear_device_step!(runtime,duration,duration,method)
        before=deepcopy(runtime.accepted_state)
        high,low,expected=setprecision(BigFloat,192) do
            alpha=BigFloat((method===:backward_euler ? 1.0 : 2.0)/duration)
            history=alpha*BigFloat(inductance)*BigFloat(incoming_current)
            method===:backward_euler || (history-=BigFloat(resistance)*BigFloat(incoming_current)-BigFloat(incoming_voltage))
            impedance=BigFloat(resistance)+alpha*BigFloat(inductance)
            voltage=-history+big"1e-6"*impedance
            high=Float64(voltage);low=Float64(voltage-BigFloat(high))
            expected=Float64((BigFloat(high)+BigFloat(low)+history)/impedance)
            high,low,expected
        end
        @test high+low==high
        current=zeros(2);jacobian=zeros(2,2)
        nonlinear_current_jacobian_compensated!(current,jacobian,runtime,[high,0.0],[low,0.0],duration)
        @test isapprox(current[1],expected;atol=1e-16,rtol=1e-12)
        @test current[2]==-current[1]
        @test sum(current)==0.0
        @test runtime.accepted_state.coil_current_a==before.coil_current_a
        @test runtime.accepted_state.coil_voltage_v==before.coil_voltage_v
        @test runtime.accepted_state.accepted_time_s==before.accepted_time_s
        @test runtime.candidate.coil_current_a[1]==current[1]
        current_correction = zeros(2)
        nonlinear_current_jacobian_compensated!(current,jacobian,current_correction,
            runtime,[high,0.0],[low,0.0],duration)
        precise_current = setprecision(BigFloat,212) do
            alpha = BigFloat((method===:backward_euler ? 1.0 : 2.0)/duration)
            history = alpha*BigFloat(inductance)*BigFloat(incoming_current)
            method===:backward_euler ||
                (history += BigFloat(incoming_voltage)-BigFloat(resistance)*BigFloat(incoming_current))
            (BigFloat(high)+BigFloat(low)+history)/(BigFloat(resistance)+alpha*BigFloat(inductance))
        end
        @test current_correction[1] != 0.0
        @test current_correction[2] == -current_correction[1]
        @test current[1]+current_correction[1] == current[1]
        @test setprecision(BigFloat,192) do
            abs(BigFloat(current[1])+BigFloat(current_correction[1])-precise_current) < big"1e-30"
        end
        first_current=copy(current);first_candidate=copy(runtime.candidate.coil_current_a)
        nonlinear_current_jacobian_compensated!(current,jacobian,current_correction,
            runtime,[high,0.0],[low,0.0],duration)
        @test current==first_current && runtime.candidate.coil_current_a==first_candidate
        nonlinear_current_jacobian_compensated!(current,jacobian,runtime,[high,0.0],[low,0.0],duration)
        @test current==first_current
        @test runtime.candidate.coil_current_a==Float64.(first_candidate)
    end
end
