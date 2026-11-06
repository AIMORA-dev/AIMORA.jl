using Test,AIMORA,LinearAlgebra
using AIMORA.NonlinearNetwork

function compensated_history_transformer_runtime(;bctran=false,duration=1e-5)
    transformer = AIMORA.TransformerApparatus
    provenance = AIMORA.StudyCore.ParameterProvenance(
        "Analytic scalar RLC accepted-history fixture","SI","explicit scalar coefficients",
        "exact synthetic coefficients","public companion-history numerical contract",
        AIMORA.StudyCore.PhysicalModelParameter)
    source = transformer.TransformerSourceRecord(:companion_history,repeat("b",64),provenance)
    connection = transformer.TransformerConnectionTopology(
        node_order=[:line,:neutral],coil_order=[:coil],winding_order=[:winding],
        phase_order=[:phase_a],coil_winding=[:winding],coil_phase=[:phase_a],
        incidence=reshape([1.0,-1.0],2,1),vector_group="I0",clock_number=0,phase_shift_rad=0.0)
    matrices = transformer.TransformerTerminalMatrices(reshape([0.25],1,1),reshape([8e-3],1,1);
        capacitance_f=reshape([1e-6],1,1))
    model = bctran ? transformer.BCTRANTransformerModel(matrices;
        positive_pair_reconstruction_residual=0.0,zero_pair_reconstruction_residual=0.0,
        inverse_reconstruction_residual=0.0) : transformer.LowFrequencyTransformerModel(matrices)
    specification = transformer.TransformerApparatusSpecification(:companion_history,
        bctran ? transformer.BCTRANTerminalTier : transformer.LowFrequencyTerminalTier,
        connection,model,transformer.TransformerRuntimeSettings(timestep_s=duration,initialization_frequency_hz=60.0);
        phase_count=1,rated_power_va=1e6,rated_voltage_v=13.8e3,rated_frequency_hz=60.0,
        sources=[source],uncertainty="exact synthetic coefficients",
        validity_domain="scalar accepted RLC history and deterministic snapshot, not equipment ratings")
    transformer.transformer_apparatus_runtime(transformer.prepare_transformer_apparatus(specification),[1,2])
end

@testset "compensated transformer acceptance retains its physical companion histories" begin
    transformer = AIMORA.TransformerApparatus
    duration = 1e-5
    for bctran in (false,true),method in (:trapezoidal,:backward_euler)
        runtime = compensated_history_transformer_runtime(;bctran,duration)
        legacy_snapshot = transformer.transformer_apparatus_runtime_snapshot(runtime)
        @test legacy_snapshot.schema_version == 1
        voltage,low = [0.1,0.0],[2.0^-60,0.0]
        prepare_nonlinear_device_step!(runtime,duration,duration,method)
        before = deepcopy(runtime.accepted_state)
        current,jacobian,current_low = zeros(2),zeros(2,2),zeros(2)
        nonlinear_current_jacobian_compensated!(current,jacobian,current_low,runtime,voltage,low,duration)
        @test runtime.accepted_state.coil_current_a == before.coil_current_a
        expected = setprecision(BigFloat,106) do
            v = BigFloat(voltage[1])+BigFloat(low[1])
            rate = BigFloat((method===:backward_euler ? 1.0 : 2.0)/duration)
            (;voltage=v,current=v/(BigFloat(0.25)+rate*BigFloat(8e-3)),
                capacitor=rate*(BigFloat(1e-6)*v))
        end
        accept_nonlinear_device_compensated_state!(runtime,voltage,low,current,jacobian,duration)
        @test isequal(only(runtime.accepted_state.coil_current_a),expected.current)
        @test isequal(only(runtime.accepted_state.coil_voltage_v),expected.voltage)
        @test isequal(only(runtime.accepted_state.capacitor_current_a),expected.capacitor)
        @test runtime.accepted_state.accepted_step_count == 1
        finish_nonlinear_device_step!(runtime)
        snapshot = transformer.transformer_apparatus_runtime_snapshot(runtime)
        @test snapshot.schema_version == 2
        restored = compensated_history_transformer_runtime(;bctran,duration)
        transformer.restore_transformer_apparatus_runtime_snapshot!(restored,snapshot)
        @test transformer.transformer_apparatus_runtime_snapshot(restored).deterministic_signature_sha256 ==
            snapshot.deterministic_signature_sha256
        for name in (:coil_current_a,:capacitor_current_a,:coil_voltage_v)
            @test isequal(getfield(restored.accepted_state,name),getfield(runtime.accepted_state,name))
            @test getfield(restored.accepted_state,name) isa Vector{BigFloat}
            @test all(value->precision(value)==106,getfield(restored.accepted_state,name))
        end
        invalid_state = deepcopy(snapshot.state)
        invalid_state.extended_history_corrections[:coil_current_a][1] = eps(invalid_state.coil_current_a[1])
        invalid = transformer.TransformerApparatusRuntimeSnapshot(2,snapshot.preparation_signature_sha256,
            invalid_state,transformer._transformer_snapshot_signature(snapshot.preparation_signature_sha256,invalid_state))
        @test_throws ArgumentError transformer.restore_transformer_apparatus_runtime_snapshot!(restored,invalid)
        @test transformer.transformer_apparatus_runtime_snapshot(restored).deterministic_signature_sha256 ==
            snapshot.deterministic_signature_sha256
        transformer.restore_transformer_apparatus_runtime_snapshot!(restored,legacy_snapshot)
        @test transformer.transformer_apparatus_runtime_snapshot(restored).deterministic_signature_sha256 ==
            legacy_snapshot.deterministic_signature_sha256
        @test all(name->getfield(restored.accepted_state,name) isa Vector{Float64},
            (:coil_current_a,:capacitor_current_a,:coil_voltage_v))
    end
end

@testset "compensated scalar snapshot pairs reject loss and preserve signed zero" begin
    core = AIMORA.StudyCore
    for value in setprecision(BigFloat,106) do
        [BigFloat(-0.0),BigFloat(0.0),BigFloat(1)+BigFloat(2)^(-100),BigFloat(-51)+BigFloat(2)^(-70)]
    end
        pair = core._encode_compensated_scalar_pair(value)
        @test pair.status === :ok
        restored = core._decode_compensated_scalar_pair(pair.high,pair.low)
        @test restored.status === :ok && isequal(restored.value,value)
        @test precision(restored.value) == 106
    end
    @test core._decode_compensated_scalar_pair(1.0,1.0).status === :noncanonical
    @test core._decode_compensated_scalar_pair(1.0,2.0^-300).status === :precision_loss
    @test core._decode_compensated_scalar_pair(Inf,0.0).status === :nonfinite
    unrepresentable = setprecision(BigFloat,256) do
        BigFloat(1)+BigFloat(2)^(-200)
    end
    @test core._encode_compensated_scalar_pair(unrepresentable).status === :precision_loss
end
