using Test

@testset "two-level bridge conduction supports production tuple owners" begin
    nonlinear = AIMORA.Nonlinear
    controlled = ntuple(_ -> nonlinear.IGBTSwitch(1, 0), 6)
    diodes = ntuple(_ -> nonlinear.DiodeValveSwitch(1, 0), 6)
    diodes[6].last_voltage = -700.0
    @test !AIMORA.EMTStudy._synchronize_two_level_bridge_conduction!(
        controlled, diodes, (0.0, 0.0, 2.0), 0.0,
    )
    @test !diodes[6].closed
    diodes[6].last_voltage = nonlinear.power_semiconductor_forward_turn_on_voltage(diodes[6])
    @test AIMORA.EMTStudy._synchronize_two_level_bridge_conduction!(
        controlled, diodes, (0.0, 0.0, 2.0), 0.0,
    )
    @test diodes[6].closed
end

@testset "two-level bridge respects incoming natural turn-on guards" begin
    nonlinear=AIMORA.Nonlinear
    emt=AIMORA.EMTStudy
    controlled=[nonlinear.IGBTSwitch(1,0) for _ in 1:6]
    diodes=[nonlinear.DiodeValveSwitch(1,0) for _ in 1:6]
    diode=diodes[6]
    diode.last_voltage=-700.0
    initial_transitions=getfield.(vcat(controlled,diodes),:topology_transition_count)
    @test nonlinear.power_semiconductor_forward_turn_on_residual(diode)<0
    @test !emt._synchronize_two_level_bridge_conduction!(controlled,diodes,(0.0,0.0,2.0),0.0)
    @test !diode.closed
    @test getfield.(vcat(controlled,diodes),:topology_transition_count)==initial_transitions

    # An admissible future trial is not permission to pre-close a reverse-
    # biased accepted owner. Junction displacement current can carry the leg
    # until the actual forward-voltage event; phase-current sign is not that event.
    predicted=deepcopy(vcat(controlled,diodes))
    nonlinear.apply_power_semiconductor_forward_turn_on!(predicted[12],5e-6)
    predicted[12].last_forward_current_a=0.25
    @test !emt._synchronize_two_level_bridge_conduction!(controlled,diodes,(0.0,0.0,2.0),0.0;
        mapped_devices=predicted)
    @test !diode.closed
    @test getfield.(vcat(controlled,diodes),:topology_transition_count)==initial_transitions

    # Preserve admission on the exact physical root and inside its forward side;
    # no extra margin or tolerance is introduced.
    diode.last_voltage=nonlinear.power_semiconductor_forward_turn_on_voltage(diode)
    @test nonlinear.power_semiconductor_forward_turn_on_residual(diode)==0.0
    @test emt._synchronize_two_level_bridge_conduction!(controlled,diodes,(0.0,0.0,2.0),0.0)
    @test diode.closed
    @test diode.topology_transition_count==1
    @test !emt._synchronize_two_level_bridge_conduction!(controlled,diodes,(0.0,0.0,2.0),0.0)
    @test diode.topology_transition_count==1
end

@testset "two-level bridge retains gate-owned bidirectional channels" begin
    nonlinear=AIMORA.Nonlinear
    controlled=[nonlinear.MOSFETSwitch(1,0) for _ in 1:6]
    diodes=[nonlinear.DiodeValveSwitch(1,0) for _ in 1:6]
    controlled[1].gate_driver.applied_on=true
    controlled[1].last_voltage=-2.0
    @test nonlinear.power_semiconductor_forward_turn_on_residual(controlled[1])===nothing
    @test AIMORA.EMTStudy._synchronize_two_level_bridge_conduction!(controlled,diodes,(1.0,0.0,0.0),0.0)
    @test controlled[1].closed
    @test controlled[1].topology_transition_count==1
    @test all(!d.closed for d in diodes)
end
