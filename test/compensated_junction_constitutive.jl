using Test,AIMORA
using AIMORA.NonlinearNetwork

@testset "compensated junction acceptance retains charge beside a large channel current" begin
    nonlinear=AIMORA.Nonlinear
    for method in (:backward_euler,:trapezoidal), low in (0.0,2.0^-70)
        junction=nonlinear.NonlinearJunctionChargeFidelity(1e-6,1.0,0.5;
            initial_voltage_v=0.1,voltage_domain_v=(-10.0,10.0))
        device=nonlinear.DiodeValveSwitch(1,0;initially_closed=true,on_conductance=1e9,
            extended_fidelity=nonlinear.PowerSemiconductorExtendedFidelity(junction_charge=junction))
        step=5e-6
        voltage=[nextfloat(0.1),0.0];correction=[low,0.0]
        prepare_nonlinear_device_step!(device,step,step,method)
        previous_charge=junction.previous_charge_c
        previous_current=junction.last_displacement_current_a
        expected=setprecision(BigFloat,256) do
            # Independent positive-voltage limit: Q=C*V, C constant.
            charge=BigFloat(junction.zero_bias_capacitance_f)*(BigFloat(voltage[1])+BigFloat(low))
            current=(method===:trapezoidal ? 2 : 1)*(charge-BigFloat(previous_charge))/BigFloat(step)
            method===:trapezoidal && (current-=BigFloat(previous_current))
            (;charge=Float64(charge),current=Float64(current),raw_current=current,
                total=BigFloat(device.on_conductance)*(BigFloat(voltage[1])+BigFloat(low))+current,
                capacitance=junction.zero_bias_capacitance_f)
        end
        current,current_low,jacobian=zeros(2),zeros(2),zeros(2,2)
        nonlinear_current_jacobian_compensated!(current,jacobian,current_low,device,voltage,correction,step)
        @test abs(expected.current)<eps(current[1])
        @test current[1]+current_low[1]==current[1]
        @test current[2]==-current[1] && current_low[2]==-current_low[1]
        @test setprecision(BigFloat,106) do
            abs(BigFloat(current[1])+BigFloat(current_low[1])-expected.total)<=8eps(BigFloat(expected.total))
        end
        @test junction.previous_charge_c==previous_charge
        accept_nonlinear_device_compensated_state!(device,voltage,correction,current,current_low,jacobian,step)
        @test junction.previous_charge_c==expected.charge
        @test junction.last_charge_c==expected.charge
        @test junction.last_displacement_current_a ≈ expected.current rtol=1e-12 atol=1e-30
        @test junction.last_capacitance_f==expected.capacitance
    end
end

@testset "compensated semiconductor rejects unrepresentable terminal output" begin
    nonlinear=AIMORA.Nonlinear
    junction=nonlinear.NonlinearJunctionChargeFidelity(1e-6,1.0,0.5;
        voltage_domain_v=(-10.0,10.0))
    device=nonlinear.DiodeValveSwitch(1,0;initially_closed=true,on_conductance=1e308,
        extended_fidelity=nonlinear.PowerSemiconductorExtendedFidelity(junction_charge=junction))
    prepare_nonlinear_device_step!(device,5e-6,5e-6,:backward_euler)
    current,low,jacobian=zeros(2),zeros(2),zeros(2,2)
    @test_throws ArgumentError nonlinear_current_jacobian_compensated!(current,jacobian,low,
        device,[2.0,0.0],zeros(2),5e-6)
    @test all(iszero,current) && all(iszero,low) && all(iszero,jacobian)
end

@testset "compensated junction constitutive charge is independent of accepted-history rounding" begin
    nonlinear=AIMORA.Nonlinear
    setprecision(BigFloat,106) do
        junction=nonlinear.NonlinearJunctionChargeFidelity(1e-6,1.0,0.5;
            initial_voltage_v=0.1,voltage_domain_v=(-10.0,10.0))
        previous=(junction.previous_voltage_v,junction.previous_charge_c,
            junction.last_capacitance_f,junction.last_charge_c)
        voltage=BigFloat(junction.previous_voltage_v)
        capacitance=BigFloat(junction.zero_bias_capacitance_f)
        charge,derivative=nonlinear._power_semiconductor_compensated_junction_charge_capacitance(junction,voltage)
        @test isequal(charge,capacitance*voltage)
        @test isequal(derivative,capacitance)
        delta=BigFloat(eps(junction.previous_voltage_v))/1024
        @test Float64(voltage+delta)==Float64(voltage-delta)==junction.previous_voltage_v
        for displacement in (-delta,delta)
            mapped,mapped_derivative=nonlinear._power_semiconductor_compensated_junction_charge_capacitance(
                junction,voltage+displacement)
            @test isequal(mapped,capacitance*(voltage+displacement))
            @test isequal(mapped_derivative,capacitance)
            # A linear capacitance has no Taylor remainder except arithmetic
            # rounding at the declared106-bit constitutive precision.
            @test abs(mapped-charge-capacitance*displacement)<=8eps(capacitance*voltage)
        end
        @test previous==(junction.previous_voltage_v,junction.previous_charge_c,
            junction.last_capacitance_f,junction.last_charge_c)
        @test_throws DomainError nonlinear._power_semiconductor_compensated_junction_charge_capacitance(
            junction,BigFloat(11))
    end
end

@testset "compensated graded junction agrees with its independent constitutive integral" begin
    nonlinear=AIMORA.Nonlinear
    for initial in (-0.1,-1.0,-4.0),grading in (0.0,0.25,0.5)
        junction=nonlinear.NonlinearJunctionChargeFidelity(1e-6,1.0,grading;
            initial_voltage_v=initial,voltage_domain_v=(-10.0,10.0))
        setprecision(BigFloat,106) do
            voltage=BigFloat(initial)
            charge,derivative=nonlinear._power_semiconductor_compensated_junction_charge_capacitance(junction,voltage)
            expected=setprecision(BigFloat,256) do
                c,vj,m=BigFloat(junction.zero_bias_capacitance_f),BigFloat(junction.junction_voltage_v),BigFloat(grading)
                normalized=1-BigFloat(voltage)/vj
                (;charge=-c*vj*(normalized^(1-m)-1)/(1-m),capacitance=c*normalized^(-m))
            end
            @test abs(charge-expected.charge)<=32eps(BigFloat(expected.charge))
            @test abs(derivative-expected.capacitance)<=32eps(BigFloat(expected.capacitance))
            @test derivative>0
        end
    end
end
