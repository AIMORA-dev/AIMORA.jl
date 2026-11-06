using Test, AIMORA

@testset "capacitor history preserves a non-grounded constant physical voltage" begin
    branches = AIMORA.Branches
    for high in ([1.0,-0.1],[1000.0,-999.1]),
        method in (branches.TrapezoidalCompanion,branches.BackwardEulerCompanion)
        low = [eps(high[1])/8,eps(high[2])/16]
        branch = branches.CapacitorBranch(1,2,0.01)
        branch.v_prev = high[1]-high[2]
        branches.update!(branch,high,low,low,1e-10,Val(method);incoming_voltage=high)
        @test branch.i_prev == branch.i_last == 0.0
        @test branch.v_prev == high[1]-high[2]
        # An explicitly different initial component voltage is not overwritten
        # by the incoming nodal state or inferred to be a dependent history.
        branch.v_prev = 2.5
        branch.i_prev = 0.125
        conductance,_ = branches.companion(branch,1e-10,Val(method))
        expected = setprecision(BigFloat,192) do
            physical = BigFloat(high[1])+BigFloat(low[1])-BigFloat(high[2])-BigFloat(low[2])
            initial = BigFloat(2.5)+BigFloat(low[1])-BigFloat(low[2])
            Float64(BigFloat(conductance)*(physical-initial)-
                (method === branches.TrapezoidalCompanion ? BigFloat(0.125) : BigFloat(0)))
        end
        branches.update!(branch,high,low,low,1e-10,Val(method);incoming_voltage=high)
        @test branch.i_prev == branch.i_last == expected
    end
end

@testset "capacitor acceptance preserves the compensated KCL current" begin
    branches=AIMORA.Branches
    for method in (branches.TrapezoidalCompanion,branches.BackwardEulerCompanion),
        grounded in (true,false), duration in (1e-10,5e-6)
        high=[1000.0,999.0]
        low=[eps(1000.0)/8,-eps(999.0)/16]
        incoming=[-eps(1000.0)/16,eps(999.0)/32]
        branch=branches.CapacitorBranch(1,grounded ? 0 : 2,0.01)
        branch.v_prev=grounded ? high[1] : high[1]-high[2]
        branch.i_prev=0.375
        conductance,_=branches.companion(branch,duration,Val(method))
        expected=setprecision(BigFloat,192) do
            endpoint=BigFloat(high[1])+BigFloat(low[1])-
                (grounded ? BigFloat(0) : BigFloat(high[2])+BigFloat(low[2]))
            previous=BigFloat(branch.v_prev)+BigFloat(incoming[1])-
                (grounded ? BigFloat(0) : BigFloat(incoming[2]))
            Float64(BigFloat(conductance)*(endpoint-previous)-
                (method===branches.TrapezoidalCompanion ? BigFloat(branch.i_prev) : BigFloat(0)))
        end
        ordinary=deepcopy(branch)
        branches.update!(ordinary,high,duration,Val(method))
        @test ordinary.i_prev!=expected
        branches.update!(branch,high,low,incoming,duration,Val(method))
        @test branch.i_prev==branch.i_last==expected
        @test branch.v_prev==(grounded ? high[1] : high[1]-high[2])
        @test low==[eps(1000.0)/8,-eps(999.0)/16]
        @test incoming==[-eps(1000.0)/16,eps(999.0)/32]
        previous_current=branch.i_prev
        branches.update!(branch,high,low,low,duration,Val(method))
        @test branch.i_prev==(method===branches.TrapezoidalCompanion ? -previous_current : 0.0)
    end
end

@testset "R-L validation retains ordinary and compensated histories" begin
    B=AIMORA.Branches
    for compensated in (false,true)
        branch=B.SeriesRLBranch(1,2,2.0,0.0)
        retained=setprecision(BigFloat,256) do
            BigFloat(3)+BigFloat(2)^(-100)
        end
        branch.i_prev=compensated ? retained : 3.0
        branch.i_last=branch.i_prev
        original=branch.i_prev
        @test B.companion(branch,1e-5) == (0.5,0.0)
        @test branch.i_prev === original && branch.i_last === original
        for value in (Inf,-Inf,NaN)
            branch.i_prev=compensated ? BigFloat(value) : value
            @test_throws ArgumentError B.companion(branch,1e-5)
        end
        branch.i_prev=original
        branch.i_last=compensated ? BigFloat(NaN) : NaN
        @test_throws ArgumentError B.companion(branch,1e-5)
    end
end
