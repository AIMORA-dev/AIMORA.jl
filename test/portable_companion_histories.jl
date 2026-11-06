using Test,AIMORA
if AIMORA.solver_available()
@testset "portable represented branch histories preserve precision and zero sign" begin
    emt=AIMORA.EMTStudy;branch=AIMORA.Branches.SeriesRLBranch(1,0,0.1,0.001)
    restored=AIMORA.Branches.SeriesRLBranch(1,0,0.1,0.001)
    values=BigFloat[BigFloat(-0.0;precision=106),BigFloat(0.0;precision=106),
        setprecision(BigFloat,106) do;BigFloat(1)+BigFloat(2)^(-90);end,
        setprecision(BigFloat,256) do;BigFloat(1)+BigFloat(2)^(-200);end]
    for value in values
        branch.i_prev=deepcopy(value);branch.i_last=deepcopy(value)
        fields=emt._portable_emt_element_state_fields(branch,1,:inductor)
        emt._restore_portable_emt_element_state!(restored,Dict(f.identity=>f for f in fields),1)
        @test isequal(restored.i_prev,value) && precision(restored.i_prev)==precision(value)
        @test isequal(restored.i_last,value) && precision(restored.i_last)==precision(value)
        encoded=emt._portable_dyadic_history_value(value)
        @test startswith(encoded,"represented_binary_v2:")
    end
    for value in (-0.0,0.0)
        branch.i_prev=branch.i_last=value
        fields=emt._portable_emt_element_state_fields(branch,1,:inductor)
        emt._restore_portable_emt_element_state!(restored,Dict(f.identity=>f for f in fields),1)
        @test restored.i_prev===value && restored.i_last===value
    end
    legacy=emt._portable_dyadic_history_decode("represented_binary_v1:106:1:2","inductor")
    @test legacy==BigFloat(0.5) && precision(legacy)==106
    legacy_zero=emt._portable_dyadic_history_decode("represented_binary_v1:106:0:1","inductor")
    @test iszero(legacy_zero) && !signbit(legacy_zero)
    for invalid in ("represented_binary_v2:106:1:1:1","represented_binary_v2:106:0:1:2",
        "represented_binary_v2:106:01:1:0","represented_binary_v2:106:1:3:0",
        "represented_binary_v2:0:0:1:1")
        @test_throws AIMORA.PortableSnapshots.PortableSnapshotFailure emt._portable_dyadic_history_decode(invalid,"inductor")
    end
end
end
