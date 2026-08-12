using AIMORA.EMTTaskPlatform
using SHA

function task_spec(
    name,
    family,
    period;
    phase = 0 // 1,
    delay = 0 // 1,
    priority = 0,
    reads = String[],
    writes = String[],
    predecessors = String[],
    effects = EMTTaskEffect[],
)
    return EMTTaskSpec(
        name,
        family,
        emt_logical_time(0),
        emt_logical_time(period),
        emt_logical_time(phase),
        emt_logical_time(delay);
        priority,
        read_resources = reads,
        write_resources = writes,
        predecessors,
        effects,
    )
end

@testset "public exact multirate EMT task plan" begin
    families = EMTTaskFamily[
        ProtectionEMTTask,
        CarrierEMTTask,
        ConverterControlEMTTask,
        MechanicalEMTTask,
        SourceEMTTask,
        ThermalEMTTask,
        InterfaceEMTTask,
        UserDefinedEMTTask,
    ]
    specifications = EMTTaskSpec[
        task_spec(
            "task_$(index)",
            family,
            index // 1_000_000;
            phase = (index - 1) // 10_000_000,
            delay = index // 100_000_000,
            priority = index - 4,
            reads = ["input_$index"],
            writes = ["output_$index"],
            effects = index == 1 ?
                [InvalidateEMTPowerHistory, InvalidateEMTOutput] : EMTTaskEffect[],
        ) for (index, family) in pairs(families)
    ]
    plan = emt_task_plan(
        specifications;
        start = emt_logical_time(0),
        stop = emt_logical_time(1 // 100),
    )
    @test plan.quantum == emt_logical_time(1 // 100_000_000)
    @test plan.horizon_ticks == 1_000_000
    @test length(plan.entries) == 8
    @test Set(entry.spec.family for entry in plan.entries) == Set(families)
    @test collect(plan.execution_order) == collect(1:8)
    @test occursin(r"^[0-9a-f]{64}$", plan.signature_sha256)
    @test plan.signature_sha256 == emt_task_plan(
        specifications;
        start = emt_logical_time(0),
        stop = emt_logical_time(1 // 100),
    ).signature_sha256
    @test_throws ArgumentError emt_logical_time(1.0e-6)
end

@testset "task dependency order and exact collision refusal" begin
    writer = task_spec(
        "writer",
        ConverterControlEMTTask,
        10 // 1_000_000;
        writes = ["shared_command"],
        priority = 4,
    )
    reader = task_spec(
        "reader",
        ProtectionEMTTask,
        20 // 1_000_000;
        reads = ["shared_command"],
        predecessors = ["writer"],
        priority = -4,
    )
    plan = emt_task_plan(
        [reader, writer];
        start = emt_logical_time(0),
        stop = emt_logical_time(100 // 1_000_000),
    )
    @test [plan.entries[index].spec for index in plan.execution_order] == [writer, reader]

    unordered = task_spec(
        "reader",
        ProtectionEMTTask,
        20 // 1_000_000;
        reads = ["shared_command"],
    )
    conflict = try
        emt_task_plan(
            [writer, unordered];
            start = emt_logical_time(0),
            stop = emt_logical_time(100 // 1_000_000),
        )
        nothing
    catch error
        error
    end
    @test conflict isa EMTTaskPlatformFailure
    @test conflict.code == :unordered_task_resource_conflict

    cyclic_writer = task_spec(
        "writer",
        ConverterControlEMTTask,
        10 // 1_000_000;
        writes = ["shared_command"],
        predecessors = ["reader"],
    )
    cycle = try
        emt_task_plan(
            [cyclic_writer, reader];
            start = emt_logical_time(0),
            stop = emt_logical_time(100 // 1_000_000),
        )
        nothing
    catch error
        error
    end
    @test cycle isa EMTTaskPlatformFailure
    @test cycle.code == :cyclic_task_dependencies
end

@testset "typed occurrence checkpoint result and task interface" begin
    specification = task_spec("source", SourceEMTTask, 1 // 1_000_000)
    plan = emt_task_plan(
        [specification];
        start = emt_logical_time(0),
        stop = emt_logical_time(10 // 1_000_000),
    )
    occurrence = EMTTaskOccurrence(
        "source",
        SourceEMTTask,
        emt_logical_time(2 // 1_000_000),
        2.0e-6,
        EMTTaskWriteStage,
        0,
        3,
        3,
        3,
        3,
    )
    state_signature = bytes2hex(SHA.sha256("state"))
    checkpoint = EMTTaskCheckpoint(
        plan.signature_sha256,
        occurrence.exact_instant,
        (activation_index = 3, held_output = 2.0),
        state_signature,
    )
    signature = emt_task_result_signature(plan, [occurrence], state_signature)
    result = EMTTaskResult(
        true,
        plan.signature_sha256,
        [occurrence];
        checkpoint,
        task_counts = ["source" => 3],
        maximum_pending_depth = 1,
        effects = [InvalidateEMTOutput],
        deterministic_signature_sha256 = signature,
    )
    @test result.accepted
    @test result.checkpoint === checkpoint
    @test result.deterministic_signature_sha256 == signature

    struct MissingTaskInterface <: AbstractEMTTask end
    task = MissingTaskInterface()
    @test_throws MethodError emt_task_name(task)
    @test_throws MethodError emt_task_checkpoint(task)
end

@testset "exact calendar metamorphisms and scale refusal" begin
    normalized = task_spec(
        "normalized",
        SourceEMTTask,
        2 // 6;
        phase = 2 // 12,
        delay = 4 // 24,
    )
    canonical = task_spec(
        "normalized",
        SourceEMTTask,
        1 // 3;
        phase = 1 // 6,
        delay = 1 // 6,
    )
    normalized_plan = emt_task_plan(
        [normalized];
        start = emt_logical_time(0),
        stop = emt_logical_time(4 // 3),
    )
    canonical_plan = emt_task_plan(
        [canonical];
        start = emt_logical_time(0),
        stop = emt_logical_time(4 // 3),
    )
    @test normalized_plan.signature_sha256 == canonical_plan.signature_sha256
    @test normalized_plan.quantum == canonical_plan.quantum == emt_logical_time(1 // 6)

    base = task_spec(
        "translated",
        ThermalEMTTask,
        7 // 1_000_000;
        phase = 2 // 1_000_000,
        delay = 3 // 1_000_000,
    )
    translated = EMTTaskSpec(
        base.name,
        base.family,
        emt_logical_time(11 // 1_000_000),
        base.period,
        base.phase,
        base.computational_delay,
    )
    base_plan = emt_task_plan(
        [base];
        start = emt_logical_time(0),
        stop = emt_logical_time(100 // 1_000_000),
    )
    translated_plan = emt_task_plan(
        [translated];
        start = emt_logical_time(11 // 1_000_000),
        stop = emt_logical_time(111 // 1_000_000),
    )
    @test base_plan.horizon_ticks == translated_plan.horizon_ticks
    @test base_plan.quantum == translated_plan.quantum
    @test base_plan.entries[1].first_activation_tick ==
        translated_plan.entries[1].first_activation_tick
    @test base_plan.entries[1].period_ticks == translated_plan.entries[1].period_ticks
    @test base_plan.entries[1].delay_ticks == translated_plan.entries[1].delay_ticks

    scaled = EMTTaskSpec(
        base.name,
        base.family,
        base.epoch,
        5 * base.period,
        5 * base.phase,
        5 * base.computational_delay,
    )
    scaled_plan = emt_task_plan(
        [scaled];
        start = emt_logical_time(0),
        stop = emt_logical_time(500 // 1_000_000),
    )
    @test scaled_plan.horizon_ticks == base_plan.horizon_ticks
    @test scaled_plan.quantum == 5 * base_plan.quantum
    @test scaled_plan.entries[1].first_activation_tick == base_plan.entries[1].first_activation_tick
    @test scaled_plan.entries[1].period_ticks == base_plan.entries[1].period_ticks
    @test scaled_plan.entries[1].delay_ticks == base_plan.entries[1].delay_ticks

    million = task_spec("million", CarrierEMTTask, 1 // 1_000_000_000)
    million_plan = emt_task_plan(
        [million];
        start = emt_logical_time(0),
        stop = emt_logical_time(999_999 // 1_000_000_000),
    )
    @test million_plan.entries[1].first_activation_tick == 0
    @test million_plan.entries[1].period_ticks == 1
    @test million_plan.entries[1].first_activation_tick +
        999_999 * million_plan.entries[1].period_ticks == million_plan.horizon_ticks

    activation_limit = try
        emt_task_plan(
            [million];
            start = emt_logical_time(0),
            stop = emt_logical_time(1_000_000 // 1_000_000_000),
        )
        nothing
    catch error
        error
    end
    @test activation_limit isa EMTTaskPlatformFailure
    @test activation_limit.code == :task_activation_limit_exceeded

    huge_denominators = (
        big(2)^120,
        big(3)^70,
    )
    overflowing_specs = EMTTaskSpec[
        task_spec(
            "overflow_$index",
            UserDefinedEMTTask,
            (denominator - 1) // denominator,
        ) for (index, denominator) in pairs(huge_denominators)
    ]
    calendar_overflow = try
        emt_task_plan(
            overflowing_specs;
            start = emt_logical_time(0),
            stop = emt_logical_time(1),
        )
        nothing
    catch error
        error
    end
    @test calendar_overflow isa EMTTaskPlatformFailure
    @test calendar_overflow.code == :logical_calendar_overflow

    horizon_limit = try
        emt_task_plan(
            [million];
            start = emt_logical_time(0),
            stop = emt_logical_time(2),
        )
        nothing
    catch error
        error
    end
    @test horizon_limit isa EMTTaskPlatformFailure
    @test horizon_limit.code == :task_horizon_limit_exceeded
end
