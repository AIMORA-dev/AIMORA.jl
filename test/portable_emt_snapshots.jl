using SHA
using Test
using AIMORA
using Printf
using AIMORA.EMTTaskPlatform

function portable_source_function_fixed_field_line(
    source_type::Integer,
    node::AbstractString,
)
    return @sprintf(
        "%2d%-6s%2d%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f",
        source_type,
        node,
        0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        2.0,
    )
end

function portable_source_function_deck()
    return AIMORA.DeckParser.parse_deck_lines(
        [
            "BEGIN NEW DATA CASE",
            "BLANK SOURCE",
            portable_source_function_fixed_field_line(1, "N1"),
            portable_source_function_fixed_field_line(2, "N2"),
            portable_source_function_fixed_field_line(3, "N3"),
            "BLANK BRANCH",
            "0,N1,0,1.0,0.0,0.0",
            "0,N2,0,1.0,0.0,0.0",
            "0,N3,0,1.0,0.0,0.0",
            "bergeron_line SIGNAL_LINE N1 LOAD 10.0 0.5 0.25 0.9",
            "conductance LOAD_SHUNT LOAD 0 1.0",
            "over16_output N1_VOLTAGE N1",
            "over16_output N2_VOLTAGE N2",
            "over16_output N3_VOLTAGE N3",
            "END",
        ];
        source = "portable source-function restart test",
    )
end

function portable_source_signal_program(; first_slot_offset::Real=0.0)
    samples = [
        [10.0 + first_slot_offset, 0.0, 30.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [15.0 + first_slot_offset, 0.0, 35.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [20.0 + first_slot_offset, 0.0, 40.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ]
    table = AIMORA.EMTStudy.TabulatedSourceSignalProvider(
        [0.0, 0.5, 1.0],
        samples;
        extrapolation = :hold,
    )
    analytic = AIMORA.Sources.AnalyticSourceSignal(
        13,
        5.0,
        0.0,
        -2.0,
        0.0,
        nextfloat(1.0),
    )
    return AIMORA.EMTStudy.SourceSignalProgram(
        table,
        [
            AIMORA.EMTStudy.AnalyticSourceSlot(
                2,
                analytic;
                assignment = :replace,
                unit = :V,
            ),
        ],
    )
end

function portable_control_network_feedback_deck()
    return AIMORA.DeckParser.parse_deck_lines(
        [
            "BEGIN NEW DATA CASE",
            " 5.0E-05 3.0E-04",
            "   20000       1       1       1       1       0       0       0",
            "TACS HYBRID",
            " 1FILTER  +SENSE",
            "       1.0       0.0",
            "       1.0     0.001",
            " 0DRIVEN  +FILTER",
            " 0MODUL   +FILTER",
            "90SENSE",
            "91OBSIN",
            "93OBSOUT",
            "23PULSE          2.0    0.0002    0.0001                         0.0003",
            "24RAMP           3.0   0.00015       0.0                         0.0003",
            "33FILTERSENSE OBSIN OBSOUTPULSE RAMP",
            "BLANK CARD ENDING TACS DATA",
            "  SENSE                       1.0",
            "  DRIVEN                      1.0",
            "  MODUL                       1.0",
            "  SWOUT                       1.0",
            "  OBSOUT                      1.0",
            "BLANK CARD ENDING BRANCHES",
            "13DRIVENSWOUT                                                         FILTER  01",
            " 0OBSIN OBSOUT      -1.0       1.0       0.0       0.0CLOSED                  00",
            "BLANK CARD ENDING SWITCHES",
            "11SENSE         10.0       0.0       0.0       0.0       0.0       0.0       0.0",
            "11OBSIN          4.0       0.0       0.0       0.0       0.0       0.0       0.0",
            "17DRIVEN 0       1.0       2.0       0.0       0.0       0.0       0.0       1.0",
            "11DRIVEN         1.0       0.0       0.0       0.0       0.0       0.0       1.0",
            "17MODUL  0       2.0       3.0       0.0       0.0       0.0       0.0       1.0",
            "11MODUL          1.0       0.0       0.0       0.0       0.0       0.0       1.0",
            "BLANK CARD ENDING SOURCES",
            "  SENSE DRIVEN MODUL SWOUT OBSINOBSOUT",
            "BLANK CARD ENDING NODE VOLTAGE OUTPUT",
            "BLANK CARD ENDING PLOTS",
            "BEGIN NEW DATA CASE",
            "BLANK CARD ENDING THE CASE",
        ];
        source = "portable control-network feedback restart test",
    )
end

function portable_snapshot_test_metadata(; profile = :portable_full)
    return AIMORA.PortableSnapshots.PortableSnapshotMetadata(
        profile,
        bytes2hex(sha256("project")),
        bytes2hex(sha256("model")),
        bytes2hex(sha256("topology")),
        bytes2hex(sha256("settings")),
        Int128(7) // Int128(2_000_000),
        7,
        ["emt.fixed_step", "emt.events", "emt.tasks"],
        "AIMORA-authored synthetic portable snapshot test";
        writer_version = "AIMORA.jl/test",
        creator_platform = "portable-test-platform",
    )
end

@testset "portable EMT registered state inventory" begin
    state = AIMORA.PortableSnapshots.PortableSnapshotStateInventory([
        AIMORA.PortableSnapshots.PortableSnapshotStateField(
            "network.node_voltage",
            "emt.network",
            :algebraic,
            "V",
            ["node"],
            AIMORA.PortableSnapshots.portable_snapshot_array(
                Float64[1.0, -0.5];
                unit = "V",
                axes = ["node"],
            ),
        ),
        AIMORA.PortableSnapshots.PortableSnapshotStateField(
            "execution.accepted_step",
            "emt.execution",
            :discrete,
            "1",
            String[],
            7,
        ),
    ])
    @test state.scalar_count == 3
    @test getfield.(state.fields, :identity) ==
        ["execution.accepted_step", "network.node_voltage"]
    encoded = AIMORA.PortableSnapshots.portable_state_inventory_record(state)
    decoded = AIMORA.PortableSnapshots.portable_state_inventory(encoded)
    @test decoded.signature_sha256 == state.signature_sha256
    @test decoded.scalar_count == state.scalar_count
    @test AIMORA.PortableSnapshots._portable_value_bytes(
        AIMORA.PortableSnapshots.portable_state_inventory_record(decoded),
    ) == AIMORA.PortableSnapshots._portable_value_bytes(encoded)
    @test_throws AIMORA.PortableSnapshots.PortableSnapshotFailure AIMORA.PortableSnapshots.PortableSnapshotStateInventory([
        state.fields[1],
        state.fields[1],
    ])
end

function portable_snapshot_test_record()
    voltage = AIMORA.PortableSnapshots.portable_snapshot_array(
        reshape(Float64[1.0, -0.5, 0.25, 0.0], 2, 2);
        unit = "V",
        axes = ["node", "sample"],
    )
    return AIMORA.PortableSnapshots.PortableSnapshotRecord(
        "aimora.test.network_state.v1",
        Pair{String,Any}[
            "time_s" => (Int128(7) // Int128(2_000_000)),
            "voltage" => voltage,
            "step" => 7,
            "valid" => true,
            "mode" => "accepted",
        ],
    )
end

@testset "portable EMT snapshot canonical envelope" begin
    public_section = AIMORA.PortableSnapshots.PortableSnapshotSection(
        "network.accepted_state",
        1,
        0,
        :public,
        portable_snapshot_test_record(),
    )
    private_section = AIMORA.PortableSnapshots.PortableSnapshotSection(
        "backend.reconstruction_state",
        1,
        0,
        :private_reconstructible,
        AIMORA.PortableSnapshots.PortableSnapshotRecord(
            "aimora.test.backend_reconstruction.v1",
            Pair{String,Any}[
                "factorization_rebuild" => true,
                "history_signature" => bytes2hex(sha256("history")),
            ],
        ),
    )
    snapshot = AIMORA.PortableSnapshots.PortableEMTSnapshot(
        portable_snapshot_test_metadata(),
        [public_section, private_section],
    )
    reversed_snapshot = AIMORA.PortableSnapshots.PortableEMTSnapshot(
        portable_snapshot_test_metadata(),
        [private_section, public_section],
    )
    bytes = AIMORA.PortableSnapshots.portable_snapshot_bytes(snapshot)
    @test bytes == AIMORA.PortableSnapshots.portable_snapshot_bytes(reversed_snapshot)
    descriptor = AIMORA.PortableSnapshots.portable_snapshot_descriptor(snapshot)
    @test descriptor.schema_major == 1
    @test descriptor.schema_minor == 1
    @test descriptor.canonical_bytes == length(bytes)
    @test descriptor.content_sha256 == bytes2hex(sha256(bytes[1:(end - 32)]))
    @test getfield.(descriptor.sections, :identity) ==
        ["backend.reconstruction_state", "network.accepted_state"]
    @test getfield.(descriptor.sections, :visibility) ==
        [:private_reconstructible, :public]
    @test isempty(descriptor.migration_path)
    @test descriptor.metadata.numeric_profile ==
        "ieee754_binary64_finite_preserve_signed_zero"
    @test descriptor.metadata.compression == "none"

    mktempdir() do directory
        first_path = joinpath(directory, "first.aimora-snapshot")
        second_path = joinpath(directory, "second.aimora-snapshot")
        written = AIMORA.PortableSnapshots.write_portable_emt_snapshot(first_path, snapshot)
        AIMORA.PortableSnapshots.write_portable_emt_snapshot(second_path, snapshot)
        @test written.content_sha256 == descriptor.content_sha256
        @test read(first_path) == read(second_path) == bytes
        inspected = AIMORA.PortableSnapshots.inspect_portable_emt_snapshot(first_path)
        @test inspected.content_sha256 == descriptor.content_sha256
        @test inspected.metadata.represented_time_s == Int128(7) // Int128(2_000_000)
        @test_throws AIMORA.PortableSnapshots.PortableSnapshotFailure AIMORA.PortableSnapshots.read_portable_emt_snapshot(first_path)
        restored = AIMORA.PortableSnapshots.read_portable_emt_snapshot(first_path; allow_private = true)
        @test AIMORA.PortableSnapshots.portable_snapshot_bytes(restored) == bytes
        network_record = only(section.value for section in restored.sections if section.identity == "network.accepted_state")
        fields = Dict(network_record.fields)
        @test fields["step"] == 7
        @test fields["time_s"] == Int128(7) // Int128(2_000_000)
        @test AIMORA.PortableSnapshots.portable_snapshot_array_values(fields["voltage"]) ==
            reshape(Float64[1.0, -0.5, 0.25, 0.0], 2, 2)

        corrupted = copy(bytes)
        corrupted[end - 40] = xor(corrupted[end - 40], 0x01)
        corrupted_path = joinpath(directory, "corrupted.aimora-snapshot")
        write(corrupted_path, corrupted)
        corruption = try
            AIMORA.PortableSnapshots.inspect_portable_emt_snapshot(corrupted_path)
            nothing
        catch error
            error
        end
        @test corruption isa AIMORA.PortableSnapshots.PortableSnapshotFailure
        @test corruption.code == :integrity
        @test_throws AIMORA.PortableSnapshots.PortableSnapshotFailure AIMORA.PortableSnapshots.inspect_portable_emt_snapshot(
            first_path;
            maximum_file_bytes = 16,
        )

        legacy_path = joinpath(directory, "legacy-minor-zero.aimora-snapshot")
        legacy_bytes = AIMORA.PortableSnapshots._portable_snapshot_bytes(
            snapshot;
            schema_minor = UInt16(0),
        )
        write(legacy_path, legacy_bytes)
        legacy_descriptor = AIMORA.PortableSnapshots.inspect_portable_emt_snapshot(
            legacy_path,
        )
        @test legacy_descriptor.schema_minor == 0
        @test legacy_descriptor.migration_path == ["aimora.portable_emt.1.0_to_1.1"]
        migrated = AIMORA.PortableSnapshots.read_portable_emt_snapshot(
            legacy_path;
            allow_private = true,
        )
        @test migrated.metadata.numeric_profile ==
            "ieee754_binary64_finite_preserve_signed_zero"
        @test migrated.metadata.compression == "none"
        @test AIMORA.PortableSnapshots.portable_snapshot_bytes(migrated) != legacy_bytes
    end
end

@testset "portable EMT snapshot public and refusal boundaries" begin
    public_section = AIMORA.PortableSnapshots.PortableSnapshotSection(
        "reference.accepted_state",
        1,
        0,
        :public,
        portable_snapshot_test_record(),
    )
    public_snapshot = AIMORA.PortableSnapshots.PortableEMTSnapshot(
        portable_snapshot_test_metadata(profile = :portable_public_reference),
        [public_section],
    )
    mktempdir() do directory
        path = joinpath(directory, "public.aimora-snapshot")
        AIMORA.PortableSnapshots.write_portable_emt_snapshot(path, public_snapshot)
        restored = AIMORA.PortableSnapshots.read_portable_emt_snapshot(path)
        @test restored.metadata.profile == :portable_public_reference
        @test only(restored.sections).visibility == :public
    end
    @test_throws AIMORA.PortableSnapshots.PortableSnapshotFailure AIMORA.PortableSnapshots.PortableEMTSnapshot(
        portable_snapshot_test_metadata(profile = :portable_public_reference),
        [AIMORA.PortableSnapshots.PortableSnapshotSection(
            "backend.reconstruction_state",
            1,
            0,
            :private_reconstructible,
            portable_snapshot_test_record(),
        )],
    )
    @test_throws AIMORA.PortableSnapshots.PortableSnapshotFailure AIMORA.PortableSnapshots.PortableSnapshotRecord(
        "aimora.test.duplicate.v1",
        Pair{String,Any}["state" => 1, "state" => 2],
    )
    @test_throws AIMORA.PortableSnapshots.PortableSnapshotFailure AIMORA.PortableSnapshots.PortableSnapshotSection(
        "Invalid Section",
        1,
        0,
        :public,
        nothing,
    )
    @test_throws AIMORA.PortableSnapshots.PortableSnapshotFailure AIMORA.PortableSnapshots.portable_snapshot_bytes(
        AIMORA.PortableSnapshots.PortableEMTSnapshot(
            portable_snapshot_test_metadata(),
            [AIMORA.PortableSnapshots.PortableSnapshotSection(
                "network.invalid_state",
                1,
                0,
                :public,
                NaN,
            )],
        ),
    )
    unsigned = AIMORA.PortableSnapshots.portable_snapshot_array(
        UInt64[0, typemax(UInt64)];
        unit = "1",
        axes = ["sample"],
    )
    signed = AIMORA.PortableSnapshots.portable_snapshot_array(
        Int64[typemin(Int64), 0, typemax(Int64)];
        unit = "1",
        axes = ["sample"],
    )
    @test AIMORA.PortableSnapshots.portable_snapshot_array_values(unsigned) ==
        UInt64[0, typemax(UInt64)]
    @test AIMORA.PortableSnapshots.portable_snapshot_array_values(signed) ==
        Int64[typemin(Int64), 0, typemax(Int64)]
    nonfinite_bytes = IOBuffer()
    AIMORA.PortableSnapshots._write_checkpoint_integer(
        nonfinite_bytes,
        reinterpret(UInt64, Inf),
    )
    @test_throws AIMORA.PortableSnapshots.PortableSnapshotFailure AIMORA.PortableSnapshots.PortableSnapshotArray(
        :float64,
        [1],
        "V",
        ["sample"],
        take!(nonfinite_bytes),
    )
    invalid_utf8 = IOBuffer(UInt8[0x01, 0x00, 0x00, 0x00, 0xff])
    @test_throws AIMORA.PortableSnapshots.PortableSnapshotFailure AIMORA.PortableSnapshots._read_portable_string(
        invalid_utf8,
        "test text";
        maximum_bytes = 16,
    )
end

if AIMORA.solver_available()
    function portable_element_owner_roundtrip(source, candidate, name::Symbol)
        captured = AIMORA.PortableSnapshots.PortableSnapshotStateInventory(
            AIMORA.EMTStudy._portable_emt_element_state_fields(source, 1, name),
        )
        fields = Dict(field.identity => field for field in captured.fields)
        AIMORA.EMTStudy._restore_portable_emt_element_state!(candidate, fields, 1)
        restored = AIMORA.PortableSnapshots.PortableSnapshotStateInventory(
            AIMORA.EMTStudy._portable_emt_element_state_fields(candidate, 1, name),
        )
        return captured, restored
    end

    @testset "portable mutable branch and nonlinear owners" begin
        capacitor = AIMORA.Branches.CapacitorBranch(1, 0, 2.0e-6)
        capacitor.i_prev = 1.25
        capacitor.v_prev = -3.5
        capacitor.i_last = 0.75
        captured, restored = portable_element_owner_roundtrip(
            capacitor,
            AIMORA.Branches.CapacitorBranch(1, 0, 1.0e-6),
            :capacitor,
        )
        @test restored.signature_sha256 == captured.signature_sha256

        resistance = [1.0 0.1; 0.1 1.2]
        inductance = [2.0e-3 0.2e-3; 0.2e-3 2.5e-3]
        coupled_series = AIMORA.Branches.CoupledSeriesRLBranch(
            [1, 2],
            [0, 0],
            resistance,
            inductance,
        )
        coupled_series.previous_current .= [1.0, -0.5]
        coupled_series.previous_voltage .= [2.0, 3.0]
        coupled_series.last_current .= [0.25, -0.125]
        coupled_series.conductance_workspace .= [4.0 0.5; 0.5 5.0]
        coupled_series.history_current_workspace .= [6.0, 7.0]
        coupled_series.port_voltage_workspace .= [8.0, 9.0]
        coupled_series.current_workspace .= [10.0, 11.0]
        coupled_series.cached_dt_s = 50.0e-6
        captured, restored = portable_element_owner_roundtrip(
            coupled_series,
            AIMORA.Branches.CoupledSeriesRLBranch(
                [1, 2],
                [0, 0],
                resistance,
                inductance,
            ),
            :coupled_series,
        )
        @test restored.signature_sha256 == captured.signature_sha256

        coupled_inductive = AIMORA.Branches.CoupledInductiveBranch(
            [1, 2],
            [0, 0],
            [-0.2 0.02; 0.02 -0.3],
            2.0 * pi * 60.0,
        )
        coupled_inductive.previous_current .= [1.5, -0.75]
        coupled_inductive.previous_voltage .= [4.0, -2.0]
        coupled_inductive.last_current .= [0.5, -0.25]
        coupled_inductive.conductance_workspace .= [3.0 0.1; 0.1 4.0]
        coupled_inductive.history_current_workspace .= [5.0, 6.0]
        coupled_inductive.port_voltage_workspace .= [7.0, 8.0]
        coupled_inductive.current_workspace .= [9.0, 10.0]
        captured, restored = portable_element_owner_roundtrip(
            coupled_inductive,
            AIMORA.Branches.CoupledInductiveBranch(
                [1, 2],
                [0, 0],
                [-0.2 0.02; 0.02 -0.3],
                2.0 * pi * 60.0,
            ),
            :coupled_inductive,
        )
        @test restored.signature_sha256 == captured.signature_sha256

        saturable = AIMORA.Nonlinear.SaturableInductorBranch(
            1,
            0,
            0.1,
            2.0e-3,
            0.5e-3,
            2.0,
        )
        saturable.i_prev = 3.0
        saturable.v_prev = 4.0
        saturable.i_last = 3.5
        saturable.last_inductance = 0.5e-3
        saturable.saturated = true
        captured, restored = portable_element_owner_roundtrip(
            saturable,
            AIMORA.Nonlinear.SaturableInductorBranch(
                1,
                0,
                0.1,
                2.0e-3,
                0.5e-3,
                2.0,
            ),
            :saturable,
        )
        @test restored.signature_sha256 == captured.signature_sha256

        slope = AIMORA.Nonlinear.SaturatedTransformerNonlinearSlopeBranch(
            1,
            0,
            2.5,
        )
        captured, restored = portable_element_owner_roundtrip(
            slope,
            AIMORA.Nonlinear.SaturatedTransformerNonlinearSlopeBranch(1, 0, 1.0),
            :transformer_slope,
        )
        @test restored.signature_sha256 == captured.signature_sha256
    end

    mutable struct PortableTaskReplayOwner
        action_count::Int
        sampled_input::Float64
        held_output::Float64
        gate_high::Bool
        edge_count::Int
    end

    struct PortableTaskAction end
    function (::PortableTaskAction)(owner::PortableTaskReplayOwner, _time_s, _index)
        owner.action_count += 1
        return nothing
    end

    struct PortableTaskRead end
    (::PortableTaskRead)(owner::PortableTaskReplayOwner, _time_s, _index) =
        owner.sampled_input

    struct PortableTaskCompute end
    (::PortableTaskCompute)(
        _owner::PortableTaskReplayOwner,
        input,
        _time_s,
        sample_index,
    ) = input + sample_index

    struct PortableTaskWrite end
    function (::PortableTaskWrite)(
        owner::PortableTaskReplayOwner,
        value,
        _time_s,
        _sample_index,
    )
        owner.held_output = value
        return nothing
    end

    struct PortablePWMRead end
    (::PortablePWMRead)(::PortableTaskReplayOwner, _time_s, cycle_index) =
        isodd(cycle_index) ? 0.5 : 0.25

    struct PortablePWMWrite end
    function (::PortablePWMWrite)(
        owner::PortableTaskReplayOwner,
        high,
        _time_s,
        _edge_index,
    )
        owner.gate_high = high
        owner.edge_count += 1
        return nothing
    end

    function portable_exact_scheduler()
        tick_s = 1.0e-6
        action = AIMORA.EMTStudy.EMTExactSampledTask(
            :source_refresh,
            4.0e-6,
            PortableTaskAction();
            tick_s,
            priority = -2,
        )
        control = AIMORA.EMTStudy.EMTExactSampledControlTask(
            :sampled_controller,
            3.0e-6,
            PortableTaskRead(),
            PortableTaskCompute(),
            PortableTaskWrite();
            tick_s,
            computational_delay_s = tick_s,
            initial_output = 0.0,
            priority = -1,
            power_history_invalidating = true,
        )
        pwm = AIMORA.EMTStudy.EMTExactPWMTask(
            :carrier_gate,
            4.0e-6,
            PortablePWMRead(),
            PortablePWMWrite();
            tick_s,
            priority = 0,
        )
        return AIMORA.EMTStudy.EMTExactSampledTaskScheduler(
            tick_s;
            tasks = [action, control, pwm],
        )
    end

    mutable struct PortableGeneralTaskState
        gain::Float64
        accepted_count::Int
    end

    function AIMORA.EMTStudy.portable_emt_task_state(state::PortableGeneralTaskState)
        return AIMORA.PortableSnapshots.PortableSnapshotRecord(
            "aimora.emt.general_task_test_state.v1",
            Pair{String,Any}[
                "accepted_count" => state.accepted_count,
                "gain" => state.gain,
            ],
        )
    end

    function AIMORA.EMTStudy.restore_portable_emt_task_state(
        ::PortableGeneralTaskState,
        record::AIMORA.PortableSnapshots.PortableSnapshotRecord,
    )
        record.schema_id == "aimora.emt.general_task_test_state.v1" || error(
            "portable general task test state schema changed",
        )
        fields = Dict(record.fields)
        Set(keys(fields)) == Set(("accepted_count", "gain")) || error(
            "portable general task test state fields changed",
        )
        return PortableGeneralTaskState(
            Float64(fields["gain"]),
            Int(fields["accepted_count"]),
        )
    end

    mutable struct PortableGeneralTaskOwner
        input::Float64
        output::Float64
    end

    struct PortableGeneralTaskRead end
    (::PortableGeneralTaskRead)(
        _state,
        owner::PortableGeneralTaskOwner,
        _instant,
        _activation_index,
    ) = owner.input

    struct PortableGeneralTaskCompute end
    function (::PortableGeneralTaskCompute)(
        state::PortableGeneralTaskState,
        _owner::PortableGeneralTaskOwner,
        input,
        _instant,
        activation_index,
    )
        return state.gain * input + activation_index
    end

    struct PortableGeneralTaskWrite end
    function (::PortableGeneralTaskWrite)(
        state::PortableGeneralTaskState,
        owner::PortableGeneralTaskOwner,
        value,
        _instant,
        _activation_index,
    )
        owner.output = value
        state.accepted_count += 1
        return nothing
    end

    function portable_general_scheduler()
        specification = EMTTaskSpec(
            "portable_source_control",
            SourceEMTTask,
            emt_logical_time(0),
            emt_logical_time(3 // 1_000_000),
            emt_logical_time(0),
            emt_logical_time(1 // 1_000_000);
            read_resources = ["source_measurement"],
            write_resources = ["source_command"],
            effects = [InvalidateEMTPowerHistory, InvalidateEMTOutput],
        )
        plan = emt_task_plan(
            [specification];
            start = emt_logical_time(0),
            stop = emt_logical_time(12 // 1_000_000),
        )
        task = AIMORA.EMTStudy.GeneralEMTTask(
            only(plan.entries),
            PortableGeneralTaskRead(),
            PortableGeneralTaskCompute(),
            PortableGeneralTaskWrite();
            state = PortableGeneralTaskState(2.0, 0),
            initial_output = 0.0,
        )
        return AIMORA.EMTStudy.GeneralEMTTaskScheduler(plan, (task,))
    end

    function run_general_scheduler_through!(scheduler, owner, endpoint_s)
        while true
            instant = AIMORA.EMTStudy.next_general_task_instant(scheduler)
            instant === nothing && return scheduler
            Float64(instant) <= endpoint_s || return scheduler
            AIMORA.EMTStudy.run_due_general_tasks!(scheduler, owner, instant)
        end
    end

    struct PortableHybridTaskAction end
    function (::PortableHybridTaskAction)(owner, _time_s, execution_index)
        runtime = something(owner.runtime.context.control_system_runtime)
        runtime.state.values[:DRIVEN] = Float64(execution_index)
        return nothing
    end

    struct PortableHybridEventTransition end
    function (::PortableHybridEventTransition)(owner, time_s)
        runtime = something(owner.runtime.context.control_system_runtime)
        runtime.state.values[:MODUL] = 1.0 + time_s
        return nothing
    end

    function portable_hybrid_prepared()
        timestep_s = 50.0e-6
        return AIMORA.EMTStudy.prepare_emt_study(
            portable_control_network_feedback_deck();
            dt_s = timestep_s,
            t_end_s = 300.0e-6,
        )
    end

    function portable_hybrid_integrator(
        workspace;
        event_name::Symbol = :control_reference_change,
    )
        task = AIMORA.EMTStudy.EMTExactSampledTask(
            :control_supervisor,
            100.0e-6,
            PortableHybridTaskAction();
            tick_s = 1.0e-6,
            power_history_invalidating = true,
        )
        scheduler = AIMORA.EMTStudy.EMTExactSampledTaskScheduler(
            1.0e-6;
            tasks = [task],
        )
        event = AIMORA.EMTStudy.EMTHybridEventSurface(
            event_name,
            _owner -> nothing,
            PortableHybridEventTransition();
            priority = -1,
            repeatable = false,
            candidate_time = _owner -> 125.0e-6,
        )
        return AIMORA.EMTStudy.configure_emt_hybrid_execution(
            workspace;
            event_surfaces = [event],
            scheduler,
            include_device_events = false,
        )
    end

    function portable_hybrid_integrator(; event_name::Symbol = :control_reference_change)
        workspace = AIMORA.EMTStudy.EMTStudyWorkspace(portable_hybrid_prepared())
        return portable_hybrid_integrator(workspace; event_name)
    end

    @testset "portable exact sampled-task scheduler state" begin
        program_signature = bytes2hex(sha256("portable exact task callback program"))
        scheduler = portable_exact_scheduler()
        owner = PortableTaskReplayOwner(0, 2.0, 0.0, false, 0)
        for tick in 0:6
            AIMORA.EMTStudy.run_due_emt_sampled_tasks!(
                scheduler,
                owner,
                tick * scheduler.tick_s,
            )
        end
        inventory = AIMORA.EMTStudy.portable_emt_task_scheduler_state_inventory(
            scheduler;
            program_signature_sha256 = program_signature,
        )
        @test any(field -> field.identity == "scheduler.tasks", inventory.fields)
        @test any(field -> field.identity == "scheduler.occurrences", inventory.fields)
        @test length(scheduler.occurrences) > 0
        @test length(scheduler.tasks[2].samples) > 0
        @test length(scheduler.tasks[2].writes) > 0
        @test length(scheduler.tasks[3].occurrences) > 0

        restored = portable_exact_scheduler()
        AIMORA.EMTStudy.restore_portable_emt_task_scheduler_state_inventory!(
            restored,
            inventory;
            program_signature_sha256 = program_signature,
        )
        restored_inventory = AIMORA.EMTStudy.portable_emt_task_scheduler_state_inventory(
            restored;
            program_signature_sha256 = program_signature,
        )
        @test restored_inventory.signature_sha256 == inventory.signature_sha256
        before_mismatch = restored_inventory.signature_sha256
        mismatch = try
            AIMORA.EMTStudy.restore_portable_emt_task_scheduler_state_inventory!(
                restored,
                inventory;
                program_signature_sha256 = bytes2hex(sha256("changed callback program")),
            )
            nothing
        catch error
            error
        end
        @test mismatch isa AIMORA.PortableSnapshots.PortableSnapshotFailure
        @test mismatch.code == :task_program_mismatch
        @test AIMORA.EMTStudy.portable_emt_task_scheduler_state_inventory(
            restored;
            program_signature_sha256 = program_signature,
        ).signature_sha256 == before_mismatch

        restored_owner = deepcopy(owner)
        for tick in 7:12
            time_s = tick * scheduler.tick_s
            AIMORA.EMTStudy.run_due_emt_sampled_tasks!(scheduler, owner, time_s)
            AIMORA.EMTStudy.run_due_emt_sampled_tasks!(restored, restored_owner, time_s)
        end
        @test (
            restored_owner.action_count,
            restored_owner.sampled_input,
            restored_owner.held_output,
            restored_owner.gate_high,
            restored_owner.edge_count,
        ) == (
            owner.action_count,
            owner.sampled_input,
            owner.held_output,
            owner.gate_high,
            owner.edge_count,
        )
        @test AIMORA.EMTStudy.portable_emt_task_scheduler_state_inventory(
            restored;
            program_signature_sha256 = program_signature,
        ).signature_sha256 == AIMORA.EMTStudy.portable_emt_task_scheduler_state_inventory(
            scheduler;
            program_signature_sha256 = program_signature,
        ).signature_sha256
    end

    @testset "portable general rational-task scheduler state" begin
        program_signature = bytes2hex(sha256("portable general task callback program"))
        scheduler = portable_general_scheduler()
        owner = PortableGeneralTaskOwner(3.0, 0.0)
        run_general_scheduler_through!(scheduler, owner, 6.0e-6)
        @test only(scheduler.tasks).activation_count == 3
        @test !isempty(only(scheduler.tasks).pending)
        @test !isempty(scheduler.occurrences)

        inventory = AIMORA.EMTStudy.portable_emt_task_scheduler_state_inventory(
            scheduler;
            program_signature_sha256 = program_signature,
        )
        restored = portable_general_scheduler()
        AIMORA.EMTStudy.restore_portable_emt_task_scheduler_state_inventory!(
            restored,
            inventory;
            program_signature_sha256 = program_signature,
        )
        @test AIMORA.EMTStudy.portable_emt_task_scheduler_state_inventory(
            restored;
            program_signature_sha256 = program_signature,
        ).signature_sha256 == inventory.signature_sha256

        restored_owner = deepcopy(owner)
        run_general_scheduler_through!(scheduler, owner, Inf)
        run_general_scheduler_through!(restored, restored_owner, Inf)
        @test (restored_owner.input, restored_owner.output) ==
            (owner.input, owner.output)
        @test only(restored.tasks).state.gain == only(scheduler.tasks).state.gain
        @test only(restored.tasks).state.accepted_count ==
            only(scheduler.tasks).state.accepted_count
        @test restored.occurrences == scheduler.occurrences
        @test AIMORA.EMTStudy.portable_emt_task_scheduler_state_inventory(
            restored;
            program_signature_sha256 = program_signature,
        ).signature_sha256 == AIMORA.EMTStudy.portable_emt_task_scheduler_state_inventory(
            scheduler;
            program_signature_sha256 = program_signature,
        ).signature_sha256
    end

    @testset "portable hybrid event and task coordinator state" begin
        task_signature = bytes2hex(sha256("portable hybrid task callback program"))
        event_signature = bytes2hex(sha256("portable hybrid event callback program"))
        baseline = portable_hybrid_integrator()
        candidate = portable_hybrid_integrator()
        for _ in 1:3
            AIMORA.EMTStudy.advance_emt_hybrid_step!(baseline)
            AIMORA.EMTStudy.advance_emt_hybrid_step!(candidate)
        end
        @test only(baseline.surface_fired)
        @test length(baseline.occurrences) == 1
        @test length(baseline.scheduler.occurrences) > 0
        plain_capture_failure = try
            AIMORA.EMTStudy.portable_emt_state_inventory(baseline.workspace)
            nothing
        catch error
            error
        end
        @test plain_capture_failure isa AIMORA.PortableSnapshots.PortableSnapshotFailure
        @test plain_capture_failure.code == :hybrid_owner_required
        mktempdir() do directory
            local_checkpoint_failure = try
                AIMORA.EMTStudy.write_emt_checkpoint(
                    joinpath(directory, "incomplete-hybrid.ckpt"),
                    baseline.workspace,
                )
                nothing
            catch error
                error
            end
            @test local_checkpoint_failure isa ArgumentError
            @test occursin("coordinating integrator", sprint(showerror, local_checkpoint_failure))
            @test isempty(readdir(directory))
        end
        inventory = AIMORA.EMTStudy.portable_emt_hybrid_state_inventory(
            baseline;
            task_program_signature_sha256 = task_signature,
            event_program_signature_sha256 = event_signature,
        )
        public_inventory = AIMORA.EMTStudy._portable_emt_state_inventory(
            baseline.workspace,
        )
        plain_candidate = AIMORA.EMTStudy.EMTStudyWorkspace(
            portable_hybrid_prepared(),
        )
        missing_owner = try
            AIMORA.EMTStudy.restore_portable_emt_state_inventory!(
                plain_candidate,
                public_inventory,
            )
            nothing
        catch error
            error
        end
        @test missing_owner isa AIMORA.PortableSnapshots.PortableSnapshotFailure
        @test missing_owner.code == :hybrid_owner_required
        @test plain_candidate.execution_mode == :unselected

        monolithic_candidate = AIMORA.EMTStudy.EMTStudyWorkspace(
            portable_hybrid_prepared(),
        )
        AIMORA.EMTStudy.evaluate_emt_study!(monolithic_candidate)
        monolithic_signature = AIMORA.EMTStudy.portable_emt_state_inventory(
            monolithic_candidate,
        ).signature_sha256
        mode_mismatch = try
            AIMORA.EMTStudy.restore_portable_emt_state_inventory!(
                monolithic_candidate,
                public_inventory,
            )
            nothing
        catch error
            error
        end
        @test mode_mismatch isa AIMORA.PortableSnapshots.PortableSnapshotFailure
        @test mode_mismatch.code == :execution_mode_mismatch
        @test AIMORA.EMTStudy.portable_emt_state_inventory(
            monolithic_candidate,
        ).signature_sha256 == monolithic_signature
        candidate.surface_fired .= false
        empty!(candidate.occurrences)
        candidate.accepted_interval_count += 7
        only(candidate.scheduler.tasks).execution_count = 0
        empty!(candidate.scheduler.occurrences)
        AIMORA.EMTStudy.restore_portable_emt_hybrid_state_inventory!(
            candidate,
            inventory;
            task_program_signature_sha256 = task_signature,
            event_program_signature_sha256 = event_signature,
        )
        @test AIMORA.EMTStudy.portable_emt_hybrid_state_inventory(
            candidate;
            task_program_signature_sha256 = task_signature,
            event_program_signature_sha256 = event_signature,
        ).signature_sha256 == inventory.signature_sha256

        before_mismatch = AIMORA.EMTStudy.portable_emt_hybrid_state_inventory(
            candidate;
            task_program_signature_sha256 = task_signature,
            event_program_signature_sha256 = event_signature,
        ).signature_sha256
        mismatch = try
            AIMORA.EMTStudy.restore_portable_emt_hybrid_state_inventory!(
                candidate,
                inventory;
                task_program_signature_sha256 = task_signature,
                event_program_signature_sha256 = bytes2hex(sha256("changed event program")),
            )
            nothing
        catch error
            error
        end
        @test mismatch isa AIMORA.PortableSnapshots.PortableSnapshotFailure
        @test mismatch.code == :hybrid_event_program_mismatch
        @test AIMORA.EMTStudy.portable_emt_hybrid_state_inventory(
            candidate;
            task_program_signature_sha256 = task_signature,
            event_program_signature_sha256 = event_signature,
        ).signature_sha256 == before_mismatch

        changed_surface = portable_hybrid_integrator(
            event_name = :different_control_reference_change,
        )
        for _ in 1:3
            AIMORA.EMTStudy.advance_emt_hybrid_step!(changed_surface)
        end
        surface_failure = try
            AIMORA.EMTStudy.restore_portable_emt_hybrid_state_inventory!(
                changed_surface,
                inventory;
                task_program_signature_sha256 = task_signature,
                event_program_signature_sha256 = event_signature,
            )
            nothing
        catch error
            error
        end
        @test surface_failure isa AIMORA.PortableSnapshots.PortableSnapshotFailure
        @test surface_failure.code == :hybrid_event_program_mismatch

        project_signature = bytes2hex(sha256("portable hybrid project"))
        model_signature = bytes2hex(sha256("portable hybrid model"))
        settings_signature = bytes2hex(sha256("portable hybrid settings"))
        prepared = portable_hybrid_prepared()
        disk_restored = mktempdir() do directory
            first_path = joinpath(directory, "hybrid-first.aimora-snapshot")
            second_path = joinpath(directory, "hybrid-second.aimora-snapshot")
            first_descriptor = AIMORA.EMTStudy.write_portable_emt_hybrid_snapshot(
                first_path,
                baseline;
                project_signature_sha256 = project_signature,
                model_signature_sha256 = model_signature,
                settings_signature_sha256 = settings_signature,
                task_program_signature_sha256 = task_signature,
                event_program_signature_sha256 = event_signature,
                provenance = "AIMORA-authored portable hybrid split replay",
            )
            second_descriptor = AIMORA.EMTStudy.write_portable_emt_hybrid_snapshot(
                second_path,
                baseline;
                project_signature_sha256 = project_signature,
                model_signature_sha256 = model_signature,
                settings_signature_sha256 = settings_signature,
                task_program_signature_sha256 = task_signature,
                event_program_signature_sha256 = event_signature,
                provenance = "AIMORA-authored portable hybrid split replay",
            )
            @test read(first_path) == read(second_path)
            @test first_descriptor.content_sha256 == second_descriptor.content_sha256
            @test first_descriptor.metadata.accepted_step == 3
            @test !occursin("/home/", String(read(first_path)))
            result = AIMORA.EMTStudy.read_portable_emt_hybrid_snapshot(
                first_path,
                prepared;
                hybrid_factory = workspace -> portable_hybrid_integrator(workspace),
                project_signature_sha256 = project_signature,
                model_signature_sha256 = model_signature,
                settings_signature_sha256 = settings_signature,
                task_program_signature_sha256 = task_signature,
                event_program_signature_sha256 = event_signature,
            )
            @test result.reconstructed
            @test result.descriptor.content_sha256 == first_descriptor.content_sha256
            result.integrator
        end

        baseline_trace = AIMORA.EMTStudy.evaluate_emt_hybrid_study!(baseline)
        candidate_trace = AIMORA.EMTStudy.evaluate_emt_hybrid_study!(candidate)
        disk_trace = AIMORA.EMTStudy.evaluate_emt_hybrid_study!(disk_restored)
        @test candidate_trace.time_s == baseline_trace.time_s
        @test candidate_trace.voltage_pu == baseline_trace.voltage_pu
        @test candidate_trace.output_pu == baseline_trace.output_pu
        @test disk_trace.time_s == baseline_trace.time_s
        @test disk_trace.voltage_pu == baseline_trace.voltage_pu
        @test disk_trace.output_pu == baseline_trace.output_pu
        @test candidate.occurrences == baseline.occurrences
        @test disk_restored.occurrences == baseline.occurrences
        @test candidate.scheduler.occurrences == baseline.scheduler.occurrences
        @test disk_restored.scheduler.occurrences == baseline.scheduler.occurrences
        @test AIMORA.EMTStudy.portable_emt_hybrid_state_inventory(
            candidate;
            task_program_signature_sha256 = task_signature,
            event_program_signature_sha256 = event_signature,
        ).signature_sha256 == AIMORA.EMTStudy.portable_emt_hybrid_state_inventory(
            baseline;
            task_program_signature_sha256 = task_signature,
            event_program_signature_sha256 = event_signature,
        ).signature_sha256
    end

    @testset "portable EMT accepted workspace state inventory" begin
        timestep_s = 1.0e-4
        parsed = AIMORA.DeckParser.parse_deck_lines(
            [
                "BEGIN NEW DATA CASE",
                "source src source 1.0e9 1.0 60.0 0.0 1.0",
                "bergeron_line line source bus 2.0 $(2.0 * timestep_s) $timestep_s 0.8",
                "resistor load bus 0 2.0",
                "time_switch tie bus 0 $(8.0 * timestep_s) $(20.0 * timestep_s) false 20.0 0.0",
            ];
            source = "portable accepted-workspace state test",
        )
        prepared = AIMORA.EMTStudy.prepare_emt_study(
            parsed;
            dt_s = timestep_s,
            t_end_s = 3.0 * timestep_s,
        )
        accepted = AIMORA.EMTStudy.EMTStudyWorkspace(prepared)
        AIMORA.EMTStudy.evaluate_emt_study!(accepted)
        inventory = AIMORA.EMTStudy.portable_emt_state_inventory(accepted)
        @test inventory.scalar_count == 79
        @test length(inventory.fields) == 68
        @test any(field -> field.identity == "network.accepted_node_voltage", inventory.fields)
        @test any(field -> field.identity == "model.i00000002.from_wave_history", inventory.fields)

        candidate = AIMORA.EMTStudy.EMTStudyWorkspace(prepared)
        AIMORA.EMTStudy.restore_portable_emt_state_inventory!(candidate, inventory)
        restored = AIMORA.EMTStudy.portable_emt_state_inventory(candidate)
        @test restored.signature_sha256 == inventory.signature_sha256
        @test candidate.runtime.context.system.v == accepted.runtime.context.system.v
        @test candidate.runtime.context.recorded_step_indices ==
            accepted.runtime.context.recorded_step_indices

        project_signature = bytes2hex(sha256("portable workspace project"))
        model_signature = bytes2hex(sha256("portable workspace model"))
        settings_signature = bytes2hex(sha256("portable workspace settings"))
        mktempdir() do directory
            first_path = joinpath(directory, "first-workspace.aimora-snapshot")
            second_path = joinpath(directory, "second-workspace.aimora-snapshot")
            first_descriptor = AIMORA.EMTStudy.write_portable_emt_workspace_snapshot(
                first_path,
                accepted;
                project_signature_sha256 = project_signature,
                model_signature_sha256 = model_signature,
                settings_signature_sha256 = settings_signature,
                provenance = "AIMORA-authored synthetic portable split replay",
            )
            second_descriptor = AIMORA.EMTStudy.write_portable_emt_workspace_snapshot(
                second_path,
                accepted;
                project_signature_sha256 = project_signature,
                model_signature_sha256 = model_signature,
                settings_signature_sha256 = settings_signature,
                provenance = "AIMORA-authored synthetic portable split replay",
            )
            @test read(first_path) == read(second_path)
            @test first_descriptor.content_sha256 == second_descriptor.content_sha256
            @test filesize(first_path) < 50_000
            @test AIMORA.PortableSnapshots.inspect_portable_emt_snapshot(
                first_path,
            ).content_sha256 == first_descriptor.content_sha256
            @test_throws AIMORA.PortableSnapshots.PortableSnapshotFailure AIMORA.PortableSnapshots.read_portable_emt_snapshot(
                first_path,
            )
            @test !occursin("AIMORASolvers", String(read(first_path)))
            @test !occursin("/home/", String(read(first_path)))

            restored_result = AIMORA.EMTStudy.read_portable_emt_workspace_snapshot(
                first_path,
                prepared;
                project_signature_sha256 = project_signature,
                model_signature_sha256 = model_signature,
                settings_signature_sha256 = settings_signature,
            )
            @test restored_result.reconstructed
            @test restored_result.descriptor.content_sha256 ==
                first_descriptor.content_sha256
            @test_throws AIMORA.PortableSnapshots.PortableSnapshotFailure AIMORA.EMTStudy.read_portable_emt_workspace_snapshot(
                first_path,
                prepared;
                project_signature_sha256 = project_signature,
                model_signature_sha256 = model_signature,
                settings_signature_sha256 = bytes2hex(sha256("stale settings")),
            )

            full_prepared = AIMORA.EMTStudy.prepare_emt_study(
                parsed;
                dt_s = timestep_s,
                t_end_s = 6.0 * timestep_s,
            )
            full_workspace = AIMORA.EMTStudy.EMTStudyWorkspace(full_prepared)
            full_trace = AIMORA.EMTStudy.evaluate_emt_study!(full_workspace)
            request = AIMORA.DeckParser.parse_emt_restart_request([
                "START AGAIN",
                lpad("9999", 8),
            ])
            resumed = AIMORA.EMTStudy.restart_emt_study!(
                restored_result.workspace,
                request;
                additional_time_s = 3.0 * timestep_s,
            )
            @test resumed.trace.time_s == full_trace.time_s
            @test resumed.trace.voltage_pu == full_trace.voltage_pu
            @test resumed.trace.output_pu == full_trace.output_pu
            @test resumed.checkpoint_state_error == 0.0
        end
    end


    @testset "portable EMT scheduled parameter alteration state" begin
        timestep_s = 1.0e-4
        parsed = AIMORA.DeckParser.parse_deck_lines(
            [
                "source src BUS1 4.0 0.0 60.0 0.0 1.0",
                "BLANK BRANCH",
                "1,BUS1,0,1.0,0.001,1.0e-6",
                "END",
            ];
            source = "portable scheduled parameter alteration test",
        )
        alteration = AIMORA.EMTStudy.SeriesRLCAlteration(
            :branch_fixed_1,
            timestep_s,
            2.0,
            2.0e-3,
            2.0e-6,
        )
        prepared = AIMORA.EMTStudy.prepare_emt_study(
            parsed;
            dt_s = timestep_s,
            t_end_s = 3.0 * timestep_s,
            series_rlc_alterations = [alteration],
        )
        accepted = AIMORA.EMTStudy.EMTStudyWorkspace(prepared)
        AIMORA.EMTStudy.evaluate_emt_study!(accepted)
        inventory = AIMORA.EMTStudy.portable_emt_state_inventory(accepted)
        @test length(inventory.fields) == 82
        @test inventory.scalar_count == 68
        @test length(accepted.runtime.context.series_rlc_alteration_records) == 1

        candidate = AIMORA.EMTStudy.EMTStudyWorkspace(prepared)
        AIMORA.EMTStudy.restore_portable_emt_state_inventory!(candidate, inventory)
        restored = AIMORA.EMTStudy.portable_emt_state_inventory(candidate)
        @test restored.signature_sha256 == inventory.signature_sha256
        @test candidate.runtime.context.series_rlc_alteration_records ==
            accepted.runtime.context.series_rlc_alteration_records
        restored_branch = candidate.runtime.context.system.elements[2]
        @test (restored_branch.r, restored_branch.l, restored_branch.c) ==
            (2.0, 2.0e-3, 2.0e-6)

        project_signature = bytes2hex(sha256("portable alteration project"))
        model_signature = bytes2hex(sha256("portable alteration model"))
        settings_signature = bytes2hex(sha256("portable alteration settings"))
        mktempdir() do directory
            path = joinpath(directory, "alteration.aimora-snapshot")
            AIMORA.EMTStudy.write_portable_emt_workspace_snapshot(
                path,
                accepted;
                project_signature_sha256 = project_signature,
                model_signature_sha256 = model_signature,
                settings_signature_sha256 = settings_signature,
                provenance = "AIMORA-authored portable parameter alteration test",
            )
            restored_result = AIMORA.EMTStudy.read_portable_emt_workspace_snapshot(
                path,
                prepared;
                project_signature_sha256 = project_signature,
                model_signature_sha256 = model_signature,
                settings_signature_sha256 = settings_signature,
            )
            full_prepared = AIMORA.EMTStudy.prepare_emt_study(
                parsed;
                dt_s = timestep_s,
                t_end_s = 6.0 * timestep_s,
                series_rlc_alterations = [alteration],
            )
            full_workspace = AIMORA.EMTStudy.EMTStudyWorkspace(full_prepared)
            full_trace = AIMORA.EMTStudy.evaluate_emt_study!(full_workspace)
            request = AIMORA.DeckParser.parse_emt_restart_request([
                "START AGAIN",
                lpad("9999", 8),
            ])
            resumed = AIMORA.EMTStudy.restart_emt_study!(
                restored_result.workspace,
                request;
                additional_time_s = 3.0 * timestep_s,
            )
            @test resumed.trace.time_s == full_trace.time_s
            @test resumed.trace.voltage_pu == full_trace.voltage_pu
            @test resumed.trace.output_pu == full_trace.output_pu
            @test resumed.checkpoint_state_error == 0.0
        end
    end


    @testset "portable EMT tabulated and analytic source-function state" begin
        timestep_s = 0.25
        parsed = portable_source_function_deck()
        AIMORA.DeckParser.assert_deck_valid!(parsed)
        program = portable_source_signal_program()
        prepared = AIMORA.EMTStudy.prepare_emt_study(
            parsed;
            dt_s = timestep_s,
            t_end_s = 0.5,
            source_signal_provider = program,
        )
        accepted = AIMORA.EMTStudy.EMTStudyWorkspace(prepared)
        AIMORA.EMTStudy.evaluate_emt_study!(accepted)
        runtime = accepted.runtime.context.source_function_runtime
        @test runtime !== nothing
        @test runtime.executed_step_count > 0
        @test runtime.external_signal_count == runtime.executed_step_count
        @test runtime.analytic_execution_count == runtime.executed_step_count
        @test runtime.next_input_row_index == runtime.executed_step_count + 1
        @test length(runtime.stage_samples) == runtime.executed_step_count

        inventory = AIMORA.EMTStudy.portable_emt_state_inventory(accepted)
        identities = Set(getfield.(inventory.fields, :identity))
        @test "source_function.present" in identities
        @test "source_function.card.voltbc_values" in identities
        @test "source_function.accepted_slot_values" in identities
        @test "source_function.stage.accepted" in identities
        @test "source_function.provider.interpolation.sample_value" in identities
        @test "source_function.provider.analytic_source_type" in identities

        candidate = AIMORA.EMTStudy.EMTStudyWorkspace(prepared)
        AIMORA.EMTStudy.restore_portable_emt_state_inventory!(candidate, inventory)
        restored_runtime = candidate.runtime.context.source_function_runtime
        @test AIMORA.EMTStudy.portable_emt_state_inventory(candidate).signature_sha256 ==
            inventory.signature_sha256
        @test restored_runtime.state.voltbc_values == runtime.state.voltbc_values
        @test getindex.(restored_runtime.slot_values) == getindex.(runtime.slot_values)
        @test getindex.(restored_runtime.row_slot_values) ==
            getindex.(runtime.row_slot_values)
        @test getfield.(restored_runtime.stage_samples, :time_s) ==
            getfield.(runtime.stage_samples, :time_s)
        @test getfield.(restored_runtime.stage_samples, :accepted_values) ==
            getfield.(runtime.stage_samples, :accepted_values)

        changed_prepared = AIMORA.EMTStudy.prepare_emt_study(
            parsed;
            dt_s = timestep_s,
            t_end_s = 0.5,
            source_signal_provider = portable_source_signal_program(
                first_slot_offset = 0.5,
            ),
        )
        changed_candidate = AIMORA.EMTStudy.EMTStudyWorkspace(changed_prepared)
        changed_error = try
            AIMORA.EMTStudy.restore_portable_emt_state_inventory!(
                changed_candidate,
                inventory,
            )
            nothing
        catch error
            error
        end
        @test changed_error isa AIMORA.PortableSnapshots.PortableSnapshotFailure
        @test changed_error.code == :settings_mismatch

        project_signature = bytes2hex(sha256("portable source-function project"))
        model_signature = bytes2hex(sha256("portable source-function model"))
        settings_signature = bytes2hex(sha256("portable source-function settings"))
        mktempdir() do directory
            path = joinpath(directory, "source-function.aimora-snapshot")
            descriptor = AIMORA.EMTStudy.write_portable_emt_workspace_snapshot(
                path,
                accepted;
                project_signature_sha256 = project_signature,
                model_signature_sha256 = model_signature,
                settings_signature_sha256 = settings_signature,
                provenance = "AIMORA-authored portable source-function restart test",
            )
            @test descriptor.canonical_bytes == filesize(path)
            @test !occursin("/home/", String(read(path)))
            restored_result = AIMORA.EMTStudy.read_portable_emt_workspace_snapshot(
                path,
                prepared;
                project_signature_sha256 = project_signature,
                model_signature_sha256 = model_signature,
                settings_signature_sha256 = settings_signature,
            )

            full_prepared = AIMORA.EMTStudy.prepare_emt_study(
                parsed;
                dt_s = timestep_s,
                t_end_s = 1.0,
                source_signal_provider = program,
            )
            full_workspace = AIMORA.EMTStudy.EMTStudyWorkspace(full_prepared)
            full_trace = AIMORA.EMTStudy.evaluate_emt_study!(full_workspace)
            request = AIMORA.DeckParser.parse_emt_restart_request([
                "START AGAIN",
                lpad("9999", 8),
            ])
            resumed = AIMORA.EMTStudy.restart_emt_study!(
                restored_result.workspace,
                request;
                additional_time_s = 0.5,
            )
            @test resumed.trace.time_s == full_trace.time_s
            @test resumed.trace.voltage_pu == full_trace.voltage_pu
            @test resumed.trace.output_pu == full_trace.output_pu
            @test resumed.checkpoint_state_error == 0.0
            resumed_runtime = restored_result.workspace.runtime.context.source_function_runtime
            full_runtime = full_workspace.runtime.context.source_function_runtime
            @test getfield.(resumed_runtime.stage_samples, :time_s) ==
                getfield.(full_runtime.stage_samples, :time_s)
            @test getfield.(resumed_runtime.stage_samples, :accepted_values) ==
                getfield.(full_runtime.stage_samples, :accepted_values)
        end
    end


    @testset "portable EMT control-network state" begin
        timestep_s = 50.0e-6
        parsed = portable_control_network_feedback_deck()
        AIMORA.DeckParser.assert_deck_valid!(parsed)
        prepared = AIMORA.EMTStudy.prepare_emt_study(
            parsed;
            dt_s = timestep_s,
            t_end_s = 150.0e-6,
        )
        accepted = AIMORA.EMTStudy.EMTStudyWorkspace(prepared)
        AIMORA.EMTStudy.evaluate_emt_study!(accepted)
        runtime = accepted.runtime.context.control_system_runtime
        @test runtime !== nothing
        @test runtime.executed_step_count > 0
        @test !isempty(runtime.state.function_states)
        @test any(!isempty, getfield.(runtime.state.function_states, :history_terms))
        @test !isempty(runtime.switch_elements)

        inventory = AIMORA.EMTStudy.portable_emt_state_inventory(accepted)
        identities = Set(getfield.(inventory.fields, :identity))
        @test "control.present" in identities
        @test "control.identity" in identities
        @test "control.values.accepted" in identities
        @test "control.function.history.values" in identities
        @test "control.source.waveform.types" in identities
        controlled_switch_index = findfirst(
            element -> element isa AIMORA.TACS.TACSControlledSwitch,
            accepted.runtime.context.system.elements,
        )
        @test controlled_switch_index !== nothing
        controlled_switch_prefix =
            "model.i" * lpad(string(controlled_switch_index), 8, '0')
        @test "$controlled_switch_prefix.closed" in identities

        candidate = AIMORA.EMTStudy.EMTStudyWorkspace(prepared)
        AIMORA.EMTStudy.restore_portable_emt_state_inventory!(candidate, inventory)
        restored_runtime = candidate.runtime.context.control_system_runtime
        @test AIMORA.EMTStudy.portable_emt_state_inventory(candidate).signature_sha256 ==
            inventory.signature_sha256
        @test restored_runtime.state.values == runtime.state.values
        @test getfield.(restored_runtime.state.function_states, :history_terms) ==
            getfield.(runtime.state.function_states, :history_terms)
        @test getfield.(restored_runtime.switch_elements, :closed) ==
            getfield.(runtime.switch_elements, :closed)

        project_signature = bytes2hex(sha256("portable control-network project"))
        model_signature = bytes2hex(sha256("portable control-network model"))
        settings_signature = bytes2hex(sha256("portable control-network settings"))
        mktempdir() do directory
            path = joinpath(directory, "control-network.aimora-snapshot")
            AIMORA.EMTStudy.write_portable_emt_workspace_snapshot(
                path,
                accepted;
                project_signature_sha256 = project_signature,
                model_signature_sha256 = model_signature,
                settings_signature_sha256 = settings_signature,
                provenance = "AIMORA-authored portable control-network restart test",
            )
            restored_result = AIMORA.EMTStudy.read_portable_emt_workspace_snapshot(
                path,
                prepared;
                project_signature_sha256 = project_signature,
                model_signature_sha256 = model_signature,
                settings_signature_sha256 = settings_signature,
            )
            full_prepared = AIMORA.EMTStudy.prepare_emt_study(
                parsed;
                dt_s = timestep_s,
                t_end_s = 300.0e-6,
            )
            full_workspace = AIMORA.EMTStudy.EMTStudyWorkspace(full_prepared)
            full_trace = AIMORA.EMTStudy.evaluate_emt_study!(full_workspace)
            request = AIMORA.DeckParser.parse_emt_restart_request([
                "START AGAIN",
                lpad("9999", 8),
            ])
            resumed = AIMORA.EMTStudy.restart_emt_study!(
                restored_result.workspace,
                request;
                additional_time_s = 150.0e-6,
            )
            @test resumed.trace.time_s == full_trace.time_s
            @test resumed.trace.voltage_pu == full_trace.voltage_pu
            @test resumed.trace.output_pu == full_trace.output_pu
            @test resumed.checkpoint_state_error == 0.0
            resumed_runtime =
                restored_result.workspace.runtime.context.control_system_runtime
            full_runtime = full_workspace.runtime.context.control_system_runtime
            @test resumed_runtime.state.values == full_runtime.state.values
            @test getfield.(resumed_runtime.state.function_states, :history_terms) ==
                getfield.(full_runtime.state.function_states, :history_terms)
            @test getfield.(resumed_runtime.switch_elements, :closed) ==
                getfield.(full_runtime.switch_elements, :closed)
        end
    end
end
