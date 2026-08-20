const MeasurementChainTest = AIMORA.MeasurementChains
using LinearAlgebra

function measurement_test_provenance(unit::AbstractString)
    return AIMORA.StudyCore.ParameterProvenance(
        "AIMORA synthetic measurement test",
        unit,
        "No transformation",
        "Exact synthetic input within the declared test domain",
        "Generic nonvendor measurement-chain test",
        AIMORA.StudyCore.PhysicalModelParameter,
    )
end

function measurement_test_acquisition(;
    tick_s=1.0e-3,
    sample_period_ticks=1,
    first_sample_tick=0,
    delay_ticks=0,
    window_count=20,
    nominal_frequency_hz=50.0,
    positive_sequence_threshold=1.0e-9,
    frequency_update_separation=1,
    maximum_retained_samples=1,
)
    return MeasurementChainTest.MeasurementAcquisitionSettings(
        tick_s=tick_s,
        sample_period_ticks=sample_period_ticks,
        first_sample_tick=first_sample_tick,
        delay_ticks=delay_ticks,
        window_weights_newest_first=ones(window_count),
        nominal_frequency_hz=nominal_frequency_hz,
        positive_sequence_threshold=positive_sequence_threshold,
        frequency_update_separation=frequency_update_separation,
        maximum_retained_samples=maximum_retained_samples,
    )
end

function measurement_test_specification(;
    id=:measurement_test,
    family=MeasurementChainTest.ElectronicVoltageSensorMeasurement,
    channel_names=[:value],
    quantity=:voltage,
    unit="V",
    orientation="positive from the declared observed terminal to reference",
    phase_order=Symbol[],
    conditioning=MeasurementChainTest.AnalogMeasurementStateSpace(),
    quantizer=MeasurementChainTest.UnquantizedMeasurement(-1.0e6, 1.0e6),
    acquisition=measurement_test_acquisition(),
    minimum_input=-1.0e6,
    maximum_input=1.0e6,
    maximum_spectral_frequency_hz=nothing,
)
    spectral_frequency = maximum_spectral_frequency_hz === nothing ?
        0.35 * MeasurementChainTest.measurement_sample_rate_hz(acquisition) :
        maximum_spectral_frequency_hz
    return MeasurementChainTest.MeasurementChainSpecification(
        id,
        family;
        channel_names=channel_names,
        quantity=quantity,
        unit=unit,
        orientation=orientation,
        phase_order=phase_order,
        conditioning=conditioning,
        quantizer=quantizer,
        acquisition=acquisition,
        minimum_input=minimum_input,
        maximum_input=maximum_input,
        maximum_spectral_frequency_hz=spectral_frequency,
        maximum_timestep_s=acquisition.tick_s,
        provenance=measurement_test_provenance(unit),
    )
end

function accept_measurement_test_value!(runtime, values, tick)
    MeasurementChainTest.prepare_measurement_analog_step!(
        runtime,
        values,
        runtime.specification.acquisition.tick_s,
    )
    return MeasurementChainTest.accept_measurement_analog_step!(runtime, tick)
end

@testset "measurement-chain contracts and refusals" begin
    @test length(instances(MeasurementChainTest.MeasurementProductFamily)) == 7
    quantizer = MeasurementChainTest.UniformMeasurementQuantizer(
        lower_limit=-1.0,
        upper_limit=1.0,
        engineering_step=0.1,
        minimum_code=-10,
        maximum_code=10,
    )
    @test quantizer.minimum_code == -10
    @test quantizer.maximum_code == 10
    @test_throws ArgumentError MeasurementChainTest.UniformMeasurementQuantizer(
        lower_limit=-2.0,
        upper_limit=2.0,
        engineering_step=0.1,
        minimum_code=-10,
        maximum_code=10,
    )
    @test_throws ArgumentError MeasurementChainTest.AnalogMeasurementStateSpace(
        reshape([1.0], 1, 1),
        [1.0],
        [1.0],
        0.0,
    )
    @test_throws ArgumentError measurement_test_specification(
        family=MeasurementChainTest.ThreePhaseSampledMeasurement,
    )

    specification = measurement_test_specification()
    contract = MeasurementChainTest.measurement_chain_contract(specification)
    @test contract.maturity == :implemented
    @test contract.owner == "AIMORA.MeasurementChains"
    @test length(MeasurementChainTest.measurement_chain_signature(specification)) == 64
    @test MeasurementChainTest.measurement_sample_rate_hz(
        specification.acquisition,
    ) == 1.0e3
    @test MeasurementChainTest.measurement_due_at_tick(
        specification.acquisition,
        0,
    )
    @test !MeasurementChainTest.measurement_due_at_tick(
        specification.acquisition,
        -1,
    )
end

@testset "measurement clipping and deterministic quantizer ties" begin
    expected_codes = Dict(
        MeasurementChainTest.MeasurementTiesToEven => (0, 0),
        MeasurementChainTest.MeasurementTiesAwayFromZero => (1, -1),
        MeasurementChainTest.MeasurementTiesTowardZero => (0, 0),
    )
    for (tie_rule, (positive_code, negative_code)) in expected_codes
        quantizer = MeasurementChainTest.UniformMeasurementQuantizer(
            lower_limit=-1.0,
            upper_limit=1.0,
            engineering_step=0.1,
            minimum_code=-10,
            maximum_code=10,
            tie_rule=tie_rule,
        )
        positive_runtime = MeasurementChainTest.MeasurementChainRuntime(
            measurement_test_specification(quantizer=quantizer),
        )
        positive_sample = only(accept_measurement_test_value!(
            positive_runtime,
            [0.05],
            0,
        ))
        @test only(positive_sample.codes) == positive_code

        negative_runtime = MeasurementChainTest.MeasurementChainRuntime(
            measurement_test_specification(quantizer=quantizer),
        )
        negative_sample = only(accept_measurement_test_value!(
            negative_runtime,
            [-0.05],
            0,
        ))
        @test only(negative_sample.codes) == negative_code

        clipped_runtime = MeasurementChainTest.MeasurementChainRuntime(
            measurement_test_specification(
                quantizer=quantizer,
                minimum_input=-2.0,
                maximum_input=2.0,
            ),
        )
        clipped_sample = only(accept_measurement_test_value!(
            clipped_runtime,
            [1.5],
            0,
        ))
        @test only(clipped_sample.instantaneous) == 1.0
        @test only(clipped_sample.codes) == 10
        @test only(clipped_sample.clipped)
    end
end

@testset "measurement analog trial, exact acquisition, and delayed release" begin
    conditioning = MeasurementChainTest.AnalogMeasurementStateSpace(
        reshape([-200.0], 1, 1),
        [200.0],
        [1.0],
        0.0,
    )
    acquisition = measurement_test_acquisition(
        sample_period_ticks=2,
        delay_ticks=2,
    )
    runtime = MeasurementChainTest.MeasurementChainRuntime(
        measurement_test_specification(
            conditioning=conditioning,
            acquisition=acquisition,
        ),
    )
    MeasurementChainTest.prepare_measurement_analog_step!(runtime, [1.0], 1.0e-3)
    candidate_state = only(runtime.candidate.state)
    @test candidate_state ≈ 1.0 / 11.0 atol=1.0e-15
    MeasurementChainTest.discard_measurement_analog_step!(runtime)
    @test only(runtime.analog_state) == 0.0
    @test !runtime.candidate_active

    @test isempty(accept_measurement_test_value!(runtime, [1.0], 0))
    sampled_value = only(runtime.analog_output)
    @test isempty(accept_measurement_test_value!(runtime, [1.0], 1))
    released = accept_measurement_test_value!(runtime, [1.0], 2)
    @test length(released) == 1
    @test only(released).source_tick == 0
    @test only(released).release_tick == 2
    @test only(only(released).instantaneous) == sampled_value

    missed_release_runtime = MeasurementChainTest.MeasurementChainRuntime(
        measurement_test_specification(
            acquisition=measurement_test_acquisition(
                sample_period_ticks=4,
                delay_ticks=1,
            ),
        ),
    )
    @test isempty(accept_measurement_test_value!(missed_release_runtime, [0.2], 0))
    prior_output = copy(missed_release_runtime.analog_output)
    MeasurementChainTest.prepare_measurement_analog_step!(
        missed_release_runtime,
        [0.3],
        1.0e-3,
    )
    refusal = try
        MeasurementChainTest.accept_measurement_analog_step!(missed_release_runtime, 2)
        nothing
    catch error
        error
    end
    @test refusal isa MeasurementChainTest.MeasurementChainRefusal
    @test refusal.code == :missed_measurement_release
    @test missed_release_runtime.analog_output == prior_output
    @test missed_release_runtime.last_accepted_tick == 0
    @test missed_release_runtime.candidate_active
    MeasurementChainTest.discard_measurement_analog_step!(missed_release_runtime)

    buffered_runtime = MeasurementChainTest.MeasurementChainRuntime(
        measurement_test_specification(),
    )
    released_buffer = MeasurementChainTest.MeasurementSample[]
    MeasurementChainTest.prepare_measurement_analog_step!(
        buffered_runtime,
        0.25,
        1.0e-3,
    )
    buffered_release = MeasurementChainTest.accept_measurement_analog_step!(
        released_buffer,
        buffered_runtime,
        0,
    )
    @test buffered_release === released_buffer
    @test length(buffered_release) == 1

    missed_acquisition_runtime = MeasurementChainTest.MeasurementChainRuntime(
        measurement_test_specification(
            acquisition=measurement_test_acquisition(sample_period_ticks=2),
        ),
    )
    only(accept_measurement_test_value!(missed_acquisition_runtime, [0.2], 0))
    MeasurementChainTest.prepare_measurement_analog_step!(
        missed_acquisition_runtime,
        [0.3],
        1.0e-3,
    )
    refusal = try
        MeasurementChainTest.accept_measurement_analog_step!(
            missed_acquisition_runtime,
            3,
        )
        nothing
    catch error
        error
    end
    @test refusal isa MeasurementChainTest.MeasurementChainRefusal
    @test refusal.code == :missed_measurement_acquisition
    @test missed_acquisition_runtime.last_accepted_tick == 0
end

@testset "causal three-phase measurement estimators" begin
    acquisition = measurement_test_acquisition(window_count=20)
    specification = measurement_test_specification(
        family=MeasurementChainTest.ThreePhaseSampledMeasurement,
        channel_names=[:voltage_a, :voltage_b, :voltage_c],
        phase_order=[:a, :b, :c],
        acquisition=acquisition,
    )
    runtime = MeasurementChainTest.MeasurementChainRuntime(specification)
    amplitude = 100.0 * sqrt(2.0)
    latest_sample = nothing
    for tick in 0:20
        angle = 2.0 * pi * 50.0 * tick * acquisition.tick_s
        values = amplitude .* [
            cos(angle),
            cos(angle - 2.0 * pi / 3.0),
            cos(angle + 2.0 * pi / 3.0),
        ]
        latest_sample = only(accept_measurement_test_value!(runtime, values, tick))
        if tick < 19
            @test latest_sample.quality == :window_incomplete
        end
    end
    @test latest_sample.quality == :valid
    @test latest_sample.sliding_rms ≈ fill(100.0, 3) atol=2.0e-13
    @test latest_sample.fundamental_rms_phasors ≈ ComplexF64[
        100.0,
        100.0 * cis(-2.0 * pi / 3.0),
        100.0 * cis(2.0 * pi / 3.0),
    ] atol=3.0e-13
    @test abs(latest_sample.sequence_phasors.zero) <= 2.0e-13
    @test latest_sample.sequence_phasors.positive ≈ 100.0 atol=3.0e-13
    @test abs(latest_sample.sequence_phasors.negative) <= 2.0e-13
    @test latest_sample.frequency_hz ≈ 50.0 atol=1.0e-12
end

@testset "measurement snapshot continuation" begin
    acquisition = measurement_test_acquisition(
        sample_period_ticks=2,
        delay_ticks=1,
        window_count=10,
    )
    specification = measurement_test_specification(acquisition=acquisition)
    uninterrupted = MeasurementChainTest.MeasurementChainRuntime(specification)
    for tick in 0:8
        accept_measurement_test_value!(uninterrupted, [sin(0.2 * tick)], tick)
    end
    snapshot = MeasurementChainTest.measurement_chain_snapshot(uninterrupted)
    restored = MeasurementChainTest.MeasurementChainRuntime(specification)
    MeasurementChainTest.restore_measurement_chain_snapshot!(restored, snapshot)
    @test MeasurementChainTest.measurement_chain_result_signature(restored) ==
        MeasurementChainTest.measurement_chain_result_signature(uninterrupted)

    for tick in 9:24
        values = [sin(0.2 * tick)]
        uninterrupted_samples = accept_measurement_test_value!(uninterrupted, values, tick)
        restored_samples = accept_measurement_test_value!(restored, values, tick)
        @test getfield.(uninterrupted_samples, :deterministic_signature_sha256) ==
            getfield.(restored_samples, :deterministic_signature_sha256)
    end
    @test MeasurementChainTest.measurement_chain_result_signature(restored) ==
        MeasurementChainTest.measurement_chain_result_signature(uninterrupted)

    stale_specification = measurement_test_specification(id=:stale_measurement_test)
    stale_runtime = MeasurementChainTest.MeasurementChainRuntime(stale_specification)
    @test_throws ArgumentError MeasurementChainTest.restore_measurement_chain_snapshot!(
        stale_runtime,
        snapshot,
    )
    @test stale_runtime.last_accepted_tick == -1

    invalid_state = merge(snapshot.state, (previous_input=[0.0, 1.0],))
    invalid_snapshot = MeasurementChainTest.MeasurementChainSnapshot(
        snapshot.schema_version,
        snapshot.specification_signature_sha256,
        invalid_state,
        MeasurementChainTest._measurement_snapshot_signature(
            snapshot.specification_signature_sha256,
            invalid_state,
        ),
    )
    unchanged_signature = MeasurementChainTest.measurement_chain_result_signature(restored)
    @test_throws DimensionMismatch MeasurementChainTest.restore_measurement_chain_snapshot!(
        restored,
        invalid_snapshot,
    )
    @test MeasurementChainTest.measurement_chain_result_signature(restored) ==
        unchanged_signature
end

@testset "measurement output retention and delayed-queue compaction" begin
    acquisition = measurement_test_acquisition(
        delay_ticks=4,
        maximum_retained_samples=2,
    )
    runtime = MeasurementChainTest.MeasurementChainRuntime(
        measurement_test_specification(acquisition=acquisition),
    )
    for tick in 0:2_100
        accept_measurement_test_value!(runtime, [sin(0.01 * tick)], tick)
    end
    @test length(runtime.samples) == 2
    @test getfield.(runtime.samples, :source_tick) == [2_095, 2_096]
    @test length(runtime.delayed_samples) <= 1_030
    @test runtime.delayed_sample_head <= length(runtime.delayed_samples) + 1
    @test_throws ArgumentError measurement_test_acquisition(
        maximum_retained_samples=0,
    )
end

function measurement_test_comtrade_record(; digital_count=3)
    analog_channels = MeasurementChainTest.ComtradeAnalogChannel[
        MeasurementChainTest.ComtradeAnalogChannel(
            1,
            "voltage_a";
            phase="A",
            circuit="generic_bus",
            unit="V",
            scale=0.5,
            offset=1.0,
            skew_us=0.0,
            minimum_raw=-100,
            maximum_raw=100,
            primary_ratio=120.0,
            secondary_ratio=1.0,
            primary_secondary=:secondary,
        ),
        MeasurementChainTest.ComtradeAnalogChannel(
            2,
            "current_a";
            phase="A",
            circuit="generic_branch",
            unit="A",
            scale=0.25,
            offset=-0.5,
            skew_us=2.0,
            minimum_raw=-100,
            maximum_raw=100,
            primary_ratio=600.0,
            secondary_ratio=1.0,
            primary_secondary=:primary,
        ),
    ]
    digital_channels = MeasurementChainTest.ComtradeDigitalChannel[
        MeasurementChainTest.ComtradeDigitalChannel(
            index,
            "status_$(index)";
            phase=index <= 3 ? string(Char('A' + index - 1)) : "",
            circuit="generic_status",
            normal_state=isodd(index),
        ) for index in 1:digital_count
    ]
    configuration = MeasurementChainTest.ComtradeConfiguration(
        "AIMORA synthetic station",
        "AIMORA measurement recorder",
        MeasurementChainTest.Comtrade2013;
        analog_channels=analog_channels,
        digital_channels=digital_channels,
        nominal_frequency_hz=50.0,
        sample_rates=[
            MeasurementChainTest.ComtradeSampleRate(1.0e3, 2),
            MeasurementChainTest.ComtradeSampleRate(2.0e3, 4),
        ],
        start_timestamp=MeasurementChainTest.ComtradeTimestamp(
            2026,
            8,
            14,
            12,
            30,
            0,
            123_456_000;
            fractional_digits=6,
        ),
        trigger_timestamp=MeasurementChainTest.ComtradeTimestamp(
            2026,
            8,
            14,
            12,
            30,
            0,
            125_000_000;
            fractional_digits=6,
        ),
        encoding=MeasurementChainTest.ComtradeASCII,
        time_multiplier=1.0,
    )
    raw_values = [
        -10.0 20.0
        0.0 10.0
        10.0 0.0
        20.0 -10.0
    ]
    digital_values = falses(4, digital_count)
    for row in axes(digital_values, 1), column in axes(digital_values, 2)
        digital_values[row, column] = isodd(row + column)
    end
    return MeasurementChainTest.ComtradeRecord(
        configuration,
        1:4,
        [0, 1_000, 2_000, 2_500],
        raw_values,
        digital_values,
    )
end

function measurement_test_binary_data(record, encoding)
    io = IOBuffer()
    digital_count = size(record.digital_values, 2)
    for row in eachindex(record.sample_numbers)
        write(io, htol(UInt32(record.sample_numbers[row])))
        write(io, htol(UInt32(record.timestamp_counts[row])))
        for value in view(record.raw_analog_values, row, :)
            if encoding === MeasurementChainTest.ComtradeBinary16
                write(io, htol(reinterpret(UInt16, Int16(value))))
            elseif encoding === MeasurementChainTest.ComtradeFloat32
                write(io, htol(reinterpret(UInt32, Float32(value))))
            else
                error("unsupported test encoding")
            end
        end
        for word_number in 1:cld(digital_count, 16)
            word = UInt16(0)
            for bit_index in 0:15
                channel = 16 * (word_number - 1) + bit_index + 1
                channel <= digital_count || break
                record.digital_values[row, channel] &&
                    (word |= UInt16(1) << bit_index)
            end
            write(io, htol(word))
        end
    end
    return take!(io)
end

@testset "bounded COMTRADE deterministic round trip" begin
    record = measurement_test_comtrade_record(digital_count=17)
    @test record.analog_values[:, 1] == 0.5 .* record.raw_analog_values[:, 1] .+ 1.0
    @test record.analog_values[:, 2] == 0.25 .* record.raw_analog_values[:, 2] .- 0.5
    @test record.time_s == [0.0, 1.0e-3, 2.0e-3, 2.5e-3]

    ascii = MeasurementChainTest.write_comtrade_record(
        record;
        encoding=MeasurementChainTest.ComtradeASCII,
    )
    repeated_ascii = MeasurementChainTest.write_comtrade_record(
        record;
        encoding=MeasurementChainTest.ComtradeASCII,
    )
    @test repeated_ascii.configuration_text == ascii.configuration_text
    @test repeated_ascii.data_bytes == ascii.data_bytes
    @test repeated_ascii.deterministic_signature_sha256 ==
        ascii.deterministic_signature_sha256
    parsed_ascii = MeasurementChainTest.read_comtrade_record(
        ascii.configuration_text,
        ascii.data_bytes,
    )
    @test parsed_ascii.configuration.revision === MeasurementChainTest.Comtrade2013
    @test parsed_ascii.configuration.encoding === MeasurementChainTest.ComtradeASCII
    @test parsed_ascii.configuration.station_name == record.configuration.station_name
    @test parsed_ascii.configuration.recording_device_id ==
        record.configuration.recording_device_id
    @test parsed_ascii.configuration.analog_channels == record.configuration.analog_channels
    @test parsed_ascii.configuration.digital_channels == record.configuration.digital_channels
    @test parsed_ascii.configuration.sample_rates == record.configuration.sample_rates
    @test parsed_ascii.configuration.start_timestamp == record.configuration.start_timestamp
    @test parsed_ascii.configuration.trigger_timestamp == record.configuration.trigger_timestamp
    @test parsed_ascii.raw_analog_values == record.raw_analog_values
    @test parsed_ascii.analog_values == record.analog_values
    @test parsed_ascii.digital_values == record.digital_values
    @test MeasurementChainTest.comtrade_record_signature(parsed_ascii) ==
        MeasurementChainTest.comtrade_record_signature(record)

    binary32 = MeasurementChainTest.write_comtrade_record(
        record;
        encoding=MeasurementChainTest.ComtradeBinary32,
    )
    parsed_binary32 = MeasurementChainTest.read_comtrade_record(
        binary32.configuration_text,
        binary32.data_bytes,
    )
    @test parsed_binary32.configuration.encoding === MeasurementChainTest.ComtradeBinary32
    @test parsed_binary32.sample_numbers == record.sample_numbers
    @test parsed_binary32.timestamp_counts == record.timestamp_counts
    @test parsed_binary32.raw_analog_values == record.raw_analog_values
    @test parsed_binary32.digital_values == record.digital_values

    mktempdir() do directory
        configuration_path = joinpath(directory, "synthetic_measurement.cfg")
        data_path = joinpath(directory, "synthetic_measurement.dat")
        file_serialization = MeasurementChainTest.write_comtrade_files(
            record,
            configuration_path,
            data_path;
            encoding=MeasurementChainTest.ComtradeBinary32,
        )
        @test read(configuration_path, String) == file_serialization.configuration_text
        @test read(data_path) == file_serialization.data_bytes
        from_files = MeasurementChainTest.read_comtrade_files(
            configuration_path,
            data_path,
        )
        @test from_files.raw_analog_values == record.raw_analog_values
        @test from_files.digital_values == record.digital_values
        @test_throws MeasurementChainTest.ComtradeRefusal MeasurementChainTest.write_comtrade_files(
            record,
            configuration_path,
            data_path,
        )
        @test_throws MeasurementChainTest.ComtradeRefusal MeasurementChainTest.read_comtrade_files(
            replace(configuration_path, ".cfg" => ".txt"),
            data_path,
        )
    end
end

@testset "bounded COMTRADE registered imports and refusals" begin
    record = measurement_test_comtrade_record()
    authored = MeasurementChainTest.write_comtrade_record(record)
    binary32_marker = "\nBINARY32\n"

    @test_throws ArgumentError MeasurementChainTest.ComtradeRecord(
        record.configuration,
        record.sample_numbers,
        [0, 1_000, 2_000, 3_000],
        record.raw_analog_values,
        record.digital_values,
    )

    for (encoding, marker) in (
        (MeasurementChainTest.ComtradeBinary16, "BINARY"),
        (MeasurementChainTest.ComtradeFloat32, "FLOAT32"),
    )
        configuration_text = replace(
            authored.configuration_text,
            "\nASCII\n" => "\n$(marker)\n",
        )
        parsed = MeasurementChainTest.read_comtrade_record(
            configuration_text,
            measurement_test_binary_data(record, encoding),
        )
        @test parsed.configuration.encoding === encoding
        @test parsed.raw_analog_values == record.raw_analog_values
        @test parsed.digital_values == record.digital_values
    end

    lines = split(chomp(authored.configuration_text), '\n')
    configuration_1999 = join(vcat(
        [replace(lines[1], ",2013" => ",1999")],
        lines[2:(end - 2)],
    ), '\n') * "\n"
    parsed_1999 = MeasurementChainTest.read_comtrade_record(
        configuration_1999,
        authored.data_bytes,
    )
    @test parsed_1999.configuration.revision === MeasurementChainTest.Comtrade1999

    lines_1991 = copy(lines[1:(end - 3)])
    lines_1991[1] = replace(lines_1991[1], ",2013" => "")
    timestamp_line = 2 + length(record.configuration.analog_channels) +
        length(record.configuration.digital_channels) + 2 +
        length(record.configuration.sample_rates) + 1
    lines_1991[timestamp_line] = replace(
        lines_1991[timestamp_line],
        "14/08/2026" => "08/14/2026",
    )
    lines_1991[timestamp_line + 1] = replace(
        lines_1991[timestamp_line + 1],
        "14/08/2026" => "08/14/2026",
    )
    configuration_1991 = join(lines_1991, '\n') * "\n"
    parsed_1991 = MeasurementChainTest.read_comtrade_record(
        configuration_1991,
        authored.data_bytes,
    )
    @test parsed_1991.configuration.revision === MeasurementChainTest.Comtrade1991
    @test parsed_1991.configuration.time_multiplier == 1.0

    unsupported_revision = replace(authored.configuration_text, ",2013" => ",2001")
    refusal = try
        MeasurementChainTest.read_comtrade_record(
            unsupported_revision,
            authored.data_bytes,
        )
        nothing
    catch error
        error
    end
    @test refusal isa MeasurementChainTest.ComtradeRefusal
    @test refusal.code == :unsupported_comtrade_revision

    unsupported_time_code = replace(
        authored.configuration_text,
        "\n,\n,\n" => "\nUTC,UTC\n,\n",
    )
    refusal = try
        MeasurementChainTest.read_comtrade_record(
            unsupported_time_code,
            authored.data_bytes,
        )
        nothing
    catch error
        error
    end
    @test refusal isa MeasurementChainTest.ComtradeRefusal
    @test refusal.code == :unsupported_time_code_metadata

    binary32 = MeasurementChainTest.write_comtrade_record(
        record;
        encoding=MeasurementChainTest.ComtradeBinary32,
    )
    @test_throws MeasurementChainTest.ComtradeRefusal MeasurementChainTest.read_comtrade_record(
        binary32.configuration_text,
        binary32.data_bytes[1:(end - 1)],
    )
    @test_throws MeasurementChainTest.ComtradeRefusal MeasurementChainTest.read_comtrade_record(
        authored.configuration_text,
        vcat(authored.data_bytes, collect(codeunits("5,3000,0,0,0,0,0\n"))),
    )
    @test_throws ArgumentError MeasurementChainTest.write_comtrade_record(
        record;
        encoding=MeasurementChainTest.ComtradeFloat32,
    )
    @test !occursin(binary32_marker, authored.configuration_text)
end

function measurement_test_transformer_source()
    return AIMORA.TransformerApparatus.TransformerSourceRecord(
        :synthetic_instrument_transformer,
        repeat("b", 64),
        measurement_test_provenance("SI"),
    )
end

function measurement_test_linear_transformer_specification()
    transformer = AIMORA.TransformerApparatus
    connection = transformer.TransformerConnectionTopology(
        node_order=[:primary_terminal, :secondary_terminal],
        coil_order=[:primary_coil, :secondary_coil],
        winding_order=[:primary_winding, :secondary_winding],
        phase_order=[:phase_a],
        coil_winding=[:primary_winding, :secondary_winding],
        coil_phase=[:phase_a, :phase_a],
        incidence=Matrix{Float64}(I, 2, 2),
        vector_group="Ii0",
        clock_number=0,
        phase_shift_rad=0.0,
    )
    matrices = transformer.TransformerTerminalMatrices(
        [0.1 0.0; 0.0 0.2],
        [0.2 0.18; 0.18 0.2],
    )
    return transformer.TransformerApparatusSpecification(
        :linear_instrument_transformer,
        transformer.LowFrequencyTerminalTier,
        connection,
        transformer.LowFrequencyTransformerModel(matrices),
        transformer.TransformerRuntimeSettings(
            timestep_s=10.0e-6,
            initialization_frequency_hz=50.0,
        );
        phase_count=1,
        rated_power_va=1.0e3,
        rated_voltage_v=100.0,
        rated_frequency_hz=50.0,
        sources=[measurement_test_transformer_source()],
        uncertainty="exact AIMORA synthetic instrument-transformer fixture",
        validity_domain="generic one-phase linear instrument-transformer test",
    )
end

function measurement_test_burden(; connected=true)
    return MeasurementChainTest.MeasurementBurden(
        connected=connected,
        series_resistance_ohm=connected ? 2.0 : 0.0,
        series_inductance_h=connected ? 1.0e-3 : 0.0,
        shunt_capacitance_f=connected ? 1.0e-7 : 0.0,
        cable_resistance_ohm=connected ? 0.1 : 0.0,
        cable_inductance_h=connected ? 10.0e-6 : 0.0,
        cable_capacitance_f=connected ? 10.0e-9 : 0.0,
        provenance=measurement_test_provenance("ohm,H,F"),
    )
end

function measurement_test_instrument_definition(
    family=MeasurementChainTest.LinearCurrentTransformerMeasurement;
    burden=measurement_test_burden(),
)
    return MeasurementChainTest.InstrumentTransformerMeasurementDefinition(
        :generic_instrument_transformer,
        family,
        measurement_test_linear_transformer_specification();
        primary_coil_index=1,
        secondary_coil_index=2,
        primary_turns=1.0,
        secondary_turns=1.0,
        burden=burden,
        maximum_secondary_impedance_ohm=100.0,
        secondary_output_sign=family in (
            MeasurementChainTest.InductiveVoltageTransformerMeasurement,
            MeasurementChainTest.CouplingCapacitorVoltageTransformerMeasurement,
        ) ? 1.0 : -1.0,
    )
end

@testset "instrument-transformer burden and coupled measurement runtime" begin
    definition = measurement_test_instrument_definition()
    readiness = MeasurementChainTest.instrument_transformer_measurement_readiness(
        definition,
    )
    @test readiness.ready
    @test isempty(readiness.reasons)
    @test isfinite(readiness.burden_impedance_ohm)
    burden_branches = MeasurementChainTest.measurement_burden_branches(
        definition.burden,
        2,
        0,
    )
    @test getfield.(burden_branches, :kind) == [:series_rl, :shunt_capacitance]
    @test MeasurementChainTest.measurement_burden_impedance(
        definition.burden,
        0.0,
    ) ≈ 2.1

    open_definition = measurement_test_instrument_definition(
        burden=measurement_test_burden(connected=false),
    )
    open_readiness = MeasurementChainTest.instrument_transformer_measurement_readiness(
        open_definition,
    )
    @test !open_readiness.ready
    @test :current_transformer_secondary_open in open_readiness.reasons
    @test isinf(abs(open_readiness.burden_impedance_ohm))

    transformer = AIMORA.TransformerApparatus
    preparation = transformer.prepare_transformer_apparatus(
        definition.apparatus;
        initialization_mode=transformer.DeenergizedTransformerInitialization,
    )
    acquisition = measurement_test_acquisition(
        tick_s=10.0e-6,
        window_count=20,
    )
    measurement_specification = measurement_test_specification(
        family=MeasurementChainTest.LinearCurrentTransformerMeasurement,
        channel_names=[:secondary_current],
        quantity=:current,
        unit="A",
        orientation="positive from the secondary dotted terminal into the burden",
        acquisition=acquisition,
    )
    runtime = MeasurementChainTest.instrument_transformer_measurement_runtime(
        definition,
        preparation,
        [1, 2],
        measurement_specification,
    )
    @test length(runtime.measurement_runtime.samples) == 1
    @test only(runtime.measurement_runtime.samples).source_tick == 0
    @test AIMORA.NonlinearNetwork.nonlinear_terminal_nodes(runtime) == [1, 2]

    current = zeros(2)
    jacobian = zeros(2, 2)
    AIMORA.NonlinearNetwork.prepare_nonlinear_device_step!(
        runtime,
        10.0e-6,
        10.0e-6,
        :trapezoidal,
    )
    AIMORA.NonlinearNetwork.nonlinear_current_jacobian!(
        current,
        jacobian,
        runtime,
        [10.0, 0.0],
        10.0e-6,
    )
    AIMORA.NonlinearNetwork.accept_nonlinear_device_state!(
        runtime,
        [10.0, 0.0],
        current,
        jacobian,
        10.0e-6,
    )
    AIMORA.NonlinearNetwork.finish_nonlinear_device_step!(runtime)
    output = MeasurementChainTest.instrument_transformer_measurement_output(runtime)
    @test output.accepted_time_s == 10.0e-6
    @test output.held_measurement == -output.secondary_current_a
    @test output.primary_referred_measurement == -output.secondary_current_a
    @test output.latest_sample.source_tick == 1
    @test length(output.deterministic_signature_sha256) == 64
    @test all(isfinite, current)
    @test all(isfinite, jacobian)

    snapshot = MeasurementChainTest.instrument_transformer_measurement_snapshot(runtime)
    restored = MeasurementChainTest.instrument_transformer_measurement_runtime(
        definition,
        preparation,
        [1, 2],
        measurement_specification,
    )
    MeasurementChainTest.restore_instrument_transformer_measurement_snapshot!(
        restored,
        snapshot,
    )
    @test MeasurementChainTest.instrument_transformer_measurement_signature(restored) ==
        MeasurementChainTest.instrument_transformer_measurement_signature(runtime)

    voltage_definition = measurement_test_instrument_definition(
        MeasurementChainTest.InductiveVoltageTransformerMeasurement,
    )
    @test MeasurementChainTest.instrument_transformer_measurement_readiness(
        voltage_definition,
    ).ready
    @test_throws ArgumentError MeasurementChainTest.InstrumentTransformerMeasurementDefinition(
        :invalid_magnetic_ct,
        MeasurementChainTest.MagneticCurrentTransformerMeasurement,
        definition.apparatus;
        primary_coil_index=1,
        secondary_coil_index=2,
        primary_turns=1.0,
        secondary_turns=1.0,
        burden=definition.burden,
        maximum_secondary_impedance_ohm=100.0,
        secondary_output_sign=-1.0,
    )
end

function measurement_test_electronic_sensor_definition(
    family=MeasurementChainTest.ElectronicVoltageSensorMeasurement;
    coupling=family === MeasurementChainTest.ElectronicVoltageSensorMeasurement ?
        MeasurementChainTest.VoltageShuntElectronicLoading :
        MeasurementChainTest.NonloadingElectronicObservation,
)
    transducer = MeasurementChainTest.AnalogMeasurementStateSpace(
        reshape([-100.0], 1, 1),
        [100.0],
        [1.0],
        0.0,
    )
    return MeasurementChainTest.ElectronicSensorDefinition(
        family === MeasurementChainTest.ElectronicVoltageSensorMeasurement ?
            :generic_electronic_voltage_sensor : :generic_electronic_current_sensor,
        family;
        coupling=coupling,
        transducer=transducer,
        input_scale=family === MeasurementChainTest.ElectronicVoltageSensorMeasurement ?
            0.01 : 1.0,
        input_offset=0.0,
        minimum_observed_input=-1.0e4,
        maximum_observed_input=1.0e4,
        loading_resistance_ohm=coupling ===
            MeasurementChainTest.VoltageShuntElectronicLoading ? 1.0e6 : 0.0,
        loading_capacitance_f=coupling ===
            MeasurementChainTest.VoltageShuntElectronicLoading ? 10.0e-12 : 0.0,
        provenance=measurement_test_provenance("V or A"),
    )
end

@testset "electronic sensor transducer state, loading, and restart" begin
    definition = measurement_test_electronic_sensor_definition()
    loading = MeasurementChainTest.electronic_sensor_loading_branches(
        definition,
        1,
        0,
    )
    @test getfield.(loading, :kind) == [:resistance, :shunt_capacitance]
    acquisition = measurement_test_acquisition()
    specification = measurement_test_specification(
        family=MeasurementChainTest.ElectronicVoltageSensorMeasurement,
        channel_names=[:sensor_voltage],
        quantity=:voltage,
        unit="V",
        conditioning=definition.transducer,
        acquisition=acquisition,
    )
    runtime = MeasurementChainTest.electronic_sensor_runtime(
        definition,
        specification;
        initial_observed_input=0.0,
    )
    @test only(runtime.measurement_runtime.samples).source_tick == 0
    MeasurementChainTest.prepare_electronic_sensor_step!(runtime, 100.0, 1.0e-3)
    @test only(runtime.measurement_runtime.candidate.output) ≈ 1.0 / 21.0 atol=1.0e-15
    MeasurementChainTest.discard_electronic_sensor_step!(runtime)
    @test runtime.accepted_observed_input == 0.0
    @test runtime.candidate_observed_input === nothing

    MeasurementChainTest.prepare_electronic_sensor_step!(runtime, 100.0, 1.0e-3)
    sample = only(MeasurementChainTest.accept_electronic_sensor_step!(runtime, 1))
    output = MeasurementChainTest.electronic_sensor_output(runtime)
    @test output.observed_input == 100.0
    @test output.transducer_input == 1.0
    @test output.conditioned_output ≈ 1.0 / 21.0 atol=1.0e-15
    @test output.held_measurement == sample.instantaneous[1]
    @test output.latest_sample.source_tick == 1
    @test length(output.deterministic_signature_sha256) == 64

    snapshot = MeasurementChainTest.electronic_sensor_snapshot(runtime)
    restored = MeasurementChainTest.electronic_sensor_runtime(
        definition,
        specification;
        initial_observed_input=0.0,
    )
    MeasurementChainTest.restore_electronic_sensor_snapshot!(restored, snapshot)
    @test MeasurementChainTest.electronic_sensor_signature(restored) ==
        MeasurementChainTest.electronic_sensor_signature(runtime)

    current_definition = measurement_test_electronic_sensor_definition(
        MeasurementChainTest.ElectronicCurrentSensorMeasurement,
    )
    @test isempty(MeasurementChainTest.electronic_sensor_loading_branches(
        current_definition,
        1,
        0,
    ))
    @test_throws ArgumentError measurement_test_electronic_sensor_definition(
        MeasurementChainTest.ElectronicCurrentSensorMeasurement;
        coupling=MeasurementChainTest.VoltageShuntElectronicLoading,
    )
    @test_throws MeasurementChainTest.MeasurementChainRefusal MeasurementChainTest.prepare_electronic_sensor_step!(
        restored,
        2.0e4,
        1.0e-3,
    )
end

function measurement_test_cvt_definition(; maximum_spectral_frequency_hz=500.0)
    electromagnetic_unit = measurement_test_instrument_definition(
        MeasurementChainTest.CouplingCapacitorVoltageTransformerMeasurement,
    )
    high_voltage_capacitance = 1.0e-9
    intermediate_voltage_capacitance = 10.0e-9
    equivalent_capacitance = inv(
        inv(high_voltage_capacitance) + inv(intermediate_voltage_capacitance),
    )
    compensation_inductance = inv((2.0 * pi * 50.0)^2 * equivalent_capacitance)
    return MeasurementChainTest.CouplingCapacitorVoltageTransformerDefinition(
        :generic_coupling_capacitor_voltage_transformer,
        electromagnetic_unit;
        high_voltage_capacitance_f=high_voltage_capacitance,
        intermediate_voltage_capacitance_f=intermediate_voltage_capacitance,
        compensation_resistance_ohm=25.0,
        compensation_inductance_h=compensation_inductance,
        suppression_resistance_ohm=1.0e5,
        suppression_capacitance_f=1.0e-9,
        maximum_spectral_frequency_hz=maximum_spectral_frequency_hz,
        maximum_timestep_s=10.0e-6,
        provenance=measurement_test_provenance("F,ohm,H"),
    )
end

@testset "explicit CVT divider, compensation, suppression, and measurement" begin
    definition = measurement_test_cvt_definition()
    readiness = MeasurementChainTest.cvt_measurement_readiness(definition)
    @test readiness.ready
    @test readiness.divider_ratio ≈ 1.0 / 11.0 atol=1.0e-16
    @test readiness.series_resonance_hz ≈ 50.0 atol=1.0e-12
    branches = MeasurementChainTest.cvt_network_branches(
        definition;
        line_node=3,
        divider_node=4,
        electromagnetic_primary_node=1,
        secondary_node=2,
    )
    @test getfield.(branches, :id) == [
        :high_voltage_coupling_capacitor,
        :intermediate_voltage_divider_capacitor,
        :compensation_reactor,
        :ferroresonance_suppression_resistance,
        :ferroresonance_suppression_capacitance,
        :secondary_burden_1,
        :secondary_burden_2,
    ]
    @test all(branch -> branch.positive_node != branch.negative_node, branches)
    @test_throws ArgumentError measurement_test_cvt_definition(
        maximum_spectral_frequency_hz=40.0,
    )

    transformer = AIMORA.TransformerApparatus
    preparation = transformer.prepare_transformer_apparatus(
        definition.electromagnetic_unit.apparatus;
        initialization_mode=transformer.DeenergizedTransformerInitialization,
    )
    acquisition = measurement_test_acquisition(
        tick_s=10.0e-6,
        window_count=20,
    )
    specification = measurement_test_specification(
        family=MeasurementChainTest.CouplingCapacitorVoltageTransformerMeasurement,
        channel_names=[:secondary_voltage],
        quantity=:voltage,
        unit="V",
        orientation="secondary dotted terminal to reference",
        acquisition=acquisition,
    )
    runtime = MeasurementChainTest.cvt_measurement_runtime(
        definition,
        preparation,
        [1, 2],
        specification,
    )
    current = zeros(2)
    jacobian = zeros(2, 2)
    AIMORA.NonlinearNetwork.prepare_nonlinear_device_step!(
        runtime,
        10.0e-6,
        10.0e-6,
        :trapezoidal,
    )
    AIMORA.NonlinearNetwork.nonlinear_current_jacobian!(
        current,
        jacobian,
        runtime,
        [20.0, 2.0],
        10.0e-6,
    )
    AIMORA.NonlinearNetwork.accept_nonlinear_device_state!(
        runtime,
        [20.0, 2.0],
        current,
        jacobian,
        10.0e-6,
    )
    AIMORA.NonlinearNetwork.finish_nonlinear_device_step!(runtime)
    output = MeasurementChainTest.cvt_measurement_output(
        definition,
        runtime;
        line_voltage_v=220.0,
        divider_voltage_v=20.0,
        electromagnetic_primary_voltage_v=20.0,
        compensation_current_a=0.1,
    )
    @test output.measurement.family ===
        MeasurementChainTest.CouplingCapacitorVoltageTransformerMeasurement
    @test output.stored_energy_j > 0.0
    @test output.divider_ratio == readiness.divider_ratio
    @test output.series_resonance_hz == readiness.series_resonance_hz
    @test length(output.deterministic_signature_sha256) == 64
end

function measurement_test_magnetic_transformer_specification()
    transformer = AIMORA.TransformerApparatus
    source = measurement_test_transformer_source()
    lower_curve = transformer.TellinenLimitingCurve(
        [-1_000.0, -500.0, 0.0, 500.0, 1_000.0],
        [-1.5, -1.0, -0.2, 0.5, 1.0],
    )
    upper_curve = transformer.TellinenLimitingCurve(
        [-1_000.0, -500.0, 0.0, 500.0, 1_000.0],
        [-1.0, -0.5, 0.2, 1.0, 1.5],
    )
    material = transformer.TellinenTransformerMagneticMaterial(
        lower_curve,
        upper_curve,
        source;
        integration_field_increment_a_per_m=2.0,
    )
    graph = transformer.TransformerMagneticGraph(
        node_order=[:magnetic_node],
        branches=[
            transformer.MagneticBranchGeometry(
                :measurement_core_limb;
                length_m=0.8,
                cross_section_m2=0.01,
            ),
            transformer.MagneticBranchGeometry(
                :measurement_return_limb;
                length_m=0.8,
                cross_section_m2=0.01,
            ),
        ],
        incidence=reshape([1.0, -1.0], 1, 2),
        winding_turns=[10.0 0.0; 0.0 100.0],
        materials=[material],
    )
    connection = transformer.TransformerConnectionTopology(
        node_order=[:primary_terminal, :secondary_terminal],
        coil_order=[:primary_coil, :secondary_coil],
        winding_order=[:primary_winding, :secondary_winding],
        phase_order=[:phase_a],
        coil_winding=[:primary_winding, :secondary_winding],
        coil_phase=[:phase_a, :phase_a],
        incidence=Matrix{Float64}(I, 2, 2),
        vector_group="Ii0",
        clock_number=0,
        phase_shift_rad=0.0,
    )
    model = transformer.MagneticEquivalentCircuitModel(
        [0.01 0.0; 0.0 1.0],
        [1.0e-4 0.0; 0.0 1.0e-3],
        graph,
    )
    return transformer.TransformerApparatusSpecification(
        :magnetic_current_transformer,
        transformer.MagneticEquivalentCircuitTier,
        connection,
        model,
        transformer.TransformerRuntimeSettings(
            timestep_s=10.0e-6,
            initialization_frequency_hz=50.0,
        );
        phase_count=1,
        rated_power_va=1.0e3,
        rated_voltage_v=100.0,
        rated_frequency_hz=50.0,
        sources=[source],
        uncertainty="exact AIMORA synthetic Tellinen current-transformer fixture",
        validity_domain="generic one-phase saturating/remanent CT test",
    )
end

@testset "saturating and remanent CT through canonical magnetic apparatus" begin
    apparatus = measurement_test_magnetic_transformer_specification()
    definition = MeasurementChainTest.InstrumentTransformerMeasurementDefinition(
        :generic_saturating_current_transformer,
        MeasurementChainTest.MagneticCurrentTransformerMeasurement,
        apparatus;
        primary_coil_index=1,
        secondary_coil_index=2,
        primary_turns=10.0,
        secondary_turns=100.0,
        burden=measurement_test_burden(),
        maximum_secondary_impedance_ohm=100.0,
        secondary_output_sign=-1.0,
    )
    @test MeasurementChainTest.instrument_transformer_measurement_readiness(
        definition,
    ).ready
    transformer = AIMORA.TransformerApparatus
    residual_flux = [1.0e-3, 1.0e-3]
    preparation = transformer.prepare_transformer_apparatus(
        apparatus;
        initialization_mode=transformer.DeenergizedTransformerInitialization,
        initial_branch_flux_wb=residual_flux,
    )
    @test preparation.initial_branch_flux_wb == residual_flux
    specification = measurement_test_specification(
        family=MeasurementChainTest.MagneticCurrentTransformerMeasurement,
        channel_names=[:secondary_current],
        quantity=:current,
        unit="A",
        orientation="positive from the secondary dotted terminal into the burden",
        acquisition=measurement_test_acquisition(
            tick_s=10.0e-6,
            window_count=20,
        ),
    )
    runtime = MeasurementChainTest.instrument_transformer_measurement_runtime(
        definition,
        preparation,
        [1, 2],
        specification,
    )
    current = zeros(2)
    jacobian = zeros(2, 2)
    AIMORA.NonlinearNetwork.prepare_nonlinear_device_step!(
        runtime,
        10.0e-6,
        10.0e-6,
        :trapezoidal,
    )
    AIMORA.NonlinearNetwork.nonlinear_current_jacobian!(
        current,
        jacobian,
        runtime,
        [2.0, 0.0],
        10.0e-6,
    )
    AIMORA.NonlinearNetwork.accept_nonlinear_device_state!(
        runtime,
        [2.0, 0.0],
        current,
        jacobian,
        10.0e-6,
    )
    AIMORA.NonlinearNetwork.finish_nonlinear_device_step!(runtime)
    output = MeasurementChainTest.instrument_transformer_measurement_output(runtime)
    @test output.family === MeasurementChainTest.MagneticCurrentTransformerMeasurement
    @test all(isfinite, current)
    @test all(isfinite, jacobian)
    @test isfinite(output.ampere_turn_residual_at)
    @test runtime.apparatus_runtime.accepted_state.branch_flux_wb != residual_flux
    @test all(
        state -> state !== nothing,
        runtime.apparatus_runtime.accepted_state.tellinen_state,
    )
    @test length(MeasurementChainTest.instrument_transformer_measurement_signature(runtime)) == 64
end
