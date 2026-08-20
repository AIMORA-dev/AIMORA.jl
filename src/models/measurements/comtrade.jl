export ComtradeRevision,
       Comtrade1991,
       Comtrade1999,
       Comtrade2013,
       ComtradeEncoding,
       ComtradeASCII,
       ComtradeBinary16,
       ComtradeBinary32,
       ComtradeFloat32,
       ComtradeTimestamp,
       ComtradeAnalogChannel,
       ComtradeDigitalChannel,
       ComtradeSampleRate,
       ComtradeConfiguration,
       ComtradeRecord,
       ComtradeReadLimits,
       ComtradeSerializedRecord,
       ComtradeRefusal,
       read_comtrade_record,
       read_comtrade_files,
       write_comtrade_record,
       write_comtrade_files,
       comtrade_record_signature

@enum ComtradeRevision begin
    Comtrade1991
    Comtrade1999
    Comtrade2013
end

@enum ComtradeEncoding begin
    ComtradeASCII
    ComtradeBinary16
    ComtradeBinary32
    ComtradeFloat32
end

const _COMTRADE_REVISION_TEXT = Dict(
    Comtrade1991 => "1991",
    Comtrade1999 => "1999",
    Comtrade2013 => "2013",
)

const _COMTRADE_ENCODING_TEXT = Dict(
    ComtradeASCII => "ASCII",
    ComtradeBinary16 => "BINARY",
    ComtradeBinary32 => "BINARY32",
    ComtradeFloat32 => "FLOAT32",
)

struct ComtradeRefusal <: Exception
    code::Symbol
    message::String
    line::Int
    field::Int
    byte_offset::Int
end

function ComtradeRefusal(
    code::Symbol,
    message::AbstractString;
    line::Integer=0,
    field::Integer=0,
    byte_offset::Integer=-1,
)
    return ComtradeRefusal(
        code,
        String(message),
        Int(line),
        Int(field),
        Int(byte_offset),
    )
end

function Base.showerror(io::IO, refusal::ComtradeRefusal)
    print(io, String(refusal.code), ": ", refusal.message)
    refusal.line > 0 && print(io, " at CFG/DAT line ", refusal.line)
    refusal.field > 0 && print(io, ", field ", refusal.field)
    refusal.byte_offset >= 0 && print(io, ", byte offset ", refusal.byte_offset)
end

function _validate_comtrade_text(value::AbstractString, role::AbstractString; empty=false)
    text = String(value)
    empty || !isempty(strip(text)) || throw(ArgumentError("$role must not be empty"))
    occursin(',', text) && throw(ArgumentError("$role must not contain a comma"))
    occursin('\n', text) && throw(ArgumentError("$role must not contain a newline"))
    occursin('\r', text) && throw(ArgumentError("$role must not contain a newline"))
    return text
end

struct ComtradeTimestamp
    year::Int
    month::Int
    day::Int
    hour::Int
    minute::Int
    second::Int
    nanosecond::Int
    fractional_digits::Int

    function ComtradeTimestamp(
        year::Integer,
        month::Integer,
        day::Integer,
        hour::Integer,
        minute::Integer,
        second::Integer,
        nanosecond::Integer=0;
        fractional_digits::Integer=6,
    )
        date = Date(Int(year), Int(month), Int(day))
        0 <= hour <= 23 || throw(ArgumentError("COMTRADE hour must be in 0 through 23"))
        0 <= minute <= 59 || throw(ArgumentError("COMTRADE minute must be in 0 through 59"))
        0 <= second <= 59 || throw(ArgumentError("COMTRADE second must be in 0 through 59"))
        0 <= nanosecond < 1_000_000_000 || throw(ArgumentError(
            "COMTRADE fractional second must be below one second",
        ))
        0 <= fractional_digits <= 9 || throw(ArgumentError(
            "COMTRADE fractional-second precision must be 0 through 9 digits",
        ))
        precision = Int(fractional_digits)
        quantum = 10^(9 - precision)
        Int(nanosecond) % quantum == 0 || throw(ArgumentError(
            "COMTRADE fractional second is not representable at its declared precision",
        ))
        return new(
            Dates.year(date),
            Dates.month(date),
            Dates.day(date),
            Int(hour),
            Int(minute),
            Int(second),
            Int(nanosecond),
            precision,
        )
    end
end

function _comtrade_timestamp_nanoseconds(timestamp::ComtradeTimestamp)
    day_offset = Int128(Dates.value(
        Date(timestamp.year, timestamp.month, timestamp.day) - Date(1, 1, 1),
    ))
    second_of_day =
        Int128(3600 * timestamp.hour + 60 * timestamp.minute + timestamp.second)
    return (day_offset * 86_400 + second_of_day) * 1_000_000_000 +
        timestamp.nanosecond
end

struct ComtradeAnalogChannel
    index::Int
    id::String
    phase::String
    circuit::String
    unit::String
    scale::Float64
    offset::Float64
    skew_us::Float64
    minimum_raw::Int64
    maximum_raw::Int64
    primary_ratio::Float64
    secondary_ratio::Float64
    primary_secondary::Symbol

    function ComtradeAnalogChannel(
        index::Integer,
        id::AbstractString;
        phase::AbstractString="",
        circuit::AbstractString="",
        unit::AbstractString,
        scale::Real,
        offset::Real=0.0,
        skew_us::Real=0.0,
        minimum_raw::Integer,
        maximum_raw::Integer,
        primary_ratio::Real=1.0,
        secondary_ratio::Real=1.0,
        primary_secondary::Symbol=:secondary,
    )
        index > 0 || throw(ArgumentError("COMTRADE analog index must be positive"))
        identifier = _validate_comtrade_text(id, "COMTRADE analog channel id")
        phase_text = _validate_comtrade_text(phase, "COMTRADE analog phase"; empty=true)
        circuit_text = _validate_comtrade_text(
            circuit,
            "COMTRADE analog circuit";
            empty=true,
        )
        unit_text = _validate_comtrade_text(unit, "COMTRADE analog unit")
        numeric = Float64.((scale, offset, skew_us, primary_ratio, secondary_ratio))
        all(isfinite, numeric) || throw(ArgumentError(
            "COMTRADE analog scaling, skew, and ratios must be finite",
        ))
        numeric[1] != 0.0 || throw(ArgumentError("COMTRADE analog scale must be nonzero"))
        numeric[4] > 0.0 && numeric[5] > 0.0 || throw(ArgumentError(
            "COMTRADE primary and secondary ratios must be positive",
        ))
        minimum = Int64(minimum_raw)
        maximum = Int64(maximum_raw)
        minimum < maximum || throw(ArgumentError(
            "COMTRADE raw analog limits must be strictly ordered",
        ))
        primary_secondary in (:primary, :secondary) || throw(ArgumentError(
            "COMTRADE analog primary/secondary flag must be :primary or :secondary",
        ))
        return new(
            Int(index),
            identifier,
            phase_text,
            circuit_text,
            unit_text,
            numeric[1],
            numeric[2],
            numeric[3],
            minimum,
            maximum,
            numeric[4],
            numeric[5],
            primary_secondary,
        )
    end
end

struct ComtradeDigitalChannel
    index::Int
    id::String
    phase::String
    circuit::String
    normal_state::Bool

    function ComtradeDigitalChannel(
        index::Integer,
        id::AbstractString;
        phase::AbstractString="",
        circuit::AbstractString="",
        normal_state::Bool=false,
    )
        index > 0 || throw(ArgumentError("COMTRADE digital index must be positive"))
        return new(
            Int(index),
            _validate_comtrade_text(id, "COMTRADE digital channel id"),
            _validate_comtrade_text(phase, "COMTRADE digital phase"; empty=true),
            _validate_comtrade_text(circuit, "COMTRADE digital circuit"; empty=true),
            normal_state,
        )
    end
end

struct ComtradeSampleRate
    rate_hz::Float64
    final_sample::Int

    function ComtradeSampleRate(rate_hz::Real, final_sample::Integer)
        rate = Float64(rate_hz)
        isfinite(rate) && rate > 0.0 || throw(ArgumentError(
            "COMTRADE sample rate must be finite and positive",
        ))
        final_sample > 0 || throw(ArgumentError(
            "COMTRADE sample-rate final sample must be positive",
        ))
        return new(rate, Int(final_sample))
    end
end

struct ComtradeConfiguration
    station_name::String
    recording_device_id::String
    revision::ComtradeRevision
    analog_channels::Vector{ComtradeAnalogChannel}
    digital_channels::Vector{ComtradeDigitalChannel}
    nominal_frequency_hz::Float64
    sample_rates::Vector{ComtradeSampleRate}
    start_timestamp::ComtradeTimestamp
    trigger_timestamp::ComtradeTimestamp
    encoding::ComtradeEncoding
    time_multiplier::Float64

    function ComtradeConfiguration(
        station_name::AbstractString,
        recording_device_id::AbstractString,
        revision::ComtradeRevision;
        analog_channels::AbstractVector{ComtradeAnalogChannel},
        digital_channels::AbstractVector{ComtradeDigitalChannel},
        nominal_frequency_hz::Real,
        sample_rates::AbstractVector{ComtradeSampleRate},
        start_timestamp::ComtradeTimestamp,
        trigger_timestamp::ComtradeTimestamp,
        encoding::ComtradeEncoding,
        time_multiplier::Real=1.0,
    )
        station = _validate_comtrade_text(station_name, "COMTRADE station name")
        device = _validate_comtrade_text(
            recording_device_id,
            "COMTRADE recording-device id",
        )
        analog = collect(analog_channels)
        digital = collect(digital_channels)
        !isempty(analog) || !isempty(digital) || throw(ArgumentError(
            "COMTRADE configuration must declare at least one channel",
        ))
        getfield.(analog, :index) == collect(eachindex(analog)) || throw(ArgumentError(
            "COMTRADE analog channel indices must be contiguous and ordered",
        ))
        getfield.(digital, :index) == collect(eachindex(digital)) || throw(ArgumentError(
            "COMTRADE digital channel indices must be contiguous and ordered",
        ))
        channel_identities = vcat(
            [(kind=:analog, id=channel.id) for channel in analog],
            [(kind=:digital, id=channel.id) for channel in digital],
        )
        length(channel_identities) == length(unique(channel_identities)) || throw(
            ArgumentError("COMTRADE channel identities must be unique within each kind"),
        )
        frequency = Float64(nominal_frequency_hz)
        isfinite(frequency) && frequency > 0.0 || throw(ArgumentError(
            "COMTRADE nominal frequency must be finite and positive",
        ))
        rates = collect(sample_rates)
        !isempty(rates) || throw(ArgumentError(
            "COMTRADE configuration must declare at least one sample-rate section",
        ))
        final_samples = getfield.(rates, :final_sample)
        issorted(final_samples; lt=<) || throw(ArgumentError(
            "COMTRADE sample-rate section endpoints must increase strictly",
        ))
        multiplier = Float64(time_multiplier)
        isfinite(multiplier) && multiplier > 0.0 || throw(ArgumentError(
            "COMTRADE time multiplier must be finite and positive",
        ))
        _comtrade_timestamp_nanoseconds(trigger_timestamp) >=
            _comtrade_timestamp_nanoseconds(start_timestamp) || throw(ArgumentError(
            "COMTRADE trigger timestamp must not precede the record start",
        ))
        return new(
            station,
            device,
            revision,
            analog,
            digital,
            frequency,
            rates,
            start_timestamp,
            trigger_timestamp,
            encoding,
            multiplier,
        )
    end
end

struct ComtradeRecord
    configuration::ComtradeConfiguration
    sample_numbers::Vector{Int64}
    timestamp_counts::Vector{Int64}
    time_s::Vector{Float64}
    raw_analog_values::Matrix{Float64}
    analog_values::Matrix{Float64}
    digital_values::BitMatrix

    function ComtradeRecord(
        configuration::ComtradeConfiguration,
        sample_numbers::AbstractVector{<:Integer},
        timestamp_counts::AbstractVector{<:Integer},
        raw_analog_values::AbstractMatrix{<:Real},
        digital_values::AbstractMatrix{Bool},
    )
        samples = Int64.(sample_numbers)
        timestamps = Int64.(timestamp_counts)
        raw = Matrix{Float64}(raw_analog_values)
        digital = BitMatrix(digital_values)
        sample_count = configuration.sample_rates[end].final_sample
        length(samples) == sample_count || throw(DimensionMismatch(
            "COMTRADE sample-number count does not match its final sample-rate endpoint",
        ))
        length(timestamps) == sample_count || throw(DimensionMismatch(
            "COMTRADE timestamp count does not match its final sample-rate endpoint",
        ))
        size(raw) == (sample_count, length(configuration.analog_channels)) ||
            throw(DimensionMismatch("COMTRADE raw analog matrix has incompatible dimensions"))
        size(digital) == (sample_count, length(configuration.digital_channels)) ||
            throw(DimensionMismatch("COMTRADE digital matrix has incompatible dimensions"))
        samples == collect(Int64, 1:sample_count) || throw(ArgumentError(
            "the bounded COMTRADE subset requires contiguous one-based sample numbers",
        ))
        all(>=(0), timestamps) &&
            (length(timestamps) <= 1 || all(diff(timestamps) .> 0)) || throw(ArgumentError(
            "COMTRADE timestamp counts must be nonnegative and increase strictly",
        ))
        all(isfinite, raw) || throw(ArgumentError(
            "the bounded COMTRADE subset does not accept missing or nonfinite analog samples",
        ))
        engineering = similar(raw)
        for (column, channel) in enumerate(configuration.analog_channels)
            all(value -> channel.minimum_raw <= value <= channel.maximum_raw, view(raw, :, column)) ||
                throw(ArgumentError(
                    "COMTRADE raw analog sample is outside its declared channel limits",
                ))
            @views engineering[:, column] .= channel.scale .* raw[:, column] .+ channel.offset
        end
        all(isfinite, engineering) || throw(ArgumentError(
            "COMTRADE analog scaling produced a nonfinite engineering value",
        ))
        timestamp_tick_s = 1.0e-6 * configuration.time_multiplier
        previous_final_sample = 0
        for section in configuration.sample_rates
            first_sample = previous_final_sample + 1
            for sample_index in (first_sample + 1):section.final_sample
                observed_interval_s = timestamp_tick_s *
                    (timestamps[sample_index] - timestamps[sample_index - 1])
                expected_interval_s = inv(section.rate_hz)
                tolerance_s = 0.5 * timestamp_tick_s + 64.0 * eps(Float64) *
                    max(observed_interval_s, expected_interval_s, 1.0)
                abs(observed_interval_s - expected_interval_s) <= tolerance_s ||
                    throw(ArgumentError(
                        "COMTRADE timestamp interval contradicts its sample-rate section",
                    ))
            end
            previous_final_sample = section.final_sample
        end
        times = timestamp_tick_s .* Float64.(timestamps)
        trigger_offset_s = Float64(
            _comtrade_timestamp_nanoseconds(configuration.trigger_timestamp) -
                _comtrade_timestamp_nanoseconds(configuration.start_timestamp),
        ) * 1.0e-9
        trigger_offset_s <= last(times) + timestamp_tick_s || throw(ArgumentError(
            "COMTRADE trigger timestamp lies after the recorded samples",
        ))
        return new(configuration, samples, timestamps, times, raw, engineering, digital)
    end
end

struct ComtradeReadLimits
    maximum_channels::Int
    maximum_samples::Int
    maximum_data_bytes::Int

    function ComtradeReadLimits(;
        maximum_channels::Integer=16_384,
        maximum_samples::Integer=10_000_000,
        maximum_data_bytes::Integer=2^30,
    )
        maximum_channels > 0 || throw(ArgumentError("COMTRADE channel limit must be positive"))
        maximum_samples > 0 || throw(ArgumentError("COMTRADE sample limit must be positive"))
        maximum_data_bytes > 0 || throw(ArgumentError("COMTRADE byte limit must be positive"))
        return new(Int(maximum_channels), Int(maximum_samples), Int(maximum_data_bytes))
    end
end

struct ComtradeSerializedRecord
    configuration_text::String
    data_bytes::Vector{UInt8}
    encoding::ComtradeEncoding
    deterministic_signature_sha256::String
end

mutable struct _ComtradeLineCursor
    lines::Vector{String}
    index::Int
end

function _next_comtrade_line!(cursor::_ComtradeLineCursor, role::AbstractString)
    cursor.index < length(cursor.lines) || throw(ComtradeRefusal(
        :truncated_configuration,
        "missing $role",
        line=cursor.index + 1,
    ))
    cursor.index += 1
    return strip(replace(cursor.lines[cursor.index], '\x1a' => ""))
end

_comtrade_fields(line::AbstractString) = String.(strip.(split(line, ','; keepempty=true)))

function _require_comtrade_field_count(fields, expected, line, role)
    length(fields) == expected || throw(ComtradeRefusal(
        :invalid_configuration_field_count,
        "$role requires $expected comma-separated fields, observed $(length(fields))",
        line=line,
    ))
    return fields
end

function _parse_comtrade_integer(value, line, field, role; type=Int64)
    parsed = tryparse(type, strip(value))
    parsed === nothing && throw(ComtradeRefusal(
        :invalid_configuration_integer,
        "$role is not a representable integer",
        line=line,
        field=field,
    ))
    return parsed
end

function _parse_comtrade_float(value, line, field, role)
    parsed = tryparse(Float64, strip(value))
    parsed !== nothing && isfinite(parsed) || throw(ComtradeRefusal(
        :invalid_configuration_number,
        "$role is not a finite decimal number",
        line=line,
        field=field,
    ))
    return parsed
end

function _parse_comtrade_revision(value, line)
    for (revision, text) in _COMTRADE_REVISION_TEXT
        strip(value) == text && return revision
    end
    throw(ComtradeRefusal(
        :unsupported_comtrade_revision,
        "only COMTRADE 1991, 1999, and 2013 are registered",
        line=line,
        field=3,
    ))
end

function _parse_comtrade_encoding(value, line)
    normalized = uppercase(strip(value))
    for (encoding, text) in _COMTRADE_ENCODING_TEXT
        normalized == text && return encoding
    end
    throw(ComtradeRefusal(
        :unsupported_comtrade_encoding,
        "only ASCII, BINARY, BINARY32, and FLOAT32 DAT encodings are registered",
        line=line,
    ))
end

function _parse_comtrade_timestamp(value, revision, line)
    fields = _require_comtrade_field_count(
        _comtrade_fields(value),
        2,
        line,
        "COMTRADE timestamp",
    )
    date_fields = split(fields[1], '/'; keepempty=true)
    time_fields = split(fields[2], ':'; keepempty=true)
    length(date_fields) == 3 && length(time_fields) == 3 || throw(ComtradeRefusal(
        :invalid_comtrade_timestamp,
        "COMTRADE timestamp must contain a three-field date and time",
        line=line,
    ))
    date_values = Int[
        _parse_comtrade_integer(item, line, index, "COMTRADE date field"; type=Int)
        for (index, item) in enumerate(date_fields)
    ]
    day, month, year = revision === Comtrade1991 ?
        (date_values[2], date_values[1], date_values[3]) : Tuple(date_values)
    second_parts = split(time_fields[3], '.'; keepempty=true)
    1 <= length(second_parts) <= 2 || throw(ComtradeRefusal(
        :invalid_comtrade_timestamp,
        "COMTRADE timestamp seconds contain more than one decimal point",
        line=line,
    ))
    hour = _parse_comtrade_integer(time_fields[1], line, 2, "COMTRADE hour"; type=Int)
    minute = _parse_comtrade_integer(time_fields[2], line, 2, "COMTRADE minute"; type=Int)
    second = _parse_comtrade_integer(second_parts[1], line, 2, "COMTRADE second"; type=Int)
    fractional_text = length(second_parts) == 2 ? second_parts[2] : ""
    all(isdigit, fractional_text) || throw(ComtradeRefusal(
        :invalid_comtrade_timestamp,
        "COMTRADE fractional second must contain decimal digits only",
        line=line,
    ))
    length(fractional_text) <= 9 || throw(ComtradeRefusal(
        :unsupported_timestamp_precision,
        "COMTRADE timestamps above nanosecond precision are unsupported",
        line=line,
    ))
    nanosecond = isempty(fractional_text) ? 0 :
        parse(Int, rpad(fractional_text, 9, '0'))
    try
        return ComtradeTimestamp(
            year,
            month,
            day,
            hour,
            minute,
            second,
            nanosecond;
            fractional_digits=length(fractional_text),
        )
    catch error
        throw(ComtradeRefusal(
            :invalid_comtrade_timestamp,
            sprint(showerror, error),
            line=line,
        ))
    end
end

function _parse_comtrade_configuration(configuration_text::AbstractString, limits)
    isvalid(configuration_text) || throw(ComtradeRefusal(
        :invalid_configuration_encoding,
        "COMTRADE CFG text is not valid UTF-8",
    ))
    normalized = replace(String(configuration_text), "\r\n" => "\n", '\r' => '\n')
    cursor = _ComtradeLineCursor(String.(split(normalized, '\n'; keepempty=true)), 0)
    identity_fields = _comtrade_fields(_next_comtrade_line!(cursor, "station identity"))
    revision = if length(identity_fields) == 2
        Comtrade1991
    elseif length(identity_fields) == 3
        _parse_comtrade_revision(identity_fields[3], cursor.index)
    else
        throw(ComtradeRefusal(
            :invalid_configuration_field_count,
            "COMTRADE station identity requires two or three fields",
            line=cursor.index,
        ))
    end
    channel_fields = _require_comtrade_field_count(
        _comtrade_fields(_next_comtrade_line!(cursor, "channel counts")),
        3,
        cursor.index,
        "COMTRADE channel counts",
    )
    total_count = _parse_comtrade_integer(
        channel_fields[1],
        cursor.index,
        1,
        "COMTRADE total channel count";
        type=Int,
    )
    analog_match = match(r"^(\d+)[Aa]$", channel_fields[2])
    digital_match = match(r"^(\d+)[Dd]$", channel_fields[3])
    analog_match !== nothing && digital_match !== nothing || throw(ComtradeRefusal(
        :invalid_channel_count_kind,
        "COMTRADE typed channel counts must end in A and D",
        line=cursor.index,
    ))
    analog_count = parse(Int, only(analog_match.captures))
    digital_count = parse(Int, only(digital_match.captures))
    total_count == analog_count + digital_count || throw(ComtradeRefusal(
        :inconsistent_channel_count,
        "COMTRADE total channel count does not equal analog plus digital counts",
        line=cursor.index,
    ))
    0 < total_count <= limits.maximum_channels || throw(ComtradeRefusal(
        :channel_limit_exceeded,
        "COMTRADE channel count is outside the configured resource limit",
        line=cursor.index,
    ))

    analog_channels = ComtradeAnalogChannel[]
    for expected_index in 1:analog_count
        fields = _comtrade_fields(_next_comtrade_line!(cursor, "analog channel"))
        allowed_count = revision === Comtrade1991 ? (10, 13) : (13,)
        length(fields) in allowed_count || throw(ComtradeRefusal(
            :invalid_configuration_field_count,
            "COMTRADE analog channel has an unsupported field count",
            line=cursor.index,
        ))
        index = _parse_comtrade_integer(
            fields[1],
            cursor.index,
            1,
            "COMTRADE analog channel index";
            type=Int,
        )
        index == expected_index || throw(ComtradeRefusal(
            :unordered_channel_index,
            "COMTRADE analog channel indices must be contiguous and ordered",
            line=cursor.index,
            field=1,
        ))
        primary = length(fields) == 13 ?
            _parse_comtrade_float(fields[11], cursor.index, 11, "primary ratio") : 1.0
        secondary = length(fields) == 13 ?
            _parse_comtrade_float(fields[12], cursor.index, 12, "secondary ratio") : 1.0
        flag = length(fields) == 13 ? uppercase(fields[13]) : "S"
        flag in ("P", "S") || throw(ComtradeRefusal(
            :invalid_primary_secondary_flag,
            "COMTRADE analog channel flag must be P or S",
            line=cursor.index,
            field=13,
        ))
        try
            push!(
                analog_channels,
                ComtradeAnalogChannel(
                    index,
                    fields[2];
                    phase=fields[3],
                    circuit=fields[4],
                    unit=fields[5],
                    scale=_parse_comtrade_float(fields[6], cursor.index, 6, "analog scale"),
                    offset=_parse_comtrade_float(fields[7], cursor.index, 7, "analog offset"),
                    skew_us=_parse_comtrade_float(fields[8], cursor.index, 8, "analog skew"),
                    minimum_raw=_parse_comtrade_integer(
                        fields[9], cursor.index, 9, "analog minimum raw value",
                    ),
                    maximum_raw=_parse_comtrade_integer(
                        fields[10], cursor.index, 10, "analog maximum raw value",
                    ),
                    primary_ratio=primary,
                    secondary_ratio=secondary,
                    primary_secondary=flag == "P" ? :primary : :secondary,
                ),
            )
        catch error
            error isa ComtradeRefusal && rethrow()
            throw(ComtradeRefusal(
                :invalid_analog_channel,
                sprint(showerror, error),
                line=cursor.index,
            ))
        end
    end

    digital_channels = ComtradeDigitalChannel[]
    for expected_index in 1:digital_count
        fields = _comtrade_fields(_next_comtrade_line!(cursor, "digital channel"))
        allowed_count = revision === Comtrade1991 ? (3, 5) : (5,)
        length(fields) in allowed_count || throw(ComtradeRefusal(
            :invalid_configuration_field_count,
            "COMTRADE digital channel has an unsupported field count",
            line=cursor.index,
        ))
        index = _parse_comtrade_integer(
            fields[1], cursor.index, 1, "COMTRADE digital channel index"; type=Int,
        )
        index == expected_index || throw(ComtradeRefusal(
            :unordered_channel_index,
            "COMTRADE digital channel indices must be contiguous and ordered",
            line=cursor.index,
            field=1,
        ))
        normal_text = length(fields) == 5 ? fields[5] : "0"
        normal_text in ("0", "1") || throw(ComtradeRefusal(
            :invalid_digital_normal_state,
            "COMTRADE digital normal state must be zero or one",
            line=cursor.index,
            field=5,
        ))
        try
            push!(
                digital_channels,
                ComtradeDigitalChannel(
                    index,
                    fields[2];
                    phase=fields[3],
                    circuit=length(fields) == 5 ? fields[4] : "",
                    normal_state=normal_text == "1",
                ),
            )
        catch error
            throw(ComtradeRefusal(
                :invalid_digital_channel,
                sprint(showerror, error),
                line=cursor.index,
            ))
        end
    end

    nominal_frequency = _parse_comtrade_float(
        _next_comtrade_line!(cursor, "nominal frequency"),
        cursor.index,
        1,
        "COMTRADE nominal frequency",
    )
    rate_count = _parse_comtrade_integer(
        _next_comtrade_line!(cursor, "sample-rate count"),
        cursor.index,
        1,
        "COMTRADE sample-rate count";
        type=Int,
    )
    rate_count > 0 || throw(ComtradeRefusal(
        :timestamp_only_rate_section_unsupported,
        "the bounded COMTRADE subset requires one or more explicit sample rates",
        line=cursor.index,
    ))
    sample_rates = ComtradeSampleRate[]
    for _ in 1:rate_count
        fields = _require_comtrade_field_count(
            _comtrade_fields(_next_comtrade_line!(cursor, "sample-rate section")),
            2,
            cursor.index,
            "COMTRADE sample-rate section",
        )
        try
            push!(
                sample_rates,
                ComtradeSampleRate(
                    _parse_comtrade_float(fields[1], cursor.index, 1, "sample rate"),
                    _parse_comtrade_integer(
                        fields[2], cursor.index, 2, "sample-rate final sample"; type=Int,
                    ),
                ),
            )
        catch error
            error isa ComtradeRefusal && rethrow()
            throw(ComtradeRefusal(
                :invalid_sample_rate_section,
                sprint(showerror, error),
                line=cursor.index,
            ))
        end
    end
    sample_rates[end].final_sample <= limits.maximum_samples || throw(ComtradeRefusal(
        :sample_limit_exceeded,
        "COMTRADE final sample exceeds the configured resource limit",
        line=cursor.index,
    ))
    start_timestamp = _parse_comtrade_timestamp(
        _next_comtrade_line!(cursor, "start timestamp"),
        revision,
        cursor.index,
    )
    trigger_timestamp = _parse_comtrade_timestamp(
        _next_comtrade_line!(cursor, "trigger timestamp"),
        revision,
        cursor.index,
    )
    encoding = _parse_comtrade_encoding(
        _next_comtrade_line!(cursor, "DAT encoding"),
        cursor.index,
    )
    time_multiplier = revision === Comtrade1991 ? 1.0 : _parse_comtrade_float(
        _next_comtrade_line!(cursor, "time multiplier"),
        cursor.index,
        1,
        "COMTRADE time multiplier",
    )
    if revision === Comtrade2013
        time_code = _next_comtrade_line!(cursor, "time-code metadata")
        leap_second = _next_comtrade_line!(cursor, "time-quality metadata")
        all(isempty, _comtrade_fields(time_code)) || throw(ComtradeRefusal(
            :unsupported_time_code_metadata,
            "time-code and local-code interpretation is outside the registered subset",
            line=cursor.index - 1,
        ))
        all(isempty, _comtrade_fields(leap_second)) || throw(ComtradeRefusal(
            :unsupported_leap_second_metadata,
            "time-quality and leap-second interpretation is outside the registered subset",
            line=cursor.index,
        ))
    end
    for line_index in (cursor.index + 1):length(cursor.lines)
        isempty(strip(replace(cursor.lines[line_index], '\x1a' => ""))) || throw(
            ComtradeRefusal(
                :unexpected_configuration_content,
                "unexpected content follows the registered COMTRADE configuration",
                line=line_index,
            ),
        )
    end
    try
        return ComtradeConfiguration(
            identity_fields[1],
            identity_fields[2],
            revision;
            analog_channels=analog_channels,
            digital_channels=digital_channels,
            nominal_frequency_hz=nominal_frequency,
            sample_rates=sample_rates,
            start_timestamp=start_timestamp,
            trigger_timestamp=trigger_timestamp,
            encoding=encoding,
            time_multiplier=time_multiplier,
        )
    catch error
        throw(ComtradeRefusal(:invalid_configuration, sprint(showerror, error)))
    end
end

function _parse_ascii_comtrade_data(configuration, data_text, limits)
    isvalid(data_text) || throw(ComtradeRefusal(
        :invalid_data_encoding,
        "COMTRADE ASCII DAT content is not valid UTF-8",
    ))
    normalized = replace(String(data_text), "\r\n" => "\n", '\r' => '\n')
    lines = String.(split(normalized, '\n'; keepempty=true))
    while !isempty(lines) && isempty(strip(last(lines)))
        pop!(lines)
    end
    sample_count = configuration.sample_rates[end].final_sample
    length(lines) == sample_count || throw(ComtradeRefusal(
        :inconsistent_data_sample_count,
        "COMTRADE DAT row count does not match the final sample-rate endpoint",
    ))
    analog_count = length(configuration.analog_channels)
    digital_count = length(configuration.digital_channels)
    samples = Vector{Int64}(undef, sample_count)
    timestamps = Vector{Int64}(undef, sample_count)
    raw = Matrix{Float64}(undef, sample_count, analog_count)
    digital = falses(sample_count, digital_count)
    expected_fields = 2 + analog_count + digital_count
    for row in 1:sample_count
        fields = _require_comtrade_field_count(
            _comtrade_fields(lines[row]),
            expected_fields,
            row,
            "COMTRADE ASCII data row",
        )
        samples[row] = _parse_comtrade_integer(
            fields[1], row, 1, "COMTRADE sample number",
        )
        timestamps[row] = _parse_comtrade_integer(
            fields[2], row, 2, "COMTRADE timestamp count",
        )
        for column in 1:analog_count
            raw[row, column] = _parse_comtrade_float(
                fields[2 + column], row, 2 + column, "COMTRADE raw analog sample",
            )
        end
        for column in 1:digital_count
            field_index = 2 + analog_count + column
            fields[field_index] in ("0", "1") || throw(ComtradeRefusal(
                :invalid_digital_sample,
                "COMTRADE digital samples must be zero or one",
                line=row,
                field=field_index,
            ))
            digital[row, column] = fields[field_index] == "1"
        end
    end
    return ComtradeRecord(configuration, samples, timestamps, raw, digital)
end

_read_little_endian(io::IO, ::Type{UInt16}) = ltoh(read(io, UInt16))
_read_little_endian(io::IO, ::Type{UInt32}) = ltoh(read(io, UInt32))

function _parse_binary_comtrade_data(configuration, data_bytes, limits)
    length(data_bytes) <= limits.maximum_data_bytes || throw(ComtradeRefusal(
        :data_byte_limit_exceeded,
        "COMTRADE DAT content exceeds the configured byte limit",
    ))
    analog_count = length(configuration.analog_channels)
    digital_count = length(configuration.digital_channels)
    digital_word_count = cld(digital_count, 16)
    analog_bytes = configuration.encoding === ComtradeBinary16 ? 2 : 4
    row_bytes = 8 + analog_bytes * analog_count + 2 * digital_word_count
    sample_count = configuration.sample_rates[end].final_sample
    expected_bytes = Base.checked_mul(row_bytes, sample_count)
    length(data_bytes) == expected_bytes || throw(ComtradeRefusal(
        :inconsistent_binary_data_size,
        "COMTRADE binary DAT byte count does not match its configured rows",
        byte_offset=min(length(data_bytes), expected_bytes),
    ))
    samples = Vector{Int64}(undef, sample_count)
    timestamps = Vector{Int64}(undef, sample_count)
    raw = Matrix{Float64}(undef, sample_count, analog_count)
    digital = falses(sample_count, digital_count)
    io = IOBuffer(data_bytes)
    for row in 1:sample_count
        row_offset = position(io)
        samples[row] = Int64(_read_little_endian(io, UInt32))
        timestamps[row] = Int64(_read_little_endian(io, UInt32))
        for column in 1:analog_count
            if configuration.encoding === ComtradeBinary16
                value = reinterpret(Int16, _read_little_endian(io, UInt16))
                value == typemin(Int16) && throw(ComtradeRefusal(
                    :missing_analog_sample_unsupported,
                    "the bounded COMTRADE subset refuses binary missing-value sentinels",
                    line=row,
                    field=2 + column,
                    byte_offset=row_offset + 8 + 2 * (column - 1),
                ))
                raw[row, column] = Float64(value)
            elseif configuration.encoding === ComtradeBinary32
                value = reinterpret(Int32, _read_little_endian(io, UInt32))
                value == typemin(Int32) && throw(ComtradeRefusal(
                    :missing_analog_sample_unsupported,
                    "the bounded COMTRADE subset refuses binary32 missing-value sentinels",
                    line=row,
                    field=2 + column,
                    byte_offset=row_offset + 8 + 4 * (column - 1),
                ))
                raw[row, column] = Float64(value)
            else
                value = reinterpret(Float32, _read_little_endian(io, UInt32))
                isfinite(value) || throw(ComtradeRefusal(
                    :nonfinite_analog_sample,
                    "COMTRADE FLOAT32 analog sample must be finite",
                    line=row,
                    field=2 + column,
                    byte_offset=row_offset + 8 + 4 * (column - 1),
                ))
                raw[row, column] = Float64(value)
            end
        end
        for word_number in 1:digital_word_count
            word_index = word_number - 1
            word = _read_little_endian(io, UInt16)
            for bit_index in 0:15
                channel = 16 * word_index + bit_index + 1
                channel <= digital_count || break
                digital[row, channel] = !iszero(word & (UInt16(1) << bit_index))
            end
        end
    end
    return ComtradeRecord(configuration, samples, timestamps, raw, digital)
end

function read_comtrade_record(
    configuration_text::AbstractString,
    data::Union{AbstractString,AbstractVector{UInt8}};
    limits::ComtradeReadLimits=ComtradeReadLimits(),
)
    configuration = _parse_comtrade_configuration(configuration_text, limits)
    data_size = data isa AbstractString ? ncodeunits(data) : length(data)
    data_size <= limits.maximum_data_bytes || throw(ComtradeRefusal(
        :data_byte_limit_exceeded,
        "COMTRADE DAT content exceeds the configured byte limit",
    ))
    if configuration.encoding === ComtradeASCII
        data_text = data isa AbstractString ? data : try
            String(copy(data))
        catch error
            throw(ComtradeRefusal(:invalid_data_encoding, sprint(showerror, error)))
        end
        return _parse_ascii_comtrade_data(configuration, data_text, limits)
    end
    data isa AbstractVector{UInt8} || throw(ComtradeRefusal(
        :binary_data_requires_bytes,
        "binary COMTRADE DAT content must be supplied as bytes",
    ))
    return _parse_binary_comtrade_data(configuration, data, limits)
end

function _validate_comtrade_pair_paths(configuration_path, data_path)
    lowercase(splitext(configuration_path)[2]) == ".cfg" || throw(ComtradeRefusal(
        :unknown_configuration_extension,
        "COMTRADE configuration path must use the .cfg extension",
    ))
    lowercase(splitext(data_path)[2]) == ".dat" || throw(ComtradeRefusal(
        :unknown_data_extension,
        "COMTRADE data path must use the .dat extension",
    ))
    splitext(basename(configuration_path))[1] == splitext(basename(data_path))[1] ||
        throw(ComtradeRefusal(
            :mismatched_file_pair,
            "COMTRADE CFG and DAT paths must share one basename",
        ))
    return nothing
end

function _read_bounded_comtrade_file(path, role, maximum_bytes)
    isfile(path) || throw(ComtradeRefusal(
        Symbol("missing_", role, "_file"),
        "COMTRADE $(uppercase(role)) file does not exist",
    ))
    filesize(path) <= maximum_bytes || throw(ComtradeRefusal(
        Symbol(role, "_byte_limit_exceeded"),
        "COMTRADE $(uppercase(role)) file exceeds the configured byte limit",
    ))
    try
        return read(path)
    catch error
        throw(ComtradeRefusal(
            Symbol(role, "_file_read_failed"),
            sprint(showerror, error),
        ))
    end
end

function read_comtrade_files(
    configuration_path::AbstractString,
    data_path::AbstractString;
    limits::ComtradeReadLimits=ComtradeReadLimits(),
)
    configuration_file = abspath(configuration_path)
    data_file = abspath(data_path)
    _validate_comtrade_pair_paths(configuration_file, data_file)
    configuration_bytes = _read_bounded_comtrade_file(
        configuration_file,
        "configuration",
        limits.maximum_data_bytes,
    )
    configuration_text = try
        String(configuration_bytes)
    catch error
        throw(ComtradeRefusal(:invalid_configuration_encoding, sprint(showerror, error)))
    end
    data_bytes = _read_bounded_comtrade_file(
        data_file,
        "data",
        limits.maximum_data_bytes,
    )
    return read_comtrade_record(configuration_text, data_bytes; limits=limits)
end

function _format_comtrade_number(value::Real)
    isfinite(value) || throw(ArgumentError("COMTRADE output values must be finite"))
    return @sprintf("%.17g", Float64(value))
end

function _format_comtrade_timestamp(timestamp::ComtradeTimestamp)
    date = string(
        lpad(timestamp.day, 2, '0'),
        '/',
        lpad(timestamp.month, 2, '0'),
        '/',
        lpad(timestamp.year, 4, '0'),
    )
    time = string(
        lpad(timestamp.hour, 2, '0'),
        ':',
        lpad(timestamp.minute, 2, '0'),
        ':',
        lpad(timestamp.second, 2, '0'),
    )
    if timestamp.fractional_digits > 0
        fractional = lpad(timestamp.nanosecond, 9, '0')[1:timestamp.fractional_digits]
        time *= "." * fractional
    end
    return date * "," * time
end

function _write_comtrade_configuration(record, encoding)
    configuration = record.configuration
    io = IOBuffer()
    println(io, configuration.station_name, ',', configuration.recording_device_id, ",2013")
    analog_count = length(configuration.analog_channels)
    digital_count = length(configuration.digital_channels)
    println(io, analog_count + digital_count, ',', analog_count, "A,", digital_count, 'D')
    for channel in configuration.analog_channels
        println(
            io,
            channel.index,
            ',',
            channel.id,
            ',',
            channel.phase,
            ',',
            channel.circuit,
            ',',
            channel.unit,
            ',',
            _format_comtrade_number(channel.scale),
            ',',
            _format_comtrade_number(channel.offset),
            ',',
            _format_comtrade_number(channel.skew_us),
            ',',
            channel.minimum_raw,
            ',',
            channel.maximum_raw,
            ',',
            _format_comtrade_number(channel.primary_ratio),
            ',',
            _format_comtrade_number(channel.secondary_ratio),
            ',',
            channel.primary_secondary === :primary ? 'P' : 'S',
        )
    end
    for channel in configuration.digital_channels
        println(
            io,
            channel.index,
            ',',
            channel.id,
            ',',
            channel.phase,
            ',',
            channel.circuit,
            ',',
            channel.normal_state ? 1 : 0,
        )
    end
    println(io, _format_comtrade_number(configuration.nominal_frequency_hz))
    println(io, length(configuration.sample_rates))
    for rate in configuration.sample_rates
        println(io, _format_comtrade_number(rate.rate_hz), ',', rate.final_sample)
    end
    println(io, _format_comtrade_timestamp(configuration.start_timestamp))
    println(io, _format_comtrade_timestamp(configuration.trigger_timestamp))
    println(io, _COMTRADE_ENCODING_TEXT[encoding])
    println(io, _format_comtrade_number(configuration.time_multiplier))
    println(io, ',')
    println(io, ',')
    return String(take!(io))
end

function _write_little_endian(io::IO, value::UInt16)
    write(io, htol(value))
end

function _write_little_endian(io::IO, value::UInt32)
    write(io, htol(value))
end

function _write_comtrade_ascii_data(record)
    io = IOBuffer()
    analog_count = size(record.raw_analog_values, 2)
    digital_count = size(record.digital_values, 2)
    for row in eachindex(record.sample_numbers)
        print(io, record.sample_numbers[row], ',', record.timestamp_counts[row])
        for column in 1:analog_count
            print(io, ',', _format_comtrade_number(record.raw_analog_values[row, column]))
        end
        for column in 1:digital_count
            print(io, ',', record.digital_values[row, column] ? 1 : 0)
        end
        print(io, '\n')
    end
    return take!(io)
end

function _write_comtrade_binary32_data(record)
    io = IOBuffer()
    digital_count = size(record.digital_values, 2)
    digital_word_count = cld(digital_count, 16)
    for row in eachindex(record.sample_numbers)
        sample = record.sample_numbers[row]
        timestamp = record.timestamp_counts[row]
        0 <= sample <= typemax(UInt32) || throw(ArgumentError(
            "COMTRADE sample number is outside the BINARY32 unsigned range",
        ))
        0 <= timestamp <= typemax(UInt32) || throw(ArgumentError(
            "COMTRADE timestamp count is outside the BINARY32 unsigned range",
        ))
        _write_little_endian(io, UInt32(sample))
        _write_little_endian(io, UInt32(timestamp))
        for value in view(record.raw_analog_values, row, :)
            rounded = round(Int64, value)
            value == rounded && typemin(Int32) < rounded <= typemax(Int32) || throw(
                ArgumentError(
                    "COMTRADE BINARY32 export requires exact non-sentinel Int32 raw values",
                ),
            )
            _write_little_endian(io, reinterpret(UInt32, Int32(rounded)))
        end
        for word_number in 1:digital_word_count
            word_index = word_number - 1
            word = UInt16(0)
            for bit_index in 0:15
                channel = 16 * word_index + bit_index + 1
                channel <= digital_count || break
                record.digital_values[row, channel] &&
                    (word |= UInt16(1) << bit_index)
            end
            _write_little_endian(io, word)
        end
    end
    return take!(io)
end

function write_comtrade_record(
    record::ComtradeRecord;
    encoding::ComtradeEncoding=ComtradeASCII,
)
    encoding in (ComtradeASCII, ComtradeBinary32) || throw(ArgumentError(
        "deterministic AIMORA COMTRADE export supports only 2013 ASCII and BINARY32",
    ))
    configuration_text = _write_comtrade_configuration(record, encoding)
    data_bytes = encoding === ComtradeASCII ?
        _write_comtrade_ascii_data(record) : _write_comtrade_binary32_data(record)
    signature = bytes2hex(sha256(vcat(
        collect(codeunits(configuration_text)),
        UInt8(0),
        data_bytes,
    )))
    return ComtradeSerializedRecord(
        configuration_text,
        data_bytes,
        encoding,
        signature,
    )
end

function write_comtrade_files(
    record::ComtradeRecord,
    configuration_path::AbstractString,
    data_path::AbstractString;
    encoding::ComtradeEncoding=ComtradeASCII,
)
    configuration_file = abspath(configuration_path)
    data_file = abspath(data_path)
    _validate_comtrade_pair_paths(configuration_file, data_file)
    ispath(configuration_file) && throw(ComtradeRefusal(
        :configuration_file_exists,
        "COMTRADE writer refuses to replace an existing CFG file",
    ))
    ispath(data_file) && throw(ComtradeRefusal(
        :data_file_exists,
        "COMTRADE writer refuses to replace an existing DAT file",
    ))
    isdir(dirname(configuration_file)) || throw(ComtradeRefusal(
        :missing_configuration_directory,
        "COMTRADE CFG parent directory does not exist",
    ))
    isdir(dirname(data_file)) || throw(ComtradeRefusal(
        :missing_data_directory,
        "COMTRADE DAT parent directory does not exist",
    ))
    serialized = write_comtrade_record(record; encoding=encoding)
    configuration_temporary, configuration_io =
        mktemp(dirname(configuration_file); cleanup=false)
    data_temporary, data_io = mktemp(dirname(data_file); cleanup=false)
    data_published = false
    try
        write(configuration_io, serialized.configuration_text)
        write(data_io, serialized.data_bytes)
        close(configuration_io)
        close(data_io)
        mv(data_temporary, data_file)
        data_published = true
        mv(configuration_temporary, configuration_file)
        return serialized
    catch error
        isopen(configuration_io) && close(configuration_io)
        isopen(data_io) && close(data_io)
        isfile(configuration_temporary) && rm(configuration_temporary)
        isfile(data_temporary) && rm(data_temporary)
        data_published && !isfile(configuration_file) && isfile(data_file) && rm(data_file)
        error isa ComtradeRefusal && rethrow()
        throw(ComtradeRefusal(:file_write_failed, sprint(showerror, error)))
    end
end

function comtrade_record_signature(record::ComtradeRecord)
    serialized = write_comtrade_record(record; encoding=ComtradeASCII)
    return serialized.deterministic_signature_sha256
end
