module Sources

using ..Branches: CurrentInjection, TheveninSource

export constant_current_injection,
       constant_thevenin_source,
       AnalyticSourceSignal,
       analytic_current_injection_source,
       analytic_source_value,
       analytic_thevenin_source,
       sinusoidal_thevenin_source,
       sinusoidal_value

mutable struct AnalyticSourceSignal
    source_type::Int
    crest::Float64
    time1_s::Float64
    angular_frequency_or_rate::Float64
    start_time_s::Float64
    stop_time_s::Float64

    function AnalyticSourceSignal(
        source_type::Integer,
        crest::Real,
        time1_s::Real,
        angular_frequency_or_rate::Real,
        start_time_s::Real,
        stop_time_s::Real,
    )
        source_crest = Float64(crest)
        source_time1 = Float64(time1_s)
        source_rate = Float64(angular_frequency_or_rate)
        source_start = Float64(start_time_s)
        source_stop = Float64(stop_time_s)
        return new(
            Int(source_type),
            source_crest,
            source_time1,
            source_rate,
            source_start,
            source_stop,
        )
    end
end

function (signal::AnalyticSourceSignal)(time_s::Real)
    return analytic_source_value(
        signal.source_type,
        signal.crest,
        signal.time1_s,
        signal.angular_frequency_or_rate,
        signal.start_time_s,
        signal.stop_time_s,
        time_s,
    )
end

"""
    sinusoidal_value(t, amplitude_pu, frequency_hz; phase_rad=0.0, offset_pu=0.0)

Return a scalar per-unit sinusoid used by the small EMT source builders.
"""
function sinusoidal_value(t::Real, amplitude_pu::Real, frequency_hz::Real;
                          phase_rad::Real=0.0, offset_pu::Real=0.0)::Float64
    return Float64(offset_pu) +
           Float64(amplitude_pu) *
           sin(2.0 * pi * Float64(frequency_hz) * Float64(t) + Float64(phase_rad))
end

function constant_thevenin_source(node::Int, conductance::Real, voltage_pu::Real)
    value = Float64(voltage_pu)
    return TheveninSource(node, Float64(conductance), _t -> value)
end

function sinusoidal_thevenin_source(node::Int, conductance::Real, amplitude_pu::Real,
                                    frequency_hz::Real; phase_rad::Real=0.0,
                                    offset_pu::Real=0.0)
    amplitude = Float64(amplitude_pu)
    frequency = Float64(frequency_hz)
    phase = Float64(phase_rad)
    offset = Float64(offset_pu)
    return TheveninSource(
        node,
        Float64(conductance),
        t -> sinusoidal_value(t, amplitude, frequency; phase_rad=phase, offset_pu=offset),
    )
end

function analytic_source_value(source_type::Integer, crest::Real, time1::Real,
                               sfreq::Real, tstart::Real, tstop::Real,
                               t::Real)::Float64
    row_type = Int(source_type)
    start_time = Float64(tstart)
    stop_time = Float64(tstop)
    current_time = Float64(t)
    if current_time < start_time || current_time >= stop_time
        return 0.0
    end

    ts = current_time - start_time
    source_crest = Float64(crest)
    source_time1 = Float64(time1)
    source_sfreq = Float64(sfreq)
    if row_type == 11
        return source_crest
    elseif row_type == 12
        source_time1 <= 0.0 && return source_crest
        return ts < source_time1 ? ts / source_time1 * source_crest : source_crest
    elseif row_type == 13
        if source_time1 > 0.0 && ts < source_time1
            return ts / source_time1 * source_crest
        end
        return source_crest + source_sfreq * (ts - max(source_time1, 0.0))
    elseif row_type == 14
        return source_crest * cos(source_sfreq * ts + source_time1)
    elseif row_type == 15
        return source_crest * (exp(source_sfreq * ts) - exp(source_time1 * ts))
    end
    throw(ArgumentError("unsupported OVER5A source type $row_type"))
end

function analytic_thevenin_source(node::Int, conductance::Real, source_type::Integer,
                                  crest::Real, time1::Real, sfreq::Real,
                                  tstart::Real, tstop::Real)
    signal = AnalyticSourceSignal(
        source_type,
        crest,
        time1,
        sfreq,
        tstart,
        tstop,
    )
    return TheveninSource(
        node,
        Float64(conductance),
        signal,
    )
end

function analytic_current_injection_source(node::Int, source_type::Integer, crest::Real,
                                           time1::Real, sfreq::Real,
                                           tstart::Real, tstop::Real)
    signal = AnalyticSourceSignal(
        source_type,
        crest,
        time1,
        sfreq,
        tstart,
        tstop,
    )
    return CurrentInjection(
        node,
        signal,
    )
end

function constant_current_injection(node::Int, current_pu::Real)
    value = Float64(current_pu)
    return CurrentInjection(node, _t -> value)
end

end
