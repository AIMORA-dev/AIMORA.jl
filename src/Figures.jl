module Figures

using Printf

export write_inverter_power_svg,
       write_inverter_current_svg,
       write_hybrid_summary_svg,
       write_unified_voltage_svg,
       write_unified_power_svg

function ensure_dir(path::AbstractString)
    isdir(path) || mkpath(path)
end

function svg_escape(text)
    s = string(text)
    s = replace(s, "&" => "&amp;")
    s = replace(s, "<" => "&lt;")
    s = replace(s, ">" => "&gt;")
    return s
end

function points(xs, ys, xscale, yscale)
    return join((@sprintf("%.2f,%.2f", xscale(x), yscale(y)) for (x, y) in zip(xs, ys)), " ")
end

function nice_limits(values; pad = 0.05)
    low = minimum(values)
    high = maximum(values)
    span = max(high - low, 1.0e-9)
    return low - pad * span, high + pad * span
end

function write_line_svg(path::AbstractString; title, ylabel, series, vline_ms = nothing, ylim = nothing)
    ensure_dir(dirname(path))
    width, height = 1200.0, 680.0
    left, right, top, bottom = 84.0, 34.0, 76.0, 86.0
    plot_w = width - left - right
    plot_h = height - top - bottom

    xs = reduce(vcat, [s.x for s in series])
    ys = reduce(vcat, [s.y for s in series])
    xmin, xmax = minimum(xs), maximum(xs)
    ymin, ymax = ylim === nothing ? nice_limits(ys) : ylim
    xscale(x) = left + (x - xmin) / max(xmax - xmin, 1.0e-9) * plot_w
    yscale(y) = top + (ymax - y) / max(ymax - ymin, 1.0e-9) * plot_h

    open(path, "w") do io
        println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$(Int(width))" height="$(Int(height))" viewBox="0 0 $(Int(width)) $(Int(height))">""")
        println(io, """<rect width="100%" height="100%" fill="#f7f7f4"/>""")
        println(io, """<rect x="$left" y="$top" width="$plot_w" height="$plot_h" fill="#ffffff" stroke="#c8c8c0"/>""")

        for i in 0:5
            y = top + i * plot_h / 5
            value = ymax - i * (ymax - ymin) / 5
            println(io, """<line x1="$left" y1="$y" x2="$(left + plot_w)" y2="$y" stroke="#ddddd5" stroke-width="1"/>""")
            @printf(io, "<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"end\" font-family=\"Arial\" font-size=\"13\" fill=\"#444\">%.2f</text>\n", left - 10, y + 4, value)
        end
        for i in 0:6
            x = left + i * plot_w / 6
            value = xmin + i * (xmax - xmin) / 6
            println(io, """<line x1="$x" y1="$top" x2="$x" y2="$(top + plot_h)" stroke="#eeeeea" stroke-width="1"/>""")
            @printf(io, "<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"middle\" font-family=\"Arial\" font-size=\"13\" fill=\"#444\">%.0f</text>\n", x, top + plot_h + 28, value)
        end

        if vline_ms !== nothing
            x = xscale(vline_ms)
            println(io, """<line x1="$x" y1="$top" x2="$x" y2="$(top + plot_h)" stroke="#333" stroke-width="1.5" stroke-dasharray="5 5"/>""")
        end

        legend_x = left + plot_w - 210
        legend_y = top + 26
        for (i, s) in enumerate(series)
            y = legend_y + (i - 1) * 24
            dash = get(s, :dash, false) ? " stroke-dasharray=\"8 5\"" : ""
            println(io, """<line x1="$legend_x" y1="$y" x2="$(legend_x + 34)" y2="$y" stroke="$(s.color)" stroke-width="3"$dash/>""")
            println(io, """<text x="$(legend_x + 44)" y="$(y + 5)" font-family="Arial" font-size="14" fill="#222">$(svg_escape(s.label))</text>""")
        end

        for s in series
            dash = get(s, :dash, false) ? " stroke-dasharray=\"8 5\"" : ""
            println(io, """<polyline fill="none" stroke="$(s.color)" stroke-width="2.5"$dash points="$(points(s.x, s.y, xscale, yscale))"/>""")
        end

        println(io, """<text x="$(width / 2)" y="36" text-anchor="middle" font-family="Arial" font-size="22" font-weight="700" fill="#1f1f1d">$(svg_escape(title))</text>""")
        println(io, """<text x="$(width / 2)" y="$(height - 24)" text-anchor="middle" font-family="Arial" font-size="16" fill="#333">Time (ms)</text>""")
        println(io, """<text transform="translate(24 $(height / 2)) rotate(-90)" text-anchor="middle" font-family="Arial" font-size="16" fill="#333">$(svg_escape(ylabel))</text>""")
        println(io, "</svg>")
    end
end

function write_inverter_power_svg(path::AbstractString, rows)
    t = [r[1] * 1000.0 for r in rows]
    write_line_svg(
        path;
        title = "Julia Inverter Power Response",
        ylabel = "Power (pu)",
        vline_ms = 50.0,
        ylim = (-0.02, 0.88),
        series = [
            (label = "P reference", x = t, y = [r[8] for r in rows], color = "#6f6f69", dash = true),
            (label = "P output", x = t, y = [r[6] for r in rows], color = "#176b87", dash = false),
            (label = "Q reference", x = t, y = [r[9] for r in rows], color = "#8a817c", dash = true),
            (label = "Q output", x = t, y = [r[7] for r in rows], color = "#9a4d9e", dash = false),
        ],
    )
end

function write_inverter_current_svg(path::AbstractString, rows)
    t = [r[1] * 1000.0 for r in rows]
    write_line_svg(
        path;
        title = "Julia Inverter dq Current Tracking",
        ylabel = "Current (pu)",
        vline_ms = 50.0,
        ylim = (-0.12, 0.88),
        series = [
            (label = "id reference", x = t, y = [r[4] for r in rows], color = "#6f6f69", dash = true),
            (label = "id", x = t, y = [r[2] for r in rows], color = "#176b87", dash = false),
            (label = "iq reference", x = t, y = [r[5] for r in rows], color = "#8a817c", dash = true),
            (label = "iq", x = t, y = [r[3] for r in rows], color = "#9a4d9e", dash = false),
        ],
    )
end

function write_hybrid_summary_svg(path::AbstractString; network_engine, inverter_engine, inverter_bus, coupling_voltage_pu, final_p_pu, final_q_pu, samples, fortran_process_ok, report_stop_detected)
    ensure_dir(dirname(path))
    open(path, "w") do io
        println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="680" viewBox="0 0 1200 680">""")
        println(io, """<rect width="100%" height="100%" fill="#f7f7f4"/>""")
        println(io, """<rect x="70" y="78" width="610" height="500" fill="#fff" stroke="#c8c8c0"/>""")
        println(io, """<rect x="730" y="78" width="380" height="500" fill="#fff" stroke="#c8c8c0"/>""")
        println(io, """<text x="600" y="38" text-anchor="middle" font-family="Arial" font-size="22" font-weight="700" fill="#1f1f1d">Hybrid Run Summary</text>""")

        lines = [
            ("Network engine", network_engine),
            ("Inverter engine", inverter_engine),
            ("Inverter bus", inverter_bus),
            ("Coupling voltage", @sprintf("%.4f pu", coupling_voltage_pu)),
            ("Fortran process OK", fortran_process_ok ? "yes" : "no"),
            ("Report stop flag", report_stop_detected ? "yes" : "no"),
            ("Samples", string(samples)),
            ("Final P/Q", @sprintf("%.4f / %.4f pu", final_p_pu, final_q_pu)),
        ]
        y = 130
        for (label, value) in lines
            println(io, """<text x="105" y="$y" font-family="Arial" font-size="17" fill="#5a5a54">$(svg_escape(label))</text>""")
            println(io, """<text x="340" y="$y" font-family="Arial" font-size="17" fill="#1f1f1d">$(svg_escape(value))</text>""")
            y += 48
        end

        bars = [("Final P", final_p_pu, "#176b87"), ("Final Q", final_q_pu, "#9a4d9e"), ("Bus V", coupling_voltage_pu, "#5f7f52")]
        for (i, (label, value, color)) in enumerate(bars)
            x = 785 + (i - 1) * 105
            h = 390 * value / 1.08
            y0 = 520 - h
            println(io, """<rect x="$x" y="$y0" width="58" height="$h" fill="$color"/>""")
            println(io, """<text x="$(x + 29)" y="$(y0 - 12)" text-anchor="middle" font-family="Arial" font-size="15" fill="#222">$(@sprintf("%.3f", value))</text>""")
            println(io, """<text x="$(x + 29)" y="552" text-anchor="middle" font-family="Arial" font-size="15" fill="#333">$(svg_escape(label))</text>""")
        end
        println(io, "</svg>")
    end
end

function write_unified_voltage_svg(path::AbstractString, rows)
    t = [r[1] * 1000.0 for r in rows]
    write_line_svg(
        path;
        title = "Reduced Feeder EMT Bus Voltage",
        ylabel = "Voltage (pu)",
        vline_ms = 50.0,
        ylim = (0.94, 1.08),
        series = [
            (label = "Load/inverter bus", x = t, y = [r[2] for r in rows], color = "#176b87", dash = false),
            (label = "Feeder source node", x = t, y = [r[3] for r in rows], color = "#8a817c", dash = true),
            (label = "Source", x = t, y = [r[4] for r in rows], color = "#5f7f52", dash = true),
        ],
    )
end

function write_unified_power_svg(path::AbstractString, rows)
    t = [r[1] * 1000.0 for r in rows]
    write_line_svg(
        path;
        title = "Reduced Feeder Inverter Power",
        ylabel = "Power (pu)",
        vline_ms = 50.0,
        ylim = (-0.02, 0.90),
        series = [
            (label = "P reference", x = t, y = [r[12] for r in rows], color = "#6f6f69", dash = true),
            (label = "P output", x = t, y = [r[6] for r in rows], color = "#176b87", dash = false),
            (label = "Q reference", x = t, y = [r[13] for r in rows], color = "#8a817c", dash = true),
            (label = "Q output", x = t, y = [r[7] for r in rows], color = "#9a4d9e", dash = false),
        ],
    )
end

end
