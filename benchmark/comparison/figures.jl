#=
    julia --project=benchmark benchmark/comparison/figures.jl

Render the comparison figures of the documentation from the committed results
in `results/`. `SQA_FIGURES_OUTPUT` selects another output directory, e.g. for
previews. `make_table.jl` reuses the kernel and scaling renderers on its run.
=#
module ComparisonFigures

    import JSON
    import CairoMakie as Makie
    using CairoMakie: rich, superscript, subscript

    const RESULTS = joinpath(@__DIR__, "results")
    const DOCS_ASSETS = normpath(joinpath(@__DIR__, "..", "..", "docs", "src", "assets"))

    # 624 units at CairoMakie's 0.75 pt per unit give a 6.5 in wide vector figure.
    const WIDTH = 624
    const INK = "#222222"
    const MUTED = "#6e6e6e"
    const RULE = "#e4e4e4"
    const SHADE = "#f0f0f0"

    const PACKAGES = [
        (key = "sqa", label = "SQA", color = "#0072B2", marker = :circle, size = 9),
        (key = "quantumalgebra", label = "QuantumAlgebra.jl", color = "#E69F00", marker = :rect, size = 8),
        (key = "sympy", label = "SymPy", color = "#CC79A7", marker = :utriangle, size = 10),
        (key = "openfermion", label = "OpenFermion", color = "#009E73", marker = :diamond, size = 10),
        (key = "sneg", label = "sneg", color = "#D55E00", marker = :cross, size = 10),
    ]
    const SQA = first(PACKAGES)
    const OTHERS = PACKAGES[2:end]
    const HIGH_ORDER_KEYS = Dict("SQA" => "sqa", "QuantumAlgebra" => "quantumalgebra")

    const THEME = Makie.Theme(;
        fontsize = 14, textcolor = INK, backgroundcolor = :white, figure_padding = (6, 14, 6, 8),
        Axis = (;
            backgroundcolor = :white, spinecolor = INK, spinewidth = 0.8,
            topspinevisible = false, rightspinevisible = false,
            xgridvisible = false, ygridcolor = RULE, ygridwidth = 0.8,
            xtickcolor = INK, ytickcolor = INK, xtickwidth = 0.8, ytickwidth = 0.8, xticksize = 4, yticksize = 4,
            xticklabelsize = 13, yticklabelsize = 13, xlabelsize = 14, ylabelsize = 14,
            xlabelpadding = 4, ylabelpadding = 6,
            titlesize = 14, titlefont = :bold, titlealign = :left, titlegap = 3,
            subtitlesize = 12, subtitlecolor = MUTED, subtitlegap = 1,
        ),
        Legend = (;
            framevisible = false, labelsize = 14, orientation = :horizontal, padding = (0, 0, 0, 0),
            colgap = 16, patchlabelgap = 5, patchsize = (22, 12), tellwidth = false,
        ),
    )

    const KERNEL_GROUPS = [
        "Jaynes–Cummings" => [
            "jc_build" => rich("build H"), "jc_H2" => rich("H", superscript("2")), "jc_heisenberg" => rich("[H, a]"),
        ],
        "Bose–Hubbard, M = 8" => ["bose_hubbard_build" => rich("build H"), "bose_hubbard_H2" => rich("H", superscript("2"))],
        "Tavis–Cummings, indexed" => [
            "tavis_cummings_build" => rich("build H"),
            "tavis_cummings_comm" => rich("[H, σ", superscript("+"), subscript("j"), "σ", superscript("−"), subscript("j"), "]"),
        ],
    ]

    const SCALING_PANELS = [
        (title = rich("Nested commutator"), xlabel = "depth d", points = [d => "jc_nested_d$d" for d in 1:7]),
        (title = rich("Power H", superscript("n")), xlabel = "power n", points = [n => "jc_pow_n$n" for n in 2:7]),
        (
            title = rich("Normal ordering (a a", superscript("†"), ")", superscript("n")),
            xlabel = "power n", points = [n => "fock_reorder_n$n" for n in 2:2:14],
        ),
        (
            title = rich("Bose–Hubbard chain"), xlabel = "sites M",
            points = [2 => "bh_chain_M2", 4 => "bh_chain_M4", 8 => "bose_hubbard_build", 16 => "bh_chain_M16", 32 => "bh_chain_M32"],
        ),
    ]

    const RATIO_ORDERS = (4, 12, 20)

    runtime(data::AbstractDict, package::String, key::String) =
        haskey(data, package) && haskey(data[package]["results"], key) ?
        Float64(data[package]["results"][key]["time_ns"]) : NaN

    function significant(x::Real, digits::Int)
        value = round(x; sigdigits = digits)
        text = isinteger(value) ? string(Int(value)) : string(value)
        return value >= 10_000 ? replace(text, r"(?<=\d)(?=(\d{3})+$)" => "\u2009") : text
    end

    function format_time(ns::Real, digits::Int, reference::Real = ns)
        for (scale, unit) in ((1.0e9, "s"), (1.0e6, "ms"), (1.0e3, "µs"))
            reference >= 0.9995scale && return significant(ns / scale, digits) * " " * unit
        end
        return significant(ns, digits) * " ns"
    end

    format_ratio(x::Real, digits::Int) = significant(x, digits) * "×"

    decades(low::Real, high::Real) = [10.0^k for k in ceil(Int, log10(low)):floor(Int, log10(high))]

    function log_limits(values::AbstractVector{<:Real}, below::Real, above::Real)
        finite = filter(isfinite, values)
        isempty(finite) && return (0.5, 2.0)
        return (min(minimum(finite), 1.0) / below, max(maximum(finite), 1.0) * above)
    end

    marker_element(style) = Makie.MarkerElement(; marker = style.marker, color = style.color, markersize = style.size + 1)

    line_element(style) = Makie.LineElement(; color = style.color, linewidth = style === SQA ? 3 : 1.6)

    function package_legend(position, entries)
        return Makie.Legend(position, [entry for (entry, _) in entries], [label for (_, label) in entries])
    end

    function save(fig::Makie.Figure, path::String)
        Makie.save(path, fig; pt_per_unit = 0.75)
        return path
    end

    """
        kernels_figure(data, path)

    Dot plot of each package's runtime relative to SQA on the fixed-size kernels, grouped by model.
    """
    function kernels_figure(data::AbstractDict, path::String)
        heads = String[]; head_ys = Float64[]
        rows = empty(last(first(KERNEL_GROUPS))); ys = Float64[]
        y = 0.0
        for (group, members) in KERNEL_GROUPS
            measured = filter(member -> isfinite(runtime(data, SQA.key, first(member))), members)
            isempty(measured) && continue
            y -= isempty(heads) ? 0.0 : 1.4
            push!(heads, group); push!(head_ys, y)
            for member in measured
                y -= 1.0
                push!(rows, member); push!(ys, y)
            end
        end
        keys = first.(rows); names = last.(rows)
        top = isempty(head_ys) ? 0.0 : first(head_ys)
        bottom = isempty(ys) ? -1.0 : last(ys)
        baselines = Float64[runtime(data, SQA.key, key) for key in keys]
        ratios = Dict(style.key => Float64[runtime(data, style.key, key) for key in keys] ./ baselines for style in OTHERS)
        low, high = log_limits(reduce(vcat, values(ratios); init = Float64[]), 1.0, 2.2)
        low = 10.0^floor(log10(low))
        offsets = Dict(zip((style.key for style in OTHERS), (0.22, 0.075, -0.075, -0.22)))
        present = [style for style in OTHERS if any(isfinite, ratios[style.key])]

        return Makie.with_theme(THEME) do
            fig = Makie.Figure(; size = (WIDTH, 360))
            entries = [
                Makie.LineElement(; color = SQA.color, linewidth = 3, points = Makie.Point2f[(0.5, -0.2), (0.5, 1.2)]) => "SQA";
                [marker_element(style) => style.label for style in present]
            ]
            package_legend(fig[1, 1:3], entries)
            labels = Makie.Axis(fig[2, 1]; width = 174)
            ax = Makie.Axis(
                fig[2, 2]; xscale = log10, xlabel = "runtime relative to SQA",
                xticks = (decades(low, high), [format_ratio(x, 2) for x in decades(low, high)]),
                leftspinevisible = false, yticksvisible = false, yticklabelsvisible = false, ygridvisible = false,
            )
            times = Makie.Axis(fig[2, 3]; width = 64)
            for side in (labels, times)
                Makie.hidedecorations!(side)
                Makie.hidespines!(side)
                Makie.xlims!(side, 0, 1)
            end
            Makie.linkyaxes!(labels, ax, times)
            Makie.xlims!(ax, low, high)
            Makie.ylims!(ax, bottom - 0.7, top + 0.7)
            Makie.colgap!(fig.layout, 1, 0)
            Makie.colgap!(fig.layout, 2, 0)
            Makie.rowgap!(fig.layout, 1, 8)

            Makie.vspan!(ax, low, 1.0; color = SHADE)
            Makie.hlines!(ax, ys; color = RULE, linewidth = 0.8)
            for side in (labels, ax, times)
                Makie.hlines!(side, head_ys[2:end] .+ 0.7; color = "#c8c8c8", linewidth = 0.8)
            end
            Makie.text!(
                ax, sqrt(low), top; text = "faster\nthan SQA", align = (:center, :center),
                justification = :center, fontsize = 12, color = MUTED, lineheight = 0.95,
            )
            Makie.text!(labels, zeros(length(heads)), head_ys; text = heads, font = :bold, fontsize = 13, align = (:left, :center))
            for (name, row) in zip(names, ys)
                Makie.text!(labels, 0.07, row; text = name, fontsize = 13, align = (:left, :center))
            end
            Makie.text!(times, 1.0, top; text = "SQA time", fontsize = 12, color = MUTED, align = (:right, :center))
            Makie.text!(times, ones(length(ys)), ys; text = [format_time(t, 3) for t in baselines], fontsize = 13, color = MUTED, align = (:right, :center))

            Makie.vlines!(ax, 1.0; color = SQA.color, linewidth = 2.5)
            for style in OTHERS
                Makie.scatter!(
                    ax, ratios[style.key], ys .+ offsets[style.key];
                    color = style.color, marker = style.marker, markersize = style.size,
                )
            end
            for (i, row) in enumerate(ys), style in OTHERS
                ratio = ratios[style.key][i]
                isfinite(ratio) && ratio < 1 || continue
                Makie.text!(
                    ax, 1.3, row; text = style.label * " " * format_ratio(1 / ratio, 2) * " faster",
                    fontsize = 12, align = (:left, :center),
                )
            end
            save(fig, path)
        end
    end

    """
        scaling_figure(data, path)

    Small multiples of each package's runtime relative to SQA against expression size.
    """
    function scaling_figure(data::AbstractDict, path::String)
        series = map(SCALING_PANELS) do panel
            xs = Float64.(first.(panel.points))
            baselines = [runtime(data, SQA.key, key) for key in last.(panel.points)]
            ratios = Dict(style.key => [runtime(data, style.key, key) for key in last.(panel.points)] ./ baselines for style in OTHERS)
            (; xs, baselines, ratios)
        end
        all_ratios = reduce(vcat, (reduce(vcat, values(s.ratios)) for s in series); init = Float64[])
        low, high = log_limits(all_ratios, 1.6, 2.5)
        present = [style for style in OTHERS if any(s -> any(isfinite, s.ratios[style.key]), series)]

        return Makie.with_theme(THEME) do
            fig = Makie.Figure(; size = (WIDTH, 480))
            entries = [
                [line_element(SQA), marker_element(SQA)] => SQA.label;
                [[line_element(style), marker_element(style)] => style.label for style in present]
            ]
            package_legend(fig[1, 2:3], entries)
            axes = Makie.Axis[]
            for (index, (panel, s)) in enumerate(zip(SCALING_PANELS, series))
                row, col = divrem(index - 1, 2) .+ (2, 2)
                measured = filter(isfinite, s.baselines)
                subtitle = isempty(measured) ? "" :
                    "SQA " * format_time(minimum(measured), 3) * " to " * format_time(maximum(measured), 3)
                ax = Makie.Axis(
                    fig[row, col]; title = panel.title, subtitle, xlabel = panel.xlabel,
                    xscale = panel === last(SCALING_PANELS) ? log2 : identity,
                    xticks = (s.xs, [string(Int(x)) for x in s.xs]), yscale = log10,
                    yticks = (decades(low, high), [format_ratio(y, 2) for y in decades(low, high)]),
                    yticklabelsvisible = col == 2,
                )
                push!(axes, ax)
                pad = panel === last(SCALING_PANELS) ? (first(s.xs) / 1.25, last(s.xs) * 1.25) :
                    (first(s.xs) - 0.06(last(s.xs) - first(s.xs)), last(s.xs) + 0.06(last(s.xs) - first(s.xs)))
                Makie.xlims!(ax, pad...)
                Makie.scatterlines!(
                    ax, s.xs, [isfinite(b) ? 1.0 : NaN for b in s.baselines]; color = SQA.color, linewidth = 3,
                    marker = SQA.marker, markersize = SQA.size,
                )
                for style in OTHERS
                    any(isfinite, s.ratios[style.key]) || continue
                    Makie.scatterlines!(
                        ax, s.xs, s.ratios[style.key]; color = style.color, linewidth = 1.6,
                        marker = style.marker, markersize = style.size,
                    )
                end
            end
            Makie.linkyaxes!(axes...)
            Makie.ylims!(first(axes), low, high)
            Makie.Label(fig[2:3, 1], "runtime relative to SQA"; rotation = π / 2, fontsize = 14)
            Makie.colgap!(fig.layout, 1, 4)
            Makie.colgap!(fig.layout, 2, 18)
            Makie.rowgap!(fig.layout, 1, 6)
            Makie.rowgap!(fig.layout, 2, 14)
            save(fig, path)
        end
    end

    """
        high_order_figure(high_order, path)

    Runtime of the verified Schrieffer–Wolff expansion against perturbative order, with hollow markers
    for points timed with fewer than three samples.
    """
    function high_order_figure(high_order::AbstractDict, path::String)
        rows = high_order["results"]
        series = map(PACKAGES[1:2]) do style
            selected = sort(filter(row -> HIGH_ORDER_KEYS[row["package"]] == style.key, rows); by = row -> row["order"])
            orders = [Int(row["order"]) for row in selected]
            times = [Float64(row["time_ns"]) for row in selected]
            sparse = [length(row["samples_ns"]) < 3 for row in selected]
            (; style, orders, times, sparse)
        end
        orders = sort(unique(reduce(vcat, (s.orders for s in series))))
        low, high = extrema(reduce(vcat, (s.times for s in series)))

        return Makie.with_theme(THEME) do
            fig = Makie.Figure(; size = (WIDTH, 340))
            entries = [[line_element(s.style), marker_element(s.style)] => s.style.label for s in series]
            package_legend(fig[1, 1], entries)
            ax = Makie.Axis(
                fig[2, 1]; xlabel = "perturbative order", ylabel = "runtime", yscale = log10,
                xticks = orders, yticks = (decades(low / 2, high * 2), [format_time(t, 3) for t in decades(low / 2, high * 2)]),
            )
            Makie.rowgap!(fig.layout, 1, 8)
            Makie.xlims!(ax, first(orders) - 0.8, last(orders) + 2.2)
            Makie.ylims!(ax, low / 2, high * 2)
            sqa, qa = series
            for order in RATIO_ORDERS
                order in sqa.orders && order in qa.orders || continue
                i = findfirst(==(order), sqa.orders)
                j = findfirst(==(order), qa.orders)
                Makie.linesegments!(
                    ax, [Makie.Point2(order, sqa.times[i] * 1.5), Makie.Point2(order, qa.times[j] / 1.5)];
                    color = MUTED, linewidth = 1,
                )
                Makie.text!(
                    ax, order - 0.25, sqrt(sqa.times[i] * qa.times[j]); text = format_ratio(qa.times[j] / sqa.times[i], 2),
                    fontsize = 12, align = (:right, :center),
                )
            end
            for s in reverse(series)
                linewidth = s.style === SQA ? 3 : 1.6
                Makie.lines!(ax, s.orders, s.times; color = s.style.color, linewidth)
                Makie.scatter!(
                    ax, s.orders, s.times; marker = s.style.marker, markersize = s.style.size + 1,
                    color = [sparse ? Makie.RGBf(1, 1, 1) : Makie.to_color(s.style.color) for sparse in s.sparse], strokecolor = s.style.color, strokewidth = 1.5,
                )
                Makie.text!(
                    ax, last(s.orders) + 0.5, last(s.times); text = format_time(last(s.times), 2, maximum(last, (t.times for t in series))),
                    fontsize = 12, color = s.style.color, align = (:left, :center),
                )
            end
            save(fig, path)
        end
    end

    function main()
        output = get(ENV, "SQA_FIGURES_OUTPUT", DOCS_ASSETS)
        mkpath(output)
        run = JSON.parsefile(joinpath(RESULTS, "run.json"))
        data = Dict{String, Any}(name => JSON.parsefile(joinpath(RESULTS, "$name.json")) for name in run["packages"])
        high_order = JSON.parsefile(joinpath(RESULTS, "high_order.json"))
        println(kernels_figure(data, joinpath(output, "comparison_kernels.svg")))
        println(scaling_figure(data, joinpath(output, "comparison_scaling.svg")))
        println(high_order_figure(high_order, joinpath(output, "comparison_high_order.svg")))
        return nothing
    end

end

if abspath(PROGRAM_FILE) == @__FILE__
    ComparisonFigures.main()
end
