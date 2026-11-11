export leaderboard, entry_polar, channels, show_potential

const _NBINS = 24

function _histogram(θ, w; nbins = _NBINS)
    edges = range(-π, π, nbins + 1)
    h = zeros(nbins)
    for (a, wa) in zip(θ, w)
        h[clamp(searchsortedlast(edges, a), 1, nbins)] += wa
    end
    W = sum(w)
    n_eff = W^2 / sum(abs2, w) # kish ESS
    p = h ./ W
    width = step(edges)
    se = sqrt.(p .* (1 .- p) ./ n_eff) ./ width
    return collect(edges), p ./ width, se
end

_closed(v) = [v; v[1]]

function entry_polar(result)
    ref = reference()
    ϕ = range(0, 2π, 200)
    panels = map(((-1, "entry into A", ref.angles_A), (1, "entry into B", ref.angles_B))) do (to, name, ref_θ)
        θ, w = entry_angles(result, to)
        edges, d_ref, _ = _histogram(ref_θ, fill(1.0, length(ref_θ)))
        mids = _closed((edges[1:(end - 1)] .+ edges[2:end]) ./ 2)
        scale = 1 / maximum(d_ref)
        radial(d) = (r = 1 .+ scale .* _closed(d); (r .* cos.(mids), r .* sin.(mids)))
        p = plot(cos.(ϕ), sin.(ϕ); c = :gray, lw = 1.5, label = "core-set boundary",
                 aspect_ratio = 1, framestyle = :none, title = name, legend = :outerbottom)
        if !isempty(θ)
            _, d, se = _histogram(θ, w)
            if result isa AMSResult
                # paths within one AMS simulation are correlated
                # Each AMS simulation m gives one sample, the fraction h[m][b] of its reactive paths in bin b
                # These samples are combined with weights p̂ₘ
                width = step(range(-π, π, _NBINS + 1))
                ps = [s.p for s in result.simulations]
                hs = map(result.simulations) do s
                    a = to == 1 ? s.angles : [rem2pi(π - x, RoundNearest) for x in s.angles]
                    isempty(a) ? zeros(_NBINS) : _histogram(a, fill(1.0, length(a)))[2] .* width
                end
                M = length(ps)
                se = [sqrt(sum(abs2, ps .* (getindex.(hs, b) .- d[b] * width)) / (M * (M - 1))) /
                      (sum(ps) / M) / width for b in 1:_NBINS]
            end
            if !(result isa AMSResult && length(result.simulations) < 2)
                lo, hi = radial(max.(d .- 2se, 0)), radial(d .+ 2se)
                plot!(p, Shape([hi[1]; reverse(lo[1])], [hi[2]; reverse(lo[2])]); c = :red,
                      alpha = 0.25, lw = 0, label = "± 2 s.e.")
            end
            name = result isa BruteForceResult ? "this simulation" : label(result.setup)
            plot!(p, radial(d)...; lw = 2, c = :red, label = name)
        end
        plot!(p, radial(d_ref)...; lw = 2, c = :black, label = "reference")
    end
    return plot(panels...; layout = (1, 2), size = (900, 480))
end

function leaderboard(method::AbstractString)
    recs = leaderboard_records(method)
    isempty(recs) && return plot(; framestyle = :none,
                                 annotations = (0.5, 0.5, text("no experiments yet", 14)))
    method == "parrep" && return _leaderboard_parrep(recs)
    method == "ams" && return _leaderboard_ams(recs)
    throw(ArgumentError("unknown method \"$method\"; use \"parrep\" or \"ams\""))
end

function _leaderboard_parrep(recs)
    xmax = maximum(r -> last(r.trace_wall), recs)
    p = plot([0, xmax], [0, xmax]; c = :black, ls = :dash, label = "sequential",
             xlabel = "wall-time (force evaluations)", ylabel = "simulation time (steps)",
             legend = :outerright, size = (1100, 550), left_margin = 8Plots.mm,
             bottom_margin = 6Plots.mm)
    for (i, r) in enumerate(recs)
        name = "$i. $(label(r.setup))"
        plot!(p, [0; r.trace_wall], [0; r.trace_sim]; lw = 2, label = name)
    end
    return p
end

function _leaderboard_ams(recs)
    relvar(r) = (k = r.k_simulations; m = sum(k) / length(k); sum(abs2, k .- m) / (length(k) - 1) / m^2)
    p = plot(; xscale = :log10, yscale = :log10,
             xlabel = "total cost per simulation (force evaluations)",
             ylabel = "relative variance of k̂ per simulation", legend = :outerright, size = (1100, 550),
             left_margin = 8Plots.mm, bottom_margin = 6Plots.mm)
    costs = [r.total for r in recs]
    cmin, cmax = minimum(costs) / 3, maximum(costs) * 3
    for (i, r) in enumerate(recs)
        v = relvar(r)
        # relative variance × cost is constant along this line
        plot!(p, [cmin, cmax], v * r.total ./ [cmin, cmax]; c = i, ls = :dash, lw = 1, label = "")
        scatter!(p, [r.total], [v]; m = :circle, ms = 7, c = i, label = "$i. $(label(r.setup))")
    end
    plot!(p; xlims = (cmin, cmax))
    return p
end

function show_potential()
    xs, ys = range(-2.5, 2.5, 200), range(-1.6, 2.8, 200)
    p = contourf(xs, ys, (a, b) -> min(V(Point(a, b)), 0.0); levels = 30, c = :viridis,
                 aspect_ratio = 1, xlims = (-2.5, 2.5), ylims = (-1.6, 2.8), size = (700, 620))
    circle(b) = (b.c[1] .+ b.r .* cos.(range(0, 2π, 100)), b.c[2] .+ b.r .* sin.(range(0, 2π, 100)))
    plot!(p, circle(core_a)...; c = :white, lw = 2, label = "A")
    plot!(p, circle(core_b)...; c = :white, lw = 2, ls = :dash, label = "B")
    plot!(p, circle(core_top)...; c = :gray, lw = 1, label = "top core-set (unphysical)")
    plot!(p, [0, 0], [y_threshold, 2.8]; c = :orange, lw = 2, label = "top pathway")
    plot!(p, [0, 0], [-1.6, y_threshold]; c = :cyan, lw = 2, label = "bottom pathway")
    return p
end

export show_states, show_levels

const _XS = range(-2.5, 2.5, 300)
const _YS = range(-1.6, 2.8, 300)

_potential_background(title) =
    contourf(_XS, _YS, (a, b) -> min(V(Point(a, b)), 0.0); levels = 30, c = :viridis,
             aspect_ratio = 1, xlims = extrema(_XS), ylims = extrema(_YS), colorbar = false,
             title = title)

function show_states()
    colors = (:red, :cyan, :orange)
    panels = map(("saddle hyperplane", "optimized", "basin", "core complement")) do name
        sys = state_system(name)
        p = _potential_background(name)
        for (i, Ω) in enumerate(sys.states)
            contour!(p, _XS, _YS, (a, b) -> Float64(inside(Ω, Point(a, b))); levels = [0.5],
                     c = colors[i], lw = 2.5)
        end
        p
    end
    return plot(panels...; layout = (2, 2), size = (1000, 1000))
end

function show_levels()
    styles = (("low", :solid), ("medium", :dash), ("high", :dot))
    panels = map(("x", "twisted", "committor")) do cv
        p = _potential_background("ξ = $cv")
        ξ = collective_variable(AMSDesign(cv = cv))
        for (z, ls) in styles
            zmin = level_min(AMSDesign(cv = cv, zmin = z))
            contour!(p, _XS, _YS, (a, b) -> ξ(Point(a, b)); levels = [zmin], c = :white,
                     lw = 2, ls = ls)
        end
        p
    end
    return plot(panels...; layout = (1, 3), size = (1500, 520),
                plot_title = "Σ: low (solid), medium (dashed), high (dotted)")
end
