export Reference, reference, rate, entry_angles, top_fraction, hill_rate, simulation_rates

struct Reference
    angles_A::Vector{Float64}
    angles_B::Vector{Float64}
    top_geom::Vector{Bool}
    p_top::Float64
    p_top_se::Float64
    k_AB::Float64
    k_AB_se::Float64
    n_events::Int
end

const _REFERENCE = Ref{Reference}()

function _load_reference!(path)
    isfile(path) || return
    time_A = parse(Float64, match(r"time_since_A=(\S+)", readline(path))[1])
    data = first(readdlm(path, ','; skipstart = 1, header = true))
    events = [ReactiveEvent(row[2], row[3], row[4], row[5] == true) for row in eachrow(data)]
    chains = data[:, 1]
    τ = reduce(vcat, [_A_to_B_durations(events[chains .== c]) for c in unique(chains)])
    tg = [e.top_geom for e in events]
    p_top = sum(tg) / length(tg)
    _REFERENCE[] = Reference([e.angle for e in events if e.to == -1],
                             [e.angle for e in events if e.to == 1], tg, p_top,
                             sqrt(p_top * (1 - p_top) / length(tg)),
                             count(e -> e.to == 1, events) / time_A, last(_rate_se(τ)),
                             length(events))
end

reference() = _REFERENCE[]

function _rate_se(τ)
    n = length(τ)
    m = sum(τ) / n
    s² = sum(abs2, τ .- m) / (n - 1)
    return (1 / m, sqrt(s² / n) / m^2)
end

function _A_to_B_durations(events)
    τ = Float64[]
    t_A = 0.0
    for e in events
        e.to == 1 ? push!(τ, e.time - t_A) : (t_A = e.time)
    end
    return τ
end

rate(r::Union{ParRepResult, BruteForceResult}) = _rate_se(_A_to_B_durations(r.events))

function entry_angles(r::Union{ParRepResult, BruteForceResult}, to::Int)
    θ = [e.angle for e in r.events if e.to == to]
    return θ, fill(1.0, length(θ))
end

function entry_angles(res::AMSResult, to::Int)
    θ, w = Float64[], Float64[]
    for simulation in res.simulations
        isempty(simulation.angles) && continue
        append!(θ, to == 1 ? simulation.angles :
                   [rem2pi(π - a, RoundNearest) for a in simulation.angles])
        append!(w, fill(simulation.p / length(simulation.angles),
                        length(simulation.angles)))
    end
    return θ, w
end

function top_fraction(r::ParRepResult)
    v = collect(skipmissing(r.top_coarse))
    isempty(v) && return (missing, missing)
    f = sum(v) / length(v)
    return (f, sqrt(f * (1 - f) / length(v)))
end

function top_fraction(res::AMSResult)
    ps = [simulation.p for simulation in res.simulations]
    fs = [isempty(simulation.top_geom) ? 0.0 :
          sum(simulation.top_geom) / length(simulation.top_geom)
          for simulation in res.simulations]
    f = sum(ps .* fs) / sum(ps)
    M = length(ps)
    pm = sum(ps) / M
    se = sqrt(sum(abs2, ps .* (fs .- f)) / (M * (M - 1))) / pm
    return (f, se)
end

function hill_rate(res::AMSResult)
    M = length(res.simulations)
    X = [[simulation.p, simulation.T_loop,
          simulation.p > 0 ? simulation.p * simulation.T_reac : 0.0]
         for simulation in res.simulations]
    m = sum(X) / M
    p, T_loop, p_T_react = m
    Σ = sum((x - m) * (x - m)' for x in X) / (M - 1)
    D = (1 - p) * T_loop + p_T_react # Hill rate denominator
    ∇k = [T_loop + c, -p * (1 - p), -p] / D^2
    return (p / D, sqrt(∇k' * Σ * ∇k / M)) # standard error using delta method
end

simulation_rates(res::AMSResult) = [simulation.p == 0 ? 0.0 :
                                    1 / ((1 - simulation.p) / simulation.p * simulation.T_loop +
                                         simulation.T_reac)
                                    for simulation in res.simulations]

total_cost(res::AMSResult) = sum(r.cost_loops.total + r.cost_splitting.total
                                 for r in res.simulations)
wall_cost(res::AMSResult) = sum(r.cost_loops.wall + r.cost_splitting.wall
                                for r in res.simulations)