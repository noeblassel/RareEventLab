export parallel_replica

function parallel_replica(setup::ParRepDesign; n_transitions::Int)
    X = core_a.c
    coarse = CoarseTrajectory(state_system(setup.states), X)
    Ω = current_state(coarse)
    observer = MilestoneObserver(X)
    x, previous, alive = fill(X, setup.N), fill(X, setup.N), zeros(Bool, setup.N)
    stats = ParRepStats()

    while n_reactions(observer) < n_transitions
        while !in_core(setup, Ω, X)
            X′ = langevin_step(X, Random.default_rng())
            charge!(stats, 1)
            observe_step!(observer, X, X′)
            if X′ ∉ Ω
                Ω = next_state!(coarse, X′, clock(observer))
            end
            X = X′
        end

        fill!(x, X)
        fill!(alive, true)
        for _ in 1:n_steps(decorrelation_time(setup, Ω))
            parallel_step!(x, previous, alive, Ω, stats)
            observe_step!(observer, X, x[1])
            X = x[1]
            X ∈ Ω || break
        end
        if X ∉ Ω
            Ω = next_state!(coarse, X, clock(observer))
            stats.n_aborts += 1
            trace!(stats, observer)
            continue
        end

        survivors = findall(alive)
        N_eff = length(survivors)
        n = 0
        while count(alive) == N_eff
            parallel_step!(x, previous, alive, Ω, stats)
            n += 1
        end
        j = findfirst(i -> !alive[i], survivors)
        T_par = (n - 1 + (j - 1) / N_eff) * Δt
        advance_clock!(observer, N_eff * T_par)
        observe_step!(observer, previous[survivors[j]], x[survivors[j]])
        X = x[survivors[j]]
        Ω = next_state!(coarse, X, clock(observer))
        push!(stats.n_eff, N_eff)
        trace!(stats, observer)
    end

    return ParRepResult(setup, observer.events, coarse_pathways(coarse, observer.events),
                        clock(observer), stats.cost, stats.trace_wall, stats.trace_sim,
                        stats.n_eff, stats.n_aborts)
end

n_steps(t) = round(Int, t / Δt)

mutable struct ParRepStats
    cost::StepCounter
    trace_wall::Vector{Int}
    trace_sim::Vector{Float64}
    n_eff::Vector{Int}
    n_aborts::Int
end
ParRepStats() = ParRepStats(StepCounter(), Int[], Float64[], Int[], 0)

charge!(s::ParRepStats, n_force_calls) = (s.cost.wall += 1; s.cost.total += n_force_calls; s)
trace!(s::ParRepStats, obs) = (push!(s.trace_wall, s.cost.wall); push!(s.trace_sim, clock(obs)); s)

function parallel_step!(x, previous, alive, Ω::MetastableState, s::ParRepStats)
    charge!(s, count(alive))
    @threads for i in eachindex(x)
        if alive[i]
            previous[i] = x[i]
            x[i] = langevin_step(x[i], Random.default_rng())
            alive[i] = x[i] ∈ Ω
        end
    end
    return x
end
