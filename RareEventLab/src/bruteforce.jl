export brute_force, BruteForceResult

struct BruteForceResult
    systems::Vector{String}
    events::Vector{ReactiveEvent}
    top_coarse::Matrix{Union{Bool, Missing}} # event × state definition
    time_since::NTuple{2, Float64}
    sim_time::Float64
    n_steps::Int
end

function brute_force(; n_transitions::Int,
                     systems = ("saddle hyperplane", "optimized", "basin", "core complement"))
    x = core_a.c
    return _brute_force_loop(x, Tuple(CoarseTrajectory(state_system(s), x) for s in systems),
                             n_transitions)
end

function _brute_force_loop(x::Point, coarse::Tuple, n_transitions::Int)
    observer = MilestoneObserver(x)
    rng = Random.default_rng()
    n_steps = 0
    while n_reactions(observer) < n_transitions
        x′ = langevin_step(x, rng)
        observe_step!(observer, x, x′)
        foreach(c -> x′ ∈ current_state(c) || next_state!(c, x′, clock(observer)), coarse)
        x = x′
        n_steps += 1
    end
    return BruteForceResult([c.states.name for c in coarse], observer.events,
                            stack([coarse_pathways(c, observer.events) for c in coarse]),
                            observer.time_since, clock(observer), n_steps)
end
