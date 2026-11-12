export run_experiment, leaderboard_records, clear_leaderboard

_leaderboard_path(method) = joinpath(pwd(), "leaderboard_$(method).jls")

function leaderboard_records(method::AbstractString)
    path = _leaderboard_path(method)
    isfile(path) || return Any[]
    try
        return deserialize(path)
    catch
        return Any[]
    end
end

clear_leaderboard(method::AbstractString) = rm(_leaderboard_path(method); force = true)

function _record!(method, setup, summary)
    recs = leaderboard_records(method)
    push!(recs, (setup = setup, summary...))
    serialize(_leaderboard_path(method), recs)
    return nothing
end

label(s::ParRepDesign) = "$(s.states) | $(s.coresets) | $(s.tcorr) | $(s.N)"
label(s::AMSDesign) = "$(s.cv) | $(s.zmin) | $(s.zmax) | $(s.N_loop) | $(s.L) | $(s.n_corr) | $(s.loop_init) | $(s.N) | $(s.m)"

function run_experiment(setup::ParRepDesign; n_transitions::Int)
    r = parallel_replica(setup; n_transitions)
    _record!("parrep", setup,
             (label = label(setup), trace_wall = r.trace_wall,
              trace_sim = round.(Int, r.trace_sim ./ Δt), k_AB = rate(r),
              top = top_fraction(r)))
    return r
end

function run_experiment(setup::AMSDesign; M::Int)
    r = adaptive_multilevel_splitting(setup; M)
    _record!("ams", setup,
             (label = label(setup), total = total_cost(r) / M, wall = wall_cost(r) / M,
              k_simulations = simulation_rates(r), k_AB = hill_rate(r), top = top_fraction(r)))
    return r
end
