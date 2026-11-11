export adaptive_multilevel_splitting, ams_simulation, loop_phase, run_to_A_or_B, prefix_above

adaptive_multilevel_splitting(setup::AMSDesign; M::Int) =
    AMSResult(setup, [ams_simulation(setup) for _ in 1:M])

function ams_simulation(setup::AMSDesign)
    ξ = collective_variable(setup)
    z_min, z_max = level_min(setup), level_max(setup)
    N, m = setup.N, setup.m

    ν_Q_samples, T_loop, cost_loops = loop_phase(setup)

    cost = StepCounter()
    replicas = run_all_to_A_or_B!([[Y₀] for Y₀ in rand(ν_Q_samples, N)], ξ, z_min, cost)

    p̂ = 1.0
    k = 1
    while true
        Z_k = partialsort([r.level for r in replicas], m)
        killed = findall(r -> r.level <= Z_k, replicas)
        M_k = length(killed)
        (Z_k >= z_max || M_k == N) && break

        p̂ *= 1 - M_k / N
        survivors = findall(r -> r.level > Z_k, replicas)
        prefixes = [prefix_above(replicas[rand(survivors)], Z_k, ξ) for _ in killed]
        replicas[killed] = run_all_to_A_or_B!(prefixes, ξ, z_min, cost)
        k += 1
    end

    reactive = filter(r -> r.reached_B, replicas)
    N_B = length(reactive)
    p̂ *= N_B / N
    return AMSSimulation(p̂, T_loop, reactive, k - 1, cost_loops, cost)
end

struct Replica
    path::Vector{Point{Float64}}
    level::Float64 # max of ξ along the path, Inf if it reached B
    reached_B::Bool
    n_new::Int
end

function run_to_A_or_B(path, ξ, z_min, rng)
    x = last(path)
    level = maximum(ξ, path)
    n = 0
    while !(x in core_b) && !(level >= z_min && x in core_a)
        x = langevin_step(x, rng)
        push!(path, x)
        level = max(level, ξ(x))
        n += 1
    end
    return Replica(path, x in core_b ? Inf : level, x in core_b, n)
end

function run_all_to_A_or_B!(prefixes, ξ, z_min, cost::StepCounter)
    replicas = Vector{Replica}(undef, length(prefixes))
    @threads for i in eachindex(prefixes)
        replicas[i] = run_to_A_or_B(prefixes[i], ξ, z_min, Random.default_rng())
    end
    cost.wall += maximum(r -> r.n_new, replicas)
    cost.total += sum(r -> r.n_new, replicas)
    return replicas
end

prefix_above(r::Replica, Z, ξ) = r.path[1:findfirst(x -> ξ(x) > Z, r.path)]

function AMSSimulation(p̂, T_loop, reactive::Vector{Replica}, n_iterations, cost_loops, cost)
    reactions = [only(observe_path(r.path).events) for r in reactive]
    durations = [(length(r.path) - 1) * Δt for r in reactive]
    T_reac = isempty(reactive) ? NaN : sum(durations) / length(durations)
    return AMSSimulation(p̂, T_loop, T_reac, [e.angle for e in reactions],
                         [e.top_geom for e in reactions], n_iterations, cost_loops, cost)
end

function observe_path(path)
    obs = MilestoneObserver(-1, NaN, 0.0, (0.0, 0.0), ReactiveEvent[])
    foreach(t -> observe_step!(obs, path[t - 1], path[t]), 2:length(path))
    return obs
end

function _gibbs_in_A(rng)
    Vmin = V(core_a.c)
    while true
        x = core_a.c + core_a.r * Point(2rand(rng) - 1, 2rand(rng) - 1)
        x in core_a && rand(rng) < exp(-β * (V(x) - Vmin)) && return x
    end
end

function loop_phase(setup::AMSDesign; N_loop::Int = setup.N_loop, max_attempt = 10000)
    ξ, z_min = collective_variable(setup), level_min(setup)
    L, n_corr = setup.L, setup.n_corr
    paths = Vector{Tuple{Vector{Point{Float64}}, Float64, Int}}(undef, N_loop)
    @threads for i in 1:N_loop
        rng = Random.default_rng()
        steps = 0
        for attempt in 1:max_attempt
            Y = [setup.loop_init == "dirac" ? dirac_point(setup) : _gibbs_in_A(rng)]
            Δ = Float64[]
            while length(Y) <= L
                step = run_to_A_or_B([last(Y)], ξ, z_min, rng)
                steps += step.n_new
                step.reached_B && break
                push!(Y, last(step.path))
                push!(Δ, step.n_new * Δt)
            end
            length(Y) == L + 1 || continue

            paths[i] = (Y[(n_corr + 1):end], sum(Δ[(n_corr + 1):L]) / (L - n_corr), steps)
            break
        end
        isassigned(paths, i) ||
            error("$max_attempt loop paths in a row reached B before L = $L loops (Σ is likely too far from A)")
    end
    ν_Q_samples = reduce(vcat, [path[1] for path in paths])
    T_loop = sum(path[2] for path in paths) / N_loop
    steps = [path[3] for path in paths]
    return ν_Q_samples, T_loop, StepCounter(maximum(steps), sum(steps))
end
