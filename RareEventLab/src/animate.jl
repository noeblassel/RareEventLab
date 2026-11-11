export animate_algorithm

function animate_algorithm(setup::ParRepDesign; n_cycles::Int = 3, fps::Int = 15)
    sys = state_system(setup.states)
    N = setup.N
    rng = Random.default_rng()
    θs = range(0, 2π, 100)
    bgs = map(sys.states) do Ω
        bg = _potential_background("")
        contour!(bg, _XS, _YS, (a, b) -> Float64(inside(Ω, Point(a, b))); levels = [0.5],
                 c = :white, lw = 1.5, label = "")
        plot!(bg, [NaN], [NaN]; c = :white, lw = 1.5, label = "current state boundary")
        if setup.coresets == "physical"
            plot!(bg, Ω.core.c[1] .+ Ω.core.r .* cos.(θs), Ω.core.c[2] .+ Ω.core.r .* sin.(θs);
                  c = :yellow, lw = 1.5, ls = :dash, label = "current core set")
        end
        plot!(bg; legend = :bottomright, legendfontsize = 7, background_color_legend = :gray20,
              foreground_color_legend = nothing, legendfontcolor = :white)
    end
    x = core_a.c
    coarse = CoarseTrajectory(sys, x)
    wall, sim = 0, 0
    trail = [x]
    frames = Any[]
    push_frame!(title, pts, colors) = push!(frames, (title, copy(trail), copy(pts), copy(colors),
                                                     wall, sim, last(coarse.visits).state))
    for _ in 1:n_cycles
        k = 0
        while !in_core(setup, current_state(coarse), x)
            xn = langevin_step(x, rng)
            xn ∈ current_state(coarse) || next_state!(coarse, xn, sim * Δt)
            x = xn
            wall += 1; sim += 1; k += 1
            push!(trail, x)
            k % 25 == 0 && push_frame!("A) waiting phase (sequential)", [x], [:black])
        end
        Ω = current_state(coarse)
        n_corr = round(Int, decorrelation_time(setup, Ω) / Δt)
        reps = fill(x, N)
        alive = trues(N)
        trail = [x]
        aborted = false
        for n in 1:n_corr
            for i in 1:N
                alive[i] || continue
                reps[i] = langevin_step(reps[i], rng)
                alive[i] = inside(Ω, reps[i])
            end
            wall += 1; sim += 1
            push!(trail, reps[1])
            push_frame!("B) decorrelation: $(count(alive))/$N replicas alive, " *
                        "t = $(round(n * Δt; digits = 2)) / $(round(n_corr * Δt; digits = 2))",
                        reps, [a ? (i == 1 ? :black : :red) : :gray for (i, a) in enumerate(alive)])
            if !alive[1]
                push_frame!("B) reference walker exited: abort", reps, fill(:gray, N))
                next_state!(coarse, reps[1], sim * Δt)
                x = reps[1]
                aborted = true
                break
            end
        end
        aborted && continue
        surv = findall(alive)
        n = 0
        exited = 0
        while exited == 0
            n += 1
            for (j, i) in enumerate(surv)
                reps[i] = langevin_step(reps[i], rng)
                exited == 0 && !inside(Ω, reps[i]) && (exited = j)
            end
            wall += 1
            sim += length(surv)
            n % 25 == 0 && push_frame!("C) parallel phase: $(length(surv)) replicas", reps[surv],
                                       fill(:red, length(surv)))
        end
        sim -= length(surv) - exited
        x = reps[surv[exited]]
        trail = [x]
        for _ in 1:8
            push_frame!("C) replica $exited exits", [x], [:orange])
        end
        next_state!(coarse, x, sim * Δt)
    end
    anim = Animation()
    for (title, tr, pts, cols, w, s, cur) in frames
        p = deepcopy(bgs[cur])
        plot!(p, first.(tr), last.(tr); c = :black, lw = 0.6, alpha = 0.6, label = "")
        scatter!(p, first.(pts), last.(pts); c = cols, ms = 4, msw = 0, label = "",
                 title = title, titlefontsize = 10)
        annotate!(p, -2.4, 2.6, text("wall-time: $w   simulation time: $s   speedup: " *
                                     "$(round(s / max(w, 1); digits = 1))", 9, :left, :white))
        frame(anim, p)
    end
    return gif(anim, joinpath(tempdir(), "parrep.gif"); fps, show_msg = false)
end

function animate_algorithm(setup::AMSDesign; fps::Int = 5, every::Int = 10)
    ξ = collective_variable(setup)
    z_min, z_max = level_min(setup), level_max(setup)
    rng = Random.default_rng()
    N, m = setup.N, setup.m
    bg = _potential_background("")
    plot!(bg; size = (750, 650))
    contour!(bg, _XS, _YS, (a, b) -> ξ(Point(a, b)); levels = [z_min], c = :white, lw = 2)
    θs = range(0, 2π, 100)
    for (b, name) in ((core_a, "A"), (core_b, "B"))
        plot!(bg, b.c[1] .+ b.r .* cos.(θs), b.c[2] .+ b.r .* sin.(θs); c = :white, lw = 1.5,
              label = "")
        annotate!(bg, b.c[1], b.c[2] - b.r - 0.15, text(name, 10, :white))
    end
    anim = Animation()
    p_hat = 1.0
    iteration = 0
    show!(p, title, n_frames) = (plot!(p; title, titlefontsize = 10);
                                 foreach(_ -> frame(anim, p), 1:n_frames))
    status!(p) = annotate!(p, -2.4, 2.6, text("iteration $iteration   p̂ = " *
                                              "$(round(p_hat; sigdigits = 3))", 9, :left, :white))
    cols = cgrad(:hsv)[range(0, 1, length = N + 1)][1:N]
    α_running, α_reached_B, α_in_A, α_killed = 1.0, 0.8, 0.45, 0.05
    function draw(paths, alphas, Z)
        p = deepcopy(bg)
        isnothing(Z) || contour!(p, _XS, _YS, (a, b) -> ξ(Point(a, b)); levels = [Z],
                                 c = :white, ls = :dash, lw = 2)
        for i in sortperm(alphas)
            plot!(p, first.(paths[i]), last.(paths[i]); c = cols[i], alpha = alphas[i], lw = 1.5,
                  label = "")
        end
        status!(p)
        return p
    end
    settled_alpha(r) = r.reached_B ? α_reached_B : α_in_A

    starts = rand(first(loop_phase(setup; N_loop = 20)), N)

    paths = [[s] for s in starts]
    levels = ξ.(starts)
    running = trues(N)
    n_steps = 0
    while any(running)
        for i in findall(running)
            y = langevin_step(last(paths[i]), rng)
            push!(paths[i], y)
            levels[i] = max(levels[i], ξ(y))
            running[i] = !(y in core_b) && !(levels[i] >= z_min && y in core_a)
        end
        n_steps += 1
        if n_steps % every == 0 || !any(running)
            alphas = [running[i] ? α_running : last(paths[i]) in core_b ? α_reached_B : α_in_A
                      for i in 1:N]
            show!(draw(paths, alphas, nothing),
                  "B) initialization: $(count(running))/$N Y₀→Y₁ replicas running",
                  !any(running) ? 20 : 1)
        end
    end
    reps = [run_to_A_or_B(path, ξ, z_min, rng) for path in paths]

    while iteration < 200
        levels = [r.level for r in reps]
        Z = partialsort(levels, m)
        Z >= z_max && break
        killed = findall(≤(Z), levels)
        survivors = findall(>(Z), levels)
        isempty(survivors) && (p_hat = 0.0; break)
        alphas = [i in killed ? α_killed : settled_alpha(reps[i]) for i in 1:N]
        show!(draw([r.path for r in reps], alphas, Z),
              "C) $(length(killed)) replica killed at level Z=$(round(Z,sigdigits=2))", 10)
        p_hat *= 1 - length(killed) / N
        new_paths = [prefix_above(reps[rand(rng, survivors)], Z, ξ) for _ in killed]
        branch_points = last.(new_paths)
        running = trues(length(killed))
        n_steps = 0
        while any(running)
            for j in findall(running)
                y = langevin_step(last(new_paths[j]), rng)
                push!(new_paths[j], y)
                running[j] = !(y in core_a) && !(y in core_b)
            end
            n_steps += 1
            if n_steps % every == 0 || !any(running)
                paths = [r.path for r in reps]
                alphas = [settled_alpha(r) for r in reps]
                for (j, i) in enumerate(killed)
                    paths[i] = new_paths[j]
                    alphas[i] = running[j] ? α_running :
                                last(new_paths[j]) in core_b ? α_reached_B : α_in_A
                end
                p = draw(paths, alphas, Z)
                scatter!(p, first.(branch_points), last.(branch_points); c = cols[killed], ms = 6,
                         msw = 0, label = "")
                show!(p, "C) branching: " *
                         "$(count(running))/$(length(killed)) replicas running to Y₁",
                      !any(running) ? 10 : 1)
            end
        end
        for (i, path) in zip(killed, new_paths)
            reps[i] = run_to_A_or_B(path, ξ, z_min, rng)
        end
        iteration += 1
    end
    p_hat *= count(r -> r.reached_B, reps) / N
    show!(draw([r.path for r in reps], settled_alpha.(reps), nothing),
          "$(count(r -> r.reached_B, reps))/$N replicas reached B", 30)
    return gif(anim, joinpath(tempdir(), "ams.gif"); fps, show_msg = false)
end
