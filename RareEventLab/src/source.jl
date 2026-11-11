export show_source, sample_trajectory

struct _Source
    md::String
end
Base.show(io::IO, ::MIME"text/markdown", s::_Source) = print(io, s.md)
Base.show(io::IO, s::_Source) = print(io, s.md)

function show_source(f)
    file, line = functionloc(first(methods(f)))
    ls = readlines(file)
    return _Source("```julia\n" * join(ls[line:findnext(==("end"), ls, line)], "\n") * "\n```")
end

function sample_trajectory(T; x0::Point = core_a.c)
    rng = Random.default_rng()
    path = [x0]
    for _ in 1:n_steps(T)
        push!(path, langevin_step(last(path), rng))
    end
    return path
end
