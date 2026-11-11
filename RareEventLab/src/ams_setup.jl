export AMSDesign, XCV, TwistedCV, CommittorCV, collective_variable, level_min, level_max, dirac_point, AMSSimulation,
       AMSResult

struct XCV end
(::XCV)(x::Point) = x[1]

struct TwistedCV end
(::TwistedCV)(x::Point) = (x[1] - x[2]) / sqrt(2.0)

struct CommittorCV{F}
    x0::F
    y0::F
    hx::F
    hy::F
    q::Matrix{F}
end

function (ξ::CommittorCV{F})(x::Point) where {F}
    x in core_b && return one(F)
    x in core_a && return zero(F)
    nx, ny = size(ξ.q)
    u = clamp((x[1] - ξ.x0) / ξ.hx, 0.0, nx - 1.000001)
    v = clamp((x[2] - ξ.y0) / ξ.hy, 0.0, ny - 1.000001)
    i, j = floor(Int, u), floor(Int, v)
    a, b = u - i, v - j
    q = ξ.q
    @inbounds return (1 - a) * (1 - b) * q[i + 1, j + 1] + a * (1 - b) * q[i + 2, j + 1] +
                     (1 - a) * b * q[i + 1, j + 2] + a * b * q[i + 2, j + 2]
end

const _committor = Ref{CommittorCV{Float64}}()

function _load_committor!(path)
    isfile(path) || return
    nx, ny, x0, x1, y0, y1 = parse.(Float64, split(readline(path)))
    q = permutedims(readdlm(path, Float64; skipstart = 1))
    _committor[] = CommittorCV(x0, y0, (x1 - x0) / (nx - 1), (y1 - y0) / (ny - 1), q)
end

Base.@kwdef struct AMSDesign
    cv::String = "x"
    zmin::Union{String, Float64} = "medium"
    zmax::Union{String, Float64} = "B"
    N_loop::Int = 100
    L::Int = 10
    n_corr::Int = 1
    loop_init::String = "gibbs"
    N::Int = 100
    m::Int = 1
    function AMSDesign(cv, zmin, zmax, N_loop, L, n_corr, loop_init, N, m)
        cv in ("x", "twisted", "committor") ||
            throw(ArgumentError("unknown cv \"$cv\"; use \"x\", \"twisted\" or \"committor\""))
        zmin isa Real || zmin in ("low", "medium", "high") ||
            throw(ArgumentError("unknown zmin \"$zmin\"; use \"low\", \"medium\", \"high\" or a number"))
        zmax isa Real || zmax in ("B", "low", "medium", "high") ||
            throw(ArgumentError("unknown zmax \"$zmax\"; use \"B\", \"low\", \"medium\", \"high\" or a number"))
        zmax isa Real && zmax > _min_B(cv) &&
            throw(ArgumentError("need zmax ≤ min_B ξ = $(_min_B(cv)) for cv = \"$cv\", got $zmax"))
        loop_init in ("gibbs", "dirac") ||
            throw(ArgumentError("unknown loop_init \"$loop_init\"; use \"gibbs\" or \"dirac\""))
        N_loop >= 1 || throw(ArgumentError("N_loop must be ≥ 1, got $N_loop"))
        0 <= n_corr < L ||
            throw(ArgumentError("need 0 ≤ n_corr < L, got n_corr = $n_corr, L = $L"))
        1 <= m < N || throw(ArgumentError("need 1 ≤ m < N, got m = $m, N = $N"))
        s = new(cv, zmin isa Real ? Float64(zmin) : zmin, zmax isa Real ? Float64(zmax) : zmax,
                N_loop, L, n_corr, loop_init, N, m)
        level_min(s) < level_max(s) || throw(ArgumentError(
            "need z_min < z_max, got z_min = $(level_min(s)), z_max = $(level_max(s))"))
        return s
    end
end

collective_variable(s::AMSDesign) =
    s.cv == "x" ? XCV() : s.cv == "twisted" ? TwistedCV() : _committor[]

const _zmin = Dict(("x", "low") => -0.72, ("x", "medium") => -0.60,
                   ("x", "high") => -0.44, ("twisted", "low") => -0.40,
                   ("twisted", "medium") => -0.307, ("twisted", "high") => -0.145,
                   ("committor", "low") => 0.005,
                   ("committor", "medium") => 0.02, ("committor", "high") => 0.1)

const _zmax = Dict(("x", "low") => 0.448, ("x", "medium") => 0.603, ("x", "high") => 0.728,
                   ("twisted", "low") => 0.177, ("twisted", "medium") => 0.306,
                   ("twisted", "high") => 0.427,
                   ("committor", "low") => 0.9, ("committor", "medium") => 0.98,
                   ("committor", "high") => 0.995)

_min_B(cv) = cv == "x" ? core_b.c[1] - core_b.r :
             cv == "twisted" ? TwistedCV()(core_b.c) - core_b.r : 1.0

level_min(s::AMSDesign) = s.zmin isa Float64 ? s.zmin : _zmin[(s.cv, s.zmin)]
level_max(s::AMSDesign) = s.zmax isa Float64 ? s.zmax :
                          s.zmax == "B" ? _min_B(s.cv) : _zmax[(s.cv, s.zmax)]

function dirac_point(::AMSDesign)
    u = saddle_left_right - core_a.c
    return core_a.c + core_a.r * u / sqrt(sum(abs2, u))
end

struct AMSSimulation{F}
    p::F
    T_loop::F
    T_reac::F
    angles::Vector{F}
    top_geom::Vector{Bool}
    n_iterations::Int
    cost_loops::StepCounter
    cost_splitting::StepCounter
end

struct AMSResult
    setup::AMSDesign
    simulations::Vector{AMSSimulation{Float64}}
end
