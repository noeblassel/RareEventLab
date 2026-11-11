export Ball, core_a, core_b, core_top, HalfSpace, MetastableState, MetastableStateDefinitions, state_system,
       inside, assign_state, angle_on, y_threshold, Point

const Point{F} = SVector{2, F}

struct Ball{F}
    c::Point{F}
    r::F
end
Base.in(x::Point, b::Ball) = sum(abs2, x - b.c) <= b.r^2

const left_minimum = Point(-1.0480549928242202, -0.042093666306677734)
const top_minimum = Point(0.0, 1.5370820044494633)
const right_minimum = Point(1.0480549928242202, -0.042093666306677734)
const core_radius = 0.3

const core_a = Ball(left_minimum, core_radius)
const core_b = Ball(right_minimum, core_radius)
const core_top = Ball(top_minimum, core_radius)

const y_threshold = 0.6 # height separating paths going through top state from paths going through high-energy saddle point
angle_on(b::Ball, x::Point) = atan(x[2] - b.c[2], x[1] - b.c[1])

struct HalfSpace{F}
    z::Point{F}
    u::Point{F}
    c::F
end
inside(h::HalfSpace, x::Point) = (x - h.z)' * h.u < h.c

struct MetastableState{N,F}
    planes::NTuple{N, HalfSpace{F}}
    hole::Ball{F}
    core::Ball{F}
    side::Int
end

inside(Ω::MetastableState, x::Point) = all(h -> inside(h, x), Ω.planes) && !(x in Ω.hole)
Base.in(x::Point, Ω::MetastableState) = inside(Ω, x)

struct MetastableStateDefinitions{N,F}
    name::String
    states::Vector{MetastableState{N,F}}
end

const _nohole = Ball(Point(0.0, 0.0), 0.0)

# saddle points and their eigenvectors
const saddle_left_top = Point(-0.6172723078764598, 1.1027345175080963)
const saddle_right_top = Point(0.6172723078764598, 1.1027345175080963)
const saddle_left_right = Point(0.0, -0.19999999999972246)
const unstable_vec_left_top = Point(0.6080988038706289, 0.793861351075306)
const unstable_vec_right_top = Point(-0.6080988038706289, 0.793861351075306)
const unstable_vec_right_left = Point(-1.0, 0.0)

const opt_shifts = (0.3, 0.5, 0.5)

function _three_states(name, (a_top, a_low, b))
    left = MetastableState((HalfSpace(saddle_left_top, unstable_vec_left_top, a_top),
                            HalfSpace(saddle_left_right, -unstable_vec_right_left, a_low)),
                           _nohole, core_a, -1)
    top = MetastableState((HalfSpace(saddle_left_top, -unstable_vec_left_top, b),
                           HalfSpace(saddle_right_top, -unstable_vec_right_top, b)),
                          _nohole, core_top, 0)
    right = MetastableState((HalfSpace(saddle_right_top, unstable_vec_right_top, a_top),
                             HalfSpace(saddle_left_right, unstable_vec_right_left, a_low)),
                            _nohole, core_b, 1)
    return MetastableStateDefinitions(String(name), [left, top, right])
end

function state_system(name::AbstractString)
    name == "optimized" && return _three_states(name, opt_shifts)
    name == "basin" && return _three_states(name, (0.0, 0.0, 0.0))
    if name == "saddle hyperplane"
        e1 = Point(1.0, 0.0)
        return MetastableStateDefinitions(String(name),
                                          [MetastableState((HalfSpace(saddle_left_right, e1, 0.0),), _nohole, core_a, -1),
                                           MetastableState((HalfSpace(saddle_left_right, -e1, 0.0),), _nohole, core_b, 1)])
    elseif name == "core complement"
        return MetastableStateDefinitions(String(name),
                                          [MetastableState((), core_b, core_a, -1), MetastableState((), core_a, core_b, 1)])
    end
    throw(ArgumentError("unknown metastable state definition \"$name\"; use \"saddle hyperplane\", " *
                        "\"optimized\", \"basin\" or \"core complement\""))
end

function _basin(x::Point)
    l12 = (x - saddle_left_top)' * unstable_vec_left_top
    l23 = (x - saddle_right_top)' * unstable_vec_right_top
    l12 >= 0 && l23 >= 0 && return 0
    return x[1] < 0 ? -1 : 1
end

function assign_state(sys::MetastableStateDefinitions, x::Point)
    b = _basin(x)
    first_in = 0
    for (i, Ω) in enumerate(sys.states)
        inside(Ω, x) || continue
        Ω.side == b && return i
        first_in == 0 && (first_in = i)
    end

    return first_in # 0 = no man's land
end
