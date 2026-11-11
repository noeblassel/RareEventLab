export StepCounter, ReactiveEvent, MilestoneObserver, observe_step!, advance_clock!, clock,
       n_reactions

mutable struct StepCounter
    wall::Int # parallel force calls
    total::Int # total force calls
end
StepCounter() = StepCounter(0, 0)

struct ReactiveEvent
    time::Float64
    to::Int # +1: A→B, -1: B→A
    angle::Float64 # entry point around the target milestone
    top_geom::Bool # last crossing of {x = 0} is above y_threshold
end

mutable struct MilestoneObserver
    last_milestone::Int # -1: A, +1: B, 0: none yet
    last_cross_y::Float64
    clock::Float64
    time_since::NTuple{2, Float64} # since last visit of A, of B
    events::Vector{ReactiveEvent}
end

MilestoneObserver(x::Point) =
    MilestoneObserver(x in core_a ? -1 : (x in core_b ? 1 : 0), NaN, 0.0, (0.0, 0.0),
                      ReactiveEvent[])

clock(obs::MilestoneObserver) = obs.clock
n_reactions(obs::MilestoneObserver) = length(obs.events)

function advance_clock!(obs::MilestoneObserver, t)
    obs.clock += t
    obs.last_milestone == -1 && (obs.time_since = obs.time_since .+ (t, 0.0))
    obs.last_milestone == 1 && (obs.time_since = obs.time_since .+ (0.0, t))
    return obs
end

function observe_step!(obs::MilestoneObserver, x::Point, x′::Point)
    advance_clock!(obs, Δt)
    if signbit(x[1]) != signbit(x′[1])
        obs.last_cross_y = x[2] + x[1] / (x[1] - x′[1]) * (x′[2] - x[2])
    end
    for (side, milestone) in ((-1, core_a), (1, core_b))
        x′ in milestone && obs.last_milestone != side || continue
        if obs.last_milestone == -side
            push!(obs.events, ReactiveEvent(obs.clock, side, angle_on(milestone, x′),
                                            obs.last_cross_y > y_threshold))
        end
        obs.last_milestone = side
    end
    return obs
end
