export CoarseTrajectory, current_state, next_state!, coarse_pathways

struct StateVisit
    time::Float64
    state::Int
    entry::Point{Float64}
end

struct CoarseTrajectory{N, F}
    states::MetastableStateDefinitions{N, F}
    visits::Vector{StateVisit}
end

CoarseTrajectory(states::MetastableStateDefinitions, x::Point) =
    CoarseTrajectory(states, [StateVisit(0.0, assign_state(states, x), x)])

current_state(c::CoarseTrajectory) = c.states.states[last(c.visits).state]

function next_state!(c::CoarseTrajectory, x::Point, t)
    push!(c.visits, StateVisit(t, assign_state(c.states, x), x))
    return current_state(c)
end

function coarse_pathways(c::CoarseTrajectory, events::Vector{ReactiveEvent})
    pathways = Vector{Union{Bool, Missing}}(missing, length(events))
    c.states.name == "core complement" && return pathways
    side(v) = c.states.states[v.state].side
    has_top_state = any(Ω -> Ω.side == 0, c.states.states)
    via_top = false
    k = 1
    for (i, e) in enumerate(events)
        while k < length(c.visits) && c.visits[k + 1].time <= e.time
            previous, visit = c.visits[k], c.visits[k + 1]
            if side(visit) != 0 && side(visit) != side(previous)
                via_top = has_top_state ? side(previous) == 0 : visit.entry[2] > y_threshold
            end
            k += 1
        end
        pathways[i] = via_top
    end
    return pathways
end
