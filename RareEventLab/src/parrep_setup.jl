export ParRepDesign, decorrelation_time, in_core, ParRepResult

Base.@kwdef struct ParRepDesign
    states::String = "optimized" # definition of metastable states
    coresets::String = "physical" # definition of core sets
    tcorr::Union{String, Float64} = "freefem" # definitions of decorrelation times
    N::Int = 100 # number of replicas
    function ParRepDesign(states, coresets, tcorr, N)
        state_system(states)
        coresets in ("physical", "whole state") ||
            throw(ArgumentError("unknown coresets \"$coresets\"; use \"physical\" or \"whole state\""))
        tcorr isa Real || tcorr in ("heuristic", "harmonic", "freefem") ||
            throw(ArgumentError("unknown tcorr \"$tcorr\"; use \"heuristic\", \"harmonic\", \"freefem\" or a number"))
        N >= 1 || throw(ArgumentError("N must be at least 1, got $N"))
        return new(states, coresets, tcorr isa Real ? Float64(tcorr) : tcorr, N)
    end
end

const tv_tol = 0.05 # total variation threshold for mixing times

# smallest Hessian eigenvalues at the minima
const ν1_left_right, ν1_top = 8.472, 3.798

const _TCORR_FF = Dict{Tuple{String, Int}, Float64}() # for freefem decorrelation times

function _load_tcorr!(path)
    isfile(path) || return
    for row in eachrow(first(readdlm(path, ','; header = true)))
        sides = row[2] == "top" ? (0,) : (-1, 1)
        foreach(side -> _TCORR_FF[(String(row[1]), side)] = Float64(row[6]), sides)
    end
end

function decorrelation_time(setup::ParRepDesign, Ω::MetastableState)
    setup.tcorr isa Float64 && return setup.tcorr
    setup.tcorr == "harmonic" && return log(1 / tv_tol) / (Ω.side == 0 ? ν1_top : ν1_left_right)
    setup.tcorr == "freefem" && return _TCORR_FF[(setup.states, Ω.side)]
    return 0.3 * minimum(values(_TCORR_FF))
end

in_core(setup::ParRepDesign, Ω::MetastableState, x::Point) =
    setup.coresets == "whole state" || x in Ω.core

struct ParRepResult
    setup::ParRepDesign
    events::Vector{ReactiveEvent}
    top_coarse::Vector{Union{Bool, Missing}}
    sim_time::Float64
    cost::StepCounter
    trace_wall::Vector{Int}
    trace_sim::Vector{Float64}
    n_eff::Vector{Int}
    n_aborts::Int
end
