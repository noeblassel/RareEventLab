module RareEventLab

using StaticArrays, Random, Base.Threads, Serialization, DelimitedFiles, Plots

include("geometry.jl")
include("dynamics.jl")
include("tracking.jl")
include("coarse.jl")
include("parrep_setup.jl")
include("parrep.jl")
include("bruteforce.jl")
include("ams_setup.jl")
include("ams.jl")
include("analysis.jl")
include("leaderboard.jl")
include("plotting.jl")
include("source.jl")
include("animate.jl")

function __init__()
    _load_tcorr!(joinpath(@__DIR__, "..", "data", "tcorr.csv"))
    _load_committor!(joinpath(@__DIR__, "..", "data", "committor.grid"))
    _load_reference!(joinpath(@__DIR__, "..", "data", "reference.csv"))
end

end
