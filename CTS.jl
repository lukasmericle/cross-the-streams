module CTS
    using StatsBase: sample, Weights
    using Statistics: mean
    include("cellarray.jl")
    include("counter.jl")
    include("puzzle.jl")
    include("gen_solution.jl")
    include("gen_puzzle.jl")
    include("solve.jl")
end
