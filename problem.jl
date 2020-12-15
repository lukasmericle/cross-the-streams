export ProblemCol
export Problem
export rows, cols
import Base: size, string
export       size, string

using Lazy: @forward

ProblemCol = Array{ProblemNumeral, 1}

Base.string(pcol::ProblemCol) = join(Base.string.(pcol), " ")

struct Problem
    rows::Array{ProblemCol}
    cols::Array{ProblemCol}
end

rows(problem::Problem) = problem.rows
cols(problem::Problem) = problem.cols

function Base.size(problem::Problem, dim::Int)
    if dim == 1
        return length(rows(problem))
    elseif dim == 2
        return length(cols(problem))
    end
end

Base.size(problem::Problem) = (Base.size(problem, 1), Base.size(problem, 2))

function Base.string(problem::Problem)
    nr, nc = size(problem)
    s = string(nr) * " " * string(nc) * "\n"
    s *= join(string.(rows(problem)), "\n") * "\n"
    s *= join(string.(cols(problem)), "\n")
    s
end
