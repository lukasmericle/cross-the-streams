export ClueVector
export Puzzle
export rows, cols
import Base: size, string
export       size, string

using Lazy: @forward

ClueVector = Array{Clue, 1}
ClueVectorView = SubArray{Clue, 1, ClueVector}
AbstractClueVector = Union{ClueVector, ClueVectorView}

Base.string(cvec::T) where T <: AbstractClueVector = join(Base.string.(cvec), " ")

struct Puzzle
    rows::Array{ClueVector}
    cols::Array{ClueVector}
end

rows(pzl::Puzzle) = pzl.rows
cols(pzl::Puzzle) = pzl.cols

function Base.size(pzl::Puzzle, dim::Int)
    if dim == 1
        return length(rows(pzl))
    elseif dim == 2
        return length(cols(pzl))
    end
end

Base.size(pzl::Puzzle) = (Base.size(pzl, 1), Base.size(pzl, 2))

function Base.string(pzl::Puzzle)
    nr, nc = size(pzl)
    s = string(nr) * " " * string(nc) * "\n"
    s *= join(string.(rows(pzl)), "\n") * "\n"
    s *= join(string.(cols(pzl)), "\n")
    s
end
