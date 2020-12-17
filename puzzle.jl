export QuestionMark, Asterisk, Clue, ClueVector
export Puzzle
export rows, cols
export complexity
import Base: size, string
export       size, string

struct QuestionMark end
struct Asterisk end

Clue = Union{T, QuestionMark, Asterisk} where T <: Int
ClueRepr = Union{T, String} where T <: Int

ClueVector = Array{Clue, 1}
ClueVectorView = SubArray{Clue, 1, ClueVector}
AbstractClueVector = Union{ClueVector, ClueVectorView}

Base.string(n::QuestionMark) = "?"
Base.string(n::Asterisk) = "*"
Base.string(cvec::T) where T <: AbstractClueVector = join(Base.string.(cvec), " ")

struct Puzzle
    rows::Array{ClueVector}
    cols::Array{ClueVector}
end

rows(pzl::Puzzle) = pzl.rows
cols(pzl::Puzzle) = pzl.cols

Base.size(pzl::Puzzle, dim::Int) = (dim == 1) ? length(rows(pzl)) : ((dim == 2) ? length(cols(pzl)) : error("puzzle has only two dimensions"))
Base.size(pzl::Puzzle) = (Base.size(pzl, 1), Base.size(pzl, 2))

function Base.string(pzl::Puzzle)
    nr, nc = size(pzl)
    s = string(nr) * " " * string(nc) * "\n"
    s *= join(string.(rows(pzl)), "\n") * "\n"
    s *= join(string.(cols(pzl)), "\n")
    s
end

complexity(pzl::Puzzle) = complexity(MatrixStateCounter(pzl))
