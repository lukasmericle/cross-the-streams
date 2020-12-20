export QuestionMark, Asterisk, Clue, ClueVector
export Puzzle, rows, cols
export complexity, entropy

struct QuestionMark end
struct Asterisk end

Clue = Union{T, QuestionMark, Asterisk} where T <: Int
ClueRepr = Union{T, String} where T <: Int

ClueVector = Vector{Clue}
ClueVectorView = SubArray{Clue, 1, ClueVector}
AbstractClueVector = Union{ClueVector, ClueVectorView}

struct Puzzle
    rows::Vector{ClueVector}
    cols::Vector{ClueVector}
end

rows(pzl::Puzzle) = pzl.rows
cols(pzl::Puzzle) = pzl.cols

Base.size(pzl::Puzzle, dim::Int) = (dim === 1) ? length(rows(pzl)) : ((dim === 2) ? length(cols(pzl)) : error("puzzle has only two dimensions"))
Base.size(pzl::Puzzle) = (Base.size(pzl, 1), Base.size(pzl, 2))

Base.string(n::QuestionMark) = "?"
Base.string(n::Asterisk) = "*"
Base.string(cvec::T) where T <: AbstractClueVector = join(Base.string.(cvec), " ")
function Base.string(pzl::Puzzle)
    nr, nc = size(pzl)
    s = string(nr) * " " * string(nc) * "\n"
    s *= join(string.(rows(pzl)), "\n") * "\n"
    s *= join(string.(cols(pzl)), "\n")
    s
end

@views inbounds(cvec::T, i::OneDCoord) where T <: AbstractClueVector = (1 <= i[1]) && (i[1] <= length(cvec))
