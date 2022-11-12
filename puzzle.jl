export QuestionMark, Asterisk, Clue, ClueVector, MatrixStateCounter
export NonogramPuzzle, CrossTheStreamsPuzzle, rows, cols
export complexity, entropy


# types for wildcard clues, and vectors and views of clues.
struct QuestionMark end
struct Asterisk end

Clue = Union{T,QuestionMark,Asterisk} where {T<:Int}
ClueRepr = Union{T,String} where {T<:Int}

ClueVector = Vector{Clue}
ClueVectorView = SubArray{Clue,1,ClueVector}
AbstractClueVector = Union{ClueVector,ClueVectorView}


# representing a puzzle.
abstract type Puzzle end

# distinguish between puzzle types for solver
# (nonogram solution conditions are a subset of CTS solution conditions).
struct NonogramPuzzle <: Puzzle
    rows::Vector{Vector{Int}}
    cols::Vector{Vector{Int}}
end

struct CrossTheStreamsPuzzle <: Puzzle
    rows::Vector{ClueVector}
    cols::Vector{ClueVector}
end

rows(pzl::Puzzle) = pzl.rows
cols(pzl::Puzzle) = pzl.cols

Base.size(pzl::Puzzle, dim::Int) = (dim === 1) ? length(rows(pzl)) : ((dim === 2) ? length(cols(pzl)) : error("puzzle has only two dimensions"))
Base.size(pzl::Puzzle) = (length(rows(pzl)), length(cols(pzl)))

# get puzzle clue for column
@views function get_puzzle_col(svec::T) where {T<:OneDSolutionCellArray}
    clues = Clue[]
    clue = 0
    for cell = svec
        if cell
            clue += 1
        elseif (clue > 0)
            push!(clues, clue)
            clue = 0
        end
    end
    (clue > 0) && push!(clues, clue)
    clues
end

@views get_puzzle_cols(smat::T) where {T<:TwoDSolutionCellArray} = map(get_puzzle_col, eachcol(smat))
@views get_puzzle_rows(smat::T) where {T<:TwoDSolutionCellArray} = map(get_puzzle_col, eachrow(smat))


# generate a puzzle from a solution
NonogramPuzzle(smat::T) where {T<:TwoDSolutionCellArray} = NonogramPuzzle(get_puzzle_rows(smat), get_puzzle_cols(smat))

function CrossTheStreamsPuzzle(smat::T) where {T<:TwoDSolutionCellArray}
    (!iscontiguous(smat) || overcrowded(smat)) && error("solution is not valid for Cross The Streams puzzle")
    CrossTheStreamsPuzzle(get_puzzle_rows(smat), get_puzzle_cols(smat))
end


# string methods
Base.string(n::QuestionMark) = "?"
Base.string(n::Asterisk) = "*"
Base.string(cvec::T) where {T<:AbstractClueVector} = join(Base.string.(cvec), " ")
function Base.string(pzl::Puzzle)
    nr, nc = size(pzl)
    s = string(nr) * " " * string(nc) * "\n"
    s *= join(string.(rows(pzl)), "\n") * "\n"
    s *= join(string.(cols(pzl)), "\n")
    s
end


# helpers.
@views inbounds(cvec::T, i::OneDCoord) where {T<:AbstractClueVector} = (1 <= i[1]) && (i[1] <= length(cvec))


# Initialize cell matrix using Puzzle directly.
init_cmat(puzzle::Puzzle) = init_cmat(size(puzzle)...)


# Initialize MatrixStateCounter using Puzzle directly.
@views function MatrixStateCounter(puzzle::Puzzle)
    """Initialize a MatrixStateCounter from a set of clues."""
    dummy_row = init_cvec(size(puzzle, 2))
    dummy_col = init_cvec(size(puzzle, 1))
    MatrixStateCounter(
        map(row -> count_states(dummy_row, row), rows(puzzle)),
        map(col -> count_states(dummy_col, col), cols(puzzle))
    )
end

@views function MatrixStateCounter(cmat::T, puzzle::Puzzle) where {T<:TwoDAbstractCellArray}
    """Initialize a MatrixStateCounter from a set of clues and a solution state."""
    MatrixStateCounter(
        map(i -> count_states(cmat[i, :], rows(puzzle)[i]), 1:size(cmat, 1)),
        map(j -> count_states(cmat[:, j], cols(puzzle)[j]), 1:size(cmat, 2))
    )
end


# Measure complexity using Puzzle directly.
complexity(puzzle::Puzzle) = complexity(MatrixStateCounter(puzzle))
