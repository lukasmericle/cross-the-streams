export CrossTheStreamsPuzzle

# generating puzzles
## judging difficulty by:
## * number of iterations it takes solver
## * some way to compute or approximate the "information content" of the problem's description
## * some way to compute or approximate the "entropy" of the solution path
## or just tune a stochastic algorithm to produce decent results, and curate by hand


# mutate puzzle clues.
@views function mutate!(cvec::T; ps::Vector{Float64}=ones(2)) where {T<:AbstractClueVector}
    """Randomize (mutate) puzzle row."""
    fn = sample([int2qmk!, nums2ask!], Weights(ps))
    fn(cvec)
    clean!(cvec)
end

@views function clean!(cvec::T) where {T<:AbstractClueVector}
    """Eliminate duplicate Asterisks."""
    all_asks = findall(isa.(cvec, Asterisk))
    (length(all_asks) < 2) && return
    dist_to_next_ask = all_asks[2:end] .- all_asks[1:end-1]
    duplicates = findall(dist_to_next_ask .=== 1)
    (length(duplicates) === 0) && return
    deleteat!(cvec, duplicates)
end

@views function int2qmk!(cvec::T) where {T<:AbstractClueVector}
    """Replaces an int with a question mark."""
    all_ints = findall(isa.(cvec, Int))
    (length(all_ints) === 0) && return
    cvec[rand(all_ints)] = QuestionMark()
end

@views function nums2ask!(cvec::T) where {T<:AbstractClueVector}
    """Replaces anywhere from 0 to all of the clues with an asterisk."""
    if ((length(cvec) === 0) && rand() < 0.5)
        push!(cvec, Asterisk())
        return
    end
    replace_len = rand(0:length(cvec))
    replace_loc = rand(1:(length(cvec)-(replace_len-1)))
    deleteat!(cvec, replace_loc:(replace_loc+(replace_len-1)))
    insert!(cvec, replace_loc, Asterisk())
end

function mutate!(puzzle::P) where {P<:Puzzle}
    """Mutate a random row or col in a puzzle."""
    if rand(Bool)
        i = rand(1:size(puzzle, 1))
        mutate!(rows(puzzle)[i])
        return ('R', i)
    else
        j = rand(1:size(puzzle, 2))
        mutate!(cols(puzzle)[j])
        return ('C', j)
    end
end


# generate puzzles.
function CrossTheStreamsPuzzle(smat::T, difficulty::S) where {T<:TwoDSolutionCellArray,S<:AbstractFloat}
    """Generate a puzzle from an existing solution with optional mutation."""
    puzzle = CrossTheStreamsPuzzle(smat)
    (difficulty <= 0) && return puzzle
    counter = MatrixStateCounter(puzzle)
    initn = convert(BigFloat, n(counter))
    dummy_row = convert(CellVector, fill(missing, size(puzzle, 2)))
    dummy_col = convert(CellVector, fill(missing, size(puzzle, 1)))
    while ((convert(BigFloat, n(counter)) / initn) < exp(difficulty))
        rc = mutate!(puzzle)
        if rc[1] === 'R'
            i = rc[2]
            recount!(rows(counter)[i], dummy_row, rows(puzzle)[i])
        elseif rc[1] === 'C'
            j = rc[2]
            recount!(cols(counter)[j], dummy_col, cols(puzzle)[j])
        else
            error("rowcol must either refer to row ('R') or col ('C')")
        end
    end
    puzzle
end
