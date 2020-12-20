export generate_puzzle

# generating puzzles
## judging difficulty by:
## * number of iterations it takes solver
## * some way to compute or approximate the "information content" of the problem's description
## * some way to compute or approximate the "entropy" of the solution path
## or just tune a stochastic algorithm to produce decent results, and curate by hand

@views function mutate!(cvec::T; ps::Vector{Float64}=ones(4)) where T <: AbstractClueVector
    fn = sample([int2qmk!, num2ask!, merge_num_ask!, insert_ask!], Weights(ps))
    fn(cvec)
    clean!(cvec)
end

@views function clean!(cvec::T) where T <: AbstractClueVector
    """
    Eliminate duplicate Asterisks.
    """
    all_asks = findall(isa.(cvec, Asterisk))
    (length(all_asks) < 2) && return
    dist_to_next_ask = all_asks[2:end] .- all_asks[1:end-1]
    duplicates = findall(dist_to_next_ask .=== 1)
    (length(duplicates) === 0) && return
    deleteat!(cvec, duplicates)
end

@views function int2qmk!(cvec::T) where T <: AbstractClueVector
    all_ints = findall(isa.(cvec, Int))
    (length(all_ints) === 0) && return
    cvec[rand(all_ints)] = QuestionMark()
end

@views function num2ask!(cvec::T) where T <: AbstractClueVector
    (length(cvec) === 0) && return
    cvec[rand(1:length(cvec))] = Asterisk()
end

@views function merge_num_ask!(cvec::T) where T <: AbstractClueVector
    (length(cvec) < 2) && return
    all_asks = findall(isa.(cvec, Asterisk))
    (length(all_asks) === 0) && return
    merge_candidates = OneDCoord[]
    while (length(merge_candidates) === 0)
        this_ask = OneDCoord(rand(all_asks))
        merge_candidates = filter(x -> inbounds(cvec, x),
                                  broadcast(vnn -> vnn + this_ask, VON_NEUMANN_NEIGHBORHOOD_1D))
    end
    deleteat!(cvec, rand(merge_candidates))
end

@views insert_ask!(cvec::T) where T <: AbstractClueVector = insert!(cvec, rand(1:(length(cvec)+1)), Asterisk())

@views function get_puzzle_col(svec::T) where T <: OneDSolutionCellArray
    clues = Clue[]
    clue = 0
    for cell=svec
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

@views get_puzzle_cols(smat::T) where T <: TwoDSolutionCellArray = map(get_puzzle_col, eachcol(smat))
@views get_puzzle_rows(smat::T) where T <: TwoDSolutionCellArray = map(get_puzzle_col, eachrow(smat))

Puzzle(smat::T) where T <: TwoDSolutionCellArray = Puzzle(get_puzzle_rows(smat), get_puzzle_cols(smat))

function generate_puzzle(smat::T; difficulty::S=0.0) where {T <: TwoDSolutionCellArray, S <: AbstractFloat}
    puzzle = Puzzle(smat)
    counter = MatrixStateCounter(puzzle)
    initial_num_states = convert(BigFloat, num_states(counter))
    dummy_row = convert(CellVector, fill(missing, size(puzzle, 2)))
    dummy_col = convert(CellVector, fill(missing, size(puzzle, 1)))
    while ((convert(BigFloat, num_states(counter)) / initial_num_states) < exp(difficulty))
        if rand(Bool)
            i = rand(1:size(puzzle, 1))
            mutate!(rows(puzzle)[i])
            recount!(rows(counter)[i], dummy_row, rows(puzzle)[i])
        else
            j = rand(1:size(puzzle, 2))
            mutate!(cols(puzzle)[j])
            recount!(cols(counter)[j], dummy_col, cols(puzzle)[j])
        end
    end
    puzzle
end
generate_puzzle(n::Int, m::Int; difficulty::T=1.0) where T <: AbstractFloat = generate_puzzle(generate_solution(n, m), difficulty=difficulty)
generate_puzzle(n::Int; difficulty::T=1.0) where T <: AbstractFloat = generate_puzzle(n, n, difficulty=difficulty)
