export generate_puzzle

# generating puzzles
## judging difficulty by:
## * number of iterations it takes solver
## * some way to compute or approximate the "information content" of the problem's description
## * some way to compute or approximate the "entropy" of the solution path
## or just tune a stochastic algorithm to produce decent results, and curate by hand

function mutate!(cvec::T) where T <: AbstractClueVector
    # replace an integer with a question mark
    # replace a numeral with an asterisk
    # merge a numeral with an asterisk
    # insert an asterisk
end

function clean!(cvec::T) where T <: AbstractClueVector
    """
    Eliminate duplicate Asterisks.
    """
    all_asks = findall(isa.(cvec, Asterisk))
    (length(all_asks) < 2) && return
    @. dist_to_next_ask = all_asks[2:end] - all_asks[1:end-1]
    duplicates = findall(dist_to_next_ask .=== 1)
    (length(duplicates) === 0) && return
    deleteat!(cvec, duplicates)
end

function int2qmk!(cvec::T) where T <: AbstractClueVector
    all_ints = findall(isa.(cvec, Int))
    ((length(all_ints) > 0) && (cvec[rand(all_ints)] = QuestionMark()))
end

num2ask!(cvec::T) where T <: AbstractClueVector = (cvec[rand(1:length(cvec))] = Asterisk())

function merge_num_ask!(cvec::T) where T <: AbstractClueVector
    (length(cvec) < 2) && return
    all_asks = findall(isa.(cvec, Asterisk))
    (length(all_asks) === 0) && return
    merge_candidates = filter(x->inbounds(cvec, x), VON_NEUMANN_NEIGHBORHOOD_1D .+ rand(all_asks))
    (length(merge_candidates) === 0) && return
    deleteat!(cvec, rand(merge_candidates))
end

insert_ask!(cvec::T) where T <: AbstractClueVector = insert!(cvec, rand(1:(length(cvec)+1)), Asterisk())

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

function generate_puzzle(smat::T) where T <: TwoDSolutionCellArray
    prows = get_puzzle_rows(smat)
    pcols = get_puzzle_cols(smat)
    # while xxx
    #     # ...mutate rows and cols...
    # end
    Puzzle(prows, pcols)
end
generate_puzzle(n::Int, m::Int) = generate_puzzle(generate_solution(n, m))
generate_puzzle(n::Int) = generate_puzzle(n, n)
