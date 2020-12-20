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
    m = 0
    while (m === 0)
        this_ask = OneDCoord(rand(all_asks))
        merge_candidates = filter(x -> inbounds(cvec, x),
                                  broadcast(vnn -> vnn + this_ask, VON_NEUMANN_NEIGHBORHOOD_1D))
        m = length(merge_candidates)
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

function generate_puzzle(smat::T; num_mutations::Int=0) where T <: TwoDSolutionCellArray
    prows = get_puzzle_rows(smat)
    pcols = get_puzzle_cols(smat)
    for i=1:num_mutations
        rand(Bool) ? mutate!(rand(prows)) : mutate!(rand(pcols))
    end
    Puzzle(prows, pcols)
end
generate_puzzle(n::Int, m::Int; num_mutations::Int=0) = generate_puzzle(generate_solution(n, m), num_mutations=num_mutations)
generate_puzzle(n::Int; num_mutations::Int=0) = generate_puzzle(n, n, num_mutations=num_mutations)
