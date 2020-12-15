export generate_problem

# generating problems
## judging difficulty by:
## * number of iterations it takes solver
## * some way to compute or approximate the "information content" of the problem's description
## * some way to compute or approximate the "entropy" of the solution path
## or just tune a stochastic algorithm to produce decent results, and curate by hand


function mutate!(pcol::ProblemCol)
    # replace an integer with a question mark
    # replace a numeral with an asterisk
    # merge a numeral with an asterisk
    # insert an asterisk
end

function clean!(pcol::ProblemCol)
    """
    Eliminate duplicate Asterisks.
    """
    all_asks = findall(isa.(pcol, Asterisk))
    (length(all_asks) < 2) && return
    @. dist_to_next_ask = all_asks[2:end] - all_asks[1:end-1]
    duplicates = findall(dist_to_next_ask .=== 1)
    (length(duplicates) == 0) && return
    deleteat!(pcol, duplicates)
end

function int2qmk!(pcol::ProblemCol)
    all_ints = findall(isa.(pcol, Integer))
    if length(all_ints) > 0
        pcol[rand(all_ints)] = QuestionMark()
    end
end

function num2ask!(pcol::ProblemCol)
    this_num = rand(1:length(pcol))
    pcol[this_num] = Asterisk()
end

function merge_num_ask!(pcol::ProblemCol)
    (length(pcol) < 2) && return
    all_asks = findall(isa.(pcol, Asterisk))
    (length(all_asks) == 0) && return
    merge_candidates = filter(x->inbounds(pcol, x), rand(all_asks) .+ VON_NEUMANN_NEIGHBORHOOD_1D)
    (length(merge_candidates) == 0) && return
    deleteat!(pcol, rand(merge_candidates))
end

function insert_ask!(pcol::ProblemCol)
    i = rand(1:(length(pcol)+1))
    insert!(pcol, i, Asterisk())
end

function get_problem_col(scol::Union{SolutionGridCol,SubArray{Bool,1,SolutionGrid}})
    numerals = ProblemNumeral[]
    i = 1
    c = 0
    while i <= length(scol)
        if scol[i]
            c += 1
        elseif c > 0
            push!(numerals, c)
            c = 0
        end
        i += 1
    end
    (c > 0) && push!(numerals, c)
    numerals
end

get_problem_cols(solution::SolutionGrid) = map(get_problem_col, eachcol(solution))
get_problem_rows(solution::SolutionGrid) = map(get_problem_col, eachrow(solution))

function generate_problem(solution::SolutionGrid)
    prows = get_problem_rows(solution)
    pcols = get_problem_cols(solution)
    # while xxx
    #     # ...mutate rows and cols...
    # end
    Problem(prows, pcols)
end
