export solve

@views function MatrixStateCounter(puzzle::Puzzle)  # initializes based on puzzle description assuming grid is blank
    dummy_row = init_cvec(size(puzzle, 2))
    dummy_col = init_cvec(size(puzzle, 1))
    MatrixStateCounter(map(i -> count_states(dummy_row, rows(puzzle)[i]), 1:size(puzzle, 1)),
                       map(j -> count_states(dummy_col, cols(puzzle)[j]), 1:size(puzzle, 2)))
end
@views function MatrixStateCounter(cmat::T, puzzle::Puzzle) where T <: TwoDAbstractCellArray  # initializes based on puzzle description and current state of grid
    MatrixStateCounter(map(i -> count_states(cmat[i,:], rows(puzzle)[i]), 1:size(cmat, 1)),
                       map(j -> count_states(cmat[:,j], cols(puzzle)[j]), 1:size(cmat, 2)))
end

@views function Base.sort!(rcs::Vector{Tuple{Char,Int}}, counter::MatrixStateCounter)
    """
    Sort the rcs so that the easiest rows and cols
    (with the fewest possible states) are first in the list.
    """
    row_idxs, col_idxs = Int[], Int[]
    for rc=rcs
        if rc[1] === 'R'
            push!(row_idxs, rc[2])
        elseif rc[1] === 'C'
            push!(col_idxs, rc[2])
        else
            error("rowcol must either refer to row ('R') or col ('C')")
        end
    end
    sort_ord = sortperm(vcat(map(n, rows(counter)[row_idxs]),   # get the order that the rcs should be in
                             map(n, cols(counter)[col_idxs])))  # (from fewest to most states, to encourage quickly acting on good information)
    rcs[:] .= vcat(map(r -> ('R', r), row_idxs),
                   map(c -> ('C', c), col_idxs))[sort_ord]  # rewrite to the same array
end

init_cmat(puzzle::Puzzle) = init_cmat(size(puzzle)...)

complexity(puzzle::Puzzle) = complexity(MatrixStateCounter(puzzle))

function minreq_cells(cluevec_ints::T, n_qmks::Int) where T <: AbstractClueVector
    n_reqd_cells = 0
    if (length(cluevec_ints) > 0)
        n_reqd_cells += sum(cluevec_ints) + (length(cluevec_ints) - 1)  # sum + (num - 1)
    end
    if (n_qmks > 0)
        (n_reqd_cells > 0) && (n_reqd_cells += 1)  # add one extra for a gap between seq of ints and seq of qmks
        n_reqd_cells += 2 * n_qmks - 1  # sum + (num - 1), assuming all ?s represent length-1 runs
    end
    n_reqd_cells
end

function space_for_ints(cluevec_ints_after::T) where T <: AbstractClueVector
    (length(cluevec_ints_after) === 0) && return 0
    sum(cluevec_ints_after) + length(cluevec_ints_after)
end

function space_for_qmks(all_qmks::Vector{Int}, i::Int)
    num_qmks_before = count(x -> (x < i), all_qmks)
    2 * num_qmks_before, 2 * (length(all_qmks) - num_qmks_before)  # num before, num after
end

@views function count_states(cellvec::T, cluevec::S) where {T <: OneDAbstractCellArray,  S <: AbstractClueVector}
    """
    Uses dynamic programming and branch-and-bound to test and accumulate
    all possible configurations of cells based on current col state and col description.
    """
    if (length(cluevec) === 0)
        if any(skipmissing(cellvec))
            return nothing  # if any trues but no clues, bound this branch
        else
            return VectorStateCounter(cellvec; init="empty")  # if no trues, no-op
        end
    end
    all_ints = findall(isa.(cluevec, Int))
    all_qmks = findall(isa.(cluevec, QuestionMark))
    all_asks = findall(isa.(cluevec, Asterisk))
    n_reqd_cells = minreq_cells(cluevec[all_ints], length(all_qmks))
    if (n_reqd_cells > length(cellvec))
        return nothing  # if the clues won't fit into the grid, bound this branch
    elseif (length(cellvec) === 0)
        return VectorStateCounter(cellvec; init="empty")  # reach this branch if (length(cellvec) === 0) and cluevec has only Asterisk in it; no-op
    end
    counter = VectorStateCounter(cellvec; init="blank")
    if (length(all_ints) > 0)
        first_int = all_ints[1]
        cluevec_before = cluevec[1:first_int-1]
        run_length = cluevec[first_int]
        cluevec_after = cluevec[first_int+1:end]
        space_for_ints_after = space_for_ints(cluevec[all_ints[2:end]])
        space_for_qmks_before, space_for_qmks_after = space_for_qmks(all_qmks, first_int)
        last_pos_apriori = length(cellvec) - (run_length - 1)  # minus 1 because ranges are inclusive in Julia
        first_pos = 1 + space_for_qmks_before
        last_pos = last_pos_apriori - (space_for_ints_after + space_for_qmks_after)
        for pos=first_pos:last_pos   # if first_pos > last_pos, this for loop is skipped and counter is returned blank (thus this branch is bounded)
            cellvec_middle = cellvec[pos:pos+(run_length-1)]
            any(isfalse.(cellvec_middle)) && continue  # continue if filled cell(s) should be empty
            counter_middle = VectorStateCounter(cellvec_middle; init="full")
            if (pos === 1)
                counter_before = VectorStateCounter(0; init="empty")
            else
                istrue(cellvec[pos-1]) && continue  # continue if empty cell should be filled
                cellvec_before = cellvec[1:pos-2]
                counter_before = count_states(cellvec_before, cluevec_before)
                isnothing(counter_before) && continue
                counter_before = counter_before * VectorStateCounter(1; init="empty")
            end
            if (pos === last_pos_apriori)
                counter_after = VectorStateCounter(0; init="empty")
            else
                istrue(cellvec[pos+(run_length-1)+1]) && continue  # continue if empty cell should be filled
                cellvec_after = cellvec[pos+(run_length-1)+2:end]
                counter_after = count_states(cellvec_after, cluevec_after)
                isnothing(counter_after) && continue
                counter_after = VectorStateCounter(1; init="empty") * counter_after
            end
            counter += vcat(counter_before, counter_middle, counter_after)
        end
    elseif (length(all_qmks) > 0)
        first_qmk = all_qmks[1]
        space_for_qmks_after = (length(all_qmks) > 0) ? (2 * (length(all_qmks) - 1)) : 0
        max_len = length(cellvec) - space_for_qmks_after
        cluevec_before = cluevec[1:first_qmk-1]
        cluevec_middle = [1]
        cluevec_after = cluevec[first_qmk+1:end]
        while (cluevec_middle[1] <= max_len)
            counter += count_states(cellvec, vcat(cluevec_before, cluevec_middle, cluevec_after))
            cluevec_middle[1] += 1  # a question mark represents a cell run of length 1 or more
        end
    elseif (length(all_asks) > 0)
        all(ismissing.(cellvec)) && return VectorStateCounter(cellvec; init="askmissing")
        all(cellvec) && return VectorStateCounter(cellvec; init="full")
        all(.!cellvec) && return VectorStateCounter(cellvec; init="empty")
        first_ask = all_asks[1]
        cluevec_before = cluevec[1:first_ask-1]
        cluevec_middle = QuestionMark[]
        cluevec_after = cluevec[first_ask+1:end]
        while ((2 * length(cluevec_middle) - 1) <= length(cellvec))
            counter += count_states(cellvec, vcat(cluevec_before, cluevec_middle, cluevec_after))
            push!(cluevec_middle, QuestionMark())  # an asterisk represents a sequence of question marks of length 0 or more
        end
    end
    (n(counter) === 0) ? nothing : counter
end

ij2rc(ij::CartesianIndex{2}) = [('R', ij[1]), ('C', ij[2])]
ijs2rcs(ijs::Vector{CartesianIndex{2}}) = vcat(ij2rc.(ijs)...)

@views function recount!(counter::VectorStateCounter, cellvec::T, cluevec::S) where {T <: OneDAbstractCellArray,  S <: AbstractClueVector}
    new_counter = count_states(cellvec, cluevec)
    isnothing(new_counter) && return nothing  # error("No valid states found during recount")
    counter.cumul[:] .= new_counter.cumul
    counter.n = new_counter.n
    true
end
@views function recount!(counter::MatrixStateCounter, cmat::T, puzzle::Puzzle, rc::Tuple{Char, Int}) where T <: TwoDAbstractCellArray
    if (rc[1] === 'R')
        i = rc[2]
        s = recount!(rows(counter)[i], cmat[i,:], rows(puzzle)[i])
    elseif (rc[1] === 'C')
        j = rc[2]
        s = recount!(cols(counter)[j], cmat[:,j], cols(puzzle)[j])
    else
        error("rowcol must either refer to row ('R') or col ('C')")
    end
    isnothing(s) ? nothing : true
end
@views function recount!(counter::MatrixStateCounter, cmat::T, puzzle::Puzzle, rcs::Vector{Tuple{Char, Int}}) where T <: TwoDAbstractCellArray
    for rc=rcs
        s = recount!(counter, cmat, puzzle, rc)
        isnothing(s) && return nothing
    end
    true
end
@views function recount!(counter::MatrixStateCounter, cmat::T, puzzle::Puzzle, ij::CartesianIndex{2}) where T <: TwoDAbstractCellArray
    i, j = Tuple(ij)
    s1 = recount!(rows(counter)[i], cmat[i,:], rows(puzzle)[i])
    isnothing(s1) && return nothing
    s2 = recount!(cols(counter)[j], cmat[:,j], cols(puzzle)[j])
    isnothing(s2) ? nothing : true
end
@views function recount!(counter::MatrixStateCounter, cmat::T, puzzle::Puzzle, ijs::Vector{CartesianIndex{2}}) where T <: TwoDAbstractCellArray
    recount!(counter, cmat, puzzle, ijs2rcs(ijs))
end

@views function update!(cmat::T, counter::MatrixStateCounter) where T <: TwoDCellArray
    this_odds = odds(counter)
    missing_cmat = ismissing.(cmat)
    true_updates = isone.(this_odds) .& missing_cmat
    false_updates = iszero.(this_odds) .& missing_cmat
    cmat[true_updates] .= true
    cmat[false_updates] .= false
    findall(true_updates .| false_updates)  # returns ijs where updates occurred this round
end

issmat(cmat::T) where T <: TwoDAbstractCellArray = all(@. !ismissing(cmat))

function iscontiguous(cmat::T) where T <: TwoDSolutionCellArray
    start_idx = findfirst(istrue.(cmat))
    # TODO: flood fill to determine whether all trues are contiguous
end

function check_cmat(puzzle::Puzzle, cmat::T, counter::MatrixStateCounter) where T <: TwoDAbstractCellArray
    issmat(cmat) && return issolution(puzzle, convert(SolutionCellMatrix, cmat), counter)
    false
end

function issolution(puzzle::Puzzle, cmat::T, counter::MatrixStateCounter) where T <: TwoDSolutionCellArray
    !iscontiguous(cmat) && return false
    rcs = vcat(map(r -> ('R', r), findall(map(n, rows(counter)) .> 1)),
               map(c -> ('C', c), findall(map(n, cols(counter)) .> 1)))
    s = recount!(counter, cmat, puzzle, rcs)
    !isnothing(s) && (num_states(counter) === 1)
end

function solve(puzzle::Puzzle, cmat::T, counter::MatrixStateCounter, ijs::Vector{CartesianIndex{2}}) where T <: TwoDCellArray
    #TODO: add guardrails using iscrowded()
    rcs = ijs2rcs(ijs)
    sort!(rcs, counter)
    history = ijs
    while length(rcs) > 0
        rc = popfirst!(rcs)
        s = recount!(counter, cmat, puzzle, rc)
        if isnothing(s)
            cmat[history] .= missing
            recount!(counter, cmat, puzzle, history)
            println("up")
            return cmat
        end
        new_ijs = update!(cmat, counter)
        check_cmat(puzzle, cmat, counter) && return convert(SolutionCellMatrix, cmat)
        if (length(new_ijs) > 0)
            append!(history, new_ijs)
            unique!(history)
            new_rcs = ijs2rcs(new_ijs)
            append!(rcs, new_rcs)
            unique!(rcs)
            sort!(rcs, counter)
        end
    end
    # once we're out of deterministic updates, we make a best guess to continue and backtrack if it doesn't work out
    this_odds = odds(counter)
    uncertainty = entropy.(this_odds)
    @. uncertainty[iszero(uncertainty)] = 1.0
    while !all(isone.(uncertainty))
        ij = findmin(uncertainty)[2]
        println(ij)
        initial_guess = (this_odds[ij] >= 0.5)
        cmat[ij] = initial_guess
        cmat = solve(puzzle, cmat, counter, [ij])
        if ismissing(cmat[ij])
            cmat[ij] = !initial_guess
            cmat = solve(puzzle, cmat, counter, [ij])
        end
        uncertainty[ij] = 1.0  # "remove" it from the "list"
        ismissing(cmat[ij]) && continue
        check_cmat(puzzle, cmat, counter) && return convert(SolutionCellMatrix, cmat)
    end
    println("up")
    cmat  # if we've gone through the whole list and still can't find a good result, then we need to backtrack
end
function solve(puzzle::Puzzle)
    cmat = init_cmat(puzzle)
    counter = MatrixStateCounter(puzzle)
    ijs = update!(cmat, counter)
    solve(puzzle, cmat, counter, ijs)
end
