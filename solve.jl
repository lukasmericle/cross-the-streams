export solve

@views function MatrixStateCounter(pzl::Puzzle)  # initializes based on puzzle description assuming grid is blank
    dummy_row = init_cvec(size(pzl, 2))
    dummy_col = init_cvec(size(pzl, 1))
    MatrixStateCounter(map(i -> count_states(dummy_row, rows(pzl)[i]), 1:size(pzl, 1)),
                       map(j -> count_states(dummy_col, cols(pzl)[j]), 1:size(pzl, 2)))
end
@views function MatrixStateCounter(cmat::T, pzl::Puzzle) where T <: TwoDCellArray  # initializes based on puzzle description and current state of grid
    MatrixStateCounter(map(i -> count_states(cmat[i,:], rows(pzl)[i]), 1:size(cmat, 1)),
                       map(i -> count_states(cmat[:,j], cols(pzl)[j]), 1:size(cmat, 2)))
end

init_cmat(pzl::Puzzle) = init_cmat(size(pzl)...)

complexity(pzl::Puzzle) = complexity(MatrixStateCounter(pzl))

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

@views function count_states(cellvec::T, cluevec::S) where {T <: OneDCellArray,  S <: AbstractClueVector}
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

@views function recount!(counter::VectorStateCounter, cellvec::T, cluevec::S) where {T <: OneDCellArray,  S <: AbstractClueVector}
    new_counter = count_states(cellvec, cluevec)
    isnothing(new_counter) && error("No valid states found during recount")
    counter.cumul[:] .= new_counter.cumul
    counter.n = new_counter.n
end
@views function recount!(counter::MatrixStateCounter, cmat::T, pzl::Puzzle, rc::Tuple{Char, Int}) where {T <: TwoDCellArray}
    if (rc[1] === 'R')
        i = rc[2]
        recount!(rows(counter)[i], cmat[i,:], rows(pzl)[i])
    elseif (rc[1] === 'C')
        j = rc[2]
        recount!(cols(counter)[j], cmat[:,j], cols(pzl)[j])
    else
        error("rowcol must either refer to row ('R') or col ('C')")
    end
end

@views function update!(cmat::T, counter::MatrixStateCounter) where T <: TwoDCellArray
    cmat_copy = copy(cmat)
    this_odds = odds(counter)
    @. cmat[isone(this_odds)] = true
    @. cmat[iszero(this_odds)] = false
    updates = (cmat_copy .!== cmat)  # true where a missing flipped to true or false
    filter!(rc -> (((rc[1] === 'R') && (n(rows(counter)[rc[2]]) > 1))
                || ((rc[1] === 'C') && (n(cols(counter)[rc[2]]) > 1))),
            vcat(map(r -> ('R', r), findall(any.(eachrow(updates)))),  # if any updates occur in a row, that row is obsolete
                 map(c -> ('C', c), findall(any.(eachcol(updates))))))
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

function solve(pzl::Puzzle, cmat::T, counter::MatrixStateCounter) where T <: TwoDCellArray
    obsolete_rcs = update!(cmat, counter)
    while (num_states(counter) > 1)
        if (length(obsolete_rcs) > 0)
            rc = popfirst!(obsolete_rcs)
            recount!(counter, cmat, pzl, rc)
            new_obsolete_rcs = update!(cmat, counter)
            if (length(new_obsolete_rcs) > 0)
                append!(obsolete_rcs, new_obsolete_rcs)
                unique!(obsolete_rcs)
                sort!(obsolete_rcs, counter)
            end
        else
            # if we're out of ideas, we make a best guess and backtrack if it didn't work
            # TODO: implement backtracking exploration
            break
        end
    end
    convert(SolutionCellMatrix, cmat)
end
solve(pzl::Puzzle) = solve(pzl, init_cmat(pzl), MatrixStateCounter(pzl))
