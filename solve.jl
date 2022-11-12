export solve, soldiff, verify


# Sorting heuristic for solution optimization.
@views function Base.sort!(rcs::Vector{Tuple{Char,Int}}, counter::MatrixStateCounter)

    """
    Sort the rcs so that the easiest rows and cols
    (with the fewest possible states to re-check) are first in the list.
    """

    ns = []
    for rc = rcs
        if rc[1] === 'R'
            push!(ns, n(rows(counter)[rc[2]]))
        elseif rc[1] === 'C'
            push!(ns, n(cols(counter)[rc[2]]))
        else
            error("rowcol must either refer to row ('R') or col ('C')")
        end
    end

    # get the order that the rcs
    # should be in from fewest to most states,
    # to act quickly on good information.
    # overwrite this array.
    rcs[:] .= rcs[sortperm(ns)]

end


# Helper methods to get more information about valid solutions.
function minreq_cells(cluevec_ints::T, n_qmks::Int) where {T<:AbstractClueVector}

    """
    Compute the minimum number of cells
    required for a set of clues for one vector.
    """

    n_reqd_cells = 0

    # obviously we need as many cells as there are explicit black cells,
    # plus the gaps between them.
    if (length(cluevec_ints) > 0)
        n_reqd_cells += sum(cluevec_ints) + (length(cluevec_ints) - 1)  # sum + (num - 1)
    end

    # we also need at least one black cell per question mark,
    # and at least one gap between them.
    if (n_qmks > 0)
        (n_reqd_cells > 0) && (n_reqd_cells += 1)  # add one extra for a gap between seq of ints and seq of qmks
        n_reqd_cells += 2 * n_qmks - 1
    end

    # Note: number and locations of asterisks are completely unimportant
    # as we are only interested in the minimum number of cells required.

    # return
    n_reqd_cells

end

function space_for_ints(cluevec_ints_after::T) where {T<:AbstractClueVector}
    """Count space needed for ints in puzzle clue."""
    (length(cluevec_ints_after) === 0) && return 0
    sum(cluevec_ints_after) + length(cluevec_ints_after)
end

function space_for_qmks(all_qmks::Vector{Int}, i::Int)
    """Count space needed for question marks in puzzle clue."""
    num_qmks_before = count(x -> (x < i), all_qmks)
    2 * num_qmks_before, 2 * (length(all_qmks) - num_qmks_before)  # num before, num after
end


# Workhorse solver function for counting states.
@views function count_states(cellvec::T, cluevec::S) where {T<:OneDAbstractCellArray,S<:AbstractClueVector}

    """
    Uses dynamic programming and branch-and-bound to test and accumulate
    all possible configurations of cells based on current col state and col description.
    """

    # shortcut for leaf case where there are no more clues.
    if (length(cluevec) === 0)
        # if any trues but no clues, this branch is invalid.
        any(skipmissing(cellvec)) && return nothing
        # if no trues and no clues, these states are false.
        return VectorStateCounter(cellvec; init="empty")
    end

    # collect indices of ints and question marks in this clue fragment.
    all_ints = findall(isa.(cluevec, Int))
    all_qmks = findall(isa.(cluevec, QuestionMark))
    all_asks = findall(isa.(cluevec, Asterisk))

    # make sure we have the room to check multiple states.
    # shortcut for leaf case where cellvec may not fit clues.
    n_reqd_cells = minreq_cells(cluevec[all_ints], length(all_qmks))
    # if the clues won't fit into the grid, this branch is invalid.
    (n_reqd_cells > length(cellvec)) && return nothing
    # if cellvec has length = 0 and cluevec has length > 0,
    # then cluevec contains only Asterisk.
    # no-op.
    (length(cellvec) === 0) && return VectorStateCounter(0; init="askmissing")

    # passed the checks: begin solving sub-problems with recursion.
    counter = VectorStateCounter(cellvec; init="blank")

    if (length(all_ints) > 0)

        # since ints are the most explicit (and restrictive),
        # we test valid locations for those first.

        first_int = all_ints[1]
        run_length = cluevec[first_int]

        # calculate quantities for finding first and last position for run in col.
        last_pos_apriori = length(cellvec) - (run_length - 1)  # minus 1 because ranges are inclusive in Julia
        space_for_ints_after = space_for_ints(cluevec[all_ints[2:end]])
        space_before, space_for_qmks_after = space_for_qmks(all_qmks, first_int)
        space_after = space_for_ints_after + space_for_qmks_after

        # loop through all positions where this run fits.
        # if first_pos > last_pos,
        # this for loop is skipped
        # and counter is returned blank
        # (thus this branch is invalid).
        for pos = (1+space_before):(last_pos_apriori-space_after)

            # indices of the run's left and right edges
            l, r = pos, pos + run_length - 1

            # cell vector representing the run itself
            cellvec_middle = cellvec[l:r]
            # skip this position if any cell in the run is not allowed to be true.
            any(isfalse, cellvec_middle) && continue
            # else, count the run as valid.
            counter_middle = VectorStateCounter(cellvec_middle; init="full")

            # repeat, but for the cells to the left.
            if (pos === 1)
                counter_before = VectorStateCounter(0; init="empty")
            else
                # continue if mandatory empty cell is already true
                istrue(cellvec[l-1]) && continue
                # else, recurse: count states of everything to the left.
                counter_before = count_states(cellvec[1:l-2], cluevec[1:first_int-1])
                # if left requires an invalid state,
                # this run position is invalid and doesn't count.
                isnothing(counter_before) && continue
                # else, keep those counted states and record the mandatory empty cell.
                push!(counter_before, false)
            end


            if (pos === last_pos_apriori)
                counter_after = VectorStateCounter(0; init="empty")
            else
                # continue if mandatory empty cell should be filled
                istrue(cellvec[r+1]) && continue
                # else, recurse: count states of everything to the right.
                counter_after = count_states(cellvec[r+2:end], cluevec[first_int+1:end])
                # if right requires an invalid state,
                # this run position is invalid and doesn't count.
                isnothing(counter_after) && continue
                # else, keep those counted states and record the mandatory empty cell.
                pushfirst!(counter_after, false)
            end

            # collect all possible counts for cells into one counter.
            counter += vcat(counter_before, counter_middle, counter_after)

        end

    elseif (length(all_qmks) > 0)

        # if there's no ints, we test for valid states with the question marks.

        # if all cells should be true
        # and there is only one question mark,
        # that question mark perfectly explains the vector.
        # if there is more than one question mark,
        # the state is invalid.
        if all(istrue, cellvec)
            return (length(all_qmks) === 1) ? VectorStateCounter(cellvec; init="full") : nothing
        end

        # if all cells should be false
        # and we require at least one question mark,
        # the state is invalid.
        all(isfalse, cellvec) && return nothing

        first_qmk = all_qmks[1]

        # manipulate the clue vector to simulate an int using the first question mark.
        min_space_after = 2 * (length(all_qmks) - 1)
        max_qmk_len = length(cellvec) - min_space_after
        cluevec_before = cluevec[1:first_qmk-1]
        cluevec_middle = [1] # a question mark represents a cell run of length 1 or more
        cluevec_after = cluevec[first_qmk+1:end]

        # pretend the first question mark takes on a series of explicit lengths,
        # and add valid states from recursive branches.
        while (cluevec_middle[1] <= max_qmk_len)
            counter += count_states(cellvec, vcat(cluevec_before, cluevec_middle, cluevec_after))
            cluevec_middle[1] += 1
        end

    elseif (length(all_asks) > 0)

        # finally, explore all possible configurations for asterisks.

        # if this is a simple state,
        # return a shortcut VectorStateCounter.
        all(ismissing, cellvec) && return VectorStateCounter(cellvec; init="askmissing")
        all(istrue, cellvec) && return VectorStateCounter(cellvec; init="full")
        all(isfalse, cellvec) && return VectorStateCounter(cellvec; init="empty")

        # manipulate clue vectors as before.
        first_ask = all_asks[1]
        cluevec_before = cluevec[1:first_ask-1]
        cluevec_middle = QuestionMark[] # an asterisk represents a sequence of question marks of length 0 or more
        cluevec_after = cluevec[first_ask+1:end]

        # pretend the asterisk is some list of question marks,
        # and add successively more question marks find all feasible states.
        while ((2 * length(cluevec_middle) - 1) <= length(cellvec))
            counter += count_states(cellvec, vcat(cluevec_before, cluevec_middle, cluevec_after))
            push!(cluevec_middle, QuestionMark())
        end

    end

    # if our counter found any valid states, pass them up.
    # else, return `nothing`, which indicates invalid state.
    (n(counter) > 0) ? counter : nothing

end


# Recount states from scratch.
@views function recount!(counter::VectorStateCounter, cellvec::T, cluevec::S) where {T<:OneDAbstractCellArray,S<:AbstractClueVector}
    """Recount all feasible states for the given row from scratch."""
    new_counter = count_states(cellvec, cluevec)
    isnothing(new_counter) && return nothing
    counter.cumul[:] .= new_counter.cumul
    counter.n = new_counter.n
    true
end

@views function recount!(counter::MatrixStateCounter, cmat::T, puzzle::P, rc::Tuple{Char,Int}) where {P<:Puzzle,T<:TwoDAbstractCellArray}
    """Recount all feasible states for the requested row or column from scratch."""
    if (rc[1] === 'R')
        i = rc[2]
        s = recount!(rows(counter)[i], cmat[i, :], rows(puzzle)[i])
    elseif (rc[1] === 'C')
        j = rc[2]
        s = recount!(cols(counter)[j], cmat[:, j], cols(puzzle)[j])
    else
        error("rowcol must either refer to row ('R') or col ('C')")
    end
    isnothing(s) ? nothing : true
end

@views function recount!(counter::MatrixStateCounter, cmat::T, puzzle::P, rcs::Vector{Tuple{Char,Int}}) where {P<:Puzzle,T<:TwoDAbstractCellArray}
    """Recount all feasible states for the requested rows and columns from scratch."""
    # sort by ease of re-checking
    sort!(rcs, counter)
    for rc = rcs
        s = recount!(counter, cmat, puzzle, rc)
        isnothing(s) && return nothing
    end
    true
end

@views function recount!(counter::MatrixStateCounter, cmat::T, puzzle::P, ij::TwoDCoord) where {P<:Puzzle,T<:TwoDAbstractCellArray}
    """Recount all feasible states for the row and column corresponding to some index."""
    i, j = Tuple(ij)
    (
        isnothing(recount!(rows(counter)[i], cmat[i, :], rows(puzzle)[i]))
        ||
        isnothing(recount!(cols(counter)[j], cmat[:, j], cols(puzzle)[j]))
    ) && return nothing
    true
end

@views function recount!(counter::MatrixStateCounter, cmat::T, puzzle::P, ijs::Vector{TwoDCoord}) where {P<:Puzzle,T<:TwoDAbstractCellArray}
    recount!(counter, cmat, puzzle, ijs2rcs(ijs))
end


# Grid state manipulation.
@views function update!(cmat::T, counter::MatrixStateCounter, apply_crowding::Bool) where {T<:TwoDCellArray}

    """
    Update cell matrix with guaranteed-true and guaranteed-false locations,
    then return locations of the newest changes.
    """

    this_odds = odds(counter)
    missing_cmat = ismissing.(cmat)

    # change missing -> true based on counter
    true_updates = missing_cmat .& isone.(this_odds)
    cmat[true_updates] .= true

    # change missing -> false based on counter (and optional crowding)
    false_locations = iszero.(this_odds)
    if (apply_crowding)
        false_locations .|= iscrowded(cmat)
    end
    false_updates = missing_cmat .& false_locations
    cmat[false_updates] .= false

    # returns ijs where updates occurred this round.
    vcat(findall(true_updates), findall(false_updates))

end

@views function undo!(cmat::T, counter::MatrixStateCounter, puzzle::Puzzle, ijs::Vector{TwoDCoord}) where {T<:TwoDCellArray}

    """Revert changes made to grid cell locations `ijs`."""

    unique!(ijs)

    # revert changes made during this loop.
    cmat[ijs] .= missing

    # reset counter for affected rows and columns.
    recount!(counter, cmat, puzzle, ijs2rcs(ijs))

    # return original grid state.
    cmat

end


# Check solution validity.
issmat(cmat::T) where {T<:TwoDAbstractCellArray} = !any(ismissing, cmat)

@views function iscontiguous(cmat::T) where {T<:TwoDSolutionCellArray}

    """Minimal flood-fill algorithm to determine solution contiguity."""

    start_idx = findfirst(istrue.(cmat))

    # all are false or missing: no connected components.
    isnothing(start_idx) && return false

    # visit ijs in flood-fill fashion.
    ijs = [start_idx]
    visited = TwoDCoord[]
    while (length(ijs) > 0)
        ij = popfirst!(ijs)
        push!(visited, ij)
        next_ijs = filled_neighbor_sites(cmat, ij)
        if (length(next_ijs) > 0)
            filter!(ij -> !in(ij, visited), next_ijs)
            append!(ijs, next_ijs)
            unique!(ijs)
        end
    end

    # true if there is one connected component.
    (count(cmat) === length(visited))

end

@views function isvalidsolution(puzzle::P, cmat::T) where {P<:Puzzle,T<:TwoDSolutionCellArray}
    isvalidsolution(puzzle, cmat, MatrixStateCounter(cmat, puzzle))
end

@views function isvalidsolution(puzzle::P, cmat::T, counter::MatrixStateCounter) where {P<:Puzzle,T<:TwoDSolutionCellArray}
    """
    Check if grid represents solution using existing counter.
    Does a final recount before checking that counter represents solution state.
    """
    ((puzzle isa CrossTheStreamsPuzzle) && !iscontiguous(cmat)) && return false
    rcs = vcat(map(r -> ('R', r), findall(n.(rows(counter)) .> 1)),
        map(c -> ('C', c), findall(n.(cols(counter)) .> 1)))
    s = recount!(counter, cmat, puzzle, rcs)
    !isnothing(s) && (n(counter) == 1)
end


# Solve.
@views function solve(puzzle::P, cmat::T, counter::MatrixStateCounter, ijs::Vector{TwoDCoord}) where {P<:Puzzle,T<:TwoDCellArray}

    """
    Main solve function.
    """

    # While we have unchecked modified states,
    # iteratively assign guaranteed-correct values to the grid state and update counter

    history = ijs

    while (length(ijs) > 0)

        # recount feasible states
        s = recount!(counter, cmat, puzzle, ijs)

        # if current grid state is actually infeasible, undo changes
        isnothing(s) && return undo!(cmat, counter, puzzle, history)

        # update grid state based on newly counted states,
        # getting new list of updated grid cell locations.
        ijs = update!(cmat, counter, puzzle isa CrossTheStreamsPuzzle)

        # record changes in case we need to undo.
        append!(history, ijs)

    end

    # return full solution if all cells are set
    issmat(cmat) && return convert(SolutionCellMatrix, cmat)

    # Once we're out of deterministic updates,
    # we make a best guess to continue down the tree,
    # and backtrack if it doesn't work out.
    # The below implementation is
    # more "depth-first" (exploitation)
    # than "breadth-first" (exploration).
    this_odds = odds(counter)
    # uncertainty is defined as the entropy of the probability of a cell being filled.
    uncertainty = entropy.(this_odds)
    # set zero-entropy (true or false, ie, not missing) locations to ignore them.
    @. uncertainty[iszero(uncertainty)] = 1.0

    # pick the next-surest position
    # and fill it according to what it is most likely to be.
    _, ij = findmin(uncertainty)

    if ((puzzle isa CrossTheStreamsPuzzle) && iscrowded(cmat, ij))

        # enforce constraint on crowding
        cmat = solve(puzzle, cmat, counter, ij, false)
        # if enforcement fails, this is an invalid state
        ismissing(cmat[ij]) && return undo!(cmat, counter, puzzle, history)

    else

        # make an initial guess based on the probability.
        initial_guess = (this_odds[ij] >= 0.5)

        # try solving with initial guess.
        # calling solve is a no-op when the guess doesn't pan out.
        cmat = solve(puzzle, cmat, counter, ij, initial_guess)

        # if no-op indeed happened...
        if ismissing(cmat[ij])

            # try the guess the other way.
            cmat = solve(puzzle, cmat, counter, ij, !initial_guess)

            # if cmat[ij] is still missing, neither guess worked!
            # this means we are in an invalid state,
            # because it is a contradiction to say
            # that the cell value is neither true nor false.
            # undo the changes from this whole function call.
            ismissing(cmat[ij]) && return undo!(cmat, counter, puzzle, history)

        end

    end

    # return full solution if all cells are set
    issmat(cmat) && return convert(SolutionCellMatrix, cmat)

    # all acceptable exit conditions are already specified.
    # if we got here, something is wrong.
    error("should never reach the end of this function")

end

@views function solve(puzzle::P, cmat::T, counter::MatrixStateCounter, ij::TwoDCoord, guess::Bool) where {P<:Puzzle,T<:TwoDCellArray}
    """
    Solve puzzle using grid and existing counter,
    after applying a specified guess to the grid state at a specified index.
    """
    !ismissing(cmat[ij]) && error("guessing on a cell that is already assigned")
    cmat[ij] = guess
    solve(puzzle, cmat, counter, [ij])
end

@views function solve(puzzle::P, cmat::T, counter::MatrixStateCounter) where {P<:Puzzle,T<:TwoDCellArray}
    """Solve puzzle using grid and existing counter."""
    ijs = update!(cmat, counter, puzzle isa CrossTheStreamsPuzzle)
    solve(puzzle, cmat, counter, ijs)
end

@views function solve(puzzle::P) where {P<:Puzzle}
    counter = MatrixStateCounter(puzzle)
    smat = solve(puzzle, init_cmat(puzzle), counter)
    # do full recount using puzzle and solution state
    !isvalidsolution(puzzle, smat) && error("solution is invalid")
    return smat
end

@views soldiff(psol::T, sol::T) where {T<:TwoDSolutionCellArray} = convert(SolutionCellMatrix, xor.(psol, sol))
@views verify(psol::T, sol::T) where {T<:TwoDSolutionCellArray} = !any(soldiff(psol, sol))
@views verify(psol::TP, sol::TS) where {TP<:TwoDCellArray,TS<:TwoDSolutionCellArray} = false
