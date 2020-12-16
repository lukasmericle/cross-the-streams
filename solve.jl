export solve

function minreq_cells(cluevec_ints::T, n_qmks::Int) where T <: AbstractClueVector
    n_reqd_cells = 0
    if length(cluevec_ints) > 0
        n_reqd_cells += sum(cluevec_ints) + (length(cluevec_ints) - 1)  # sum + (num - 1)
    end
    if n_qmks > 0
        if n_reqd_cells > 0
            n_reqd_cells += 1  # add one extra for a gap between seq of ints and seq of qmks
        end
        n_reqd_cells += 2 * n_qmks - 1  # sum + (num - 1), assuming all ?s represent length-1 runs
    end
    n_reqd_cells
end

function space_for_ints(cluevec_ints_after::Array{Int})
    (length(cluevec_ints_after) == 0) && return 0
    sum(cluevec_ints_after) + length(cluevec_ints_after)
end

function space_for_qmks(all_qmks::Array{Int}, i::Int)
    num_qmks_before = count(x -> x < i, all_qmks)
    2 * num_qmks_before, 2 * (length(all_qmks) - num_qmks_before)  # num before, num after
end

@views function count_states(cellvec::T, cluevec::S) where {T <: OneDCellArray,  S <: AbstractClueVector}

    """
    Uses dynamic programming and branch-and-bound to test and accumulate
    all possible configurations of cells based on current col state and col description.
    """

    if (length(cluevec) == 0)
        if any(skipmissing(cellvec))
            return VectorStateCounter(cellvec; init="blank")  # if any trues but no col description, bound this branch
        else
            return VectorStateCounter(cellvec; init="empty")  # if no trues, no-op
        end
    end

    all_ints = findall(isa.(cluevec, Int))
    all_qmks = findall(isa.(cluevec, QuestionMark))
    all_asks = findall(isa.(cluevec, Asterisk))

    n_reqd_cells = minreq_cells(cluevec[all_ints], length(all_qmks))

    if (n_reqd_cells > length(cellvec))
        return VectorStateCounter(cellvec; init="blank")  # if the col description won't fit into the grid, bound this branch
    elseif (length(cellvec) === 0)
        return VectorStateCounter(cellvec; init="empty")  # reach this branch if (length(cellvec) == 0) and col has only Asterisk in it; no-op
    end

    counter = VectorStateCounter(cellvec; init="blank")

    if length(all_ints) > 0

        first_int = all_ints[1]

        cluevec_before = cluevec[1:first_int-1]
        run_length = cluevec[first_int]
        cluevec_after = cluevec[first_int+1:end]

        space_for_ints_after = space_for_ints(convert(Array{Int}, cluevec[all_ints[2:end]]))
        space_for_qmks_before, space_for_qmks_after = space_for_qmks(all_qmks, first_int)

        last_pos_apriori = length(cellvec) - (run_length - 1)

        first_pos = 1 + space_for_qmks_before
        last_pos = last_pos_apriori - space_for_ints_after - space_for_qmks_after

        # if first_pos > last_pos, this for loop is skipped and counter is returned blank (thus this branch is bounded)
        for pos=first_pos:last_pos

            cellvec_middle = cellvec[pos:pos+(run_length-1)]
            any(isfalse.(cellvec_middle)) && continue  # continue if filled cell should be empty
            counter_middle = VectorStateCounter(cellvec_middle; init="full")

            if pos == 1
                counter_before = VectorStateCounter(0; init="empty")
            else
                istrue(cellvec[pos-1]) && continue  # continue if empty cell should be filled
                cellvec_before = cellvec[1:pos-2]
                counter_before = count_states(cellvec_before, cluevec_before) * VectorStateCounter(1; init="empty")
                (n(counter_before) === 0) && continue
            end

            if pos == last_pos_apriori
                counter_after = VectorStateCounter(0; init="empty")
            else
                istrue(cellvec[pos+(run_length-1)+1]) && continue  # continue if empty cell should be filled
                cellvec_after = cellvec[pos+(run_length-1)+2:end]
                counter_after = VectorStateCounter(1; init="empty") * count_states(cellvec_after, cluevec_after)
                (n(counter_after) === 0) && continue
            end

            counter += vcat(counter_before,
                            counter_middle,
                            counter_after)

        end

    elseif length(all_qmks) > 0

        first_qmk = all_qmks[1]

        space_for_qmks_after = (length(all_qmks) > 0) ? (2 * (length(all_qmks) - 1)) : 0
        max_len = length(cellvec) - space_for_qmks_after

        cluevec_before = cluevec[1:first_qmk-1]
        cluevec_middle = [1]
        cluevec_after = cluevec[first_qmk+1:end]

        while cluevec_middle[1] <= max_len

            counter += count_states(cellvec, vcat(cluevec_before, cluevec_middle, cluevec_after))
            cluevec_middle[1] += 1

        end

    elseif length(all_asks) > 0

        all(ismissing.(cellvec)) && return VectorStateCounter(cellvec; init="askmissing")
        all(istrue.(cellvec)) && return VectorStateCounter(cellvec; init="full")
        all(isfalse.(cellvec)) && return VectorStateCounter(cellvec; init="empty")

        first_ask = all_asks[1]

        cluevec_before = cluevec[1:first_ask-1]
        cluevec_middle = QuestionMark[]
        cluevec_after = cluevec[first_ask+1:end]

        while (2 * length(cluevec_middle) - 1) <= length(cellvec)

            counter += count_states(cellvec, vcat(cluevec_before, cluevec_middle, cluevec_after))
            push!(cluevec_middle, QuestionMark())  # an asterisk represents a sequence of question marks of length 0 or more

        end

    end

    counter

end

function recount!(counter::VectorStateCounter, cellvec::T, cluevec::S) where {T <: OneDCellArray,  S <: AbstractClueVector}
    new_counter = count_states(cellvec, cluevec)
    counter.cumul[:] .= new_counter.cumul
    counter.n = new_counter.n
end

@views function recount!(counter::MatrixStateCounter, cmat::T, pzl::Puzzle, rowcol::Tuple{Char, Int}) where {T <: TwoDCellArray}
    if rowcol[1] === 'R'
        i = rowcol[2]
        recount!(rows(counter)[i], cmat[i,:], rows(pzl)[i])
    elseif rowcol[1] === 'C'
        j = rowcol[2]
        recount!(cols(counter)[j], cmat[:,j], cols(pzl)[j])
    else
        error("rowcol must either refer to row ('R') or col ('C')")
    end
end

function update!(cmat::T, counter::MatrixStateCounter) where T <: TwoDCellArray
    this_odds = odds(counter)
    cmat_copy = copy(cmat)
    @. cmat[isone(this_odds)] = true
    @. cmat[iszero(this_odds)] = false
    updates = (cmat_copy .=== cmat)  # true where a missing flipped to true or false
    obsolete_rows = findall(any.(eachrow(updates)))  # if any updates occur in a row, that row is obsolete
    obsolete_cols = findall(any.(eachcol(updates)))
    vcat(map(r -> ('R', r), obsolete_rows), map(c -> ('C', c), obsolete_cols))
end

function solve(pzl::Puzzle, cmat::T, counter::MatrixStateCounter) where T <: TwoDCellArray
    obsolete_rowcols = update!(cmat, counter)
    while any(ismissing.(cmat))
        rowcol = popfirst!(obsolete_rowcols)
        recount!(counter, cmat, pzl, rowcol)
        new_obsolete_rowcols = update!(cmat, counter)
        append!(obsolete_rowcols, new_obsolete_rowcols)
        filter(rc -> (((rc[1] === 'R') && any(ismissing.(@view cmat[rc[2],:])))
                   || ((rc[1] === 'C') && any(ismissing.(@view cmat[:,rc[2]])))),
               obsolete_rowcols)
        unique!(obsolete_rowcols)
    end
    convert(SolutionCellMatrix, cmat)
end

solve(pzl::Puzzle) = solve(pzl, init_cmat(pzl), MatrixStateCounter(pzl))
