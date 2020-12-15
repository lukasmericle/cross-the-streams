export solve, count_states, odds
import Base: vcat, size, foldl, length, iterate, setindex, setindex!, getindex, firstindex, lastindex
export       vcat, size, foldl, length, iterate, setindex, setindex!, getindex, firstindex, lastindex

using Lazy: @forward

abstract type CellRun end

struct EmptyCellRun
    width::Int

    function EmptyCellRun(width::Int)
        (width < 0) && error("Cell run cannot have negative width")
        new(width)
    end
end

struct FullCellRun
    width::Int

    function FullCellRun(width::Int)
        (width < 0) && error("Cell run cannot have negative width")
        new(width)
    end
end

width(run::T) where T <: CellRun = run.width
Base.length(run::T) where T <: CellRun = width(run)

T(gcol::GridCol) where T <: CellRun = T(length(gcol))

function init_grid(problem::Problem)
    grid = Grid(undef, size(problem)...)
    fill!(grid, missing)
    grid
end

function init_gcol(n::Int)
    gcol = GridCol(undef, n)
    fill!(gcol, missing)
    gcol
end

abstract type AbstractStateCounter end

mutable struct ColStateCounter <: AbstractStateCounter
    cumul::Array{Int, 1}
    n::Int
end

function ColStateCounter(m::Int; init::String="blank")
    if init == "blank"
        cumul = zeros(Int, m)
        n = 0
    elseif init == "empty"
        cumul = zeros(Int, m)
        n = 1
    elseif init == "full"
        cumul = ones(Int, m)
        n = 1
    elseif init == "askmissing"
        cumul = Array{Int}(undef, m)
        fill!(cumul, 2^(m-1))
        n = 2^m
    end
    ColStateCounter(cumul, n)
end

ColStateCounter() = ColStateCounter(0; init="blank")
ColStateCounter(grid::GridCol; init::String="blank") = ColStateCounter(length(grid); init=init)

cumul(counter::ColStateCounter) = counter.cumul
n(counter::ColStateCounter) = counter.n

@forward ColStateCounter.cumul size
@forward ColStateCounter.cumul length
@forward ColStateCounter.cumul iterate
@forward ColStateCounter.cumul setindex
@forward ColStateCounter.cumul setindex!
@forward ColStateCounter.cumul getindex
@forward ColStateCounter.cumul firstindex
@forward ColStateCounter.cumul lastindex

function Base.vcat(a::ColStateCounter, b::ColStateCounter)
    (length(a) == 0) && return b
    (length(b) == 0) && return a
    a.cumul = vcat(n(b) .* cumul(a), n(a) .* cumul(b))
    a.n *= n(b)
    a
end

Base.vcat(counters::Vararg{ColStateCounter}) = Base.foldl(Base.vcat, counters, EmptyCellRun(0))

function Base.vcat(a::ColStateCounter, b::FullCellRun)
    (length(a) == 0) && return b
    (length(b) == 0) && return a
    a.cumul = vcat(cumul(a), n(a) .* ones(length(b)))
    a
end

function Base.vcat(a::FullCellRun, b::ColStateCounter)
    (length(a) == 0) && return b
    (length(b) == 0) && return a
    b.cumul = vcat(n(b) .* ones(length(a)), cumul(b))
    b
end

function Base.vcat(a::ColStateCounter, b::EmptyCellRun)
    (length(a) == 0) && return b
    (length(b) == 0) && return a
    a.cumul = vcat(cumul(a), zeros(length(b)))
    a
end

function Base.vcat(a::EmptyCellRun, b::ColStateCounter)
    (length(a) == 0) && return b
    (length(b) == 0) && return a
    b.cumul = vcat(zeros(length(a)), cumul(b))
    b
end

function (Base.:+)(a::ColStateCounter, b::ColStateCounter)
    a.cumul .+= b.cumul
    a.n += b.n
    a
end

function (Base.:*)(a::Union{T, ColStateCounter}, b::Union{T, ColStateCounter}) where T <: CellRun
    Base.vcat(a, b)
end

odds(counter::ColStateCounter) = @. cumul(counter) / n(counter)

mutable struct GridStateCounter <: AbstractStateCounter
    rows::Array{ColStateCounter}
    cols::Array{ColStateCounter}
end

rows(counter::GridStateCounter) = counter.rows
cols(counter::GridStateCounter) = counter.cols

function GridStateCounter(problem::Problem)  # initializes based on problem description assuming grid is blank
    n, m = size(problem)
    dummy_row = init_gcol(m)
    dummy_col = init_gcol(n)
    GridStateCounter(
        map(i -> count_states(dummy_row, rows(problem)[i]), 1:n)
        map(j -> count_states(dummy_col, cols(problem)[j]), 1:m))  # TODO: redo with views
end

function GridStateCounter(grid::Grid, problem::Problem)  # initializes based on problem description and current state of grid
    GridStateCounter(
        map(i -> count_states(grid[i,:], rows(problem)[i]), 1:size(grid, 1))
        map(i -> count_states(grid[:,j], cols(problem)[j]), 1:size(grid, 2)))  # TODO: redo with views
end

function GridStateCounter(n::Int, m::Int; init::String="blank")
    GridStateCounter(
        [ColStateCounter(m; init=init) for _=1:n],
        [ColStateCounter(n; init=init) for _=1:m])
end

Base.size(counter::GridStateCounter, dim::Int) = (dim == 1) ? length(rows(counter)) : ((dim == 2) ? length(cols(counter)) : error("grid state counter has only two dimensions"))
Base.size(counter::GridStateCounter) = (Base.size(counter, 1), Base.size(counter, 2))

function odds_rowcol(x::T, y::T) where T <: AbstractFloat
    (isone(x) || isone(y)) && return 1.0
    (iszero(x) || iszero(y)) && return 0.0
    x * y
end

function odds(counter::GridStateCounter)
    n, m = size(counter)
    row_odds = Array{Float64}(undef, n, m)
    col_out_odds = Array{Float64}(undef, n, m)
    for i=1:n
        row_odds[i,:] .= rows(counter)[i] |> odds
    end
    for j=1:m
        col_out_odds[:,j] .= cols(counter)[j] |> odds
    end
    for (i,j)=eachindex(col_out_odds)  # TODO: vectorize
        col_out_odds[i,j] = odds_rowcol(row_odds[i,j], col_out_odds[i,j])  # reuse one of the arrays as the output
    end
    col_out_odds
end

function minreq_cells(pcol_ints::ProblemCol, n_qmks::Int)
    n_reqd_cells = 0
    if length(pcol_ints) > 0
        n_reqd_cells += sum(pcol_ints) + (length(pcol_ints) - 1)  # sum + (num - 1)
    end
    if n_qmks > 0
        if n_reqd_cells > 0
            n_reqd_cells += 1  # add one extra for a gap between seq of ints and seq of qmks
        end
        n_reqd_cells += 2 * n_qmks - 1  # sum + (num - 1), assuming all ?s represent length-1 runs
    end
    n_reqd_cells
end

function space_for_ints(pcol_ints_after::Array{Int})
    (length(pcol_ints_after) == 0) && return 0
    sum(pcol_ints_after) + length(pcol_ints_after)
end

function space_for_qmks(all_qmks::Array{Int}, i::Int)
    num_qmks_before = count(x -> x < i, all_qmks)
    2 * num_qmks_before, 2 * (length(all_qmks) - num_qmks_before)  # num before, num after
end

function count_states(gcol::GridCol, pcol::ProblemCol)

    # TODO: come up with a way to bound branches without allocating new ColStateCounters
    # with FullCellRun, EmptyCellRun, etc. as placeholders for the results that are returned at the leaves
    # of the tree. Then only when we vcat a set of counters do we accumulate the appropriate values
    # based on the types returned from the leaves. this will greatly reduce the number of allocations.
    # for bound conditions, return NaN.

    """
    Uses dynamic programming and branch-and-bound to test and accumulate
    all possible configurations of cells based on current col state and col description.
    """

    if (length(pcol) == 0)
        if any(skipmissing(gcol))
            return ColStateCounter(gcol; init="blank")  # if any trues but no col description, bound this branch
        else
            return ColStateCounter(gcol; init="empty")  # if no trues, no-op
        end
    end

    all_ints = findall(isa.(pcol, Int))
    all_qmks = findall(isa.(pcol, QuestionMark))
    all_asks = findall(isa.(pcol, Asterisk))

    n_reqd_cells = minreq_cells(pcol[all_ints], length(all_qmks))

    if n_reqd_cells > length(gcol)
        return ColStateCounter(gcol; init="blank")  # if the col description won't fit into the grid, bound this branch
    elseif length(gcol) == 0
        return ColStateCounter(gcol; init="empty")  # reach this branch if (length(gcol) == 0) and col has only Asterisk in it; no-op
    end

    counter = ColStateCounter(gcol; init="blank")

    if length(all_ints) > 0

        first_int = all_ints[1]

        pcol_before = pcol[1:first_int-1]
        run_length = pcol[first_int]
        pcol_after = pcol[first_int+1:end]

        space_for_ints_after = space_for_ints(pcol[all_ints[2:end]])
        space_for_qmks_before, space_for_qmks_after = space_for_qmks(all_qmks, first_int)

        last_pos_apriori = length(gcol) - (run_length - 1)

        if length(all_ints) > 1
            ints_after = all_ints[2:end]
            sum_ints_after = sum(pcol[ints_after])  # count black cells
            num_ints_after = length(ints_after)  # count mandatory white cells
            space_for_ints_after = sum_ints_after + num_ints_after
        else
            space_for_ints_after = 0
        end

        first_pos = 1 + space_for_qmks_before
        last_pos = last_pos_apriori - space_for_ints_after - space_for_qmks_after

        # if first_pos > last_pos, this for loop is skipped and counter is returned blank (thus this branch is bounded)
        for pos=first_pos:last_pos

            gcol_middle = gcol[pos:pos+(run_length-1)]
            any(isfalse.(gcol_middle)) && continue  # continue if filled cell should be empty
            counter_middle = ColStateCounter(gcol_middle; init="full")

            if pos == 1
                counter_before = ColStateCounter(0; init="empty")
            else
                istrue(gcol[pos-1]) && continue  # continue if empty cell should be filled
                gcol_before = gcol[1:pos-2]
                counter_before = count_states(gcol_before, pcol_before) * ColStateCounter(1; init="empty")
                (n(counter_before) == 0) && continue
            end

            if pos == last_pos_apriori
                counter_after = ColStateCounter(0; init="empty")
            else
                istrue(gcol[pos+(run_length-1)+1]) && continue  # continue if empty cell should be filled
                gcol_after = gcol[pos+(run_length-1)+2:end]
                counter_after = ColStateCounter(1; init="empty") * count_states(gcol_after, pcol_after)
                (n(counter_after) == 0) && continue
            end

            counter += vcat(counter_before, counter_middle, counter_after)

        end

    elseif length(all_qmks) > 0

        first_qmk = all_qmks[1]

        space_for_qmks_after = (length(qmks) > 0) ? (2 * (length(all_qmks) - 1)) : 0
        max_len = length(gcol) - space_for_qmks_after

        pcol_before = pcol[1:first_qmk-1]
        pcol_middle = [1]
        pcol_after = pcol[first_qmk+1:end]

        while pcol_middle[1] <= max_len

            counter += count_states(gcol, vcat(pcol_before, pcol_middle, pcol_after))
            pcol_middle[1] += 1

        end

    elseif length(all_asks) > 0

        all(ismissing.(gcol)) && return ColStateCounter(gcol; init="askmissing")
        all(istrue.(gcol)) && return ColStateCounter(gcol; init="full")
        all(isfalse.(gcol)) && return ColStateCounter(gcol; init="empty")

        first_ask = all_asks[1]

        pcol_before = pcol[1:first_ask-1]
        pcol_middle = QuestionMark[]
        pcol_after = pcol[first_ask+1:end]

        while (2 * length(pcol_middle) - 1) <= length(gcol)

            counter += count_states(gcol, vcat(pcol_before, pcol_middle, pcol_after))
            push!(pcol_middle, QuestionMark())  # an asterisk represents a sequence of question marks of length 0 or more

        end

    end

    counter

end

function recount!(counter::ColStateCounter, gcol::GridCol, pcol::ProblemCol)
    counter = count_states(gcol, pcol)
end

function update!(grid::Grid, counter::GridStateCounter)
    this_odds = odds(counter)
    @. grid[isone(this_odds)] = true
    @. grid[iszero(this_odds)] = false
    # return list of rows and cols which are now obsolete (and must be re-counted) as a result of the changes
end

function solve(problem::Problem, grid::Grid, counter::GridStateCounter)
    obsolete_rows, obsolete_cols = Int[], Int[]
    while any(ismissing.(grid))
        # ...
        new_obsolete_rows, new_obsolete_cols = update!(grid, counter)
        # ...
    end
    grid
end

solve(problem::Problem) = convert(SolutionGrid, solve(problem, init_grid(problem), GridStateCounter(problem)))
