export GridCol, Grid, SolutionGrid
export empty_neighbor_sites, iscrowded, inbounds
import Base: string
export       string

using Lazy: @forward

abstract type AbstractGridCol{T} <: AbstractVector{T} end
abstract type AbstractGrid{T} <: AbstractMatrix{T} end

GridCol = Array{CellState, 1}
SolutionGridCol = Array{Bool, 1}
OneDGrid = Union{GridCol, SolutionGridCol}

Grid = Array{CellState, 2}
SolutionGrid = Array{Bool, 2}
TwoDGrid = Union{Grid, SolutionGrid}

TwoDCoord = Tuple{Int,Int}

VON_NEUMANN_NEIGHBORHOOD_1D = [-1, 1]
VON_NEUMANN_NEIGHBORHOOD_2D = [(0,-1), (-1,0), (1,0), (0,1)]

inbounds(grid::T, i::Int, dim::Int) where T <: TwoDGrid = (1 <= i) && (i <= size(grid, dim))
inbounds(grid::T, ij::TwoDCoord) where T <: TwoDGrid = inbounds(grid, ij[1], 1) && inbounds(grid, ij[2], 2)

function empty_neighbor_sites(grid::T, ij::TwoDCoord) where T <: TwoDGrid
    filter(
        iijj -> inbounds(grid, iijj) && isempty(grid[iijj...]),
        broadcast(vnn -> ij .+ vnn, VON_NEUMANN_NEIGHBORHOOD_2D))
end

function iscrowded(grid::T, ij::TwoDCoord) where T <: TwoDGrid
    """
    Return true if any of the two-by-two's which contain (i,j)
    have three or more filled cells.
    """
    let (i,j) = ij
        for ii=(i-1):i
            !((1 <= ii) && (ii+1 <= size(grid, 1))) && continue
            for jj=(j-1):j
                !((1 <= jj) && (jj+1 <= size(grid, 2))) && continue
                (count(@view grid[ii:ii+1,jj:jj+1]) >= 3) && return true
            end
        end
    end
    false
end

function Base.string(grid::T) where T <: OneDGrid
    prod(map(x -> (ismissing(x) ? "><" : (x ? "██" : "  ")),  grid))
end
function Base.string(grid::T) where T <: TwoDGrid
    (n, m) = size(grid)
    s = "┏" * "━━"^m * "┓"
    for i=1:n
        s *= "\n" * "┃" * Base.string(grid[i,:]) * "┃"
    end
    s *= "\n" * "┗" * "━━"^m * "┛"
    s
end
