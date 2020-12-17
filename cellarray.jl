export CellState, CellVector, SolutionCellVector, CellMatrix, SolutionCellMatrix
export init_cmat, init_cvec
export empty_neighbor_sites
import Base: string
export       string

using Lazy: @forward

CellState = Union{Bool, Missing}

istrue(x::CellState) = !ismissing(x) && x
isfalse(x::CellState) = !ismissing(x) && !x
isfilled(x::CellState) = istrue(x)
isempty(x::CellState) = !isfilled(x)

CellVector = Vector{CellState}
SolutionCellVector = Vector{Bool}

CellMatrix = Matrix{CellState}
SolutionCellMatrix = Matrix{Bool}

CellVectorVectorView = SubArray{CellState, 1, CellVector}
CellMatrixVectorView = SubArray{CellState, 1, CellMatrix}
CellMatrixMatrixView = SubArray{CellState, 2, CellMatrix}

SolutionCellVectorVectorView = SubArray{Bool, 1, SolutionCellVector}
SolutionCellMatrixVectorView = SubArray{Bool, 1, SolutionCellMatrix}
SolutionCellMatrixMatrixView = SubArray{Bool, 2, SolutionCellMatrix}

OneDCellArray = Union{CellVector, CellVectorVectorView, CellMatrixVectorView}
TwoDCellArray = Union{CellMatrix, CellMatrixMatrixView}

OneDSolutionCellArray = Union{SolutionCellVector, SolutionCellVectorVectorView, SolutionCellMatrixVectorView}
TwoDSolutionCellArray = Union{SolutionCellMatrix, SolutionCellMatrixMatrixView}

OneDAbstractCellArray = Union{OneDCellArray, OneDSolutionCellArray}
TwoDAbstractCellArray = Union{TwoDCellArray, TwoDSolutionCellArray}

TwoDCoord = Tuple{Int,Int}

VON_NEUMANN_NEIGHBORHOOD_1D = [-1, 1]
VON_NEUMANN_NEIGHBORHOOD_2D = [(0,-1), (-1,0), (1,0), (0,1)]

function init_cmat(n::Int, m::Int)
    cmat = CellMatrix(undef, n, m)
    fill!(cmat, missing)
    cmat
end

function init_cvec(n::Int)
    cvec = CellVector(undef, n)
    fill!(cvec, missing)
    cvec
end

string_vec(cvec::T) where T <: OneDAbstractCellArray = prod(map(x -> (ismissing(x) ? "><" : (x ? "██" : "  ")),  cvec))
Base.string(cvec::T) where T <: OneDCellArray = string_vec(cvec)  # to avoid ambiguity, overload with custom method
Base.string(cvec::T) where T <: OneDSolutionCellArray = string_vec(cvec)

function string_mat(cmat::T) where T <: TwoDAbstractCellArray
    (n, m) = size(cmat)
    s = "┏" * "━━"^m * "┓"
    for i=1:n
        s *= "\n" * "┃" * Base.string(cmat[i,:]) * "┃"
    end
    s *= "\n" * "┗" * "━━"^m * "┛"
    s
end
Base.string(cmat::SolutionCellMatrix) = string_mat(cmat)  # to avoid ambiguity, overload with custom method
Base.string(cmat::T) where T <: TwoDAbstractCellArray = string_mat(cmat)

inbounds(cmat::T, i::Int, dim::Int) where T <: TwoDAbstractCellArray = (1 <= i) && (i <= size(cmat, dim))
inbounds(cmat::T, ij::TwoDCoord) where T <: TwoDAbstractCellArray = inbounds(cmat, ij[1], 1) && inbounds(cmat, ij[2], 2)

@views function iscrowded(cmat::T, ij::TwoDCoord) where T <: TwoDAbstractCellArray
    """
    Return true if any of the two-by-two's which contain (i,j)
    have three or more filled cells.
    """
    let (i,j) = ij
        for ii=(i-1):i
            !((1 <= ii) && (ii+1 <= size(cmat, 1))) && continue
            for jj=(j-1):j
                !((1 <= jj) && (jj+1 <= size(cmat, 2))) && continue
                (count(cmat[ii:ii+1,jj:jj+1]) >= 3) && return true
            end
        end
    end
    false
end

@views function empty_neighbor_sites(cmat::T, ij::TwoDCoord) where T <: TwoDSolutionCellArray
    filter(
        iijj -> (inbounds(cmat, iijj) && isempty(cmat[iijj...])),
        broadcast(vnn -> ij .+ vnn, VON_NEUMANN_NEIGHBORHOOD_2D))
end
