export CellMatrix, SolutionCellMatrix
export init_cmat


# cell state is `true` for black, `false` for white, and `missing` for ambiguous.
# special semantic state `isempty(x)` === `isfalse(x)` || `ismissing(x)`.
CellState = Union{Bool,Missing}

# shorthand to check states of grid cells.
# typically used e.g. `all(istrue, cellvec)`.
istrue(x::CellState) = !ismissing(x) && x
isfalse(x::CellState) = !ismissing(x) && !x
isempty(x::CellState) = ismissing(x) || !x  # == !istrue(x)


# cell vectors and matrices hold cell states.
# since a solution is by definition unambiguous,
# solution cell vectors and matrices can only take booleans.
CellVector = Vector{CellState}
CellMatrix = Matrix{CellState}
SolutionCellVector = Vector{Bool}
SolutionCellMatrix = Matrix{Bool}

# initialization shorthand for cell vectors and matrices
init_cvec(n::Int) = convert(CellVector, fill(missing, n))
init_cmat(n::Int, m::Int) = convert(CellMatrix, fill(missing, n, m))


# views avoid copies to make the program more efficient.
CellVectorVectorView = SubArray{CellState,1,CellVector}
CellMatrixVectorView = SubArray{CellState,1,CellMatrix}
CellMatrixMatrixView = SubArray{CellState,2,CellMatrix}
SolutionCellVectorVectorView = SubArray{Bool,1,SolutionCellVector}
SolutionCellMatrixVectorView = SubArray{Bool,1,SolutionCellMatrix}
SolutionCellMatrixMatrixView = SubArray{Bool,2,SolutionCellMatrix}

# A general type placeholder for either a vector or a view in one or two dimensions.
OneDCellArray = Union{CellVector,CellVectorVectorView,CellMatrixVectorView}
TwoDCellArray = Union{CellMatrix,CellMatrixMatrixView}
OneDSolutionCellArray = Union{SolutionCellVector,SolutionCellVectorVectorView,SolutionCellMatrixVectorView}
TwoDSolutionCellArray = Union{SolutionCellMatrix,SolutionCellMatrixMatrixView}

# Kind of using "abstract" types backwards here.
OneDAbstractCellArray = Union{OneDCellArray,OneDSolutionCellArray}
TwoDAbstractCellArray = Union{TwoDCellArray,TwoDSolutionCellArray}
OneOrTwoDAbstractCellArray = Union{OneDAbstractCellArray,TwoDAbstractCellArray}

# String methods
string_vec(cvec::T) where {T<:OneDAbstractCellArray} = prod(map(x -> (ismissing(x) ? "><" : (x ? "██" : "  ")), cvec))
Base.string(cvec::T) where {T<:OneDCellArray} = string_vec(cvec)  # to avoid ambiguity during dispatch, overload with custom method
Base.string(cvec::T) where {T<:OneDSolutionCellArray} = string_vec(cvec)

function string_mat(cmat::T) where {T<:TwoDAbstractCellArray}
    (n, m) = size(cmat)
    s = "┏" * "━━"^m * "┓"
    for i = 1:n
        s *= "\n" * "┃" * string(cmat[i, :]) * "┃"
    end
    s *= "\n" * "┗" * "━━"^m * "┛"
    s
end
Base.string(cmat::SolutionCellMatrix) = string_mat(cmat)  # to avoid ambiguity during dispatch, overload with custom method
Base.string(cmat::T) where {T<:TwoDAbstractCellArray} = string_mat(cmat)


# Representing coordinates in the grid.
OneDCoord = CartesianIndex{1}
TwoDCoord = CartesianIndex{2}

ij2rc(ij::TwoDCoord) = [('R', ij[1]), ('C', ij[2])]
ijs2rcs(ijs::Vector{TwoDCoord}) = unique(vcat(ij2rc.(ijs)...))


# Constants and methods to help with neighborhood checking (solution continuity)
VON_NEUMANN_NEIGHBORHOOD_1D = [OneDCoord(-1),
    OneDCoord(1)]
VON_NEUMANN_NEIGHBORHOOD_2D = [TwoDCoord(0, -1),
    TwoDCoord(-1, 0),
    TwoDCoord(1, 0),
    TwoDCoord(0, 1)]
TWO_BY_TWO_NEIGHBORHOOD_2D = [TwoDCoord(0, 0), TwoDCoord(0, 1), TwoDCoord(1, 0), TwoDCoord(1, 1)]

@views function iscrowded(cmat::T, ij::TwoDCoord) where {T<:TwoDAbstractCellArray}
    """
    Return true if any of the two-by-two's which contain (i,j)
    have three or more filled cells.
    Should only be applied to cells that are marked `missing`.
    """
    istrue(cmat[ij]) && return true  # trivially true
    # remaining applies only to false or missing values
    let (i, j) = Tuple(ij)
        for ii = (i-1):i
            !((1 <= ii) && (ii + 1 <= size(cmat, 1))) && continue
            for jj = (j-1):j
                !((1 <= jj) && (jj + 1 <= size(cmat, 2))) && continue
                (count(cmat[ii:ii+1, jj:jj+1]) === 3) && return true
            end
        end
    end
    false
end

@views function iscrowded(cmat::TC) where {TC<:TwoDAbstractCellArray}
    out = zeros(Bool, size(cmat))
    for i = 1:(size(cmat, 1)-1)
        for j = 1:(size(cmat, 2)-1)
            window = istrue.(cmat[i:i+1, j:j+1])
            cnt = count(window)
            if (cnt === 4)
                out[i:i+1, j:j+1] .= true
            elseif (cnt === 3)
                out[findfirst(.!window)+TwoDCoord(i - 1, j - 1)] = true
            end
        end
    end
    out
end

@views function overcrowded(cmat::TS) where {TS<:TwoDSolutionCellArray}
    for i = 1:(size(cmat, 1)-1)
        for j = 1:(size(cmat, 2)-1)
            all(cmat[i:i+1, j:j+1]) && return true
        end
    end
    false
end

function neighbor_sites(cmat::T, ij::TwoDCoord) where {T<:TwoDSolutionCellArray}
    filter(iijj -> inbounds(cmat, iijj),
        broadcast(vnn -> ij + vnn, VON_NEUMANN_NEIGHBORHOOD_2D))
end

function empty_neighbor_sites(cmat::T, ij::TwoDCoord) where {T<:TwoDSolutionCellArray}
    filter(x -> isempty(cmat[x]), neighbor_sites(cmat, ij))
end

function filled_neighbor_sites(cmat::T, ij::TwoDCoord) where {T<:TwoDSolutionCellArray}
    filter(x -> istrue(cmat[x]), neighbor_sites(cmat, ij))
end


# helpers.
@views inbounds(cvec::T, i::OneDCoord) where {T<:OneDAbstractCellArray} = (1 <= i[1]) && (i[1] <= length(cvec))
@views inbounds(cmat::T, i::Int, dim::Int) where {T<:TwoDAbstractCellArray} = (1 <= i) && (i <= size(cmat, dim))
@views inbounds(cmat::T, ij::TwoDCoord) where {T<:TwoDAbstractCellArray} = inbounds(cmat, ij[1], 1) && inbounds(cmat, ij[2], 2)
@views Base.count(cmat::T) where {T<:OneOrTwoDAbstractCellArray} = sum(istrue, cmat)
Base.isnothing(cmat::T) where {T<:TwoDAbstractCellArray} = false
