export CellMatrix, SolutionCellMatrix
export init_cmat

CellState = Union{Bool, Missing}

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
OneOrTwoDAbstractCellArray = Union{OneDAbstractCellArray, TwoDAbstractCellArray}

OneDCoord = CartesianIndex{1}
TwoDCoord = CartesianIndex{2}

VON_NEUMANN_NEIGHBORHOOD_1D = [OneDCoord(-1),
                               OneDCoord(1)]
VON_NEUMANN_NEIGHBORHOOD_2D = [TwoDCoord(0,-1),
                               TwoDCoord(-1,0),
                               TwoDCoord(1,0),
                               TwoDCoord(0,1)]

string_vec(cvec::T) where T <: OneDAbstractCellArray = prod(map(x -> (ismissing(x) ? "><" : (x ? "██" : "  ")),  cvec))
Base.string(cvec::T) where T <: OneDCellArray = string_vec(cvec)  # to avoid ambiguity during dispatch, overload with custom method
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
Base.string(cmat::SolutionCellMatrix) = string_mat(cmat)  # to avoid ambiguity during dispatch, overload with custom method
Base.string(cmat::T) where T <: TwoDAbstractCellArray = string_mat(cmat)

istrue(x::CellState) = !ismissing(x) && x
isfalse(x::CellState) = !ismissing(x) && !x
isempty(x::CellState) = ismissing(x) || !x  # == !istrue(x)

init_cvec(n::Int) = convert(CellVector, fill(missing, n))

init_cmat(n::Int, m::Int) = convert(CellMatrix, fill(missing, n, m))

@views inbounds(cvec::T, i::OneDCoord) where T <: OneDAbstractCellArray = (1 <= i[1]) && (i[1] <= length(cvec))
@views inbounds(cmat::T, i::Int, dim::Int) where T <: TwoDAbstractCellArray = (1 <= i) && (i <= size(cmat, dim))
@views inbounds(cmat::T, ij::TwoDCoord) where T <: TwoDAbstractCellArray = inbounds(cmat, ij[1], 1) && inbounds(cmat, ij[2], 2)

@views Base.count(cmat::T) where T <: OneOrTwoDAbstractCellArray = sum(istrue, cmat)

Base.isnothing(cmat::T) where T <: TwoDAbstractCellArray = false

@views function iscrowded(cmat::T, ij::TwoDCoord) where T <: TwoDAbstractCellArray
    """
    Return true if any of the two-by-two's which contain (i,j)
    have three or more filled cells.
    """
    istrue(cmat[ij]) && return true  # trivially true
    let (i, j) = Tuple(ij)
        for ii=(i-1):i
            !((1 <= ii) && (ii+1 <= size(cmat, 1))) && continue
            for jj=(j-1):j
                !((1 <= jj) && (jj+1 <= size(cmat, 2))) && continue
                (count(cmat[ii:ii+1,jj:jj+1]) === 3) && return true
            end
        end
    end
    false
end
@views iscrowded(cmat::TC, ijs::TA) where {TC <: TwoDAbstractCellArray, TA <: Union{AbstractArray{TwoDCoord}, CartesianIndices}} = map(ij -> iscrowded(cmat, ij), ijs)
@views iscrowded(cmat::T) where T <: TwoDAbstractCellArray = iscrowded(cmat, CartesianIndices(fill(1, size(cmat))))

function neighbor_sites(cmat::T, ij::TwoDCoord) where T <: TwoDSolutionCellArray
    filter(iijj -> inbounds(cmat, iijj),
           broadcast(vnn -> ij + vnn, VON_NEUMANN_NEIGHBORHOOD_2D))
end

function empty_neighbor_sites(cmat::T, ij::TwoDCoord) where T <: TwoDSolutionCellArray
    filter(x -> isempty(cmat[x]), neighbor_sites(cmat, ij))
end

function filled_neighbor_sites(cmat::T, ij::TwoDCoord) where T <: TwoDSolutionCellArray
    filter(x -> istrue(cmat[x]), neighbor_sites(cmat, ij))
end
