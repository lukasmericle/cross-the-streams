
export CellState

CellState = Union{Bool, Missing}

istrue(x::CellState) = !ismissing(x) && x
isfalse(x::CellState) = !ismissing(x) && !x
isfilled(x::CellState) = istrue(x)
isempty(x::CellState) = !isfilled(x)
