
export CellState

CellState = Union{Bool, Missing}

istrue(x::CellState) = !ismissing(x) && x
isfalse(x::CellState) = !ismissing(x) && !x
isfilled(x::CellState) = istrue(x)
isempty(x::CellState) = !isfilled(x)

convert(type::Type{Missing}, state::CellState) = ismissing(state) ? state : error("This cell state is a Boolean.")
convert(type::Type{Bool}, state::CellState) = !ismissing(state) ? state : error("This cell state is `missing`")
convert(type::Type{CellState}, state::Missing) = missing
convert(type::Type{CellState}, state::Bool) = state
