export VectorStateCounter, MatrixStateCounter
import Base: vcat, size, length, iterate, setindex, setindex!, getindex, firstindex, lastindex
export       vcat

using Lazy: @forward

abstract type AbstractStateCounter end

mutable struct VectorStateCounter <: AbstractStateCounter
    cumul::Array{Int, 1}
    n::Int
end

function VectorStateCounter(m::Int; init::String="blank")
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
        cumul = Vector{Int}(undef, m)
        fill!(cumul, 2^(m-1))
        n = 2^m
    end
    VectorStateCounter(cumul, n)
end

VectorStateCounter() = VectorStateCounter(0; init="blank")
VectorStateCounter(cvec::T; init::String="blank") where T <: OneDCellArray = VectorStateCounter(length(cvec); init=init)

cumul(counter::VectorStateCounter) = counter.cumul
n(counter::VectorStateCounter) = counter.n

@forward VectorStateCounter.cumul size
@forward VectorStateCounter.cumul length
@forward VectorStateCounter.cumul iterate
@forward VectorStateCounter.cumul setindex
@forward VectorStateCounter.cumul setindex!
@forward VectorStateCounter.cumul getindex
@forward VectorStateCounter.cumul firstindex
@forward VectorStateCounter.cumul lastindex

function Base.vcat(a::VectorStateCounter, b::VectorStateCounter)
    (length(a) === 0) && return b
    (length(b) === 0) && return a
    a.cumul = vcat(n(b) .* cumul(a), n(a) .* cumul(b))
    a.n *= n(b)
    a
end

Base.vcat(counters::Vararg{VectorStateCounter}) = Base.foldl(Base.vcat, counters, init=VectorStateCounter(0; init="empty"))

function (Base.:+)(a::VectorStateCounter, b::VectorStateCounter)
    a.cumul .+= b.cumul
    a.n += b.n
    a
end

function (Base.:*)(a::VectorStateCounter, b::VectorStateCounter)
    Base.vcat(a, b)
end

odds(counter::VectorStateCounter) = cumul(counter) ./ n(counter)

mutable struct MatrixStateCounter <: AbstractStateCounter
    rows::Vector{VectorStateCounter}
    cols::Vector{VectorStateCounter}
end

rows(counter::MatrixStateCounter) = counter.rows
cols(counter::MatrixStateCounter) = counter.cols

function MatrixStateCounter(pzl::Puzzle)  # initializes based on puzzle description assuming grid is blank
    n, m = size(pzl)
    dummy_row = init_cvec(m)
    dummy_col = init_cvec(n)
    MatrixStateCounter(map(i -> count_states(dummy_row, rows(pzl)[i]), 1:n),
                       map(j -> count_states(dummy_col, cols(pzl)[j]), 1:m))  # TODO: redo with views
end

function MatrixStateCounter(cmat::T, pzl::Puzzle) where T <: TwoDCellArray  # initializes based on puzzle description and current state of grid
    MatrixStateCounter(map(i -> count_states(cmat[i,:], rows(pzl)[i]), 1:size(cmat, 1)),
                       map(i -> count_states(cmat[:,j], cols(pzl)[j]), 1:size(cmat, 2)))  # TODO: redo with views
end

function MatrixStateCounter(n::Int, m::Int; init::String="blank")
    MatrixStateCounter([VectorStateCounter(m; init=init) for _=1:n],
                       [VectorStateCounter(n; init=init) for _=1:m])
end

Base.size(counter::MatrixStateCounter, dim::Int) = (dim == 1) ? length(rows(counter)) : ((dim == 2) ? length(cols(counter)) : error("grid state counter has only two dimensions"))
Base.size(counter::MatrixStateCounter) = (Base.size(counter, 1), Base.size(counter, 2))

function odds_rowcol(xy::T) where {S <: AbstractFloat, T <: AbstractArray{S,1}}
    any(isone.(xy)) && return 1.0
    any(iszero.(xy)) && return 0.0
    prod(xy)
end

function odds(counter::MatrixStateCounter)
    n, m = size(counter)
    omat = Array{Float64}(undef, 2, n, m)
    for i=1:n
        omat[1,i,:] .= odds(rows(counter)[i])
    end
    for j=1:m
        omat[2,:,j] .= odds(cols(counter)[j])
    end
    mapslices(odds_rowcol, omat, dims=[1])[1,:,:]
end

function entropy(p::T) where T <: AbstractFloat
    (isone(p) || iszero(p)) && return 0.0
    -p * log(p) - (1-p) * log(1-p)
end

entropy(counter::MatrixStateCounter) = entropy.(odds(counter))

function hint(counter::MatrixStateCounter, cmat::T) where T <: TwoDCellArray
    findmin()
end