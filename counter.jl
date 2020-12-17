export VectorStateCounter, MatrixStateCounter
import Base: vcat, isnothing, size, length, iterate, setindex, setindex!, getindex, firstindex, lastindex
export       vcat, isnothing, size, length, iterate, setindex, setindex!, getindex, firstindex, lastindex

using Lazy: @forward

abstract type AbstractStateCounter end

mutable struct VectorStateCounter <: AbstractStateCounter
    cumul::Vector{Int}
    n::Int
end
VectorStateCounter() = nothing
function VectorStateCounter(m::Int; init::String="blank")
    if (init == "blank")
        cumul = zeros(Int, m)
        n = 0
    elseif (init == "empty")
        cumul = zeros(Int, m)
        n = 1
    elseif (init == "full")
        cumul = ones(Int, m)
        n = 1
    elseif (init == "askmissing")
        cumul = fill(2^(m-1), m)
        n = 2^m
    end
    VectorStateCounter(cumul, n)
end
VectorStateCounter(cvec::T; init::String="blank") where T <: OneDCellArray = VectorStateCounter(length(cvec); init=init)

MaybeVectorStateCounter = Union{Nothing, VectorStateCounter}

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

Base.isnothing(counter::VectorStateCounter) = false

function (Base.:+)(a::VectorStateCounter, b::VectorStateCounter)
    a.cumul[:] .+= b.cumul
    a.n += b.n
    a
end
(Base.:+)(a::VectorStateCounter, b::Nothing) = a
(Base.:+)(b::Nothing, a::VectorStateCounter) = a

(Base.:*)(a::MaybeVectorStateCounter, b::MaybeVectorStateCounter) = Base.vcat(a, b)

function Base.vcat(a::VectorStateCounter, b::VectorStateCounter)
    (length(a) === 0) && return b
    (length(b) === 0) && return a
    a.cumul = vcat(n(b) .* cumul(a), n(a) .* cumul(b))
    a.n *= n(b)
    a
end
Base.vcat(a::VectorStateCounter, b::Nothing) = nothing
Base.vcat(b::Nothing, a::VectorStateCounter) = nothing
Base.vcat(counters::Vararg{MaybeVectorStateCounter}) = prod(counters)

mutable struct MatrixStateCounter <: AbstractStateCounter
    rows::Vector{VectorStateCounter}
    cols::Vector{VectorStateCounter}
end
function MatrixStateCounter(n::Int, m::Int; init::String="blank")
    MatrixStateCounter([VectorStateCounter(m; init=init) for _=1:n],
                       [VectorStateCounter(n; init=init) for _=1:m])
end

rows(counter::MatrixStateCounter) = counter.rows
cols(counter::MatrixStateCounter) = counter.cols

Base.isnothing(counter::MatrixStateCounter) = false

Base.size(counter::MatrixStateCounter, dim::Int) = (dim == 1) ? length(rows(counter)) : ((dim == 2) ? length(cols(counter)) : error("matrix state counter has only two dimensions"))
Base.size(counter::MatrixStateCounter) = (Base.size(counter, 1), Base.size(counter, 2))

function odds_rc(xy::T) where {S <: AbstractFloat, T <: AbstractArray{S,1}}
    any(isone.(xy)) && return 1.0
    any(iszero.(xy)) && return 0.0
    prod(xy)
    # sqrt(prod(xy))  # perhaps the geometric mean is better for this application?
end

odds(counter::VectorStateCounter) = cumul(counter) ./ n(counter)
function odds(counter::MatrixStateCounter)
    n, m = size(counter)
    omat = Array{Float64}(undef, 2, n, m)
    for i=1:n omat[1,i,:] .= odds(rows(counter)[i]) end
    for j=1:m omat[2,:,j] .= odds(cols(counter)[j]) end
    mapslices(odds_rc, omat, dims=[1])[1,:,:]
end

num_states(counter::VectorStateCounter) = convert(BigInt, n(counter))
max_num_states(counter::VectorStateCounter) = 2 ^ convert(BigInt, length(counter))  # this is the number of states when `cluevec = [Asterisk()]`
num_states(counter::MatrixStateCounter) = prod(num_states.(rows(counter))) * prod(num_states.(cols(counter)))
max_num_states(counter::MatrixStateCounter) = 2 ^ convert(BigInt, length(rows(counter)) + length(cols(counter)))

complexity(counter::T) where T <: AbstractStateCounter = BigFloat(num_states(counter) - 1) / BigFloat(max_num_states(counter) - 1)

entropy(p::T) where T <: AbstractFloat = (isone(p) || iszero(p)) ? 0.0 : -(p * log(p) + (1 - p) * log(1 - p)) / log(2)  # log(2) to normalize output between 0 and 1
entropy(counter::T) where T <: AbstractStateCounter = mean(entropy.(odds(counter)))
