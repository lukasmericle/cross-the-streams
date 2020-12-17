export VectorStateCounter, MatrixStateCounter

abstract type AbstractStateCounter end

mutable struct VectorStateCounter <: AbstractStateCounter
    cumul::Vector{Int}
    n::Int
end
VectorStateCounter() = nothing
function VectorStateCounter(m::Int; init::String="blank")
    if (init === "blank")
        cumul = zeros(Int, m)
        n = 0
    elseif (init === "empty")
        cumul = zeros(Int, m)
        n = 1
    elseif (init === "full")
        cumul = ones(Int, m)
        n = 1
    elseif (init === "askmissing")
        cumul = fill(2^(m-1), m)
        n = 2^m
    end
    VectorStateCounter(cumul, n)
end
VectorStateCounter(cvec::T; init::String="blank") where T <: OneDAbstractCellArray = VectorStateCounter(length(cvec); init=init)

cumul(counter::VectorStateCounter) = counter.cumul
n(counter::VectorStateCounter) = counter.n

MaybeVectorStateCounter = Union{Nothing, VectorStateCounter}

mutable struct MatrixStateCounter <: AbstractStateCounter
    rows::Vector{VectorStateCounter}
    cols::Vector{VectorStateCounter}
end
MatrixStateCounter(n::Int, m::Int; init::String="blank") = MatrixStateCounter([VectorStateCounter(m; init=init) for _=1:n],
                                                                              [VectorStateCounter(n; init=init) for _=1:m])

rows(counter::MatrixStateCounter) = counter.rows
cols(counter::MatrixStateCounter) = counter.cols

Base.copy(counter::VectorStateCounter) = VectorStateCounter(copy(cumul(counter)), copy(n(counter)))
Base.copy(counter::MatrixStateCounter) = MatrixStateCounter(copy.(rows(counter)), copy.(cols(counter)))

Base.length(counter::VectorStateCounter) = Base.length(cumul(counter))

Base.size(counter::VectorStateCounter) = (length(counter),)
Base.size(counter::MatrixStateCounter, dim::Int) = (dim === 1) ? length(rows(counter)) : ((dim === 2) ? length(cols(counter)) : error("matrix state counter has only two dimensions"))
Base.size(counter::MatrixStateCounter) = (Base.size(counter, 1), Base.size(counter, 2))

Base.isnothing(counter::VectorStateCounter) = false
Base.isnothing(counter::MatrixStateCounter) = false

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

num_states(counter::VectorStateCounter) = convert(BigInt, n(counter))
num_states(counter::MatrixStateCounter) = prod(num_states.(rows(counter))) * prod(num_states.(cols(counter)))

max_num_states(counter::VectorStateCounter) = 2 ^ convert(BigInt, length(counter))  # this is the number of states when `cluevec = [Asterisk()]`
max_num_states(counter::MatrixStateCounter) = 4 ^ convert(BigInt, prod(size(counter)))  # if every row and col have the maximum number of states, this returns the value of the prod of all of them.

complexity(counter::T) where T <: AbstractStateCounter = BigFloat(num_states(counter) - 1) / BigFloat(max_num_states(counter) - 1)

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

entropy(p::T) where T <: AbstractFloat = (isone(p) || iszero(p)) ? 0.0 : -(p * log(p) + (1 - p) * log(1 - p)) / log(2)  # log(2) to normalize output between 0 and 1
entropy(counter::T) where T <: AbstractStateCounter = mean(entropy.(odds(counter)))
