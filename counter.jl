export VectorStateCounter, MatrixStateCounter


# types for state counters.
abstract type AbstractStateCounter end

mutable struct VectorStateCounter <: AbstractStateCounter
    cumul::Vector{Int}
    n::Int
end

function VectorStateCounter(m::Int; init::String="blank")
    if (init === "blank")
        # no states counted so far.
        cumul = zeros(Int, m)
        n = 0
    elseif (init === "empty")
        # one state counted where there are no `true`s.
        cumul = zeros(Int, m)
        n = 1
    elseif (init === "full")
        # one state counted where there are only `true`s.
        cumul = ones(Int, m)
        n = 1
    elseif (init === "askmissing")
        # shortcut for the state counter that results from one asterisk.
        (m == 0) && return VectorStateCounter(Int[], 1)
        cumul = fill(2^(m - 1), m)
        n = 2^m
    end
    VectorStateCounter(cumul, n)
end

VectorStateCounter(cvec::T; init::String="blank") where {T<:OneDAbstractCellArray} = VectorStateCounter(length(cvec); init=init)

cumul(counter::VectorStateCounter) = counter.cumul
n(counter::VectorStateCounter) = convert(BigInt, counter.n)

MaybeVectorStateCounter = Union{Nothing,VectorStateCounter}

mutable struct MatrixStateCounter <: AbstractStateCounter
    rows::Vector{VectorStateCounter}
    cols::Vector{VectorStateCounter}
end

function MatrixStateCounter(n::Int, m::Int; init::String="blank")
    MatrixStateCounter(
        [VectorStateCounter(m; init=init) for _ = 1:n],
        [VectorStateCounter(n; init=init) for _ = 1:m]
    )
end

rows(counter::MatrixStateCounter) = counter.rows
cols(counter::MatrixStateCounter) = counter.cols

n(counter::MatrixStateCounter) = prod(n.(rows(counter))) * prod(n.(cols(counter)))


# helpers and common Base methods.
Base.copy(counter::VectorStateCounter) = VectorStateCounter(copy(counter.cumul), copy(counter.n))
Base.copy(counter::MatrixStateCounter) = MatrixStateCounter(copy.(rows(counter)), copy.(cols(counter)))

Base.length(counter::VectorStateCounter) = Base.length(cumul(counter))

Base.size(counter::VectorStateCounter) = (length(counter),)
Base.size(counter::MatrixStateCounter, dim::Int) = (dim === 1) ? length(rows(counter)) : ((dim === 2) ? length(cols(counter)) : error("matrix state counter has only two dimensions"))
Base.size(counter::MatrixStateCounter) = (length(rows(counter)), length(cols(counter)))

Base.isnothing(counter::VectorStateCounter) = false
Base.isnothing(counter::MatrixStateCounter) = false

maxn(counter::VectorStateCounter) = 2^convert(BigInt, length(counter))  # this is the number of states when `cluevec = [Asterisk()]`
maxn(counter::MatrixStateCounter) = 4^convert(BigInt, prod(size(counter)))  # if every row and col have the maximum number of states, this returns the value of the prod of all of them.


# overloading addition.
function (Base.:+)(a::VectorStateCounter, b::VectorStateCounter)
    a.cumul[:] .+= b.cumul
    a.n += b.n
    a
end

# no-op if trying to add `nothing`
(Base.:+)(a::VectorStateCounter, b::Nothing) = a
(Base.:+)(a::Nothing, b::VectorStateCounter) = b


# defining concatenation.
function Base.vcat(counters::Vararg{MaybeVectorStateCounter})
    any(isnothing, counters) && return nothing
    total_states = prod(n, counters)
    counters[1].cumul = vcat(((total_states / n(c)) .* cumul(c) for c = filter(cntr -> length(cntr) > 0, counters))...)
    counters[1].n = total_states
    counters[1]
end

function Base.pushfirst!(counter::VectorStateCounter, s::Bool)
    pushfirst!(counter.cumul, s ? counter.n : 0)
end

function Base.push!(counter::VectorStateCounter, s::Bool)
    push!(counter.cumul, s ? counter.n : 0)
end


# odds calculation for each location.
function odds_rc(xy::T) where {S<:AbstractFloat,T<:AbstractArray{S,1}}
    """Compute the odds from a set of independent probabilities."""
    any(isone, xy) && return 1.0
    any(iszero, xy) && return 0.0  # is a check better than just multiplying by zero?
    prod(xy)
end

odds(counter::VectorStateCounter) = cumul(counter) ./ n(counter)

function odds(counter::MatrixStateCounter)

    """
    Probability of a location being filled is
    the product of the probabilities
    contributed from the respective row and column.
    """

    omat = Array{Float64}(undef, 2, size(counter)...)
    for (i, row) = enumerate(rows(counter))
        omat[1, i, :] .= odds(row)
    end
    for (j, col) = enumerate(cols(counter))
        omat[2, :, j] .= odds(col)
    end

    mapslices(odds_rc, omat, dims=[1])[1, :, :]

end

complexity(counter::T) where {T<:AbstractStateCounter} = BigFloat(n(counter) - 1) / BigFloat(maxn(counter) - 1)

entropy(p::T) where {T<:AbstractFloat} = (isone(p) || iszero(p)) ? 0.0 : -(p * log(p) + (1 - p) * log(1 - p)) / log(2)  # log(2) to normalize output between 0 and 1
entropy(counter::T) where {T<:AbstractStateCounter} = mean(entropy.(odds(counter)))
