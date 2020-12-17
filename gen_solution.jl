export generate_solution

using StatsBase: sample, Weights
using Statistics: mean

function construct_gaussian_kernel(k::T) where T <: Real
    k = max(3, trunc(Int8, k))
    iseven(k) && (k -= 1)  # k should be odd so that the kernel is always centered on a cell
    v = float.(map(x->binomial(k-1,x), 0:(k-1)))
    kernel = v * v'  # outer product
    kernel ./= sum(kernel)
end

@views function convolve(smat::TC, kernel::Matrix{TK}, ij::TwoDCoord) where {TC <: TwoDSolutionCellArray, TK <: AbstractFloat}
    c = 0.0
    let k = size(kernel, 1), kk = trunc(Int, (k - 1) / 2), (i, j) = ij
        for u=1:k
            ii = i + (u - 1) - kk
            !inbounds(smat, ii, 1) && continue
            for v=1:k
                jj = j + (v - 1) - kk
                !inbounds(smat, jj, 2) && continue
                smat[ii, jj] && (c += kernel[u, v])
            end
        end
    end
    c
end

GOLDEN_RATIO = (sqrt(5) + 1) / 2
BASE_COEFF = 1 / GOLDEN_RATIO

KERNEL_WIDTH_COEFF = BASE_COEFF ^ (1 - log1p(666e-4))
make_kernel(smat::T) where T <: TwoDSolutionCellArray = construct_gaussian_kernel(sqrt(prod(size(smat))) ^ KERNEL_WIDTH_COEFF)

NEXT_SITE_COEFF = BASE_COEFF ^ (1 + log(666e-2))
@views function choose_next_site(smat::TC, kernel::Matrix{TK}, ijs::Vector{TwoDCoord}) where {TC <: TwoDSolutionCellArray, TK <: AbstractFloat}
    """
    Chooses the next site based on weights derived from
    the local density of true cells in the current state.
    """
    sample(1:length(ijs), Weights(1 ./ (map(ij->convolve(smat, kernel, ij), ijs) .^ NEXT_SITE_COEFF)))
end

RAND_SKIP_COEFF = BASE_COEFF ^ (1 - log1p(55555e-6))
function random_skip(smat::T, ij::TwoDCoord) where T <: TwoDSolutionCellArray
    kernel = construct_gaussian_kernel(sqrt(prod(size(smat))))
    w = convolve(smat, kernel, ij)
    (rand() < (w ^ RAND_SKIP_COEFF))
end

function generate_solution(n::Int, m::Int)
    smat = fill(false, n, m)
    ijs = TwoDCoord[(rand(1:n), rand(1:m))]
    kernel = make_kernel(smat)
    while (length(ijs) > 0)
        idx = choose_next_site(smat, kernel, ijs)
        ij = splice!(ijs, idx)
        (iscrowded(smat, ij) || random_skip(smat, ij)) && continue
        smat[ij...] = true  # if we get here, flip site to true
        next_ijs = empty_neighbor_sites(smat, ij)
        if (length(next_ijs) > 0)
            append!(ijs, next_ijs)
            unique!(ijs)
        end
    end
    smat
end
generate_solution(n::Int) = generate_solution(n, n)
