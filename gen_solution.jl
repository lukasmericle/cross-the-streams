export generate_solution
import Base: string
export       string

using StatsBase: sample, Weights
using Statistics: mean

function construct_gaussian_kernel(k::T) where T <: Real
    k = max(3, trunc(Int8, k))
    if iseven(k) k -= 1 end  # k should be odd so that the kernel is always centered on a cell
    v = float.(map(x->binomial(k-1,x), 0:(k-1)))
    kernel = v * v'  # outer product
    kernel ./= sum(kernel)
end

function convolve(grid::TG, kernel::Array{TK, 2}, ij::TwoDCoord) where {TG <: TwoDGrid, TK <: AbstractFloat}
    c = 0.0
    let k = size(kernel, 1), kk = trunc(Int, (k-1)/2), (i, j) = ij
        for u=1:k
            ii = i + (u - 1) - kk
            !inbounds(grid, ii, 1) && continue
            for v=1:k
                jj = j + (v - 1) - kk
                !inbounds(grid, jj, 2) && continue
                if istrue(grid[ii, jj])
                    c += kernel[u,v]
                end
            end
        end
    end
    c
end

GOLDEN_RATIO = (sqrt(5) + 1) / 2
BASE_COEFF = 1 / GOLDEN_RATIO

KERNEL_WIDTH_COEFF = BASE_COEFF ^ (1 - log1p(666e-4))
make_kernel(grid::T) where T<:TwoDGrid = construct_gaussian_kernel(sqrt(prod(size(grid))) ^ KERNEL_WIDTH_COEFF)

NEXT_SITE_COEFF = BASE_COEFF ^ (1 + log(666e-2))
function choose_next_site(grid::TG, kernel::Array{TK, 2}, ijs::Array{TwoDCoord}) where {TG <: TwoDGrid, TK <: AbstractFloat}
    """
    Chooses the next site based on weights derived from
    the local density of true cells in the current grid state.
    """
    inv_w = map(ij->convolve(grid, kernel, ij), ijs)  # check if `convolve.(grid, kernel, ijs)` works
    sample(1:length(ijs), Weights(@. 1 / (inv_w ^ NEXT_SITE_COEFF)))
end

RAND_SKIP_COEFF = BASE_COEFF ^ (1 - log1p(55555e-6))
function random_skip(grid::T, ij::TwoDCoord) where T <: TwoDGrid
    kernel = construct_gaussian_kernel(sqrt(prod(size(grid))))
    w = convolve(grid, kernel, ij)
    (rand() < (w ^ RAND_SKIP_COEFF))
end

function generate_solution(n::Int, m::Int)

    grid = fill(false, n, m)

    ijs = TwoDCoord[(rand(1:n), rand(1:m))]
    kernel = make_kernel(grid)

    while length(ijs) > 0

        idx = choose_next_site(grid, kernel, ijs)
        ij = splice!(ijs, idx)

        # skip (i, j) if crowded by 3 other true cells or we want to randomly skip it
        (iscrowded(grid, ij) || random_skip(grid, ij)) && continue

        # if we get here, flip site to true
        grid[ij...] = true

        # get empty neighbors
        next_ijs = empty_neighbor_sites(grid, ij)
        if length(next_ijs) > 0
            append!(ijs, next_ijs)
            unique!(ijs)
        end

    end

    grid

end

function generate_solution(n::Int)
    generate_solution(n, n)
end
