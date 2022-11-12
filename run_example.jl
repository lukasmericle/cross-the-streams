include("CTS.jl")

using .CTS
using BenchmarkTools

n = 13

sol = generate_cross_the_streams_solution(n, n)
println("Randomly generated solution:")
println(string(sol))
println()

difficulty = 6.66
puz = CrossTheStreamsPuzzle(sol, difficulty)
println("Generated puzzle from above solution (difficulty = " * string(difficulty) * "):")
println(string(puz))
println()

println("Solving this puzzle takes...")
psol = @btime begin
    solve(puz)
end
# psol = solve(puz)

println("Solved:")
println(string(psol))
println()
sols_equal = verify(psol, sol)

println("Solution equal to original? " * string(sols_equal) * "!")
if (!sols_equal)
    println("Difference is:")
    println(string(soldiff(psol, sol)))
    println("Note: this is not likely a problem with the solver, "
            *
            "but rather is a result of the ambiguity of the puzzle's clues.")
end
