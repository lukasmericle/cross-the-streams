# A generator and solver for Cross The Streams puzzles

## Puzzle

[Cross The Streams](https://logic-puzzle.fandom.com/wiki/Cross_The_Streams) is a variant of the puzzle game [Nonogram](https://en.wikipedia.org/wiki/Nonogram), sometimes also called Picross or Griddler.

Cross The Streams generalizes Nonograms by introducing wildcards into the puzzle clues. The solver code is configured to solve both Nonograms and Cross The Streams puzzles automatically based on their respective solution conditions. (Nonogram solution conditions are a subset of Cross The Streams's.)

## Implementation

This implementation uses dynamic programming to estimate grid cell occupancy probabilities for updating the solution state. It optimizes the search order by re-checking the easiest checks first, for the biggest wins early before deeper tree explorations takes place. In cases of ambiguity, a branch-and-bound search is used to explore possible solutions. Minimal heuristics are implemented: if there are no more deterministic updates to make, we make a best guess based on the entropy of the estimated distribution for each grid cell.

## Drawbacks

The solver in general is agnostic to whole-grid state, *i.e.*, it considers grid cell state probabilities independently, and considers all rows/columns independent from all other rows/columns. Notably, one thing this solver does *not* do right now (in the case of Cross The Streams) is some kind of look-ahead approach to ensure only one connected component of filled cells.

As this implements a dynamic-programming-based tree search, there are certain conditions and constraints which do not necessarily inform the search, instead only filtering infeasible solutions *ex post facto*. One way to improve this would be some kind of SAT or general constraint solver. The difficulty then is encoding the rules in an efficient manner for solving in such a way that partial grid states can be represented and reasoned over.

## Possible Improvements / Extensions

* Add a GUI for generating and solving both nonograms and cross the streams puzzles.
* Better initial behavior (avoid unnecessary recounts).
  * Estimate or calculate number of states for given clue row/col to prioritize solver's traversal.
* Improve recount speed.
  * Cache a representation of valid row/col states for faster recounts.
  * There may be helpful shortcuts for counting valid states for more exotic cases (*e.g.*. sequential `?`s, *etc.*)
  * May require switching to another solver approach (recognizing runs in row/col states and fitting clues to that).
* Improve cross the streams solution generation.
* Minimize number of validity checks (requires verifying correctness of algorithm).
