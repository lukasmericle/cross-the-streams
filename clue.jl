export Clue
export QuestionMark, Asterisk
import Base: string
export       string

struct QuestionMark end
struct Asterisk end

Clue = Union{T, QuestionMark, Asterisk} where T <: Int
ClueRepr = Union{T, String} where T <: Int

Base.string(n::QuestionMark) = "?"
Base.string(n::Asterisk) = "*"
