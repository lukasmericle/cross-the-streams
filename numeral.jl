export Numeral
export QuestionMark, Asterisk
export value
import Base: string
export       string

using Lazy: @forward

struct QuestionMark end
struct Asterisk end

ProblemNumeral = Union{T, QuestionMark, Asterisk} where T <: Int
ProblemNumeralValue = Union{T, String} where T <: Int

Base.string(n::QuestionMark) = "?"
Base.string(n::Asterisk) = "*"
