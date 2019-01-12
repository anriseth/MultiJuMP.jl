
abstract type MultiMethod end
struct WeightedSum <: MultiMethod end
struct EpsilonCons <: MultiMethod end

struct NBI <: MultiMethod
    inequality::Bool
end
NBI() = NBI(false)

abstract type ProblemLinearity end
struct LinearProblem <: ProblemLinearity end
struct NonLinearProblem <: ProblemLinearity end

"""
    `multisolve(m::Model, data::MultiData, met::MultiMethod, lin::ProblemLinearity)`

Solves the multiobjective problem of type lin with method met
"""
function multisolve end

function islinear(m::MultiData)
    m.linear ? LinearProblem() : NonLinearProblem()
end

function solvehook(m::Model; method::MultiMethod = WeightedSum(),
                   kwargs...)
    data = get_multidata(m)
    linearity = islinear(data)
    status = multisolve(m, data, method, linearity)
    return status
end
