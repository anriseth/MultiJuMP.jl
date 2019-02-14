
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


"""
Extract the Utopia point from the MultiData struct.
"""
function getutopia(m::MultiData)
    Fmin = diag(m.Phi)
    # Convert values depending on whether
    # we maximize or minimize them.
    multiplier = (obj.sense == :Max ? -1. : 1. for obj in m.objectives)
    return [mult * f for (mult, f) in zip(multiplier, Fmin)]
end

"""
Extract the Nadir point from the MultiData struct.
"""
function getnadir(m::MultiData)
    Fmax = (f for f in maximum(m.Phi, dims=2))
    # Convert values depending on whether
    # we maximize or minimize them.
    multiplier = (obj.sense == :Max ? -1. : 1. for obj in m.objectives)
    return [mult * f for (mult, f) in zip(multiplier, Fmax)]
end
