
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
    sensemap = Dict(:Min => 1.0, :Max => -1.0)
    multiplier = [sensemap[obj.sense] for obj in m.objectives]
    return Fmin .* multiplier
end

"""
Extract the Nadir point from the MultiData struct.
"""
function getnadir(m::MultiData)
    Fmax = maximum(m.Phi, dims=2)
    # Flatten array  as `maximum` returns it as a column vector.
    # TODO: Is there a better way to do this?
    Fmax = Fmax[:]

    # Convert values depending on whether
    # we maximize or minimize them.
    sensemap = Dict(:Min => 1.0, :Max => -1.0)
    multiplier = [sensemap[obj.sense] for obj in m.objectives]
    return Fmax .* multiplier
end
