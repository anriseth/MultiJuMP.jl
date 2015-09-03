# extensions to JuMP for multiobjective optimisation

module MultiJuMP

using JuMP
import Ipopt
import ReverseDiffSparse: ParametricExpression
using Compat

export MultiModel, getMultiData

# stores extension data inside JuMP Model
type MultiData
    # TODO: can this work with linear and affine expressions as well?
    #objectives::Array{ParametricExpression{0},1}
    f1
    f2

    #
    ninitial::Int
    maxlength::Float64
    #
    # utopia
    utopiavarvalues::Array{Dict, 1}
    utopia::Array{Float64,1}
    nadir::Array{Float64,1}
    paretofront
end


function MultiModel(;solver=Ipopt.IpoptSolver())
    m = Model(solver=solver)
    m.solvehook = solvehook
    m.ext[:Multi] = MultiData(ParametricExpression{0}, ParametricExpression{0},
                              2, 1.0,
                              Array{Dict}[], Float64[],
                              Float64[], Any[])
    return m
end

function getMultiData(m::Model)
    if haskey(m.ext, :Multi)
        return m.ext[:Multi]::MultiData
    else
        error("This functionality is only available for MultiJuMP models")
    end
end

function solvehook(m::Model; kwargs...)
    multim = getMultiData(m)
    f1 = multim.f1
    f2 = multim.f2
    # Stage 1
    #Normalize the objective functions in the objective space
    @setNLObjective(m, :Min, f1)
    status = solve(m, ignore_solve_hook=true)
    if status != :Optimal
        return status
    end
    push!(multim.utopiavarvalues, Dict([key => getValue(val) for (key, val) in m.varDict]))
    push!(multim.utopia, getValue(f1))

    nadir2 = getValue(f2)

    @setNLObjective(m, :Min, f2)
    status = solve(m, ignore_solve_hook=true)
    if status != :Optimal
        return status
    end
    push!(multim.utopiavarvalues, Dict([key => getValue(val) for (key, val) in m.varDict]))
    push!(multim.utopia, getValue(f2))

    # Add nadir. [with > 2 objectives, we need to use utopiavarvalues and a loop]
    push!(multim.nadir, getValue(f1), nadir2)

    # Stage 2
    # Perform multiobjective optimization using the usual weighted sum approach

    # Stage 3
    # Delete nearly overlapping solutions on the pareto front
    # Stage 4
    # Stage 5
    # Stage 6
    # Stage 7
    # Stage 8

    return status
end

end
