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
    normalf1
    normalf2

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
                              ParametricExpression{0}, ParametricExpression{0},
                              2, 1.0,
                              Dict[], Float64[],
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

function _solve_stage1(m::Model)
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
    valf1 = getValue(f1)
    valf2 = getValue(f2)
    nadir2 = valf2
    push!(multim.utopia, valf1)
    push!(multim.paretofront, [valf1, valf2])

    @setNLObjective(m, :Min, f2)
    status = solve(m, ignore_solve_hook=true)
    if status != :Optimal
        return status
    end
    push!(multim.utopiavarvalues, Dict([key => getValue(val) for (key, val) in m.varDict]))
    valf1 = getValue(f1)
    valf2 = getValue(f2)

    push!(multim.utopia, valf2)
    push!(multim.paretofront, [valf1, valf2])
    # Add nadir. [with > 2 objectives, we need to use utopiavarvalues and a loop]
    push!(multim.nadir, valf1, nadir2)

    @defNLExpr(nf1, (f1 - multim.utopia[1]) / (multim.nadir[1] - multim.utopia[1]))
    @defNLExpr(nf2, (f2 - multim.utopia[2]) / (multim.nadir[2] - multim.utopia[2]))

    multim.normalf1 = nf1
    multim.normalf2 = nf2

    return status

end

function _solve_stage2(m::Model)
    # Perform multiobjective optimization using the usual weighted sum approach
    multim = getMultiData(m)
    stepsize = 1.0 / multim.ninitial

    # We can't use multim.normalf1 in the @setNLObjective macro
    # local definitions like nf1 = multim.normalf1 works
    nf1 = multim.normalf1
    nf2 = multim.normalf2

    for weight in [stepsize:stepsize:1-stepsize]
        @defNLExpr(weightedobj, nf1 * weight + nf2 * (1.0-weight))
        @setNLObjective(m, :Min, weightedobj)
        status = solve(m, ignore_solve_hook=true)
        if status != :Optimal
            return status
        end
        valf1 = getValue(multim.f1)
        valf2 = getValue(multim.f2)

        push!(multim.paretofront, [valf1, valf2])
    end

    return :Optimal # TODO: find a better way to do this?
end

function solvehook(m::Model; kwargs...)

    # Get utopia and nadir points. Normalise functions
    status = _solve_stage1(m)
    if status != :Optimal
        return status
    end

    # Perform standard normalised weighted sum
    status = _solve_stage2(m)
    if status != :Optimal
        return status
    end


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
