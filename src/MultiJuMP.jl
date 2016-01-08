# extensions to JuMP for multiobjective optimisation

module MultiJuMP

using JuMP
import Ipopt
import JuMP: JuMPTypes, getValue

export MultiModel, SingleObjective, getMultiData

type SingleObjective
    f # JuMP-expression TODO: use JuMPTypes or something?
    sense::Symbol
    # TODO: implement bound and initial value in algorithm
    initialvalue::Dict{Symbol,Any} # Variable => Initial value
    bound::Float64 # Hard lower/upper depending on sense
end

SingleObjective() = SingleObjective(Any, :Min, Dict{Symbol,Any}(), NaN)
SingleObjective(f::JuMPTypes; sense::Symbol = :Min,
                iv::Dict{Symbol,Any} = Dict{Symbol,Any}(),
                bound::Float64 = NaN) = SingleObjective(f, sense, iv, bound)

getValue(obj::SingleObjective) =  getValue(obj.f)

function senseValue(obj::SingleObjective)
    # To get Φ and Fstar with right signs
    if obj.sense == :Max
        return -getValue(obj.f)
    else
        return getValue(obj.f)
    end
end

getValue(arr::Array{SingleObjective}) = map(getValue, arr)
senseValue(arr::Array{SingleObjective}) = map(senseValue, arr)

# stores extension data inside JuMP Model
type MultiData
    objectives::Array{SingleObjective}
    f1::SingleObjective
    f2::SingleObjective
    normalf1
    normalf2

    #
    ninitial::Int
    #
    # stored values
    utopiavarvalues::Array{Dict, 1}
    utopia::Array{Float64,1}
    nadir::Array{Float64,1}
    paretovarvalues::Array{Dict, 1}
    paretofront

    Phi::Array{Float64,2}
    Fstar::Array{Float64,1}
end


function MultiModel(;solver=Ipopt.IpoptSolver())
    m = Model(solver=solver)
    m.solvehook = solvehook
    m.ext[:Multi] = MultiData(Array(SingleObjective,0),
                              SingleObjective(), SingleObjective(),
                              Any, Any,
                              10,
                              Dict[], Float64[],
                              Float64[], Dict[], Any[],
                              Array(Float64,2,2), Array(Float64,2))
    return m
end

function getMultiData(m::Model)
    if haskey(m.ext, :Multi)
        return m.ext[:Multi]::MultiData
    else
        error("This functionality is only available for MultiJuMP models")
    end
end

function _solve_ws(m::Model)
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
    push!(multim.paretovarvalues, Dict([key => getValue(val) for (key, val) in m.varDict]))
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
    push!(multim.paretovarvalues, Dict([key => getValue(val) for (key, val) in m.varDict]))
    valf1 = getValue(f1)
    valf2 = getValue(f2)

    push!(multim.utopia, valf2)
    push!(multim.paretofront, [valf1, valf2])

    # Add nadir. [with > 2 objectives, we need to use utopiavarvalues and a loop]
    push!(multim.nadir, valf1, nadir2)

    # TODO: redo nf{1,2}: We don't need the multim.utopia[i] subtraction
    # JuMP doesn't handle divisions well, so make dummy multipliers
    multiplier1 = 1.0 / (multim.nadir[1] - multim.utopia[1])
    multiplier2 = 1.0 / (multim.nadir[2] - multim.utopia[2])
    @defNLExpr(m, nf1, (f1 - multim.utopia[1]) * multiplier1)
    @defNLExpr(m, nf2, (f2 - multim.utopia[2]) * multiplier2)

    multim.normalf1 = nf1
    multim.normalf2 = nf2

    # Perform multiobjective optimization using the usual weighted sum approach
    stepsize = 1.0 / multim.ninitial

    @defNLParam(m, w[i=1:2] == ones(2)[i])
    @setNLObjective(m, :Min, nf1*w[1] + nf2*w[2])

    weights = linspace(0,1,multim.ninitial)
    for weight in weights[2:end-1]
        setValue(w, [weight, 1-weight])
        status = solve(m, ignore_solve_hook=true)
        if status != :Optimal
            return status
        end
        valf1 = getValue(multim.f1)
        valf2 = getValue(multim.f2)

        push!(multim.paretofront, [valf1, valf2])
        push!(multim.paretovarvalues, Dict([key => getValue(val) for (key, val) in m.varDict]))
    end

    return :Optimal # TODO: find a better way to do this?
end

function solve_nbi(m::Model)
    multim = getMultiData(m)
    objectives = multim.objectives
    const sensemap = Dict(:Min => 1.0, :Max => -1.0)

    # Stage 1: Calculate Φ
    numobj = length(objectives)
    Fstar = zeros(numobj)
    Phi = zeros(numobj,numobj)

    # Individual minimisations
    for (i, objective) in enumerate(objectives)
        @setNLObjective(m, objective.sense, objective.f)
        for (key, value) in objective.initialvalue
            setValue(m.varDict[key], value)
        end
        status = solve(m, ignore_solve_hook=true);
        if status != :Optimal
            return status
        end

        push!(multim.utopiavarvalues, Dict([key => getValue(val) for (key, val) in m.varDict]))
        push!(multim.paretovarvalues, Dict([key => getValue(val) for (key, val) in m.varDict]))

        Phi[:,i] = senseValue(objectives)

        push!(multim.paretofront, getValue(objectives))
    end
    Fstar = diag(Phi)

    for j = 1:numobj
        Phi[:,j] = Phi[:,j] - Fstar
    end
    multim.Phi = Phi
    multim.Fstar = Fstar

    # Stage 2: Create NBI subproblems
    @defVar(m, t)
    # TODO: There is a bug in JuMP so it doesn't propagate all
    # the necessary information if we use @setObjective instead of NLObjective
    # TODO: test this with nlprewrite
    @setNLObjective(m, :Max, t)

    beta = zeros(numobj); beta[end] = 1.0
    @defNLParam(m, β[i=1:numobj] == beta[i])
    for (i, objective) in enumerate(objectives)
        @addNLConstraint(m, nbiconstr,
                         sum{Phi[i,j]*(β[j]-t),
                         j=1:numobj; j != i} ==
                         sensemap[objective.sense]*objective.f-Fstar[i])
    end

    # Stage 3: Solve NBI subproblems

    function traverse_chim(numpoints::Int, iter::Array{Int}, level::Int = 1)
        # Recursive function for iterating over all the β-values
        # TODO: make this clearer / faster?
        nummax = numpoints-1
        N = length(β)
        @assert 0 < level < N
        if level < N-2
            # Subtract 1 from loop as we have done the
            # individual optimisation already
            for i = 0:nummax-sum(iter[1:level-1])-1
                iter[level] = i
                traverse(numpoints, β, level+1)
            end
        end

        for i = 0:nummax-sum(iter[1:N-2])-1
            iter[N-1] = i
            iter[N] = nummax - sum(iter[1:N-1])

            if iter[N] == nummax
                # We have done the individual optimisation already
                continue
            end

            setValue(β, iter/nummax)
            status = solve(m, ignore_solve_hook=true);
            if status != :Optimal
                return status
            end

            push!(multim.paretofront, getValue(objectives))
            push!(multim.paretovarvalues, Dict([key => getValue(val) for (key, val) in m.varDict]))
        end
        return :Optimal
    end

    status = traverse_chim(multim.ninitial,
                           Array{Int}(length(objectives)))
    return status
end

function solvehook(m::Model; method = :NBI, kwargs...)
    if method == :WS
        status = _solve_ws(m)
    else
        status = solve_nbi(m)
    end

    return status
end

end
