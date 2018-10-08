# extensions to JuMP for multiobjective optimisation

module MultiJuMP

using JuMP
import Ipopt
import JuMP: JuMPTypes, getvalue
import MathProgBase
using RecipesBase, LaTeXStrings
import Base.warn
import Combinatorics.combinations
using LinearAlgebra
export MultiModel, SingleObjective, getMultiData

mutable struct SingleObjective
    f # JuMP-expression TODO: use JuMPTypes or something?
    sense::Symbol
    initialvalue::Dict{Symbol,Any} # Variable => Initial value
    # TODO: implement bound in algorithm
    bound::Float64 # Hard lower/upper depending on sense
end

SingleObjective() = SingleObjective(Any, :Min, Dict{Symbol,Any}(), NaN)
SingleObjective(f::JuMPTypes; sense::Symbol = :Min,
                iv::Dict{Symbol,Any} = Dict{Symbol,Any}(),
                bound::Float64 = NaN) = SingleObjective(f, sense, iv, bound)

getvalue(obj::SingleObjective) =  getvalue(obj.f)

function senseValue(obj::SingleObjective)
    # To get Φ and Fstar with right signs
    if obj.sense == :Max
        return -getvalue(obj.f)
    else
        return getvalue(obj.f)
    end
end

getvalue(arr::Array{SingleObjective}) = map(getvalue, arr)
senseValue(arr::Array{SingleObjective}) = map(senseValue, arr)

# stores extension data inside JuMP Model
mutable struct MultiData{Tx,To}
    objectives::Vector{SingleObjective}
    f1::SingleObjective
    f2::SingleObjective
    normalf1
    normalf2

    # Number of points in each direction
    # of Pareto submanifold
    pointsperdim::Int
    #
    # stored values
    utopiavarvalues::Vector{Vector{Tx}}
    utopia::Vector{To}
    nadir::Vector{Tx}
    paretovarvalues::Vector{Vector{To}}
    paretofront

    Phi::Array{To,2}
    Fstar::Vector{To}
end


function MultiModel(;solver=Ipopt.IpoptSolver())
    m = Model(solver=solver)
    m.solvehook = solvehook
    m.ext[:Multi] = MultiData(Array{SingleObjective}(undef, 0),
                              SingleObjective(), SingleObjective(),
                              Any, Any,
                              10,
                              Vector{Float64}[], Float64[],
                              Float64[], Vector{Float64}[], Any[],
                              Array{Float64}(undef, 2,2), Array{Float64}(undef, 2))
    return m
end

function getMultiData(m::Model)
    if haskey(m.ext, :Multi)
        return m.ext[:Multi]::MultiData
    else
        error("This functionality is only available for MultiJuMP models")
    end
end

function betas(levels,parts)
    # Sets up the tree of all possible
    # convex combinations of the objectives
    map(x->([x;levels+parts].-[0;x] .-1)/parts,
        combinations(1:(levels+parts-1),(levels-1)))
end

function _solve_ws(m::Model)
    multim = getMultiData(m)
    objectives = multim.objectives
    sensemap = Dict(:Min => 1.0, :Max => -1.0)

    vararr = [JuMP.Variable(m,i) for i in 1:MathProgBase.numvar(m)]

    numobj = length(objectives)
    Phi = zeros(numobj,numobj)

    # Individual minimisations
    for (i, objective) in enumerate(objectives)
        @NLobjective(m, objective.sense, objective.f)
        for (key, value) in objective.initialvalue
            # TODO: What is the correct way to set these values?
            setvalue(m.objDict[key], value)
        end
        status = solve(m, ignore_solve_hook=true);
        if status != :Optimal
            return status
        end

        push!(multim.utopiavarvalues, [getvalue(var) for var in vararr])
        push!(multim.paretovarvalues, [getvalue(var) for var in vararr])

        Phi[:,i] = senseValue(objectives)

        push!(multim.paretofront, getvalue(objectives))
    end
    Fmax = maximum(Phi,2)
    Fmin = minimum(Phi,2) # == diag(Phi)?
    if Fmax == Fmin
        error("The Nadir and Utopia points are equal") # I think that's what this means?
    end

    multim.Phi = Phi

    beta = zeros(numobj); beta[end] = 1.0
    @NLparameter(m, β[i=1:numobj] == beta[i])

    @NLobjective(m, :Min,
                 sum(β[i]*(sensemap[objectives[i].sense]*objectives[i].f -
                           Fmin[i])/(Fmax[i]-Fmin[i]) for i=1:numobj))

    betatree = betas(numobj, multim.pointsperdim-1)

    for betaval in betatree
        if count(t -> t != 0, betaval) == 1
            # Skip individual optimisations as
            # they are already performed
            continue
        end

        setvalue(β, betaval)

        status = solve(m, ignore_solve_hook=true);
        if status != :Optimal
            return status
        end

        push!(multim.paretofront, getvalue(objectives))
        push!(multim.paretovarvalues, [getvalue(var) for var in vararr])
    end

    return :Optimal
end

function _solve_nbi(m::Model, inequalityconstraint::Bool = false)
    multim = getMultiData(m)
    objectives = multim.objectives
    sensemap = Dict(:Min => 1.0, :Max => -1.0)

    vararr = [JuMP.Variable(m,i) for i in 1:MathProgBase.numvar(m)]

    # Stage 1: Calculate Φ
    numobj = length(objectives)
    Fstar = zeros(numobj)
    Phi = zeros(numobj,numobj)

    # Individual minimisations
    for (i, objective) in enumerate(objectives)
        @NLobjective(m, objective.sense, objective.f)
        for (key, value) in objective.initialvalue
            # TODO: What is the correct way to do this? (anonymous variables)
            setvalue(m.objDict[key], value)
        end
        status = solve(m, ignore_solve_hook=true);
        if status != :Optimal
            return status
        end

        push!(multim.utopiavarvalues, [getvalue(var) for var in vararr])
        push!(multim.paretovarvalues, [getvalue(var) for var in vararr])

        Phi[:,i] = senseValue(objectives)

        push!(multim.paretofront, getvalue(objectives))
    end
    Fstar = diag(Phi)

    for j = 1:numobj
        Phi[:,j] = Phi[:,j] - Fstar
    end
    multim.Phi = Phi
    multim.Fstar = Fstar

    # Stage 2: Create NBI subproblems
    @variable(m, t)
    # TODO: There is a bug in JuMP so it doesn't propagate all
    # the necessary information if we use @objective instead of NLObjective
    # TODO: test this with nlprewrite
    @NLobjective(m, :Max, t)

    beta = zeros(numobj); beta[end] = 1.0
    @NLparameter(m, β[i=1:numobj] == beta[i])

    if inequalityconstraint == false
        # Standard NBI
        for (i, objective) in enumerate(objectives)
            @NLconstraint(m,
                          sum(Phi[i,j]*(β[j]-t)
                              for j = 1:numobj if j != i) ==
                          sensemap[objective.sense]*objective.f-Fstar[i])
        end
    else
        # Pascoletti-Serafini extension
        for (i, objective) in enumerate(objectives)
            @NLconstraint(m,
                          sum(Phi[i,j]*(β[j]-t)
                              for j=1:numobj if j != i) >=
                          sensemap[objective.sense]*objective.f-Fstar[i])
        end
    end

    # Stage 3: Solve NBI subproblems
    betatree = betas(numobj, multim.pointsperdim-1)

    for betaval in betatree
        if count(t -> t != 0, betaval) == 1
            # Skip individual optimisations as
            # they are already performed
            continue
        end
        setvalue(β, betaval)
        status = solve(m, ignore_solve_hook=true);
        if status != :Optimal
            return status
        end

        push!(multim.paretofront, getvalue(objectives))
        push!(multim.paretovarvalues, [getvalue(var) for var in vararr])
    end

    return :Optimal
end

function _solve_eps(m::Model)
    multim = getMultiData(m)
    objectives = multim.objectives
    sensemap = Dict(:Min => 1.0, :Max => -1.0)

    vararr = [JuMP.Variable(m,i) for i in 1:MathProgBase.numvar(m)]

    numobj = length(objectives)
    Phi = zeros(numobj,numobj)

    if numobj > 2
        # TODO:
        # The logic  here becomes difficult, as the feasible region will require
        # a dependency between the constraints
        Base.error(":EPS is thought through for > 2 objectives yet")
    end

    # Individual minimisations
    for (i, objective) in enumerate(objectives)
        @NLobjective(m, objective.sense, objective.f)
        for (key, value) in objective.initialvalue
            # TODO: What is the correct way to do this? (anonymous variables)
            setvalue(m.objDict[key], value)
        end
        status = solve(m, ignore_solve_hook=true);
        if status != :Optimal
            return status
        end

        push!(multim.utopiavarvalues, [getvalue(var) for var in vararr])
        push!(multim.paretovarvalues, [getvalue(var) for var in vararr])

        Phi[:,i] = senseValue(objectives)

        push!(multim.paretofront, getvalue(objectives))
    end
    Fmax = maximum(Phi,2)
    Fmin = minimum(Phi,2) # == diag(Phi)?

    multim.Phi = Phi

    @NLobjective(m, objectives[end].sense, objectives[end].f)

    beta = zeros(numobj); beta[end] = 1.0
    @NLparameter(m, β[i=1:numobj] == beta[i])

    @NLconstraint(m, objconstr[i=1:numobj-1],
                  sensemap[objectives[i].sense]*objectives[i].f
                  <= β[i]*Fmin[i]+(1-β[i])*Fmax[i])

    betatree = betas(numobj, multim.pointsperdim-1)

    for betaval in betatree
        if count(t -> t != 0, betaval) == 1
            # Skip individual optimisations as
            # they are already performed
            continue
        end

        setvalue(β, betaval)

        status = solve(m, ignore_solve_hook=true);
        if status != :Optimal
            return status
        end

        push!(multim.paretofront, getvalue(objectives))
        push!(multim.paretovarvalues, [getvalue(var) for var in vararr])
    end

    return :Optimal
end

function solvehook(m::Model; method::Symbol = :NBI,
                   inequalityconstraint::Bool = false, kwargs...)
    if method == :WS
        # Weighted sums
        status = _solve_ws(m)
    elseif method == :EPS
        # Epsilon constraint method
        status = _solve_eps(m)
    elseif method == :NBI
        # Normal boundary intersection
        status = _solve_nbi(m, inequalityconstraint)
    else
        Base.error("Multiobjective method not recognized.")
    end

    return status
end

@recipe function f{Tx,To}(md::MultiData{Tx,To})
    seriestype --> :scatter
    xlabel --> L"$f_1$"
    ylabel --> L"$f_2$"
    zlabel --> L"$f_3$"
    label --> ""
    numobjectives = length(md.objectives)
    if numobjectives > 3 || numobjectives < 2
        Base.error("Only plotting 2d and 3d fronts")
    end

    f1arr = convert(Array{Float64},
                    [val[1] for val in md.paretofront])
    f2arr = convert(Array{Float64},
                    [val[2] for val in md.paretofront])

    title --> "Pareto front with $(length(f1arr)) points"

    if numobjectives == 2
        xyz = (f1arr, f2arr)
    else
        f3arr = convert(Array{Float64},
                        [val[3] for val in md.paretofront])
        xyz = (f1arr, f2arr, f3arr)
    end
    (xyz...,)
end

@recipe function f(m::JuMP.Model)
    plot(getMultiData(m))
end

end
