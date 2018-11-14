function _solve_ws(m::Model)
    multim = getMultiData(m)
    if multim.linear
        return _solve_ws_lin(m, multim)
    end
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
    Fmax = maximum(Phi, dims=2)
    Fmin = minimum(Phi, dims=2) # == diag(Phi)?
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
    if multim.linear
        return _solve_eps_lin(m, multim)
    end
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
    Fmax = maximum(Phi, dims=2)
    Fmin = minimum(Phi, dims=2) # == diag(Phi)?

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
