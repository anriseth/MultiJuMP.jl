function _solve_ws_lin(m::Model, multim::MultiData)
    objectives = multim.objectives
    sensemap = Dict(:Min => 1.0, :Max => -1.0)

    vararr = [JuMP.Variable(m,i) for i in 1:MathProgBase.numvar(m)]

    numobj = length(objectives)
    Phi = zeros(numobj,numobj)

    # Individual minimisations
    for (i, objective) in enumerate(objectives)
        @objective(m, objective.sense, objective.f)
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

    betatree = betas(numobj, multim.pointsperdim-1)

    for betaval in betatree
        if count(t -> t != 0, betaval) == 1
            # Skip individual optimisations as
            # they are already performed
            continue
        end
        JuMP.setobjective(m, :Min,
            sum(betaval[i]*(sensemap[objectives[i].sense]*objectives[i].f -
            Fmin[i])/(Fmax[i]-Fmin[i]) for i in Base.OneTo(numobj))
        )

        status = solve(m, ignore_solve_hook=true);
        if status != :Optimal
            return status
        end

        push!(multim.paretofront, getvalue(objectives))
        push!(multim.paretovarvalues, [getvalue(var) for var in vararr])
    end
    return :Optimal
end

function _solve_eps_lin(m::Model, multim::MultiData)
    objectives = multim.objectives
    sensemap = Dict(:Min => 1.0, :Max => -1.0)

    vararr = [JuMP.Variable(m,i) for i in Base.OneTo(MathProgBase.numvar(m))]

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
        @objective(m, objective.sense, objective.f)
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
    # maintain last objective, other as constraints
    @objective(m, objectives[end].sense, objectives[end].f)

    @constraint(m, objconstr,
                  sensemap[objectives[1].sense]*objectives[1].f
                  <= Fmax[1])
    betatree = betas(numobj, multim.pointsperdim-1)

    for betaval in 0.0:0.001:1.0

        JuMP.setRHS(objconstr, betaval*Fmin[1]+(1-betaval)*Fmax[1])

        status = solve(m, ignore_solve_hook=true);
        if status != :Optimal
            return status
        end

        push!(multim.paretofront, getvalue(objectives))
        push!(multim.paretovarvalues, [getvalue(var) for var in vararr])
    end

    return :Optimal
end