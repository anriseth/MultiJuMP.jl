
function multisolve(m::Model, mdata::MultiData, ::WeightedSum, ::LinearProblem)
    objectives = mdata.objectives
    sensemap = Dict(MOI.MIN_SENSE => 1.0, MOI.MAX_SENSE => -1.0)

    vararr = all_variables(m)

    numobj = length(objectives)
    Phi = zeros(numobj,numobj)

    # Individual minimisations
    for (i, objective) in enumerate(objectives)
        @objective(m, objective.sense, objective.f)
        for (key, value) in objective.initialvalue
            set_start_value(variable_by_name(m, key), value)
        end
        optimize!(m, ignore_optimize_hook=true);
        if !(termination_status(m) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED])
            return termination_status(m)
        end

        push!(mdata.utopiavarvalues, [value(var) for var in vararr])
        push!(mdata.paretovarvalues, [value(var) for var in vararr])

        Phi[:,i] = sensevalue.(objectives)

        push!(mdata.paretofront, value.(objectives))
    end
    Fmax = maximum(Phi, dims=2)
    Fmin = minimum(Phi, dims=2) # == diag(Phi)?
    if Fmax == Fmin
        error("The Nadir and Utopia points are equal") # I think that's what this means?
    end

    mdata.Phi = Phi

    beta = zeros(numobj); beta[end] = 1.0

    betatree = betas(numobj, mdata.pointsperdim-1)

    for betaval in betatree
        if count(t -> t != 0, betaval) == 1
            # Skip individual optimisations as
            # they are already performed
            continue
        end
        @objective(m, MOI.MIN_SENSE,
            sum(betaval[i]*(sensemap[objectives[i].sense]*objectives[i].f -
            Fmin[i])/(Fmax[i]-Fmin[i]) for i in Base.OneTo(numobj))
        )

        optimize!(m, ignore_optimize_hook=true);
        if !(termination_status(m) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED])
            return termination_status(m)
        end

        push!(mdata.paretofront, value.(objectives))
        push!(mdata.paretovarvalues, [value(var) for var in vararr])
    end
    return MOI.OPTIMAL
end

function multisolve(m::Model, mdata::MultiData, ::EpsilonCons, ::LinearProblem)
    objectives = mdata.objectives
    sensemap = Dict(MOI.MIN_SENSE => 1.0, MOI.MAX_SENSE => -1.0)

    vararr = all_variables(m)

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
            set_start_value(variable_by_name(m, key), value)
        end
        optimize!(m, ignore_optimize_hook=true);
        if !(termination_status(m) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED])
            return termination_status(m)
        end

        push!(mdata.utopiavarvalues, [value(var) for var in vararr])
        push!(mdata.paretovarvalues, [value(var) for var in vararr])

        Phi[:,i] = sensevalue.(objectives)

        push!(mdata.paretofront, value.(objectives))
    end
    Fmax = maximum(Phi, dims=2)
    Fmin = minimum(Phi, dims=2) # == diag(Phi)?

    mdata.Phi = Phi
    # maintain last objective, other as constraints
    @objective(m, objectives[end].sense, objectives[end].f)

    @constraint(m, objconstr,
                  sensemap[objectives[1].sense]*objectives[1].f
                  <= Fmax[1])
    betatree = betas(numobj, mdata.pointsperdim-1)

    for betaval in 0.0:0.001:1.0

        set_normalized_rhs(objconstr, betaval*Fmin[1]+(1-betaval)*Fmax[1])

        optimize!(m, ignore_optimize_hook=true);
        if !(termination_status(m) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED])
            return termination_status(m)
        end

        push!(mdata.paretofront, value.(objectives))
        push!(mdata.paretovarvalues, [value(var) for var in vararr])
    end

    return MOI.OPTIMAL
end
