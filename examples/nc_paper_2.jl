#==
Example taken from
Messac et. al., Amir Ismail-Yahaya, and Christopher A Mattson.
The normalized normal constraint method for generating the pareto frontier.
Structural and multidisciplinary optimization, 25(2):86–98, 2003.

This example shows that the NBI method will generate local Pareto points that
are not global and dominated points that are not locally optimal.
We can make sure all NBI points are locally Pareto optimal by using an
inequality-constrained extension.
==#

using MultiJuMP, JuMP
using Ipopt

m = multi_model(Ipopt.Optimizer)
@variable(m, 0 <= x[i=1:2] <= 5)
@NLexpression(m, f1, x[1])
@NLexpression(m, f2, x[2])
@NLconstraint(m, 5exp(-x[1])+2exp(-0.5(x[1]-3)^2) <= x[2])

obj1 = SingleObjective(f1, sense = MOI.MIN_SENSE)

# As the problem is nonconvex, we have to supply the
# initial value to get the global minimum of f2
obj2 = SingleObjective(f2, sense = MOI.MIN_SENSE,
                       iv = Dict{String, Any}("x[$i]" => j for (i,j) in enumerate([5., 0.])))

multim = get_multidata(m)
multim.objectives = [obj1, obj2]
multim.pointsperdim = 60

optimize!(m, method = NBI(true)) # inequalityconstraint = true

function pltfun(x1)
    set_start_value(x[1], x1)
    set_start_value(x[2:n], zeros(n-1))
    value(f2)
end
pltfun(x1) = 5exp(-x1)+2exp(-0.5(x1-3)^2)

plt = plot(multim, label="NBI points")
plot!(plt, pltfun, 0,5, markercolor="orange",
      label="Image boundary")
