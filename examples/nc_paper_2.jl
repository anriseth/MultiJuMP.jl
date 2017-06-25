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

m = MultiModel(solver = IpoptSolver())
@variable(m, 0 <= x[i=1:2] <= 5)
@NLexpression(m, f1, x[1])
@NLexpression(m, f2, x[2])
@NLconstraint(m, 5exp(-x[1])+2exp(-0.5(x[1]-3)^2) <= x[2])

obj1 = SingleObjective(f1, sense = :Min)

# As the problem is nonconvex, we have to supply the
# initial value to get the global minimum of f2
obj2 = SingleObjective(f2, sense = :Min,
                       iv = Dict{Symbol, Any}(:x => [5., 0.]))

multim = getMultiData(m)
multim.objectives = [obj1, obj2]
multim.pointsperdim = 60

solve(m, method = :NBI, inequalityconstraint = true)
#solve(m, method = :NBI)

function pltfun(x1)
    setvalue(x[1], x1)
    setvalue(x[2:n], zeros(n-1))
    getvalue(f2)
end
pltfun(x1) = 5exp(-x1)+2exp(-0.5(x1-3)^2)

plt = plot(multim, label="NBI points")
plot!(plt, pltfun, 0,5, markercolor="orange",
      label="Image boundary")
