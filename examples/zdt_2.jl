#==
This is the modified ZDT example taken from
Esmaile Khorram, K Khaledian, and M Khaledyan.
A numerical method for
constructing the pareto front of multi-objective optimization problems.
Journal of Computational and Applied Mathematics, 261:158–171, 2014.

This is an example of a nonconvex Pareto front, where the weighted sum method
would not work.
==#
using JuMP, MultiJuMP, Ipopt

m = MultiModel(solver=IpoptSolver())
n = 30

l = -ones(n); l[1] = 0
u = ones(n)
@variable(m, l[i] <= x[i=1:n] <= u[i])
@NLexpression(m, f1, x[1])
@NLexpression(m, g, 1 + 9*sum{x[j]^2, j=2:n}/(n-1))
@NLexpression(m, h, 1 - (f1/g)^2)
@NLexpression(m, f2, g*h)

obj1 = SingleObjective(f1, sense = :Min)

obj2 = SingleObjective(f2, sense = :Min)

multim = getMultiData(m)
multim.objectives = [obj1, obj2]
multim.pointsperdim = 30
solve(m, method = :NBI)

plotfront(multim)
