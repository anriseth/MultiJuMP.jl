#==
Example taken from
Indraneel Das and John E Dennis.
Normal-boundary intersection: A new method for generating the pareto surface in
nonlinear multicriteria optimization problems.
SIAM Journal on Optimization, 8(3):631–657, 1998.
==#

using MultiJuMP, JuMP
using Ipopt

m = multi_model(solver = IpoptSolver())
@variable(m, x[i=1:5])
@NLexpression(m, f1, sum(x[i]^2 for i=1:5))
@NLexpression(m, f2, 3x[1]+2x[2]-x[3]/3+0.01*(x[4]-x[5])^3)
@NLconstraint(m, x[1]+2x[2]-x[3]-0.5x[4]+x[5]==2)
@NLconstraint(m, 4x[1]-2x[2]+0.8x[3]+0.6x[4]+0.5x[5]^2 == 0)
@NLconstraint(m, sum(x[i]^2 for i=1:5) <= 10)

obj1 = SingleObjective(f1)
obj2 = SingleObjective(f2)

md = get_multidata(m)
md.objectives = [obj1, obj2]
md.pointsperdim = 20
solve(m, method = :NBI)

using Plots
plot(md)
