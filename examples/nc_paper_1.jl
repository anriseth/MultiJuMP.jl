#==
Example taken from
Messac et. al., Amir Ismail-Yahaya, and Christopher A Mattson.
The normalized normal constraint method for generating the pareto frontier.
Structural and multidisciplinary optimization, 25(2):86–98, 2003.

For the NC method, it is necessary to scale the objectives to get
an evenly spaced Pareto front. This example shows that this is
not necessary with the NBI method.
==#


using MultiJuMP, JuMP
using Ipopt

m = multi_model(solver = IpoptSolver())
@variable(m, x[i=1:2])
@NLexpression(m, f1, x[1])
@NLexpression(m, f2, x[2])
@NLconstraint(m, ((x[1]-20)/20)^8+(x[2]-1)^8 <= 1)

obj1 = SingleObjective(f1, sense = :Min)
obj2 = SingleObjective(f2, sense = :Min)

multim = get_multidata(m)
multim.objectives = [obj1, obj2]
multim.pointsperdim = 20
solve(m, method = :NBI)

using Plots
plot(multim)
