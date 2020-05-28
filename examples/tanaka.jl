#==
This example shows a disconnected Pareto front.

Tanaka, Masahiro, et al.
"GA-based decision support system for multicriteria optimization."
Systems, Man and Cybernetics, 1995. Intelligent Systems for the 21st Century.,
EEE International Conference on. Vol. 2. IEEE, 1995.
==#

using MultiJuMP, JuMP
using Ipopt

m = multi_model(Ipopt.Optimizer)
@variable(m, 0 <= x[i=1:2] <= π)
@NLexpression(m, f1, x[1])
@NLexpression(m, f2, x[2])
@NLconstraint(m, x[1]^2+x[2]^2 -1 - 0.1cos(16atan(x[1]/x[2])) >= 0)
@NLconstraint(m, (x[1]-0.5)^2+(x[2]-0.5)^2 <= 0.5)

obj1 = SingleObjective(f1, sense = MOI.MIN_SENSE,
                       iv = Dict{String, Any}("x[$i]" => j for (i,j) in enumerate([0., 1.])))

# As the problem is nonconvex, we have to supply the
# initial value to get the global minimum of f2
obj2 = SingleObjective(f2, sense = MOI.MIN_SENSE,
                       iv = Dict{String, Any}("x[$i]" => j for (i,j) in enumerate([1., 0.])))

multim = get_multidata(m)
multim.objectives = [obj1, obj2]
multim.pointsperdim = 40
optimize!(m, method = NBI())

using Plots
plot(multim)
