#==
This example shows a simple linear bi-objective problem.
==#

using MultiJuMP, JuMP
using Clp: ClpSolver

const mmodel = multi_model(solver = ClpSolver(), linear = true)
const y = @variable(mmodel, 0 <= y <= 10.0)
const z = @variable(mmodel, 0 <= z <= 10.0)
@constraint(mmodel, y + z <= 15.0)

# objectives
const exp_obj1 = @expression(mmodel, -y +0.05 * z)
const exp_obj2 = @expression(mmodel, 0.05 * y - z)
const obj1 = SingleObjective(exp_obj1)
const obj2 = SingleObjective(exp_obj2)

# # setting objectives in the data
const multim = get_multidata(mmodel)
multim.objectives = [obj1, obj2]

solve(mmodel, method = WeightedSum())

using Plots: plot
plot(multim)
