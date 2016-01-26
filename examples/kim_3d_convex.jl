#==
This example shows multiobjective Paretofronts with 3 objectives.
The front is convex. The NBI method does not catch the full Pareto front,
e.g. as seen by the missing arcs between f1 and f2, whilst f3 is 0.
To capture this area, we need some of the elements in β to be negative.

Taken from
Il Yong Kim and OL De Weck.
Adaptive weighted sum method for multiobjective optimization:
a new method for pareto front generation.
Structural and Multidisciplinary Optimization, 31(2):105–116, 2006.
==#

using MultiJuMP, JuMP
using Ipopt
using JLD

m = MultiModel(solver = IpoptSolver())

@defVar(m, x[i=1:3] >= 0)
@defNLExpr(m, f1, x[1])
@defNLExpr(m, f2, x[2])
@defNLExpr(m, f3, x[3])
@addNLConstraint(m, x[1]^4+2x[2]^3+5x[3]^2<=1)

obj1 = SingleObjective(f1, sense = :Max)#,
                       #iv = Dict{Symbol, Any}(:x => [0, 0., 1.2]))

obj2 = SingleObjective(f2, sense = :Max)#,
                       #iv = Dict{Symbol, Any}(:x => [0, 0., 1.2]))

obj3 = SingleObjective(f3, sense = :Max)#,
                       #iv = Dict{Symbol, Any}(:x => [0, 0., 2.]))


md = getMultiData(m)
md.objectives = [obj1, obj2, obj3]
#md.objectives = [obj1, obj2]
md.pointsperdim = 10
solve(m, method = :NBI)

numpoints = length(md.paretofront)
f1arr = convert(Array{Float64},
                [val[1] for val in md.paretofront])
f2arr = convert(Array{Float64},
                [val[2] for val in md.paretofront])
f3arr = convert(Array{Float64},
                [val[3] for val in md.paretofront])

@save("paretovals_kim_3d_convex.jld", f1arr, f2arr, f3arr)

# using PyPlot
# fig = figure()
# scatter3D(f1arr, f2arr, f3arr)
# show()
# xlabel("f1")
# ylabel("f2")
# zlabel("f3")

using Plots
backend(:pyplot)
plotfront(md)
