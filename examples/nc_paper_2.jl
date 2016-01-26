#==
Example taken from
Messac et. al., Amir Ismail-Yahaya, and Christopher A Mattson.
The normalized normal constraint method for generating the pareto frontier.
Structural and multidisciplinary optimization, 25(2):86–98, 2003.

This example shows that the NBI method will generate local Pareto points that
are not global and dominated points that are not locally optimal.
We can make sure all NBI points are locally Pareto optimal by using an
inequality-constrained extension. (TODO: implement this option and show it here)
==#

using MultiJuMP, JuMP
using Ipopt

m = MultiModel(solver = IpoptSolver())
@defVar(m, 0 <= x[i=1:2] <= 5)
@defNLExpr(m, f1, x[1])
@defNLExpr(m, f2, x[2])
@addNLConstraint(m, 5exp(-x[1])+2exp(-0.5(x[1]-3)^2) <= x[2])

obj1 = SingleObjective(f1, sense = :Min)

# As the problem is nonconvex, we have to supply the
# initial value to get the global minimum of f2
obj2 = SingleObjective(f2, sense = :Min,
                       iv = Dict{Symbol, Any}(:x => [5., 0.]))

multim = getMultiData(m)
multim.objectives = [obj1, obj2]
multim.pointsperdim = 60
solve(m, method = :NBI)

#plotfront(multim)

function plot_comparison()
    function pltfun(x1)
        setValue(x[1], x1)
        setValue(x[2:n], zeros(n-1))
        getValue(f2)
    end
    pltfun(x1) = 5exp(-x1)+2exp(-0.5(x1-3)^2)

    numpoints = length(multim.paretofront)
    f1arr = convert(Array{Float64},
                    [val[1] for val in multim.paretofront])
    f2arr = convert(Array{Float64},
                    [val[2] for val in multim.paretofront])

    discl = layer(x=f1arr, y=f2arr, Geom.point)
    ctsl = layer(x->pltfun(x), 0, 5,
                 Theme(default_color=colorant"orange"))
    plot(discl, ctsl,
         Guide.xlabel("f<sub>1</sub>"), Guide.ylabel("f<sub>2</sub>"),
         Guide.title("Pareto front with $numpoints points"))
end

frontplot = plot_comparison()
