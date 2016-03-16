#==
This is the modified ZDT example taken from
Esmaile Khorram, K Khaledian, and M Khaledyan.
A numerical method for
constructing the pareto front of multi-objective optimization problems.
Journal of Computational and Applied Mathematics, 261:158–171, 2014.

This is an example of a disconnected Pareto front.
The standard equality-based NBI method fails miserably here, and
will hit an infeasibility error somewhere along the way.
It will be interesting to see how it works with
the inequality-extended NBI method
==#

using JuMP, MultiJuMP, Ipopt

m = MultiModel(solver=IpoptSolver())
n = 30

l = -ones(n); l[1] = 0
u = ones(n)
@variable(m, l[i] <= x[i=1:n] <= u[i])
@NLexpression(m, f1, x[1])
@NLexpression(m, g, 1 + 9*sum{x[j]^2, j=2:n}/(n-1))
@NLexpression(m, h, 1 - sqrt((f1/g)) - (f1/g)*sin(10*π*f1))
@NLexpression(m, f2, g*h)

iv1 = zeros(n)
obj1 = SingleObjective(f1, sense = :Min,
                       iv = Dict{Symbol, Any}(:x => iv1))

# As the problem is nonconvex, we have to supply the
# initial value to get the global minimum of f2
iv2 = zeros(n)
iv2[1] = 0.85
obj2 = SingleObjective(f2, sense = :Min,
                       iv = Dict{Symbol, Any}(:x => iv2))

multim = getMultiData(m)
multim.objectives = [obj1, obj2]
multim.pointsperdim = 60
# NB: solve fails with an infeasibility error 65% through the algorithm,
# You can still plot the results by pasting in the plot_comparison() function below
solve(m, method = :NBI)


function plot_comparison()
    function pltfun(x1)
        setValue(x[1], x1)
        setValue(x[2:n], zeros(n-1))
        getValue(f2)
    end

    numpoints = length(multim.paretofront)
    f1arr = convert(Array{Float64},
                    [val[1] for val in multim.paretofront])
    f2arr = convert(Array{Float64},
                    [val[2] for val in multim.paretofront])

    discl = layer(x=f1arr, y=f2arr, Geom.point)
    ctsl = layer(x->pltfun(x), 0, 1,
                 Theme(default_color=colorant"orange"))
    plot(discl, ctsl,
         Guide.xlabel("f<sub>1</sub>"), Guide.ylabel("f<sub>2</sub>"),
         Guide.title("Pareto front with $numpoints points"))
end

frontplot = plot_comparison()
