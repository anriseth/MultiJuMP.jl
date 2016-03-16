# MultiJuMP

[![Build Status](https://travis-ci.org/anriseth/MultiJuMP.jl.svg?branch=master)](https://travis-ci.org/anriseth/MultiJuMP.jl)

MultiJuMP enables the user to easily run multiobjective optimisation problems
and generate Pareto fronts. The code is built as an extension of
[https://github.com/JuliaOpt/JuMP.jl](JuMP).
We have implemented three ways to trace out the Pareto front:
- Normal Boundary Intersection (`solve(m, method = :NBI)`)
- Weighted sums (`solve(m, method = :WS)`)
- Constraint methods (`solve(m, method = :EPS)`)
    * This method only works for biobjective optimisation as of now

## Installation
In Julia, call `Pkg.clone("https://github.com/anriseth/MultiJuMP.jl.git")` to install MultiJuMP.

## Usage
Have a look in the `examples/` directory for different use cases, including
tri-objective Pareto fronts.
As a usage example, we implement the test from
_Das and Dennis, 1998: Normal-boundary intersection: A new method for
generating the Pareto surface in nonlinear multicriteria optimization problems_:

```julia
using MultiJuMP, JuMP
using Ipopt

m = MultiModel(solver = IpoptSolver())
@variable(m, x[i=1:5])
@NLexpression(m, f1, sum{x[i]^2, i=1:5})
@NLexpression(m, f2, 3x[1]+2x[2]-x[3]/3+0.01*(x[4]-x[5])^3)
@NLconstraint(m, x[1]+2x[2]-x[3]-0.5x[4]+x[5]==2)
@NLconstraint(m, 4x[1]-2x[2]+0.8x[3]+0.6x[4]+0.5x[5]^2 == 0)
@NLconstraint(m, sum{x[i]^2, i=1:5} <= 10)

iv1 = [0.3, 0.5, -0.26, -0.13, 0.28] # Initial guess
obj1 = SingleObjective(f1, sense = :Min,
                       iv = Dict{Symbol,Any}(:x => iv1))
obj2 = SingleObjective(f2, sense = :Min)

md = getMultiData(m)
md.objectives = [obj1, obj2]
md.pointsperdim = 20
solve(m, method = :NBI) # method = :WS or method = :EPS
```

Plot with Plots.jl using ```plotfront(md)```, or more generally:
```julia
using Gadfly
f1arr = convert(Array{Float64},
                [val[1] for val in md.paretofront])
f2arr = convert(Array{Float64},
                [val[2] for val in md.paretofront])

nbi = plot(x=f1arr, y=f2arr, Geom.point,
           Guide.xlabel("f1"), Guide.ylabel("f2"))
```
<!-- Github bug
![Pareto front example](./img/pareto_example.svg) -->
![Pareto front example](https://cdn.rawgit.com/anriseth/MultiJuMP.jl/master/img/pareto_example.svg)


##TODO:
- Tell travis to install Ipopt
- Create 3 objective test for :EPS, :NBI and :WS
- Add bounds on the MultiObjective type (So we can ask to only search over subset of pareto front )
    * __This seems to be causing problems for NBI__
- Add objective `t` before the individual optimisations?
    * Then we can warm-start the NBI subproblems after the individual runs?
    * It seemed to cause an issue in revenue-profit optimisation,
    where including `t` caused the algorithm to find a worse, local optimum
- For 3 objectives or more: Make it possible to have different spacing in the
  components of $\beta$
- Implement Eichfelder algorithm?
