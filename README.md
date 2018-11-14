# MultiJuMP

[![Build Status](https://travis-ci.org/anriseth/MultiJuMP.jl.svg?branch=master)](https://travis-ci.org/anriseth/MultiJuMP.jl)
[![DOI](https://zenodo.org/badge/41850916.svg)](https://zenodo.org/badge/latestdoi/41850916)

MultiJuMP enables the user to easily run multiobjective optimisation problems
and generate Pareto fronts. The code is built as an extension of
[JuMP](https://github.com/JuliaOpt/JuMP.jl).
We have implemented three ways to trace out the Pareto front:
- Normal Boundary Intersection (`solve(m, method = :NBI)`)
- Weighted sums (`solve(m, method = :WS)`)
- Constraint methods (`solve(m, method = :EPS)`)
    * This method only works for biobjective optimisation as of now

**Disclaimer 1**: MultiJuMP is *not* developed or maintained by the JuMP developers.  

## Installation
In Julia, call `Pkg.add("MultiJuMP")` to install MultiJuMP.

## Usage
Have a look in the `examples/` directory for different use cases, including
tri-objective Pareto fronts.

MultiJuMP supports linear problems using the `linear=true` keyword when
calling `MultiModel(linear=true)`. Currently, only the `:EPS`
and `:WS` methods are supported for linear problems.  

As a usage example, we implement the test from
_Das and Dennis, 1998: Normal-boundary intersection: A new method for
generating the Pareto surface in nonlinear multicriteria optimization problems_:

```julia
using MultiJuMP, JuMP
using Ipopt

m = MultiModel(solver = IpoptSolver())
@variable(m, x[i=1:5])
@NLexpression(m, f1, sum(x[i]^2 for i=1:5))
@NLexpression(m, f2, 3x[1]+2x[2]-x[3]/3+0.01*(x[4]-x[5])^3)
@NLconstraint(m, x[1]+2x[2]-x[3]-0.5x[4]+x[5]==2)
@NLconstraint(m, 4x[1]-2x[2]+0.8x[3]+0.6x[4]+0.5x[5]^2 == 0)
@NLconstraint(m, sum(x[i]^2 for i=1:5) <= 10)

iv1 = [0.3, 0.5, -0.26, -0.13, 0.28] # Initial guess
obj1 = SingleObjective(f1, sense = :Min,
                       iv = Dict{Symbol,Any}(:x => iv1))
obj2 = SingleObjective(f2, sense = :Min)

md = getMultiData(m)
md.objectives = [obj1, obj2]
md.pointsperdim = 20
solve(m, method = :NBI) # method = :WS or method = :EPS
```

Plotting with `Plots.jl` is supported via recipes:
```julia
using Plots
pltnbi = plot(md)
```
<!-- Github bug
![Pareto front example](./img/pareto_example.svg) -->
![Pareto front example](https://cdn.rawgit.com/anriseth/MultiJuMP.jl/master/img/pareto_example.svg)
