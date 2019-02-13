# extensions to JuMP for multiobjective optimisation

module MultiJuMP

using JuMP
import JuMP: JuMPTypes, getvalue
import MathProgBase
using RecipesBase, LaTeXStrings
import Combinatorics
using LinearAlgebra
export multi_model, SingleObjective, get_multidata
export getutopia, getnadir
export WeightedSum, EpsilonCons, NBI

include("types.jl")
include("methods.jl")
include("linear.jl")
include("nonlinear.jl")
include("plot_recipe.jl")

end
