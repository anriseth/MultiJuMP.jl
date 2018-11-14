# extensions to JuMP for multiobjective optimisation

module MultiJuMP

using JuMP
import JuMP: JuMPTypes, getvalue
import MathProgBase
using RecipesBase, LaTeXStrings
import Combinatorics.combinations
using LinearAlgebra
export MultiModel, SingleObjective, getMultiData

include("types.jl")
include("linear.jl")
include("methods.jl")
include("plot_recipe.jl")

end
