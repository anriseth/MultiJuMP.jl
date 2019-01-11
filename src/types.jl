"""
Wrapper for single-objective function with domain type VT
"""
mutable struct SingleObjective
    f # JuMP-expression TODO: use JuMPTypes or something?
    sense::Symbol
    initialvalue::Dict{Symbol,Any} # Variable => Initial value
    # TODO: implement bound in algorithm
    bound::Float64 # Hard lower/upper depending on sense
end

SingleObjective() = SingleObjective(Any, :Min, Dict{Symbol,Any}(), NaN)
function SingleObjective(f::JuMPTypes; sense::Symbol = :Min,
                         iv::Dict{Symbol,Any} = Dict{Symbol,Any}(),
                         bound::Float64 = NaN)
    SingleObjective(f, sense, iv, bound)
end
getvalue(obj::SingleObjective) =  getvalue(obj.f)

"""
Orients the objective in the minimization sense
"""
function sensevalue(obj::SingleObjective)
    # To get Φ and Fstar with right signs
    if obj.sense == :Max
        return -getvalue(obj.f)
    else
        return getvalue(obj.f)
    end
end

"""
Stores extension data for multi-objective
optimization inside JuMP Model
"""
mutable struct MultiData{Tx,To}
    objectives::Vector{SingleObjective}
    f1::SingleObjective
    f2::SingleObjective
    normalf1
    normalf2

    linear::Bool
    # Number of points in each direction
    # of Pareto submanifold
    pointsperdim::Int
    #
    # stored values
    utopiavarvalues::Vector{Vector{Tx}}
    utopia::Vector{To}
    nadir::Vector{Tx}
    paretovarvalues::Vector{Vector{To}}
    paretofront

    Phi::Array{To,2}
    Fstar::Vector{To}
end


function multi_model(;solver=nothing, linear=false)
    m = Model(solver=solver)
    m.solvehook = solvehook

    m.ext[:Multi] = MultiData(Array{SingleObjective}(undef, 0),
                              SingleObjective(), SingleObjective(),
                              Any, Any, linear,
                              10,
                              Vector{Float64}[], Float64[],
                              Float64[], Vector{Float64}[], Any[],
                              Array{Float64}(undef, 2,2), Array{Float64}(undef, 2))
    return m
end

function multi_model(m::JuMP.AbstractModel; linear=false)
    m.solvehook = solvehook
    m.ext[:Multi] = MultiData(Array{SingleObjective}(undef, 0),
                              SingleObjective(), SingleObjective(),
                              Any, Any, linear,
                              10,
                              Vector{Float64}[], Float64[],
                              Float64[], Vector{Float64}[], Any[],
                              Array{Float64}(undef, 2,2), Array{Float64}(undef, 2))
    return m
end

"""
Gets the MultiData struct from a model if it exists
"""
function get_multidata(m::Model)
    if haskey(m.ext, :Multi)
        return m.ext[:Multi]::MultiData
    else
        error("This functionality is only available for MultiJuMP models")
    end
end

"""
Sets up the tree of all possible
convex combinations of the objectives
"""
function betas(levels, parts)
    map(Combinatorics.combinations(1:(levels+parts-1),(levels-1))) do x
        ([x;levels+parts].-[0;x] .-1)/parts
    end
end
