using MultiJuMP, JuMP
using Base.Test
using FactCheck
#using AmplNLWriter
using Ipopt

# facts("Utopia and Nadir points") do
#     # Create model
#     m = MultiModel()
#     @defVar(m, -1 <= x <= 1)
#     @defNLExpr(m, f1, x^2)
#     @defNLExpr(m, f2, (x-0.5)^2)

#     multim = getMultiData(m)
#     multim.f1 = f1
#     multim.f2 = f2

#     @fact solve(m, method = :WS) --> :Optimal
#     @fact multim.utopia --> roughly([0.0, 0.0], 1e-5)
#     @fact multim.nadir --> roughly([0.25, 0.25], 1e-5)
# end


facts("NBI optimisation") do
    # TODO: do for WS as well
    m = MultiModel(solver = IpoptSolver(print_level=0))
    @defVar(m, x[i=1:5])
    @defNLExpr(m, f1, sum{x[i]^2, i=1:5})
    @defNLExpr(m, f2, 3x[1]+2x[2]-x[3]/3+0.01*(x[4]-x[5])^3)
    @addNLConstraint(m, x[1]+2x[2]-x[3]-0.5x[4]+x[5]==2)
    @addNLConstraint(m, 4x[1]-2x[2]+0.8x[3]+0.6x[4]+0.5x[5]^2 == 0)
    @addNLConstraint(m, sum{x[i]^2, i=1:5} <= 10)

    obj1 = SingleObjective(f1)
    obj2 = SingleObjective(f2)

    multim = getMultiData(m)
    multim.objectives = [obj1, obj2]
    multim.pointsperdim = 5
    solve(m, method = :NBI)

    f1arr = convert(Array{Float64},
                    [multim.paretofront[i][1] for i in 1:multim.pointsperdim])
    f2arr = convert(Array{Float64},
                    [multim.paretofront[i][2] for i in 1:multim.pointsperdim])
    f1true = [0.555081, 10.0, 7.16982, 4.48665, 2.08008]
    f2true = [2.13057, -4.01115, -2.78066, -1.45458, 0.0513689]

    @fact f1arr --> roughly(f1true, 1e-5)
    @fact f2arr --> roughly(f2true, 1e-5)
end

facts("WS optimisation") do
    # TODO: do for EPS as well
    m = MultiModel(solver = IpoptSolver(print_level=0))
    @defVar(m, x[i=1:5])
    @defNLExpr(m, f1, sum{x[i]^2, i=1:5})
    @defNLExpr(m, f2, 3x[1]+2x[2]-x[3]/3+0.01*(x[4]-x[5])^3)
    @addNLConstraint(m, x[1]+2x[2]-x[3]-0.5x[4]+x[5]==2)
    @addNLConstraint(m, 4x[1]-2x[2]+0.8x[3]+0.6x[4]+0.5x[5]^2 == 0)
    @addNLConstraint(m, sum{x[i]^2, i=1:5} <= 10)

    obj1 = SingleObjective(f1)
    obj2 = SingleObjective(f2)

    multim = getMultiData(m)
    multim.objectives = [obj1, obj2]
    multim.pointsperdim = 5
    solve(m, method = :WS)

    f1arr = convert(Array{Float64},
                    [multim.paretofront[i][1] for i in 1:multim.pointsperdim])
    f2arr = convert(Array{Float64},
                    [multim.paretofront[i][2] for i in 1:multim.pointsperdim])
    f1true = [0.555081, 10.0, 10.0, 2.88267, 0.726877]
    f2true = [2.13057, -4.01115, -4.01115, -0.510963, 1.48613]

    @fact f1arr --> roughly(f1true, 1e-5)
    @fact f2arr --> roughly(f2true, 1e-5)
end



facts("TODO: Initial value test") do
    # TODO: Let f1 have two minima and
    # check that initial value gives correct minima
    @fact 1 --> 1
end
