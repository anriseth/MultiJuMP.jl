using MultiJuMP, JuMP
using Base.Test
using FactCheck


facts("Utopia and Nadir points") do
    # Create model
    m = MultiModel()
    @defVar(m, -1 <= x <= 1)
    @defNLExpr(f1, x^2)
    @defNLExpr(f2, (x-0.5)^2)

    multim = getMultiData(m)
    multim.f1 = f1
    multim.f2 = f2

    @fact solve(m) --> :Optimal

    @fact multim.utopia --> roughly([0.0, 0.0], 1e-5)
    @fact multim.nadir --> roughly([0.25, 0.25], 1e-5)
end


facts("Stage 2 pareto front") do
    # TODO
    @fact 1 == 1
end
