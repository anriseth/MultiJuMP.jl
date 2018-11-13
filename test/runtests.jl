using MultiJuMP, JuMP
using Test
#using AmplNLWriter
using Ipopt
using Clp: ClpSolver

# facts("Utopia and Nadir points") do
#     # Create model
#     m = MultiModel()
#     @variable(m, -1 <= x <= 1)
#     @NLexpression(m, f1, x^2)
#     @NLexpression(m, f2, (x-0.5)^2)

#     multim = getMultiData(m)
#     multim.f1 = f1
#     multim.f2 = f2

#     @fact solve(m, method = :WS) --> :Optimal
#     @fact multim.utopia --> roughly([0.0, 0.0], 1e-5)
#     @fact multim.nadir --> roughly([0.25, 0.25], 1e-5)
# end


@testset "NBI optimisation" begin
    m = MultiModel(solver = IpoptSolver(print_level=0))
    @variable(m, x[i=1:5])
    @NLexpressions m begin
        f1, sum(x[i]^2 for i=1:5)
        f2, 3x[1]+2x[2]-x[3]/3+0.01*(x[4]-x[5])^3
    end
    @NLconstraints m begin
        x[1]+2x[2]-x[3]-0.5x[4]+x[5]==2
        4x[1]-2x[2]+0.8x[3]+0.6x[4]+0.5x[5]^2 == 0
        sum(x[i]^2 for i=1:5) <= 10
    end

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

    @test isapprox(f1arr, f1true, atol=1e-2)
    @test isapprox(f2arr, f2true, atol=1e-2)
end

@testset "WS optimisation" begin
    m = MultiModel(solver = IpoptSolver(print_level=0))
    @variable(m, x[i=1:5])
    @NLexpressions m begin
        f1, sum(x[i]^2 for i=1:5)
        f2, 3x[1]+2x[2]-x[3]/3+0.01*(x[4]-x[5])^3
    end
    @NLconstraints m begin
        x[1]+2x[2]-x[3]-0.5x[4]+x[5]==2
        4x[1]-2x[2]+0.8x[3]+0.6x[4]+0.5x[5]^2 == 0
        sum(x[i]^2 for i=1:5) <= 10
    end

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

    @test isapprox(f1arr, f1true, atol=1e-2)
    @test isapprox(f2arr, f2true, atol=1e-2)
end

@testset "EPS optimisation" begin
    m = MultiModel(solver = IpoptSolver(print_level=0))
    @variable(m, x[i=1:5])
        @NLexpressions m begin
        f1, sum(x[i]^2 for i=1:5)
        f2, 3x[1]+2x[2]-x[3]/3+0.01*(x[4]-x[5])^3
    end
    @NLconstraints m begin
        x[1]+2x[2]-x[3]-0.5x[4]+x[5]==2
        4x[1]-2x[2]+0.8x[3]+0.6x[4]+0.5x[5]^2 == 0
        sum(x[i]^2 for i=1:5) <= 10
    end

    obj1 = SingleObjective(f1)
    obj2 = SingleObjective(f2)

    multim = getMultiData(m)
    multim.objectives = [obj1, obj2]
    multim.pointsperdim = 5
    solve(m, method = :EPS)

    f1arr = convert(Array{Float64},
                    [multim.paretofront[i][1] for i in 1:multim.pointsperdim])
    f2arr = convert(Array{Float64},
                    [multim.paretofront[i][2] for i in 1:multim.pointsperdim])
    f1true = [0.555081, 10.0, 7.63877, 5.27754, 2.91631]
    f2true = [2.13057, -4.01115, -2.99356, -1.86918, -0.532781]

    @test isapprox(f1arr, f1true, atol=1e-2)
    @test isapprox(f2arr, f2true, atol=1e-2)
end


@testset "WS linear" begin
    mmodel = MultiModel(solver = ClpSolver(), linear = true)
    y = @variable(mmodel, 0 <= y <= 10.0)
    z = @variable(mmodel, 0 <= z <= 10.0)
    @constraint(mmodel, y + z <= 15.0)

    # objectives
    exp_obj1 = @expression(mmodel, -y +0.05 * z)
    exp_obj2 = @expression(mmodel, 0.05 * y - z)
    obj1 = SingleObjective(exp_obj1)
    obj2 = SingleObjective(exp_obj2)

    # setting objectives in the MultiData
    multim = getMultiData(mmodel)
    multim.objectives = [obj1, obj2]

    status = solve(mmodel, method = :WS)
    @test status == :Optimal

    true_par_vals = [-10.0, -9.75, -4.5, 0.5]
    for idx in eachindex(true_par_vals)
        tf1 = true_par_vals[idx]
        tf2 = true_par_vals[4 - idx + 1]
        @test any(multim.paretofront) do v
            v[1] ≈ tf1 && v[2] ≈ tf2 
        end
    end
end

@testset "EPS linear" begin
    mmodel = MultiModel(solver = ClpSolver(), linear = true)
    y = @variable(mmodel, 0 <= y <= 10.0)
    z = @variable(mmodel, 0 <= z <= 10.0)
    @constraint(mmodel, y + z <= 15.0)

    # objectives
    exp_obj1 = @expression(mmodel, -y +0.05 * z)
    exp_obj2 = @expression(mmodel, 0.05 * y - z)
    obj1 = SingleObjective(exp_obj1)
    obj2 = SingleObjective(exp_obj2)

    # setting objectives in the MultiData
    multim = getMultiData(mmodel)
    multim.objectives = [obj1, obj2]

    status = solve(mmodel, method = :EPS)
    @test status == :Optimal
    println(multim.paretofront)
    true_par_vals = [-10.0, -9.75, -4.5, 0.5]
    for idx in eachindex(true_par_vals)
        tf1 = true_par_vals[idx]
        tf2 = true_par_vals[4 - idx + 1]
        @test any(multim.paretofront) do v
            v[1] ≈ tf1 && v[2] ≈ tf2 
        end
    end
end

@testset "TODO: Initial value test" begin
    # TODO: Let f1 have two minima and
    # check that initial value gives correct minima
    # TODO: Add anonymous variables
    @test 1 == 1
end

@testset "TODO: NBI Inequality test" begin
    # TODO: Test that the inequality-extension of NBI works
    # We can use the nc_paper_2.jl example
    @test 1 == 1
end
