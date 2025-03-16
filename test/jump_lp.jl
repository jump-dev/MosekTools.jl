function test_jump_lp(solver)
    @testset "Linear Programming" begin
        @testset "LP1" begin
            # simple 2 variable, 1 constraint problem
            # min -x
            # st   x + y <= 1   (x + y - 1 ∈ Nonpositives)
            #       x, y >= 0   (x, y ∈ Nonnegatives)

            m = Model(; solver = solver)
            @variable(m, x >= 0.0)
            @variable(m, y >= 0.0)
            @objective(m, Min, -x)

            @constraint(m, x + y <= 1)

            JuMP.solve(m)

            @test JuMP.isattached(m)
            @test JuMP.hasvariableresult(m)

            @test JuMP.terminationstatus(m) == MOI.Success
            @test JuMP.primalstatus(m) == MOI.FeasiblePoint

            @test JuMP.resultvalue(x) ≈ 1.0 atol = 1e-6
            @test JuMP.resultvalue(y) ≈ 0.0 atol = 1e-6
            @test JuMP.objectivevalue(m) ≈ -1.0 atol = 1e-6
        end
    end
end
