# To ensure that tests can run on travis we have to do a little
# hackadoodle here. The tests require a license file. We include
# a license file that is only valid for one day (the day when
# change is submitted).
# If there is no valid license file, we default to that file.

using Test

using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

const FALLBACK_URL = "mosek://solve.mosek.com:30080"

using MosekTools
const optimizer = Mosek.Optimizer(fallback = FALLBACK_URL)
MOI.set(optimizer, MOI.Silent(), true)

@testset "SolverName" begin
    @test MOI.get(optimizer, MOI.SolverName()) == "Mosek"
end

@testset "Parameters" begin
    optimizer = Mosek.Optimizer(fallback = FALLBACK_URL)
    @testset "Double Parameter" begin
        MOI.set(optimizer, MosekTools.Parameter("INTPNT_CO_TOL_DFEAS"), 1e-7)
        @test MOI.get(optimizer, MosekTools.Parameter("MSK_DPAR_INTPNT_CO_TOL_DFEAS")) == 1e-7
        MOI.set(optimizer, MosekTools.Parameter("MSK_DPAR_INTPNT_CO_TOL_DFEAS"), 1e-8)
        @test MOI.get(optimizer, MosekTools.Parameter("MSK_DPAR_INTPNT_CO_TOL_DFEAS")) == 1e-8
        @testset "with integer value" begin
            MOI.set(optimizer, MosekTools.Parameter("MSK_DPAR_INTPNT_CO_TOL_DFEAS"), 1)
            @test MOI.get(optimizer, MosekTools.Parameter("MSK_DPAR_INTPNT_CO_TOL_DFEAS")) == 1
        end
    end
    @testset "Integer Parameter" begin
        MOI.set(optimizer, MosekTools.Parameter("MSK_IPAR_INTPNT_MAX_ITERATIONS"), 100)
        @test MOI.get(optimizer, MosekTools.Parameter("MSK_IPAR_INTPNT_MAX_ITERATIONS")) == 100
        MOI.set(optimizer, MosekTools.Parameter("INTPNT_MAX_ITERATIONS"), 200)
        @test MOI.get(optimizer, MosekTools.Parameter("MSK_IPAR_INTPNT_MAX_ITERATIONS")) == 200
        @testset "with enum value" begin
            MOI.set(optimizer, MosekTools.Parameter("MSK_IPAR_OPTIMIZER"), MosekTools.Mosek.MSK_OPTIMIZER_DUAL_SIMPLEX)
            @test MOI.get(optimizer, MosekTools.Parameter("MSK_IPAR_OPTIMIZER")) == convert(Int32, MosekTools.Mosek.MSK_OPTIMIZER_DUAL_SIMPLEX)
        end
    end
    @testset "String Parameter" begin
        MOI.set(optimizer, MosekTools.Parameter("PARAM_WRITE_FILE_NAME"), "foo.txt")
        # Needs https://github.com/JuliaOpt/Mosek.jl/pull/174
        #@test MOI.get(optimizer, MosekTools.Parameter("MSK_SPAR_PARAM_WRITE_FILE_NAME")) == "foo.txt"
        MOI.set(optimizer, MosekTools.Parameter("MSK_SPAR_PARAM_WRITE_FILE_NAME"), "bar.txt")
        #@test MOI.get(optimizer, MosekTools.Parameter("MSK_SPAR_PARAM_WRITE_FILE_NAME")) == "bar.txt"
    end
end

@testset "supports_default_copy_to" begin
    @test MOIU.supports_default_copy_to(optimizer, false)
    @test MOIU.supports_default_copy_to(optimizer, true)
end

@testset "Name" begin
    MOIT.nametest(optimizer)
    MOI.empty!(optimizer)
end

const config = MOIT.TestConfig(atol=1e-3, rtol=1e-3)

@testset "Basic" begin
    @testset "Linear" begin
        MOIT.basic_constraint_tests(
            optimizer, config,
            include=[
                (MOI.SingleVariable, MOI.EqualTo{Float64}),
                (MOI.SingleVariable, MOI.LessThan{Float64}),
                (MOI.SingleVariable, MOI.GreaterThan{Float64}),
                (MOI.SingleVariable, MOI.Interval{Float64}),
                (MOI.SingleVariable, MOI.ZeroOne),
                (MOI.SingleVariable, MOI.Integer),
                (MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64}),
                (MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}),
                (MOI.ScalarAffineFunction{Float64}, MOI.GreaterThan{Float64}),
                (MOI.ScalarAffineFunction{Float64}, MOI.Interval{Float64})
        ])
    end
#    @testset "Conic" begin
#        # Doesn't work, see https://github.com/JuliaOpt/MathOptInterface.jl/pull/703
#        MOIT.basic_constraint_tests(
#            optimizer, config,
#            delete = false, # TODO
#            get_constraint_function = false, # TODO
#            get_constraint_set = false, # TODO
#            include=[
#                (MOI.VectorOfVariables, MOI.SecondOrderCone),
#                (MOI.VectorOfVariables, MOI.RotatedSecondOrderCone)
#        ])
#    end
end

const bridged = MOIB.full_bridge_optimizer(optimizer, Float64)

@testset "Copy" begin
    model = MOIB.full_bridge_optimizer(Mosek.Optimizer(), Float64)
    MOIT.copytest(bridged, model)
end

@testset "Unit" begin
    # Mosek does not support names
    MOIT.unittest(bridged,
                  config,
                  [# Does not support quadratic objective yet, needs
                   # https://github.com/JuliaOpt/MathOptInterface.jl/issues/529
                   "solve_qp_edge_cases",
                   # Find objective bound of 0.0 which is lower than 4.0
                   "solve_objbound_edge_cases",
                   # Cannot put multiple bound sets of the same type on a variable
                   "solve_integer_edge_cases"
                  ])
end

@testset "Continuous Linear" begin
    MOIT.contlineartest(bridged, config)
end

@testset "Continuous Quadratic" begin
    MOIT.qcptest(bridged, config)
end

@testset "Continuous Conic" begin
    MOIT.contconictest(bridged, config, ["exp", "rootdets", "logdet"])
end

@testset "Integer Linear" begin
    MOIT.intlineartest(optimizer, config, ["int2"])
end

# Test that objective and constraint data are copied over correctly when
# a scalar variable is transformed to a matrix one
@testset "SDP $add_before" for add_before in [true, false]
    atol = 1e-4
    rtol = 1e-3
    # Problem SDP1 - sdo1 from MOSEK docs
    # See sdp1 of MOI contconic tests for how to get the analytical solution
    δ = √(1 + (3*√2+2)*√(-116*√2+166) / 14) / 2
    ε = √((1 - 2*(√2-1)*δ^2) / (2-√2))
    y2 = 1 - ε*δ
    y1 = 1 - √2*y2
    obj = y1 + y2/2
    k = -2*δ/ε
    x2 = ((3-2obj)*(2+k^2)-4) / (4*(2+k^2)-4*√2)
    α = √(3-2obj-4x2)/2
    β = k*α

    MOI.empty!(bridged)

    X = MOI.add_variables(bridged, 6)
    x = MOI.add_variables(bridged, 3)

    cx = MOI.add_constraint(bridged, MOI.VectorOfVariables(x), MOI.SecondOrderCone(3))
    if add_before
        cX = MOI.add_constraint(bridged, MOI.VectorOfVariables(X), MOI.PositiveSemidefiniteConeTriangle(3))
    end

    c1 = MOI.add_constraint(bridged, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1., 1, 1, 1], [X[1], X[3], X[end], x[1]]), 0.), MOI.EqualTo(1.))
    c2 = MOI.add_constraint(bridged, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1., 2, 1, 2, 2, 1, 1, 1], [X; x[2]; x[3]]), 0.), MOI.EqualTo(1/2))

    objXidx = [1:3; 5:6]
    objXcoefs = 2*ones(5)
    MOI.set(bridged, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
            MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([objXcoefs; 1.0], [X[objXidx]; x[1]]), 0.0))
    MOI.set(bridged, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    if !add_before
        cX = MOI.add_constraint(bridged, MOI.VectorOfVariables(X), MOI.PositiveSemidefiniteConeTriangle(3))
    end

    @test MOI.get(bridged, MOI.TerminationStatus()) == MOI.OPTIMIZE_NOT_CALLED
    MOI.optimize!(bridged)
    @test MOI.get(bridged, MOI.TerminationStatus()) == MOI.OPTIMAL

    @test MOI.get(bridged, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
    @test MOI.get(bridged, MOI.DualStatus()) == MOI.FEASIBLE_POINT

    @test MOI.get(bridged, MOI.ObjectiveValue()) ≈ obj atol=atol rtol=rtol

    Xv = [α^2, α*β, β^2, α^2, α*β, α^2]
    xv = [√2*x2, x2, x2]
    @test MOI.get(bridged, MOI.VariablePrimal(), X) ≈ Xv atol=atol rtol=rtol
    @test MOI.get(bridged, MOI.VariablePrimal(), x) ≈ xv atol=atol rtol=rtol

    @test MOI.get(bridged, MOI.ConstraintPrimal(), cX) ≈ Xv atol=atol rtol=rtol
    @test MOI.get(bridged, MOI.ConstraintPrimal(), cx) ≈ xv atol=atol rtol=rtol
    @test MOI.get(bridged, MOI.ConstraintPrimal(), c1) ≈ 1. atol=atol rtol=rtol
    @test MOI.get(bridged, MOI.ConstraintPrimal(), c2) ≈ .5 atol=atol rtol=rtol

    cX0 = 1+(√2-1)*y2
    cX1 = 1-y2
    cX2 = -y2
    cXv = [cX0, cX1, cX0, cX2, cX1, cX0]
    @test MOI.get(bridged, MOI.ConstraintDual(), cX) ≈ cXv atol=atol rtol=rtol
    @test MOI.get(bridged, MOI.ConstraintDual(), cx) ≈ [1-y1, -y2, -y2] atol=atol rtol=rtol
    @test MOI.get(bridged, MOI.ConstraintDual(), c1) ≈ y1 atol=atol rtol=rtol
    @test MOI.get(bridged, MOI.ConstraintDual(), c2) ≈ y2 atol=atol rtol=rtol

    @test MosekTools.Mosek.getnumvar(optimizer.task) == (add_before ? 3 : 9)
end

#include("test_jump.jl")
