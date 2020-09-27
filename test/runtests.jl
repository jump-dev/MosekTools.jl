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
        MOI.set(optimizer, MOI.RawParameter("INTPNT_CO_TOL_DFEAS"), 1e-7)
        @test MOI.get(optimizer, MOI.RawParameter("MSK_DPAR_INTPNT_CO_TOL_DFEAS")) == 1e-7
        MOI.set(optimizer, MOI.RawParameter("MSK_DPAR_INTPNT_CO_TOL_DFEAS"), 1e-8)
        @test MOI.get(optimizer, MOI.RawParameter("MSK_DPAR_INTPNT_CO_TOL_DFEAS")) == 1e-8
        @testset "with integer value" begin
            MOI.set(optimizer, MOI.RawParameter("MSK_DPAR_INTPNT_CO_TOL_DFEAS"), 1)
            @test MOI.get(optimizer, MOI.RawParameter("MSK_DPAR_INTPNT_CO_TOL_DFEAS")) == 1
        end
    end
    @testset "Integer Parameter" begin
        MOI.set(optimizer, MOI.RawParameter("MSK_IPAR_INTPNT_MAX_ITERATIONS"), 100)
        @test MOI.get(optimizer, MOI.RawParameter("MSK_IPAR_INTPNT_MAX_ITERATIONS")) == 100
        MOI.set(optimizer, MOI.RawParameter("INTPNT_MAX_ITERATIONS"), 200)
        @test MOI.get(optimizer, MOI.RawParameter("MSK_IPAR_INTPNT_MAX_ITERATIONS")) == 200
        @testset "with enum value" begin
            MOI.set(optimizer, MOI.RawParameter("MSK_IPAR_OPTIMIZER"), MosekTools.Mosek.MSK_OPTIMIZER_DUAL_SIMPLEX)
            @test MOI.get(optimizer, MOI.RawParameter("MSK_IPAR_OPTIMIZER")) == convert(Int32, MosekTools.Mosek.MSK_OPTIMIZER_DUAL_SIMPLEX)
        end
    end
    @testset "String Parameter" begin
        MOI.set(optimizer, MOI.RawParameter("PARAM_WRITE_FILE_NAME"), "foo.txt")
        # Needs https://github.com/JuliaOpt/Mosek.jl/pull/174
        #@test MOI.get(optimizer, MOI.RawParameter("MSK_SPAR_PARAM_WRITE_FILE_NAME")) == "foo.txt"
        MOI.set(optimizer, MOI.RawParameter("MSK_SPAR_PARAM_WRITE_FILE_NAME"), "bar.txt")
        #@test MOI.get(optimizer, MOI.RawParameter("MSK_SPAR_PARAM_WRITE_FILE_NAME")) == "bar.txt"
    end
    @testset "TimeLimitSec" begin
        @test MOI.get(optimizer, MOI.RawParameter("MSK_DPAR_OPTIMIZER_MAX_TIME")) == -1
        @test MOI.get(optimizer, MOI.TimeLimitSec()) === nothing
        MOI.set(optimizer, MOI.TimeLimitSec(), 1.0)
        @test MOI.get(optimizer, MOI.RawParameter("MSK_DPAR_OPTIMIZER_MAX_TIME")) == 1.0
        @test MOI.get(optimizer, MOI.TimeLimitSec()) === 1.0
        MOI.set(optimizer, MOI.TimeLimitSec(), nothing)
        @test MOI.get(optimizer, MOI.RawParameter("MSK_DPAR_OPTIMIZER_MAX_TIME")) == -1
        @test MOI.get(optimizer, MOI.TimeLimitSec()) === nothing
    end
end

@testset "supports_default_copy_to" begin
    @test MOIU.supports_default_copy_to(optimizer, false)
    @test MOIU.supports_default_copy_to(optimizer, true)
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
    @testset "Conic" begin
        MOIT.basic_constraint_tests(
            optimizer, config,
            include=[
                (MOI.VectorOfVariables, MOI.ExponentialCone),
                (MOI.VectorOfVariables, MOI.DualExponentialCone),
                (MOI.VectorOfVariables, MOI.PowerCone{Float64}),
                (MOI.VectorOfVariables, MOI.DualPowerCone{Float64}),
                (MOI.VectorOfVariables, MOI.SecondOrderCone),
                (MOI.VectorOfVariables, MOI.RotatedSecondOrderCone)
        ])
    end
end

const bridged = MOIB.full_bridge_optimizer(optimizer, Float64)

# Mosek errors during `MOI.set` instead of `MOI.get` when there are duplicates.
#@testset "Name" begin
#    MOIT.nametest(bridged)
#end

@testset "Copy" begin
    model = MOIB.full_bridge_optimizer(Mosek.Optimizer(), Float64)
    MOIT.copytest(bridged, model)
end

@testset "Start" begin
    # TODO this should be checked somewhere in MOI
    @test MOI.supports(bridged, MOI.VariablePrimalStart(), MOI.VariableIndex)
end

@testset "Unit" begin
    # Mosek does not support names
    MOIT.unittest(bridged, config, [
        # TODO
        "number_threads",
        # Find objective bound of 0.0 which is lower than 4.0
        "solve_objbound_edge_cases",
        # Cannot put multiple bound sets of the same type on a variable
        "solve_integer_edge_cases",
        # Cannot mix `ZeroOne` with `GreaterThan`/`LessThan`
        "solve_zero_one_with_bounds_1",
        "solve_zero_one_with_bounds_2",
        "solve_zero_one_with_bounds_3"])
end

@testset "Continuous Linear" begin
    MOIT.contlineartest(bridged, config)
end

@testset "Continuous Quadratic" begin
    MOIT.contquadratictest(bridged, config, [
        # Non-convex
        "ncqcp",
        # QuadtoSOC does not work as the matrix is not SDP
        "socp"
    ])
end

@testset "Continuous Conic" begin
    MOIT.contconictest(bridged, config, ["rootdets", "logdet"])
end

@testset "Integer Linear" begin
    MOIT.intlineartest(bridged, config, [
        # SOS constraints not supported:
        # int2 uses SOS and indicator can be bridged to SOS
        "int2", "indicator1", "indicator2", "indicator3", "indicator4"
    ])
end
