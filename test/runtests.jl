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
const optimizer = Mosek.Optimizer()
MOI.set(optimizer, MOI.RawOptimizerAttribute("fallback"), FALLBACK_URL)
MOI.set(optimizer, MOI.Silent(), true)

@testset "SolverName" begin
    @test MOI.get(optimizer, MOI.SolverName()) == "Mosek"
end

@testset "Parameters" begin
    optimizer = Mosek.Optimizer()
    MOI.set(optimizer, MOI.RawOptimizerAttribute("fallback"), FALLBACK_URL)
    @testset "Double Parameter" begin
        MOI.set(optimizer, MOI.RawOptimizerAttribute("INTPNT_CO_TOL_DFEAS"), 1e-7)
        @test MOI.get(optimizer, MOI.RawOptimizerAttribute("MSK_DPAR_INTPNT_CO_TOL_DFEAS")) == 1e-7
        MOI.set(optimizer, MOI.RawOptimizerAttribute("MSK_DPAR_INTPNT_CO_TOL_DFEAS"), 1e-8)
        @test MOI.get(optimizer, MOI.RawOptimizerAttribute("MSK_DPAR_INTPNT_CO_TOL_DFEAS")) == 1e-8
        @testset "with integer value" begin
            MOI.set(optimizer, MOI.RawOptimizerAttribute("MSK_DPAR_INTPNT_CO_TOL_DFEAS"), 1)
            @test MOI.get(optimizer, MOI.RawOptimizerAttribute("MSK_DPAR_INTPNT_CO_TOL_DFEAS")) == 1
        end
    end
    @testset "Integer Parameter" begin
        MOI.set(optimizer, MOI.RawOptimizerAttribute("MSK_IPAR_INTPNT_MAX_ITERATIONS"), 100)
        @test MOI.get(optimizer, MOI.RawOptimizerAttribute("MSK_IPAR_INTPNT_MAX_ITERATIONS")) == 100
        MOI.set(optimizer, MOI.RawOptimizerAttribute("INTPNT_MAX_ITERATIONS"), 200)
        @test MOI.get(optimizer, MOI.RawOptimizerAttribute("MSK_IPAR_INTPNT_MAX_ITERATIONS")) == 200
        @testset "with enum value" begin
            MOI.set(optimizer, MOI.RawOptimizerAttribute("MSK_IPAR_OPTIMIZER"), MosekTools.Mosek.MSK_OPTIMIZER_DUAL_SIMPLEX)
            @test MOI.get(optimizer, MOI.RawOptimizerAttribute("MSK_IPAR_OPTIMIZER")) == convert(Int32, MosekTools.Mosek.MSK_OPTIMIZER_DUAL_SIMPLEX)
        end
    end
    @testset "String Parameter" begin
        MOI.set(optimizer, MOI.RawOptimizerAttribute("PARAM_WRITE_FILE_NAME"), "foo.txt")
        # Needs https://github.com/JuliaOpt/Mosek.jl/pull/174
        #@test MOI.get(optimizer, MOI.RawOptimizerAttribute("MSK_SPAR_PARAM_WRITE_FILE_NAME")) == "foo.txt"
        MOI.set(optimizer, MOI.RawOptimizerAttribute("MSK_SPAR_PARAM_WRITE_FILE_NAME"), "bar.txt")
        #@test MOI.get(optimizer, MOI.RawOptimizerAttribute("MSK_SPAR_PARAM_WRITE_FILE_NAME")) == "bar.txt"
    end
    @testset "TimeLimitSec" begin
        @test MOI.get(optimizer, MOI.RawOptimizerAttribute("MSK_DPAR_OPTIMIZER_MAX_TIME")) == -1
        @test MOI.get(optimizer, MOI.TimeLimitSec()) === nothing
        MOI.set(optimizer, MOI.TimeLimitSec(), 1.0)
        @test MOI.get(optimizer, MOI.RawOptimizerAttribute("MSK_DPAR_OPTIMIZER_MAX_TIME")) == 1.0
        @test MOI.get(optimizer, MOI.TimeLimitSec()) === 1.0
        MOI.set(optimizer, MOI.TimeLimitSec(), nothing)
        @test MOI.get(optimizer, MOI.RawOptimizerAttribute("MSK_DPAR_OPTIMIZER_MAX_TIME")) == -1
        @test MOI.get(optimizer, MOI.TimeLimitSec()) === nothing
    end
end

@testset "supports_default_copy_to" begin
    @test MOIU.supports_default_copy_to(optimizer, false)
    @test MOIU.supports_default_copy_to(optimizer, true)
end

const config = MOIT.Config(
    Float64, atol=1e-3, rtol=1e-3,
    exclude=Any[MOI.ConstraintName, MOI.VariableBasisStatus, MOI.ConstraintBasisStatus], # result in errors for now
)

@testset "Basic and linear tests" begin
    MOIT.runtests(optimizer, config, include=["basic", "linear", "quadratic"],
        exclude=["Indicator", "Cone", "conic"],
    )
end

@testset "Conic problems" begin
    MOIT.runtests(optimizer, config,
        include=["conic", "SecondOrderCone", "Semidefinite", "Exponential", "PowerCone", "Cone"],
        exclude=["Indicator", "basic", "linear"],
    )
end

@testset "Integer problems" begin
    MOIT.runtests(optimizer, config,
        include=["Integer", "ZeroOne"],
        exclude=[
            "Indicator", "basic", "linear", "conic", "SecondOrderCone", "Semidefinite",
            "test_variable_solve_ZeroOne_with_0_upper_bound", "test_variable_solve_ZeroOne_with_upper_bound", # issues/74
        ],
    )
end

@testset "Bridged" begin
    model = MOIB.full_bridge_optimizer(Mosek.Optimizer(), Float64)
    MOI.set(model, MOI.Silent(), true)

    # linear and basic tests
    MOIT.runtests(model, config,
        include=["basic", "linear", "quadratic"],
        exclude=[
            "Cone", "conic",
            "test_basic_ScalarQuadraticFunction_EqualTo", # non-PSD quadratic
            "test_basic_ScalarQuadraticFunction_GreaterThan",
            "test_basic_ScalarQuadraticFunction_Interval",
            "test_basic_ScalarQuadraticFunction_Semi",
            "test_basic_ScalarQuadraticFunction_ZeroOne",
            "test_basic_ScalarQuadraticFunction_Integer",
            "test_basic_VectorQuadraticFunction_Nonnegatives",
            "test_basic_VectorQuadraticFunction_Zeros",
            "nonconvex",
        ],
    )

    # conic
    MOIT.runtests(
        model, config,
        include=["conic", "SecondOrderCone", "Semidefinite", "Exponential", "PowerCone", "Cone"],
        exclude=[
            "test_basic_VectorAffineFunction_PositiveSemidefiniteConeSquare", # AssertionError: (m.x_sd[ref2id(vi)]).matrix == -1 src/variable.jl:173
            "test_basic_VectorOfVariables_PositiveSemidefiniteConeSquare",
            "test_basic_VectorAffineFunction_LogDetConeTriangle",
            "test_basic_VectorAffineFunction_NormNuclearCone",
            "test_basic_VectorOfVariables_NormNuclearCone",
            "test_basic_VectorAffineFunction_NormSpectralCone",
            "test_basic_VectorOfVariables_NormSpectralCone",
            "test_basic_VectorAffineFunction_RootDetConeTriangle",
            "test_basic_VectorOfVariables_RootDetConeTriangle",
            "test_basic_VectorAffineFunction_PositiveSemidefiniteConeTriangle", # TODO: implement get ConstraintSet for SAF
            "test_basic_VectorOfVariables_PositiveSemidefiniteConeTriangle",
            "test_basic_VectorQuadraticFunction_", # not PSD because of equality
            "test_quadratic_SecondOrderCone_basic",
            "test_basic_VectorOfVariables_LogDetConeTriangle", # Mosek.MosekError(1307, "Variable '' (1) is a member of cone '' (0).") src/msk_functions.jl:477
            "test_conic_LogDetConeTriangle_VectorOfVariables",
        ],
    )

    # integer
    MOIT.runtests(model, config,
        include=["Integer", "ZeroOne"],
        exclude=[
            "Cone", "conic",
            "test_constraint_ZeroOne_bounds", # Cannot put multiple bound sets of the same type on a variable
            "test_variable_solve_ZeroOne_with_0_upper_bound",
            "test_variable_solve_ZeroOne_with_upper_bound",
            "test_basic_ScalarQuadraticFunction_ZeroOne", # non-PSD quadratic
        ],
    )

    # other attribute tests
    MOIT.runtests(model, config,
        exclude=[
            "Cone", "conic", "Integer", "ZeroOne", "basic", "linear",
            "test_model_ListOfConstraintAttributesSet", # list not properly set
            "BoundAlreadySet", # TODO throw error if bound already set
            "test_model_ModelFilter_AbstractVariableAttribute",
            "test_model_VariableName", # Mosek currently throws when setting twice, not when getting names
            "test_model_Name_VariableName_ConstraintName",
            "test_model_duplicate_VariableName",
            "test_model_VariablePrimalStart", # able to set but not to get VariablePrimalStart
            "test_objective_qp_ObjectiveFunction_zero_ofdiag", # MOI.ListOfModelAttributesSet
            "test_objective_set_via_modify",
            "test_quadratic_nonconvex_constraint_integration",
            "test_solve_ObjectiveBound_MAX_SENSE_LP", # ObjectiveBound invalid

        ],
    )
    @test MOI.supports(model, MOI.VariablePrimalStart(), MOI.VariableIndex)
end
