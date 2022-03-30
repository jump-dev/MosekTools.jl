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

function MosekOptimizerWithFallback()
    optimizer = Mosek.Optimizer()
    MOI.set(optimizer, MOI.RawOptimizerAttribute("fallback"), FALLBACK_URL)
    MOI.set(optimizer, MOI.Silent(), true)
    return optimizer
end


@testset "SolverName" begin
    @test MOI.get(optimizer, MOI.SolverName()) == "Mosek"
end

@testset "Parameters" begin
    optimizer = MosekOptimizerWithFallback()
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

@testset "supports_incremental_interface" begin
    @test MOI.supports_incremental_interface(optimizer)
end

const config = MOIT.Config(
    Float64, atol=1e-3, rtol=1e-3,
    exclude=Any[MOI.ConstraintName, MOI.VariableBasisStatus, MOI.ConstraintBasisStatus], # result in errors for now
)

@testset "Direct optimizer tests" begin
    MOIT.runtests(optimizer, config,
        exclude=[
            # FIXME
            # Expression: MOI.add_constraint(model, x, set2)
            #   Expected: MathOptInterface.LowerBoundAlreadySet{MathOptInterface.EqualTo{Float64}, MathOptInterface.GreaterThan{Float64}}(MathOptInterface.VariableIndex(1))
            #     Thrown: ErrorException("Cannot put multiple bound sets of the same type on a variable")
            "test_model_LowerBoundAlreadySet",
            "test_model_UpperBoundAlreadySet",
            # FIXME ArgumentError: MosekTools.Optimizer does not support getting the attribute MathOptInterface.VariablePrimalStart().
            "test_model_VariablePrimalStart",
            # FIXME
            "test_model_duplicate_VariableName",
            # FIXME `MOI.ListOfConstraintAttributesSet` incorrect
            "test_model_ListOfConstraintAttributesSet",
            #  Expression: MOI.set(model, MOI.ConstraintName(), c, "c1")
            #    Expected: MathOptInterface.UnsupportedAttribute{MathOptInterface.ConstraintName}(MathOptInterface.ConstraintName(), "`ConstraintName`s are not supported for `VariableIndex` constraints.")
            #  No exception thrown
            "test_model_VariableIndex_ConstraintName",
            # FIXME Needs https://github.com/jump-dev/MathOptInterface.jl/pull/1787
            "test_variable_solve_ZeroOne_with_0_upper_bound",
        ],
    )
end

@testset "Bridged and cached" begin
    #model = MOIB.full_bridge_optimizer(Mosek.Optimizer(), Float64)
    model = MOIB.full_bridge_optimizer(MosekOptimizerWithFallback(), Float64)
    MOI.set(model, MOI.Silent(), true)

    # linear and basic tests
    MOIT.runtests(model, config,
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
            "test_conic_LogDetConeTriangle_VectorOfVariables",
            "test_variable_solve_ZeroOne_with_0_upper_bound",
            "test_variable_solve_ZeroOne_with_upper_bound",
            "test_model_ListOfConstraintAttributesSet", # list not properly set
            "BoundAlreadySet", # TODO throw error if bound already set
            "test_model_duplicate_VariableName",
            "test_model_VariablePrimalStart", # able to set but not to get VariablePrimalStart
            "test_objective_set_via_modify",
            # FIXME Mosek.MosekError(1307, "Variable '' (9) is a member of cone '' (0).")
            "test_basic_VectorQuadraticFunction_LogDetConeTriangle",
            "test_basic_VectorOfVariables_LogDetConeTriangle",
            # FIXME Needs https://github.com/jump-dev/MathOptInterface.jl/pull/1787
            "test_constraint_ZeroOne_bounds",
            "test_constraint_ZeroOne_bounds_2",
            "test_constraint_ZeroOne_bounds_3",
            "test_variable_solve_ZeroOne_with_0_upper_bound",
            "test_variable_solve_ZeroOne_with_upper_bound",
        ],
    )

    @test MOI.supports(model, MOI.VariablePrimalStart(), MOI.VariableIndex)

    model = MOI.Bridges.full_bridge_optimizer(
        MOI.Utilities.CachingOptimizer(
                MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
                #Mosek.Optimizer(),
                MosekOptimizerWithFallback(),
        ),
        Float64,
    )
    MOI.set(model, MOI.Silent(), true)

    MOI.Test.runtests(model, config,
        exclude=[
            # FIXME Mosek.MosekError(1307, "Variable '' (1) is a member of cone '' (0).") src/msk_functions.jl:477
            "test_conic_LogDetConeTriangle_VectorOfVariables",
            # FIXME Needs https://github.com/jump-dev/MathOptInterface.jl/pull/1787
            "test_constraint_ZeroOne_bounds",
            "test_variable_solve_ZeroOne_with_0_upper_bound",
            "test_variable_solve_ZeroOne_with_upper_bound",
        ],
    )
end

@testset "Matrix name" begin
    MOI.empty!(optimizer)
    x, cx = MOI.add_constrained_variables(optimizer, MOI.PositiveSemidefiniteConeTriangle(3))
    err = MOI.UnsupportedAttribute{MOI.VariableName}
    @test_throws err MOI.set(optimizer, MOI.VariableName(), x[1], "a")
    MOI.empty!(optimizer)
    model = MOI.Utilities.CachingOptimizer(MOI.Utilities.Model{Float64}(), optimizer)
    x, cx = MOI.add_constrained_variables(model, MOI.PositiveSemidefiniteConeTriangle(3))
    MOI.set(model, MOI.VariableName(), x[1], "a")
    MOI.Utilities.attach_optimizer(model) # Should drop errors silently in the copy
    @test "a" == MOI.get(model, MOI.VariableName(), x[1])
end

# See https://github.com/jump-dev/MosekTools.jl/issues/95
@testset "Mapping enums" begin
    MOI.empty!(optimizer)
    # Force variable bridging to test attribute substitution
    bridged = MOI.Bridges.Variable.Zeros{Float64}(optimizer)
    cache = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    model = MOI.Utilities.CachingOptimizer(cache, bridged)
    attr = MOI.RawOptimizerAttribute("MSK_IPAR_INTPNT_SOLVE_FORM")
    value = MosekTools.MSK_SOLVE_DUAL
    MOI.add_constrained_variables(model, MOI.Zeros(1))
    MOI.Utilities.attach_optimizer(model)
    # The function should not be used so we can use anything here as first argument
    @test MOI.Utilities.substitute_variables(x -> x^2, value) == value
    MOI.set(model, attr, value)
    # Currently Mosek returns `value.value` but we don't want the tests to fail if this gets fixed in the future
    # so we also allow `value`.
    @test MOI.get(model, attr) == value.value || MOI.get(model, attr) == value
    @test MOI.get(optimizer, attr) == value.value || MOI.get(model, attr) == value
end
