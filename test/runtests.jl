# Copyright (c) 2017: Ulf Worsøe, Mosek ApS
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module TestMosekTools

using MosekTools
using Test

import MathOptInterface as MOI

const _MOSEK_API_COUNTER = Ref{Int}(0)

if isdefined(MOI.Test, :_error_handler)
    @eval function MOI.Test._error_handler(
        err::Mosek.MosekError,
        name::String,
        warn_unsupported::Bool,
    )
        if Mosek.Rescode(err.rcode) == Mosek.MSK_RES_ERR_SERVER_STATUS
            _MOSEK_API_COUNTER[] += 1
            return  # Server returned non-ok HTTP status code
        elseif Mosek.Rescode(err.rcode) == Mosek.MSK_RES_ERR_CONE_OVERLAP_APPEND
            return  # Mosek doesn't support this problem formulation
        end
        return rethrow(err)
    end
end

function runtests()
    _MOSEK_API_COUNTER[] = 0
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$name" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    # Test that there are less than 5 API call failures during a test run
    @test _MOSEK_API_COUNTER[] < 5
    return
end

function MosekOptimizerWithFallback()
    optimizer = Mosek.Optimizer()
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("fallback"),
        "mosek://solve.mosek.com:30080",
    )
    MOI.set(optimizer, MOI.Silent(), true)
    return optimizer
end

function test_SolverName()
    @test MOI.get(Mosek.Optimizer(), MOI.SolverName()) == "Mosek"
    return
end

function test_Double_Parameter()
    optimizer = MosekOptimizerWithFallback()
    MOI.set(optimizer, MOI.RawOptimizerAttribute("INTPNT_CO_TOL_DFEAS"), 1e-7)
    @test MOI.get(
        optimizer,
        MOI.RawOptimizerAttribute("MSK_DPAR_INTPNT_CO_TOL_DFEAS"),
    ) == 1e-7
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("MSK_DPAR_INTPNT_CO_TOL_DFEAS"),
        1e-8,
    )
    @test MOI.get(
        optimizer,
        MOI.RawOptimizerAttribute("MSK_DPAR_INTPNT_CO_TOL_DFEAS"),
    ) == 1e-8
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("MSK_DPAR_INTPNT_CO_TOL_DFEAS"),
        1,
    )
    @test MOI.get(
        optimizer,
        MOI.RawOptimizerAttribute("MSK_DPAR_INTPNT_CO_TOL_DFEAS"),
    ) == 1
    return
end

function test_Integer_Parameter()
    optimizer = MosekOptimizerWithFallback()
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("MSK_IPAR_INTPNT_MAX_ITERATIONS"),
        100,
    )
    @test MOI.get(
        optimizer,
        MOI.RawOptimizerAttribute("MSK_IPAR_INTPNT_MAX_ITERATIONS"),
    ) == 100
    MOI.set(optimizer, MOI.RawOptimizerAttribute("INTPNT_MAX_ITERATIONS"), 200)
    @test MOI.get(
        optimizer,
        MOI.RawOptimizerAttribute("MSK_IPAR_INTPNT_MAX_ITERATIONS"),
    ) == 200
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("MSK_IPAR_OPTIMIZER"),
        Mosek.MSK_OPTIMIZER_DUAL_SIMPLEX,
    )
    @test MOI.get(optimizer, MOI.RawOptimizerAttribute("MSK_IPAR_OPTIMIZER")) ==
          convert(Int32, Mosek.MSK_OPTIMIZER_DUAL_SIMPLEX)
    return
end

function test_String_Parameter()
    optimizer = MosekOptimizerWithFallback()
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("PARAM_WRITE_FILE_NAME"),
        "foo.txt",
    )
    # Needs https://github.com/JuliaOpt/Mosek.jl/pull/174
    #@test MOI.get(optimizer, MOI.RawOptimizerAttribute("MSK_SPAR_PARAM_WRITE_FILE_NAME")) == "foo.txt"
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("MSK_SPAR_PARAM_WRITE_FILE_NAME"),
        "bar.txt",
    )
    #@test MOI.get(optimizer, MOI.RawOptimizerAttribute("MSK_SPAR_PARAM_WRITE_FILE_NAME")) == "bar.txt"
    return
end

function test_TimeLimitSec()
    optimizer = MosekOptimizerWithFallback()
    @test MOI.get(
        optimizer,
        MOI.RawOptimizerAttribute("MSK_DPAR_OPTIMIZER_MAX_TIME"),
    ) == -1
    @test MOI.get(optimizer, MOI.TimeLimitSec()) === nothing
    MOI.set(optimizer, MOI.TimeLimitSec(), 1.0)
    @test MOI.get(
        optimizer,
        MOI.RawOptimizerAttribute("MSK_DPAR_OPTIMIZER_MAX_TIME"),
    ) == 1.0
    @test MOI.get(optimizer, MOI.TimeLimitSec()) === 1.0
    MOI.set(optimizer, MOI.TimeLimitSec(), nothing)
    @test MOI.get(
        optimizer,
        MOI.RawOptimizerAttribute("MSK_DPAR_OPTIMIZER_MAX_TIME"),
    ) == -1
    @test MOI.get(optimizer, MOI.TimeLimitSec()) === nothing
    return
end

function test_supports_incremental_interface()
    @test MOI.supports_incremental_interface(Mosek.Optimizer())
    return
end

function test_Number_of_variables_and_deletion()
    optimizer = MosekOptimizerWithFallback()
    x = MOI.add_variable(optimizer)
    @test MOI.get(optimizer, MOI.NumberOfVariables()) == 1
    y, cy = MOI.add_constrained_variables(
        optimizer,
        MOI.PositiveSemidefiniteConeTriangle(3),
    )
    @test MOI.get(optimizer, MOI.NumberOfVariables()) == 7
    MOI.delete(optimizer, x)
    @test MOI.get(optimizer, MOI.NumberOfVariables()) == 6
    z, cz = MOI.add_constrained_variables(optimizer, MOI.SecondOrderCone(3))
    @test MOI.get(optimizer, MOI.NumberOfVariables()) == 9
    @test_throws MOI.DeleteNotAllowed MOI.delete(optimizer, y[1])
    @test_throws MOI.DeleteNotAllowed MOI.delete(optimizer, y)
    @test MOI.get(optimizer, MOI.NumberOfVariables()) == 9
    @test_throws MOI.DeleteNotAllowed MOI.delete(optimizer, z[1])
    @test MOI.get(optimizer, MOI.NumberOfVariables()) == 9
    MOI.delete(optimizer, z)
    @test MOI.get(optimizer, MOI.NumberOfVariables()) == 6
    return
end

function test_Matrix_name()
    optimizer = MosekOptimizerWithFallback()
    x, cx = MOI.add_constrained_variables(
        optimizer,
        MOI.PositiveSemidefiniteConeTriangle(3),
    )
    err = MOI.UnsupportedAttribute{MOI.VariableName}
    @test_throws err MOI.set(optimizer, MOI.VariableName(), x[1], "a")
    MOI.empty!(optimizer)
    model = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.Model{Float64}(),
        optimizer,
    )
    x, cx = MOI.add_constrained_variables(
        model,
        MOI.PositiveSemidefiniteConeTriangle(3),
    )
    MOI.set(model, MOI.VariableName(), x[1], "a")
    MOI.Utilities.attach_optimizer(model) # Should drop errors silently in the copy
    @test "a" == MOI.get(model, MOI.VariableName(), x[1])
    return
end

# See https://github.com/jump-dev/MosekTools.jl/issues/95
function test_Mapping_enums()
    optimizer = MosekOptimizerWithFallback()
    # Force variable bridging to test attribute substitution
    bridged = MOI.Bridges.Variable.Zeros{Float64}(optimizer)
    cache = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    model = MOI.Utilities.CachingOptimizer(cache, bridged)
    attr = MOI.RawOptimizerAttribute("MSK_IPAR_INTPNT_SOLVE_FORM")
    value = Mosek.MSK_SOLVE_DUAL
    MOI.add_constrained_variables(model, MOI.Zeros(1))
    MOI.Utilities.attach_optimizer(model)
    # The function should not be used so we can use anything here as first argument
    @test MOI.Utilities.substitute_variables(x -> x^2, value) == value
    MOI.set(model, attr, value)
    # Currently Mosek returns `value.value` but we don't want the tests to fail if this gets fixed in the future
    # so we also allow `value`.
    @test MOI.get(model, attr) == value.value || MOI.get(model, attr) == value
    @test MOI.get(optimizer, attr) == value.value ||
          MOI.get(model, attr) == value
    return
end

function test_LMIs()
    optimizer = MosekOptimizerWithFallback()
    @test MOI.supports_constraint(
        optimizer,
        MOI.VectorAffineFunction{Float64},
        MOI.Scaled{MOI.PositiveSemidefiniteConeTriangle},
    )
    bridged = MOI.Bridges.full_bridge_optimizer(optimizer, Float64)
    @test MOI.Bridges.bridge_type(
        bridged,
        MOI.VectorAffineFunction{Float64},
        MOI.PositiveSemidefiniteConeTriangle,
    ) === MOI.Bridges.Constraint.SetDotScalingBridge{
        Float64,
        MOI.PositiveSemidefiniteConeTriangle,
        MOI.VectorAffineFunction{Float64},
        MOI.VectorAffineFunction{Float64},
    }
    return
end

function _test_symmetric_reorder(lower, n)
    set = MOI.Scaled(MOI.PositiveSemidefiniteConeTriangle(n))
    N = MOI.dimension(set)
    @test MosekTools.reorder(
        lower,
        MOI.ScaledPositiveSemidefiniteConeTriangle,
        true,
    ) == 1:N
    @test MosekTools.reorder(
        1:N,
        MOI.ScaledPositiveSemidefiniteConeTriangle,
        false,
    ) == lower
    for (up, low) in enumerate(lower)
        @test MosekTools.reorder(up, set, true) == low
        @test MosekTools.reorder(low, set, false) == up
    end
end

function test_symmetric_reorder()
    _test_symmetric_reorder([1], 1)
    _test_symmetric_reorder([1, 2, 3], 2)
    _test_symmetric_reorder([1, 2, 4, 3, 5, 6], 3)
    _test_symmetric_reorder([1, 2, 5, 3, 6, 8, 4, 7, 9, 10], 4)
    _test_symmetric_reorder(
        [1, 2, 6, 3, 7, 10, 4, 8, 11, 13, 5, 9, 12, 14, 15],
        5,
    )
    _test_symmetric_reorder(
        [
            1,
            2,
            7,
            3,
            8,
            12,
            4,
            9,
            13,
            16,
            5,
            10,
            14,
            17,
            19,
            6,
            11,
            15,
            18,
            20,
            21,
        ],
        6,
    )
    return
end

function test_variable_basis_status()
    model = Mosek.Optimizer()
    x, cx = MOI.add_constrained_variables(
        model,
        MOI.PositiveSemidefiniteConeTriangle(2),
    )
    attr = MOI.VariableBasisStatus()
    index = MosekTools.mosek_index(model, x[1])
    err = ErrorException("$attr not supported for PSD variable $index")
    @test_throws err MOI.get(model, attr, index)
    return
end

function test_modify_psd()
    model = Mosek.Optimizer()
    x, _ = MOI.add_constrained_variables(
        model,
        MOI.PositiveSemidefiniteConeTriangle(2),
    )
    c = MOI.add_constraint(model, 1.0 * x[1], MOI.LessThan(1.0))
    change = MOI.ScalarCoefficientChange(x[1], 2.0)
    err = MOI.ModifyConstraintNotAllowed(
        c,
        change,
        "Modifying the coefficient of the variable correspond to an entry of a PSD matrix is not supported",
    )
    @test_throws err MOI.modify(model, c, change)
    attr = MOI.ObjectiveFunction{typeof(1.0x[1])}()
    MOI.set(model, attr, 1.0x[1])
    change = MOI.ScalarCoefficientChange(x[1], 2.0)
    err = MOI.ModifyObjectiveNotAllowed(
        change,
        "Modifying the coefficient of the variable correspond to an entry of a PSD matrix is not supported",
    )
    err = MOI.SetAttributeNotAllowed(
        attr,
        "Cannot set a different objective if a previous objective was set including the contribution of the entry of a PSD variable.",
    )
    @test_throws err MOI.set(model, attr, 2.0x[1])
    y = MOI.add_variable(model)
    @test_throws err MOI.set(model, attr, 1.0y)
    return
end

function test_moi_test_runtests_Mosek()
    config = MOI.Test.Config(
        Float64;
        atol = 1e-3,
        rtol = 1e-3,
        exclude = Any[
            MOI.ConstraintName,
            MOI.ConstraintBasisStatus,
            # TODO remove `MOI.delete` once it is implemented for ACC
            MOI.delete,
        ],
    )
    MOI.Test.runtests(
        MosekOptimizerWithFallback(),
        config;
        exclude = [
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
            # Expression: status in (config.optimal_status, MOI.INVALID_MODEL)
            # Evaluated: MathOptInterface.OTHER_ERROR in (MathOptInterface.OPTIMAL, MathOptInterface.INVALID_MODEL)
            "test_conic_empty_matrix",
        ],
    )
    return
end

function test_moi_test_runtests_Bridge_Mosek()
    model =
        MOI.Bridges.full_bridge_optimizer(MosekOptimizerWithFallback(), Float64)
    config = MOI.Test.Config(
        Float64;
        atol = 1e-3,
        rtol = 1e-3,
        exclude = Any[
            MOI.ConstraintName,
            MOI.ConstraintBasisStatus,
            # TODO remove `MOI.delete` once it is implemented for ACC
            MOI.delete,
        ],
    )
    MOI.Test.runtests(
        model,
        config;
        exclude = [
            "test_model_ListOfConstraintAttributesSet", # list not properly set
            "test_model_duplicate_VariableName",
            "test_model_VariablePrimalStart", # able to set but not to get VariablePrimalStart
            # Cannot put multiple bound sets of the same type on a variable
            "test_basic_VectorAffineFunction_Circuit",
            "test_basic_VectorOfVariables_Circuit",
            "test_basic_VectorQuadraticFunction_Circuit",
            "test_cpsat_Circuit",
            "test_basic_VectorOfVariables_NormCone",
            # Evaluated: MathOptInterface.OTHER_ERROR in (MathOptInterface.OPTIMAL, MathOptInterface.INVALID_MODEL)
            "test_conic_empty_matrix",
            # FIXME ConstraintPrimal incorrect, to investigate
            "test_conic_HermitianPositiveSemidefiniteConeTriangle_1",
            # FIXME ListOfConstraints incorrect
            "test_conic_SecondOrderCone_VectorAffineFunction",
        ],
    )
    @test MOI.supports(model, MOI.VariablePrimalStart(), MOI.VariableIndex)
    return
end

function test_moi_test_runtests_Bridge_Cache_Mosek()
    model = MOI.Bridges.full_bridge_optimizer(
        MOI.Utilities.CachingOptimizer(
            MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
            MosekOptimizerWithFallback(),
        ),
        Float64,
    )
    config = MOI.Test.Config(
        Float64;
        atol = 1e-3,
        rtol = 1e-3,
        exclude = Any[
            MOI.ConstraintName,
            MOI.ConstraintBasisStatus,
            # TODO remove `MOI.delete` once it is implemented for ACC
            MOI.delete,
        ],
    )
    MOI.Test.runtests(
        model,
        config;
        exclude = [
            # Evaluated: MathOptInterface.OTHER_ERROR in (MathOptInterface.OPTIMAL, MathOptInterface.INVALID_MODEL)
            "test_conic_empty_matrix",
            # FIXME ConstraintPrimal incorrect, to investigate
            "test_conic_HermitianPositiveSemidefiniteConeTriangle_1",
        ],
    )
    return
end

function test_more_SDP_tests_by_forced_bridging()
    model = MOI.Bridges.full_bridge_optimizer(
        MOI.Bridges.Constraint.RSOCtoPSD{Float64}( # Forced bridging
            MOI.Utilities.CachingOptimizer(
                MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
                MosekOptimizerWithFallback(),
            ),
        ),
        Float64,
    )
    config = MOI.Test.Config(
        Float64;
        atol = 1e-3,
        rtol = 1e-3,
        exclude = Any[
            MOI.ConstraintName,
            MOI.ConstraintBasisStatus,
            # TODO remove `MOI.delete` once it is implemented for ACC
            MOI.delete,
        ],
    )
    MOI.Test.runtests(model, config; include = ["conic_SecondOrderCone"])
    return
end

end  # module

TestMosekTools.runtests()
