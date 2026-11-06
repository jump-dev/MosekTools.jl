# Copyright (c) 2026: jump-dev contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module TestZeroOne

using Test
import MathOptInterface as MOI
import MosekTools

const Mosek = MosekTools.Mosek

function variable_bound(model, x)
    return Mosek.getvarbound(model.task, MosekTools.column(model, x).value)
end

function variable_type(model, x)
    return Mosek.getvartype(model.task, MosekTools.column(model, x).value)
end

function test_lifecycle()
    model = MosekTools.Optimizer()
    @test MOI.supports_constraint(model, MOI.VariableIndex, MOI.ZeroOne)
    x, c = MOI.add_constrained_variable(model, MOI.ZeroOne())
    integer = MOI.ConstraintIndex{MOI.VariableIndex,MOI.Integer}(x.value)
    @test MOI.is_valid(model, c)
    @test !MOI.is_valid(model, integer)
    @test MOI.get(model, MOI.ConstraintFunction(), c) == x
    @test MOI.get(model, MOI.ConstraintSet(), c) == MOI.ZeroOne()
    @test MOI.get(
        model,
        MOI.NumberOfConstraints{MOI.VariableIndex,MOI.ZeroOne}(),
    ) == 1
    @test MOI.get(
        model,
        MOI.ListOfConstraintIndices{MOI.VariableIndex,MOI.ZeroOne}(),
    ) == [c]
    @test MOI.get(model, MOI.ListOfConstraintTypesPresent()) ==
          [(MOI.VariableIndex, MOI.ZeroOne)]
    @test variable_bound(model, x) == (Mosek.MSK_BK_RA, 0.0, 1.0)
    @test variable_type(model, x) == Mosek.MSK_VAR_TYPE_INT
    # Mimic fixing a binary variable, then restoring its integrality.
    MOI.delete(model, c)
    @test !MOI.is_valid(model, c)
    @test variable_bound(model, x)[1] == Mosek.MSK_BK_FR
    @test variable_type(model, x) == Mosek.MSK_VAR_TYPE_CONT
    fixed = MOI.add_constraint(model, x, MOI.EqualTo(1.0))
    MOI.delete(model, fixed)
    c = MOI.add_constraint(model, x, MOI.ZeroOne())
    @test MOI.is_valid(model, c)
    @test !MOI.is_valid(model, integer)
    @test variable_bound(model, x) == (Mosek.MSK_BK_RA, 0.0, 1.0)
    MOI.delete(model, x)
    @test !MOI.is_valid(model, c)
    @test isempty(model.binary_bounds)
    y = MOI.add_variable(model)
    @test variable_bound(model, y)[1] == Mosek.MSK_BK_FR
    @test variable_type(model, y) == Mosek.MSK_VAR_TYPE_CONT
    MOI.add_constraint(model, y, MOI.ZeroOne())
    MOI.empty!(model)
    @test MOI.is_empty(model)
    @test isempty(model.binary_bounds)
    return
end

function test_bounds()
    for (set, changed, bounds, restored) in (
            (
                MOI.LessThan(2.0),
                MOI.LessThan(0.5),
                (0.0, 0.5),
                (Mosek.MSK_BK_UP, 0.5),
            ),
            (
                MOI.GreaterThan(-2.0),
                MOI.GreaterThan(0.5),
                (0.5, 1.0),
                (Mosek.MSK_BK_LO, 0.5),
            ),
            (
                MOI.EqualTo(1.0),
                MOI.EqualTo(0.0),
                (0.0, 0.0),
                (Mosek.MSK_BK_FX, 0.0),
            ),
            (
                MOI.Interval(-2.0, 2.0),
                MOI.Interval(0.2, 0.8),
                (0.2, 0.8),
                (Mosek.MSK_BK_RA, 0.2),
            ),
        ),
        binary_first in (false, true),
        delete_binary_first in (false, true)

        model = MosekTools.Optimizer()
        x = MOI.add_variable(model)
        if binary_first
            binary = MOI.add_constraint(model, x, MOI.ZeroOne())
            bound = MOI.add_constraint(model, x, set)
        else
            bound = MOI.add_constraint(model, x, set)
            binary = MOI.add_constraint(model, x, MOI.ZeroOne())
        end
        @test MOI.get(model, MOI.ConstraintSet(), bound) == set
        MOI.set(model, MOI.ConstraintSet(), bound, changed)
        @test MOI.get(model, MOI.ConstraintSet(), bound) == changed
        @test variable_bound(model, x) == (Mosek.MSK_BK_RA, bounds...)
        if delete_binary_first
            MOI.delete(model, binary)
            @test MOI.get(model, MOI.ConstraintSet(), bound) == changed
            @test variable_type(model, x) == Mosek.MSK_VAR_TYPE_CONT
            @test variable_bound(model, x)[1] == restored[1]
            if changed isa MOI.LessThan
                @test variable_bound(model, x)[3] == restored[2]
            else
                @test variable_bound(model, x)[2] == restored[2]
            end
            MOI.delete(model, bound)
            @test variable_bound(model, x)[1] == Mosek.MSK_BK_FR
        else
            MOI.delete(model, bound)
            @test variable_bound(model, x) == (Mosek.MSK_BK_RA, 0.0, 1.0)
            @test variable_type(model, x) == Mosek.MSK_VAR_TYPE_INT
            MOI.delete(model, binary)
            @test variable_bound(model, x)[1] == Mosek.MSK_BK_FR
        end
    end
    return
end

function test_two_bounds()
    for lower_first in (false, true), delete_lower_first in (false, true)
        model = MosekTools.Optimizer()
        x, _ = MOI.add_constrained_variable(model, MOI.ZeroOne())
        if lower_first
            lower = MOI.add_constraint(model, x, MOI.GreaterThan(0.2))
            upper = MOI.add_constraint(model, x, MOI.LessThan(0.8))
        else
            upper = MOI.add_constraint(model, x, MOI.LessThan(0.8))
            lower = MOI.add_constraint(model, x, MOI.GreaterThan(0.2))
        end
        @test variable_bound(model, x) == (Mosek.MSK_BK_RA, 0.2, 0.8)
        if delete_lower_first
            MOI.delete(model, lower)
            @test variable_bound(model, x) == (Mosek.MSK_BK_RA, 0.0, 0.8)
            MOI.delete(model, upper)
        else
            MOI.delete(model, upper)
            @test variable_bound(model, x) == (Mosek.MSK_BK_RA, 0.2, 1.0)
            MOI.delete(model, lower)
        end
        @test variable_bound(model, x) == (Mosek.MSK_BK_RA, 0.0, 1.0)
    end
    return
end

function test_integer()
    for binary_first in (false, true), delete_binary_first in (false, true)
        model = MosekTools.Optimizer()
        x = MOI.add_variable(model)
        if binary_first
            binary = MOI.add_constraint(model, x, MOI.ZeroOne())
            integer = MOI.add_constraint(model, x, MOI.Integer())
        else
            integer = MOI.add_constraint(model, x, MOI.Integer())
            binary = MOI.add_constraint(model, x, MOI.ZeroOne())
        end
        @test MOI.is_valid(model, binary)
        @test MOI.is_valid(model, integer)
        if delete_binary_first
            MOI.delete(model, binary)
            @test MOI.is_valid(model, integer)
            @test variable_bound(model, x)[1] == Mosek.MSK_BK_FR
            @test variable_type(model, x) == Mosek.MSK_VAR_TYPE_INT
            MOI.delete(model, integer)
        else
            MOI.delete(model, integer)
            @test MOI.is_valid(model, binary)
            @test variable_bound(model, x) == (Mosek.MSK_BK_RA, 0.0, 1.0)
            @test variable_type(model, x) == Mosek.MSK_VAR_TYPE_INT
            MOI.delete(model, binary)
        end
        @test variable_type(model, x) == Mosek.MSK_VAR_TYPE_CONT
    end
    return
end

function test_matrix_variable()
    model = MosekTools.Optimizer()
    x, _ = MOI.add_constrained_variables(
        model,
        MOI.PositiveSemidefiniteConeTriangle(2),
    )
    @test_throws(
        MOI.AddConstraintNotAllowed{MOI.VariableIndex,MOI.ZeroOne},
        MOI.add_constraint(model, x[1], MOI.ZeroOne()),
    )
    return
end

function test_moi_ZeroOne()
    model = MosekTools.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    MOI.Test.runtests(
        model,
        MOI.Test.Config(; atol = 1e-3, rtol = 1e-3);
        include = ["ZeroOne"],
    )
    return
end

function runtests()
    is_test(name) = startswith(string(name), "test_")
    @testset "$name" for name in filter(is_test, names(@__MODULE__; all = true))
        getfield(@__MODULE__, name)()
    end
    return
end

end

TestZeroOne.runtests()
