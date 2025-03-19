# Copyright (c) 2017: Ulf Worsøe, Mosek ApS
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# * When we say "variable" we mean an MOI variable and
# * when we say columns, we mean a Mosek variable.
# This is because Mosek variables are numbered 1:n and are columns of the `A`
# matrix while the MOI variables numbering can have holes, i.e. when variables
# are deleted.

###############################################################################
# TASK ########################################################################
###### The `task` field should not be accessed outside this section. ##########

num_columns(task::Mosek.MSKtask) = Mosek.getnumvar(task)

num_columns(m::Optimizer) = num_columns(m.task)

add_column(task::Mosek.MSKtask) = Mosek.appendvars(task, 1)

add_column(m::Optimizer) = add_column(m.task)

"""
    function init_columns(task::Mosek.MSKtask, cols::ColumnIndices)

Set the bound to free, i.e. `MSK_BK_FR`, in the internal Mosek task for the
columns cols.

See [`clear_columns`](@ref) which is kind of the reverse operation.
"""
function init_columns(task::Mosek.MSKtask, cols::ColumnIndices)
    # Set each column to a free variable
    N = length(cols.values)
    bnd = zeros(Float64, N)
    Mosek.putvarboundlist(task, cols.values, fill(Mosek.MSK_BK_FR, N), bnd, bnd)
    if DEBUG
        for col in cols.values
            Mosek.putvarname(task, col, "x$col")
        end
    end
    return
end

function init_columns(m::Optimizer, vis::Vector{MOI.VariableIndex})
    return init_columns(m.task, columns(m, vis))
end

## Name #######################################################################
###############################################################################

"""
    function clear_columns(task::Mosek.MSKtask, cols::Vector{Int32})

Delete all coefficients for columns of indices `cols`. Note that the column
is not actually deleted since it will be reused by a new variable when one
is added.

See [`init_columns`](@ref) which is kind of the reverse operation.
"""
function clear_columns(task::Mosek.MSKtask, cols::ColumnIndices)
    N = length(cols.values)
    # Objective: Clear any non-zeros in `c` vector
    Mosek.putclist(task, cols.values, zeros(Int64, N))
    # Constraints: Clear any non-zeros in columns of `A` matrix
    Mosek.putacollist(
        task,
        cols.values,
        zeros(Int64, N),
        zeros(Int64, N),
        Int32[],
        Float64[],
    )
    # Bounds: Fix the variable to zero to make it have very low footprint for
    #         mosek in case `MOI.optimize!` is called before a new variable is
    #         added to reuse this column.
    bnd = zeros(Float64, N)
    Mosek.putvarboundlist(task, cols.values, fill(Mosek.MSK_BK_FX, N), bnd, bnd)
    if DEBUG
        for col in cols.values
            # Rename deleted column to help debugging
            Mosek.putvarname(task, col, "deleted$col")
        end
    end
    return
end

function clear_columns(m::Optimizer, vis::Vector{MOI.VariableIndex})
    clear_columns(m.task, columns(m, vis))
    return
end

###############################################################################
## INDEXING ###################################################################
###############################################################################

function column(m::Optimizer, vi::MOI.VariableIndex)
    @assert iszero(m.x_sd[vi.value].matrix)
    return ColumnIndex(getindex(m.x_block, ref2id(vi)))
end

function columns(m::Optimizer, vis::Vector{MOI.VariableIndex})
    return ColumnIndices(
        Int32[(mosek_index(m, vi)::ColumnIndex).value for vi in vis],
    )
end

function is_scalar(m::Optimizer, vi::MOI.VariableIndex)
    if 1 <= vi.value <= length(m.x_sd)
        matrix_index = m.x_sd[vi.value]
        if matrix_index.matrix >= 0
            return iszero(matrix_index.matrix)
        end
    end
    return throw(MOI.InvalidIndex(vi))
end

is_matrix(m::Optimizer, vi::MOI.VariableIndex) = !is_scalar(m, vi)

function mosek_index(m::Optimizer, vi::MOI.VariableIndex)
    if is_scalar(m, vi)
        return column(m, vi)
    end
    return m.x_sd[vi.value]
end

function index_of_column(m::Optimizer, col::Int32)
    id = m.x_block.back[col]
    if iszero(id)
        return nothing
    end
    return MOI.VariableIndex(id)
end

"""
    function allocate_variable(m::Optimizer)

Allocate one variable. If there is a free column (because one variable was
deleted), we use it and don't create any new column, i.e. don't call
`appendvars`. Otherwise, we create new column.

See [`clear_variable`](@ref) which is kind of the reverse operation.
"""
function allocate_variable(m::Optimizer, vi::MOI.VariableIndex)
    numvar = num_columns(m)
    alloced = ensurefree(m.x_block, 1)
    if !iszero(alloced)
        @assert isone(alloced)
        @assert length(m.x_block) == num_columns(m) + 1
        add_column(m)
        @assert length(m.x_block) == num_columns(m)
    end
    return allocate_block(m.x_block, 1, vi.value)
end

## Delete #####################################################################
######### Delete variables by clearing the column. The column is reused when ##
######### new variables is added. While there exists a function in the Mosek ##
######### API to delete a column, it is costly and is best avoided.          ##

function throw_if_cannot_delete(m::Optimizer, vi::MOI.VariableIndex)
    id = ref2id(vi)
    if !allocated(m.x_block, ref2id(vi))
        # The matrix variables are created but not allocated so
        # if it can either be a scalar variable that was deleted
        # or a matrix variable
        if m.x_sd[ref2id(vi)].matrix == -1
            # scalar variable that was deleted
            throw(MOI.InvalidIndex(vi))
        else
            # matrix variable
            throw(
                MOI.DeleteNotAllowed(
                    vi,
                    "It is part of a positive semidefinite matrix block.",
                ),
            )
        end
    end
    if is_scalar(m, vi)
        col = column(m, vi)
        for S in [
            MOI.LessThan{Float64},
            MOI.GreaterThan{Float64},
            MOI.EqualTo{Float64},
            MOI.Interval{Float64},
            MOI.Integer,
        ]
            if has_flag(m, vi, S)
                MOI.delete(
                    m,
                    MOI.ConstraintIndex{MOI.VariableIndex,S}(vi.value),
                )
            end
        end
    end
    return
end

"""
    clear_variable(m::Optimizer, vi::MOI.VariableIndex)

Delete the `vi` variable (without actually doing any change to the internal
Mosek solver) to free the corresponding column for reuse by another variable.

See [`allocate_variable`](@ref) which is kind of the reverse operation.
"""
function clear_variable(m::Optimizer, vi::MOI.VariableIndex)
    deleteblock(m.x_block, ref2id(vi))
    return
end

###############################################################################
# MOI #########################################################################
###############################################################################

## Add ########################################################################
###############################################################################

function new_variable_index(m::Optimizer, matrix_index::MatrixIndex)
    id = create_block(m.x_block, 1)
    push!(m.x_constraints, 0x0)
    push!(m.x_sd, matrix_index)
    push!(m.variable_to_vector_constraint_id, 0)
    return MOI.VariableIndex(id)
end

function MOI.add_variable(m::Optimizer)
    vi = new_variable_index(m, MatrixIndex(0, 0, 0))
    allocate_variable(m, vi)
    init_columns(m, [vi])
    return vi
end

## Delete #####################################################################
###############################################################################

function MOI.is_valid(model::Optimizer, vi::MOI.VariableIndex)
    return allocated(model.x_block, ref2id(vi))
end

function delete_vector_of_variables_constraint(
    m::Optimizer,
    vis::Vector{MOI.VariableIndex},
)
    i = first(vis).value
    id = m.variable_to_vector_constraint_id[i]
    if id > 0
        S = type_cone(Mosek.getconeinfo(m.task, id)[1])
        ci = MOI.ConstraintIndex{MOI.VectorOfVariables,S}(i)
        if MOI.is_valid(m, ci) &&
           vis == MOI.get(m, MOI.ConstraintFunction(), ci).variables
            MOI.delete(m, ci)
        end
    end
    return
end

function MOI.delete(m::Optimizer, vis::Vector{MOI.VariableIndex})
    for vi in vis
        throw_if_cannot_delete(m, vi)
    end
    delete_vector_of_variables_constraint(m, vis)
    for vi in vis
        MOI.delete(m, vi)
    end
    return
end

function MOI.delete(m::Optimizer, vi::MOI.VariableIndex)
    throw_if_cannot_delete(m, vi)
    delete_vector_of_variables_constraint(m, [vi])
    if !iszero(m.variable_to_vector_constraint_id[vi.value])
        MOI.Utilities.throw_delete_variable_in_vov(vi)
    end
    clear_columns(m, [vi])
    clear_variable(m, vi)
    m.x_sd[vi.value] = MatrixIndex(-1, 0, 0)
    return
end

## List #######################################################################
###############################################################################

MOI.get(m::Optimizer, ::MOI.NumberOfVariables) = sum(m.x_block.size)

function MOI.get(m::Optimizer, ::MOI.ListOfVariableIndices)
    ids = allocatedlist(m.x_block)
    return MOI.VariableIndex[MOI.VariableIndex(vid) for vid in ids]
end

function MOI.get(m::Optimizer, ::MOI.ListOfVariableAttributesSet)
    ret = MOI.AbstractVariableAttribute[]
    if m.has_variable_names
        push!(ret, MOI.VariableName())
    end
    if !isempty(m.variable_primal_start)
        push!(ret, MOI.VariablePrimalStart())
    end
    return ret
end

## Name #######################################################################
###############################################################################

# We leave `supports` to `false` because it's not supported by matrix indices
# See https://github.com/jump-dev/MosekTools.jl/issues/80
# We still implement the methods for people that knowingly set names of scalar
# variables though.
MOI.supports(::Optimizer, ::MOI.VariableName, ::Type{MOI.VariableIndex}) = false

_throw_if_matrix(::ColumnIndex) = nothing

function _throw_if_matrix(::MatrixIndex)
    msg = "Mosek does not support names for positive semidefinite variables."
    return throw(MOI.UnsupportedAttribute(MOI.VariableName(), msg))
end

function MOI.set(
    m::Optimizer,
    ::MOI.VariableName,
    vi::MOI.VariableIndex,
    name::String,
)
    m.has_variable_names = true
    col = mosek_index(m, vi)
    _throw_if_matrix(col)
    Mosek.putvarname(m.task, col.value, name)
    return
end

function MOI.get(m::Optimizer, ::MOI.VariableName, vi::MOI.VariableIndex)
    col = mosek_index(m, vi)
    _throw_if_matrix(col)
    return Mosek.getvarname(m.task, col.value)
end

function MOI.get(m::Optimizer, ::Type{MOI.VariableIndex}, name::String)
    asgn, index = Mosek.getvarnameindex(m.task, name)
    if iszero(asgn)
        return nothing
    end
    return index_of_column(m, index)
end
