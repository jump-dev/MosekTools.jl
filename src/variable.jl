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
## INDEXING ###################################################################
###############################################################################

function column(m::Optimizer, vi::MOI.VariableIndex)
    @assert iszero(m.x_sd[vi.value].matrix)
    return ColumnIndex(getindex(m.x_block, vi.value))
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

function _col_to_index(m::Optimizer, col::Int32)
    @assert !iszero(m.x_block.back[col])
    return MOI.VariableIndex(m.x_block.back[col])
end

## Delete #####################################################################
######### Delete variables by clearing the column. The column is reused when ##
######### new variables is added. While there exists a function in the Mosek ##
######### API to delete a column, it is costly and is best avoided.          ##

function _throw_if_cannot_delete(m::Optimizer, vi::MOI.VariableIndex)
    if !allocated(m.x_block, vi.value)
        # The matrix variables are created but not allocated so
        # if it can either be a scalar variable that was deleted
        # or a matrix variable
        if m.x_sd[vi.value].matrix == -1
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
    # Allocate one variable. If there is a free column (because one variable was
    # deleted), we use it and don't create any new column, i.e. don't call
    # `appendvars`. Otherwise, we create new column.
    numvar = Mosek.getnumvar(m.task)
    alloced = ensurefree(m.x_block, 1)
    if !iszero(alloced)
        @assert isone(alloced)
        @assert length(m.x_block) == Mosek.getnumvar(m.task) + 1
        Mosek.appendvars(m.task, 1)
        @assert length(m.x_block) == Mosek.getnumvar(m.task)
    end
    allocate_block(m.x_block, 1, vi.value)
    Mosek.putvarbound(m.task, column(m, vi).value, Mosek.MSK_BK_FR, 0.0, 0.0)
    return vi
end

## Delete #####################################################################
###############################################################################

function MOI.is_valid(model::Optimizer, vi::MOI.VariableIndex)
    return get(model.x_block.size, vi.value, 0) == 1
end

function _delete_vector_of_variables_constraint(
    m::Optimizer,
    vis::Vector{MOI.VariableIndex},
)
    i = first(vis).value
    id = m.variable_to_vector_constraint_id[i]
    if id <= 0
        return
    end
    S = _type_cone(Mosek.getconeinfo(m.task, id)[1])
    ci = MOI.ConstraintIndex{MOI.VectorOfVariables,S}(i)
    if MOI.is_valid(m, ci) &&
       vis == MOI.get(m, MOI.ConstraintFunction(), ci).variables
        MOI.delete(m, ci)
    end
    return
end

function MOI.delete(m::Optimizer, vis::Vector{MOI.VariableIndex})
    for vi in vis
        _throw_if_cannot_delete(m, vi)
    end
    _delete_vector_of_variables_constraint(m, vis)
    for vi in vis
        MOI.delete(m, vi)
    end
    return
end

function MOI.delete(m::Optimizer, vi::MOI.VariableIndex)
    _throw_if_cannot_delete(m, vi)
    _delete_vector_of_variables_constraint(m, [vi])
    if !iszero(m.variable_to_vector_constraint_id[vi.value])
        MOI.Utilities.throw_delete_variable_in_vov(vi)
    end
    col = column(m, vi)
    # Objective: Clear any non-zero in `c` vector
    Mosek.putcj(m.task, col.value, 0.0)
    # Constraints: Clear any non-zeros in columns of `A` matrix
    Mosek.putacol(m.task, col.value, Int32[], Float64[])
    # Bounds: Fix the variable to zero to make it have very low footprint for
    #         mosek in case `MOI.optimize!` is called before a new variable is
    #         added to reuse this column.
    Mosek.putvarbound(m.task, col.value, Mosek.MSK_BK_FX, 0.0, 0.0)
    deleteblock(m.x_block, vi.value)
    m.x_sd[vi.value] = MatrixIndex(-1, 0, 0)
    delete!(m.variable_to_name, vi)
    m.name_to_variable = nothing
    return
end

## List #######################################################################
###############################################################################

MOI.get(m::Optimizer, ::MOI.NumberOfVariables) = sum(m.x_block.size)

function MOI.get(m::Optimizer, ::MOI.ListOfVariableIndices)
    return MOI.VariableIndex[
        MOI.VariableIndex(i) for (i, x) in enumerate(m.x_sd) if x.matrix >= 0
    ]
end

function MOI.get(m::Optimizer, ::MOI.ListOfVariableAttributesSet)
    ret = MOI.AbstractVariableAttribute[]
    if !isempty(m.variable_to_name)
        push!(ret, MOI.VariableName())
    end
    if !isempty(m.variable_primal_start)
        push!(ret, MOI.VariablePrimalStart())
    end
    return ret
end

## Name #######################################################################
###############################################################################

MOI.supports(::Optimizer, ::MOI.VariableName, ::Type{MOI.VariableIndex}) = true

function _putvarname(m::Optimizer, col::ColumnIndex, name::AbstractString)
    Mosek.putvarname(m.task, col.value, name)
    return
end

_putvarname(::Optimizer, ::MatrixIndex, ::AbstractString) = nothing

function MOI.set(
    m::Optimizer,
    ::MOI.VariableName,
    vi::MOI.VariableIndex,
    name::String,
)
    m.variable_to_name[vi] = name
    m.name_to_variable = nothing
    _putvarname(m, mosek_index(m, vi), name)
    return
end

function MOI.get(m::Optimizer, ::MOI.VariableName, vi::MOI.VariableIndex)
    return get(m.variable_to_name, vi, "")
end

function _rebuild_name_to_variable_index(m::Optimizer)
    m.name_to_variable = Dict{String,Union{Nothing,MOI.VariableIndex}}()
    for (x, name) in m.variable_to_name
        if isempty(name)
            continue
        end
        if haskey(m.name_to_variable, name)
            m.name_to_variable[name] = nothing
        else
            m.name_to_variable[name] = x
        end
    end
    return
end

function MOI.get(m::Optimizer, ::Type{MOI.VariableIndex}, name::String)
    if m.name_to_variable === nothing
        _rebuild_name_to_variable_index(m)
    end
    index = get(m.name_to_variable, name, missing)
    if ismissing(index)
        return nothing  # No name exists
    elseif isnothing(index)
        error("Multiple variables have the name: $name")
    end
    return index
end
