# * When we say "variable" we mean an MOI variable and
# * when we say columns, we mean a Mosek variable.
# This is because Mosek variables are numbered 1:n and are columns of the `A`
# matrix while the MOI variables numbering can have holes, i.e. when variables
# are deleted.

###############################################################################
# TASK ########################################################################
###### The `task` field should not be accessed outside this section. ##########

num_columns(task::Mosek.MSKtask) = getnumvar(task)
num_columns(m::MosekModel) = num_columns(m.task)

add_column(task::Mosek.MSKtask) = appendvars(task, 1)
add_column(m::MosekModel) = add_column(m.task)

"""
    function init_columns(task::Mosek.MSKtask, cols::Vector{Int32})

Set the bound to free, i.e. `MSK_BK_FR`, in the internal Mosek task for the
columns cols.

See [`clear_columns`](@ref) which is kind of the reverse operation.
"""
function init_columns(task::Mosek.MSKtask, cols::Vector{Int32})
    # Set each column to a free variable
    N = length(cols)
    bnd = zeros(Float64, N)
    putvarboundlist(task, cols, fill(MSK_BK_FR, N), bnd, bnd)

    if DEBUG
        for col in cols
            putvarname(task, col, "x$col")
        end
    end
    return
end
function init_columns(m::MosekModel, refs::Vector{MOI.VariableIndex})
    init_columns(m.task, columns(m, refs))
end

## Name #######################################################################
###############################################################################
function set_column_name(task::Mosek.MSKtask, col::Int32, name::String)
    putvarname(task, col, name)
end
function set_column_name(m::MosekModel, ref::MOI.VariableIndex, name::String)
    set_column_name(m.task, column(m, ref), name)
end
column_name(task::Mosek.MSKtask, col::Int32) = getvarname(task, col)
function column_name(m::MosekModel, ref::MOI.VariableIndex)
    column_name(m.task, column(m, ref))
end
function column_with_name(task::Mosek.MSKtask, name::String)
    asgn, col = getvarnameindex(task, name)
    if iszero(asgn)
        return nothing
    else
        return col
    end
end
column_with_name(m::MosekModel, name::String) = column_with_name(m.task, name)

"""
    function clear_columns(task::Mosek.MSKtask, cols::Vector{Int32})

Delete all coefficients for columns of indices `cols`. Note that the column
is not actually deleted since it will be reused by a new variable when one
is added.

See [`init_columns`](@ref) which is kind of the reverse operation.
"""
function clear_columns(task::Mosek.MSKtask, cols::Vector{Int32})
    N = length(cols)
    # Objective: Clear any non-zeros in `c` vector
    putclist(task, cols, zeros(Int64, N))

    # Constraints: Clear any non-zeros in columns of `A` matrix
    putacollist(task,
                cols,
                zeros(Int64, N),
                zeros(Int64, N),
                Int32[],
                Float64[])

    # Bounds: Fix the variable to zero to make it have very low footprint for
    #         mosek in case `MOI.optimize!` is called before a new variable is
    #         added to reuse this column.
    bnd = zeros(Float64, N)
    putvarboundlist(task, cols, fill(MSK_BK_FX, N), bnd, bnd)

    if DEBUG
        for col in cols
            # Rename deleted column to help debugging
            putvarname(task, col, "deleted$col")
        end
    end
    return
end
function clear_columns(m::MosekModel, refs::Vector{MOI.VariableIndex})
    clear_columns(m.task, columns(m, refs))
end

###############################################################################
## INDEXING ###################################################################
###############################################################################

function column(m::MosekModel, ref::MOI.VariableIndex)::Int32
    return getindex(m.x_block, ref2id(ref))
end
function columns(m::MosekModel, refs::Vector{MOI.VariableIndex})
    return Int32[column(m, ref) for ref in refs]
end
function index_of_column(m::MosekModel, col::Int32)
    id = m.x_block.back[col]
    if iszero(id)
        return nothing
    else
        return MOI.VariableIndex(id)
    end
end

"""
    function allocate_variable(m::MosekModel)

Allocate one variable. If there is a free column (because one variable was
deleted), we use it and don't create any new column, i.e. don't call
`appendvars`. Otherwise, we create new column.

See [`clear_variable`](@ref) which is kind of the reverse operation.
"""
function allocate_variable(m::MosekModel)
    @assert length(m.x_boundflags) == length(m.x_block)
    numvar = num_columns(m)
    alloced = ensurefree(m.x_block, 1)
    if !iszero(alloced)
        @assert isone(alloced)
        @assert length(m.x_block) == num_columns(m) + 1
        add_column(m)
        @assert length(m.x_block) == num_columns(m)
        push!(m.x_boundflags, 0)
        push!(m.x_numxc, 0)
    end
    m.publicnumvar += 1
    return MOI.VariableIndex(newblock(m.x_block, 1))
end

## Delete #####################################################################
######### Delete variables by clearing the column. The column is reused when ##
######### new variables is added. While there exists a function in the Mosek ##
######### API to delete a column, it is costly and is best avoided.          ##

function throw_if_cannot_delete(m::MosekModel, ref::MOI.VariableIndex)
    if !allocated(m.x_block, ref2id(ref))
        throw(MOI.InvalidIndex(ref))
    end
    if !iszero(m.x_numxc[ref2id(ref)])
        throw(CannotDelete("Cannot delete variable $ref while a bound constraint is defined on it"))
    end
end

"""
    clear_variable(m::MosekModel, ref::MOI.VariableIndex)

Delete the `ref` variable (without actually doing any change to the internal
Mosek solver) to free the corresponding column for reuse by another variable.

See [`allocate_variable`](@ref) which is kind of the reverse operation.
"""
function clear_variable(m::MosekModel, ref::MOI.VariableIndex)
    m.publicnumvar -= 1
    deleteblock(m.x_block, ref2id(ref))
end

###############################################################################
# MOI #########################################################################
###############################################################################

## Add ########################################################################
###############################################################################

function MOI.add_variables(m::MosekModel, N::Integer)
    @assert N ≥ 0
    refs = MOI.VariableIndex[allocate_variable(m) for i in 1:N]
    init_columns(m, refs)
    return refs
end
function MOI.add_variable(m::MosekModel)
    ref = allocate_variable(m)
    init_columns(m, [ref])
    return ref
end

## Delete #####################################################################
###############################################################################
function MOI.is_valid(model::MosekModel, ref::MOI.VariableIndex)
    return allocated(model.x_block, ref2id(ref))
end

function MOI.delete(m::MosekModel, refs::Vector{MOI.VariableIndex})
    for ref in refs
        throw_if_cannot_delete(m, ref)
    end
    clear_columns(m, refs)
    for ref in refs
        clear_variable(m, ref)
    end
end
function MOI.delete(m::MosekModel, ref::MOI.VariableIndex)
    throw_if_cannot_delete(m, ref)
    clear_columns(m, [ref])
    clear_variable(m, ref)
end

## List #######################################################################
###############################################################################
MOI.get(m::MosekModel, attr::MOI.NumberOfVariables) = m.publicnumvar
function MOI.get(m::MosekModel, attr::MOI.ListOfVariableIndices)
    ids = allocatedlist(m.x_block)
    return MOI.VariableIndex[MOI.VariableIndex(vid) for vid in ids]
end

## Name #######################################################################
###############################################################################
MOI.supports(::MosekModel, ::MOI.VariableName, ::Type{MOI.VariableIndex}) = true
function MOI.set(m::MosekModel, ::MOI.VariableName, ref::MOI.VariableIndex,
                 name::String)
    set_column_name(m, ref, name)
end
function MOI.get(m::MosekModel, ::MOI.VariableName, ref::MOI.VariableIndex)
    return column_name(m, ref)
end
function MOI.get(m::MosekModel, ::Type{MOI.VariableIndex}, name::String)
    col = column_with_name(m, name)
    if col === nothing
        return nothing
    else
        return index_of_column(m, col)
    end
end
