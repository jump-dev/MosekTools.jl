# * When we say "variable" we mean an MOI variable and
# * when we say columns, we mean a Mosek variable.
# This is because Mosek variables are numbered 1:n and are columns of the `A`
# matrix while the MOI variables numbering can have holes, i.e. when variables
# are deleted.

###############################################################################
## ADD

"""
    function allocate_variable(m::MosekModel)

Allocate one variable. If there is a free column (because one variable was
deleted), we use it and don't create any new column, i.e. don't call
`appendvars`. Otherwise, we create new column.

See [`clear_variable`](@ref) which is kind of the reverse operation.
"""
function allocate_variable(m::MosekModel)
    @assert length(m.x_boundflags) == length(m.x_block)
    numvar = getnumvar(m.task)
    alloced = ensurefree(m.x_block, 1)
    if !iszero(alloced)
        @assert isone(alloced)
        @assert length(m.x_block) == getnumvar(m.task) + 1
        appendvars(m.task, 1)
        @assert length(m.x_block) == getnumvar(m.task)
        push!(m.x_boundflags, 0)
        push!(m.x_numxc, 0)
    end
    m.publicnumvar += 1
    return MOI.VariableIndex(newblock(m.x_block, 1))
end

function MOI.is_valid(model::MosekModel, ref::MOI.VariableIndex)
    return allocated(model.x_block, ref2id(ref))
end

"""
    function init_columns(m::MosekModel, refs::Vector{MOI.VariableIndex})

Set the bound to free, i.e. `MSK_BK_FR`, in the internal Mosek solver for the
column corresponding to each variable in `ref`.

See [`clear_columns`](@ref) which is kind of the reverse operation.
"""
function init_columns(m::MosekModel, refs::Vector{MOI.VariableIndex})
    cols = columns(m, refs)

    # Set each column to a free variable
    N = length(cols)
    bnd = zeros(Float64, N)
    putvarboundlist(m.task, cols, fill(MSK_BK_FR, N), bnd, bnd)

    if DEBUG
        for col in cols
            putvarname(m.task, col, "x$col")
        end
    end
    return
end

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

###############################################################################
## DELETE
## Delete variables by clearing the column. The column is reused when a new
## variables is added. While there exists a function in the Mosek API to delete
## a column, it is costly and is best avoided.

function throw_if_cannot_delete(m::MosekModel, ref::MOI.VariableIndex)
    if !allocated(m.x_block, ref2id(ref))
        throw(MOI.InvalidIndex(ref))
    end
    if !iszero(m.x_numxc[ref2id(ref)])
        throw(CannotDelete("Cannot delete variable $ref while a bound constraint is defined on it"))
    end
end

function column(m::MosekModel, ref::MOI.VariableIndex)::Int32
    return getindex(m.x_block, ref2id(ref))
end
function columns(m::MosekModel, refs::Vector{MOI.VariableIndex})
    return Int32[column(m, ref) for ref in refs]
end

"""
    function clear_columns(m::MosekModel, refs::Vector{MOI.VariableIndex})

Delete all coefficients for columns of indices `cols`. Note that the column
is not actually deleted since it will be reused by a new variable when one
is added.

See [`init_columns`](@ref) which is kind of the reverse operation.
"""
function clear_columns(m::MosekModel, refs::Vector{MOI.VariableIndex})
    cols = columns(m, refs)
    N = length(cols)
    # Objective: Clear any non-zeros in `c` vector
    putclist(m.task, cols, zeros(Int64, N))

    # Constraints: Clear any non-zeros in columns of `A` matrix
    putacollist(m.task,
                cols,
                zeros(Int64, N),
                zeros(Int64, N),
                Int32[],
                Float64[])

    # Bounds: Fix the variable to zero to make it have very low footprint for
    #         mosek in case `MOI.optimize!` is called before a new variable is
    #         added to reuse this column.
    bnd = zeros(Float64,N)
    putvarboundlist(m.task, cols, fill(MSK_BK_FX, N), bnd, bnd)

    if DEBUG
        for col in cols
            # Rename deleted column to help debugging
            putvarname(m.task, col, "deleted$col")
        end
    end
    return
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

function MOI.delete(m::MosekModel, refs::Vector{MOI.VariableIndex})
    for ref in refs
        throw_if_cannot_delete(m, ref)
    end

    ids = Int[ ref2id(ref) for ref in refs ]

    sizes = Int[blocksize(m.x_block,id) for id in ids]
    N = sum(sizes)
    indexes = Array{Int}(undef,N)
    offset = 1
    for i in 1:length(ids)
        getindexes(m.x_block,ids[i],indexes,offset)
        offset += sizes[i]
    end

    clear_columns(m, refs)

    for ref in refs
        clear_variable(m, ref)
    end
end

function MOI.delete(m::MosekModel, ref::MOI.VariableIndex)
    throw_if_cannot_delete(m, ref)

    id = ref2id(ref)

    clear_columns(m, [ref])

    clear_variable(m, ref)
end

###############################################################################
## ATTRIBUTES

MOI.get(m::MosekModel, attr::MOI.NumberOfVariables) = m.publicnumvar
function MOI.get(m::MosekModel, attr::MOI.ListOfVariableIndices)
    ids = allocatedlist(m.x_block)
    return MOI.VariableIndex[MOI.VariableIndex(vid) for vid in ids]
end

function MOI.set(m::MosekModel, ::MOI.VariableName, ref::MOI.VariableIndex,
                 value::String)
    putvarname(m.task, column(m, ref), value)
end
