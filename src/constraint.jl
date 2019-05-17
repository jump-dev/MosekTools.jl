###############################################################################
# TASK ########################################################################
###### The `task` field should not be accessed outside this section. ##########

## Affine Constraints #########################################################
##################### lc ≤ Ax ≤ uc ############################################

function allocateconstraints(m::MosekModel, N::Int)
    numcon = getnumcon(m.task)
    alloced = ensurefree(m.c_block,N)
    id = newblock(m.c_block, N)

    if alloced > 0
        appendcons(m.task, alloced)
    end
    return id
end

function getconboundlist(t::Mosek.Task, subj::Vector{Int32})
    n = length(subj)
    bk = Vector{Boundkey}(undef,n)
    bl = Vector{Float64}(undef,n)
    bu = Vector{Float64}(undef,n)
    for i in 1:n
        bki,bli,bui = getconbound(t,subj[i])
        bk[i] = bki
        bl[i] = bli
        bu[i] = bui
    end
    bk,bl,bu
end

# `putbaraij` and `putbarcj` need the whole matrix as a sum of sparse mat at once
function split_scalar_matrix(m::MosekModel, terms::Vector{MOI.ScalarAffineTerm{Float64}},
                             set_sd::Function)
    cols = Int32[]
    values = Float64[]
    sd_row  = Vector{Int32}[Int32[] for i in 1:length(m.sd_dim)]
    sd_col  = Vector{Int32}[Int32[] for i in 1:length(m.sd_dim)]
    sd_coef = Vector{Float64}[Float64[] for i in 1:length(m.sd_dim)]
	function add(col::ColumnIndex, coefficient::Float64)
		push!(cols, col.value)
		push!(values, coefficient)
	end
	function add(mat::MatrixIndex, coefficient::Float64)
        coef = mat.row == mat.column ? coefficient : coefficient / 2
        push!(sd_row[mat.matrix], mat.row)
        push!(sd_col[mat.matrix], mat.column)
        push!(sd_coef[mat.matrix], coef)
    end
	for term in terms
        add(mosek_index(m, term.variable_index), term.coefficient)
    end
    for j in 1:length(m.sd_dim)
        if !isempty(sd_row[j])
            id = appendsparsesymmat(m.task, m.sd_dim[j], sd_row[j],
                                    sd_col[j], sd_coef[j])
            set_sd(j, [id], [1.0])
        end
    end
    return cols, values
end

function set_row(task::Mosek.MSKtask, row::Int32, cols::ColumnIndices,
                 values::Vector{Float64})
    putarow(task, row, cols.values, values)
end
function set_row(m::MosekModel, row::Int32,
                 terms::Vector{MOI.ScalarAffineTerm{Float64}})
    cols, values = split_scalar_matrix(m, terms,
        (j, ids, coefs) -> putbaraij(m.task, row, j, ids, coefs))
    set_row(m.task, row, ColumnIndices(cols), values)
end


function set_coefficients(task::Mosek.MSKtask, rows::Vector{Int32},
                          cols::ColumnIndices, values::Vector{Float64})
    putaijlist(task, rows, cols.values, values)
end

function set_coefficients(task::Mosek.MSKtask, rows::Vector{Int32},
                          col::ColumnIndex, values::Vector{Float64})
    n = length(rows)
    @assert n == length(values)
    set_coefficient(task, rows, ColumnIndices(fill(col.value, n)), values)
end
function set_coefficients(m::MosekModel, rows::Vector{Int32},
                          vi::MOI.VariableIndex, values::Vector{Float64})
    set_coefficient(m.task, rows, mosek_index(m, vi), values)
end

function set_coefficient(task::Mosek.MSKtask, row::Int32, col::ColumnIndex,
                         value::Float64)
    putaij(task, row, col.value, value)
end
function set_coefficient(m::MosekModel, row::Int32, vi::MOI.VariableIndex,
                         value::Float64)
    set_coefficient(m.task, row, mosek_index(m, vi), value)
end

bound_key(::Type{MOI.GreaterThan{Float64}}) = MSK_BK_LO
bound_key(::Type{MOI.LessThan{Float64}})    = MSK_BK_UP
bound_key(::Type{MOI.EqualTo{Float64}})     = MSK_BK_FX
bound_key(::Type{MOI.Interval{Float64}})    = MSK_BK_RA
add_bound(m::MosekModel, row::Int32, dom::MOI.GreaterThan{Float64}) = putconbound(m.task, row, bound_key(typeof(dom)), dom.lower, dom.lower)
add_bound(m::MosekModel, row::Int32, dom::MOI.LessThan{Float64})    = putconbound(m.task, row, bound_key(typeof(dom)), dom.upper, dom.upper)
add_bound(m::MosekModel, row::Int32, dom::MOI.EqualTo{Float64})     = putconbound(m.task, row, bound_key(typeof(dom)), dom.value, dom.value)
add_bound(m::MosekModel, row::Int32, dom::MOI.Interval{Float64})    = putconbound(m.task, row, bound_key(typeof(dom)), dom.lower, dom.upper)
function bounds_to_set(::Type{S}, bk, bl, bu) where S
    if S == MOI.GreaterThan{Float64}
        return S(bl)
    elseif S == MOI.LessThan{Float64}
        return S(bu)
    elseif S == MOI.EqualTo{Float64}
        @assert bl == bu
        return S(bu)
    else
        @assert S == MOI.Interval{Float64}
        return S(bl, bu)
    end
end
function get_bound(m::MosekModel,
                   ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, S}) where {S}
    bounds_to_set(S, getconbound(m.task, row(m, ci))...)
end

## Variable Constraints #######################################################
####################### lx ≤ x ≤ u ############################################
#######################      x ∈ K ############################################

function getvarboundlist(t::Mosek.Task, subj::Vector{Int32})
    n = length(subj)
    bk = Vector{Boundkey}(undef,n)
    bl = Vector{Float64}(undef,n)
    bu = Vector{Float64}(undef,n)
    for i in 1:n
        bki,bli,bui = getvarbound(t, subj[i])
        bk[i] = bki
        bl[i] = bli
        bu[i] = bui
    end
    bk,bl,bu
end


function delete_variable_constraint(m::MosekModel, col::ColumnIndex,
                                    ::Type{<:Union{MOI.Interval, MOI.EqualTo}})
    putvarbound(m.task, col.value, MSK_BK_FR, 0.0, 0.0)
end
function delete_variable_constraint(m::MosekModel, col::ColumnIndex,
                                    ::Type{MOI.Integer})
    putvartype(m.task, col.value, MSK_VAR_TYPE_CONT)
end
function delete_variable_constraint(m::MosekModel, col::ColumnIndex,
                                    ::Type{MOI.ZeroOne})
    putvartype(m.task, col.value, MSK_VAR_TYPE_CONT)
    putvarbound(m.task, col.value, MSK_BK_FR, 0.0, 0.0)
end
function delete_variable_constraint(m::MosekModel, col::ColumnIndex,
                                    ::Type{MOI.LessThan{Float64}})
    bk, lo, up = getvarbound(m.task, col.value)
    if bk == MSK_BK_UP
        bk = MSK_BK_FR
    else
        @assert bk == MSK_BK_RA
        bk = MSK_BK_LO
    end
    putvarbound(m.task, col.value, bk, lo, 0.0)
end
function delete_variable_constraint(m::MosekModel, col::ColumnIndex,
                                    ::Type{MOI.GreaterThan{Float64}})
    bk, lo, up = getvarbound(m.task, col.value)
    if bk == MSK_BK_LO
        bk = MSK_BK_FR
    else
        @assert bk == MSK_BK_RA
        bk = MSK_BK_UP
    end
    putvarbound(m.task, col.value, bk, 0.0, up)
end
function add_variable_constraint(m::MosekModel, col::ColumnIndex, dom::MOI.Interval)
    putvarbound(m.task, col.value, MSK_BK_RA, dom.lower, dom.upper)
end
function add_variable_constraint(m::MosekModel, col::ColumnIndex, dom::MOI.EqualTo)
    putvarbound(m.task, col.value, MSK_BK_FX, dom.value, dom.value)
end
function add_variable_constraint(m::MosekModel, col::ColumnIndex, ::MOI.Integer)
    putvartype(m.task, col.value, MSK_VAR_TYPE_INT)
end
function add_variable_constraint(m::MosekModel, col::ColumnIndex, ::MOI.ZeroOne)
    putvartype(m.task, col.value, MSK_VAR_TYPE_INT)
    putvarbound(m.task, col.value, MSK_BK_RA, 0.0, 1.0)
end
function add_variable_constraint(m::MosekModel, col::ColumnIndex, dom::MOI.LessThan)
    bk, lo, up = getvarbound(m.task, col.value)
    if bk == MSK_BK_FR
        bk = MSK_BK_UP
    else
        @assert bk == MSK_BK_LO
        bk = MSK_BK_RA
    end
    putvarbound(m.task, col.value, bk, lo, dom.upper)
end
function add_variable_constraint(m::MosekModel, col::ColumnIndex,
                                 dom::MOI.GreaterThan)
    bk, lo, up = getvarbound(m.task, col.value)
    if bk == MSK_BK_FR
        bk = MSK_BK_LO
    else
        @assert bk == MSK_BK_UP
        bk = MSK_BK_RA
    end
    putvarbound(m.task, col.value, bk, dom.lower, up)
end
function get_variable_constraint(m::MosekModel,
                                 col::ColumnIndex,
                                 ci::MOI.ConstraintIndex{MOI.SingleVariable, S}) where S
    return bounds_to_set(S, getvarbound(m.task, col.value)...)
end
function get_variable_constraint(m::MosekModel, vi::MOI.VariableIndex,
                                 ci::MOI.ConstraintIndex)
    return get_variable_constraint(m, mosek_index(m, vi), ci)
end

cone_parameter(dom :: MOI.PowerCone{Float64})     = dom.exponent
cone_parameter(dom :: MOI.DualPowerCone{Float64}) = dom.exponent
cone_parameter(dom :: C) where C <: MOI.AbstractSet = 0.0

function add_cone(m::MosekModel, cols::ColumnIndices, set)
    appendcone(m.task, cone_type(typeof(set)), cone_parameter(set), cols.values)
    id = getnumcone(m.task)
    if DEBUG
        putconename(m.task, id, "$id")
    end
    return id
end

## Name #######################################################################
###############################################################################

function set_row_name(task::Mosek.MSKtask, row::Int32, name::String)
    putconname(task, row, name)
end

function set_row_name(m::MosekModel,
                      c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}},
                      name::AbstractString)
    set_row_name(m.task, row(m, c), name)
end
function set_row_name(m::MosekModel, c::MOI.ConstraintIndex,
                      name::AbstractString)
    # Fallback for `SingleVariable` and `VectorOfVariables`.
    m.con_to_name[c] = name
end

function delete_name(m::MosekModel, ci::MOI.ConstraintIndex)
    name = MOI.get(m, MOI.ConstraintName(), ci)
    if !isempty(name)
        cis = m.constrnames[name]
        deleteat!(cis, findfirst(isequal(ci), cis))
    end
end

###############################################################################
# INDEXING ####################################################################
###############################################################################

function row(m::MosekModel,
             c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}})::Int32
    return getindex(m.c_block, c.value)
end
function columns(m::MosekModel, ci::MOI.ConstraintIndex{MOI.VectorOfVariables})
    return ColumnIndices(getcone(m.task, ci.value)[4])
end

const VectorCone = Union{MOI.SecondOrderCone,
                         MOI.RotatedSecondOrderCone,
                         MOI.PowerCone,
                         MOI.DualPowerCone,
                         MOI.ExponentialCone,
                         MOI.DualExponentialCone}

# Two `SingleVariable`-in-`S` cannot be set to the same variable if
# the two constraints
# * both set a lower bound, or
# * both set an upper bound, or
# * both set it to integer.
# The `incompatible_mask` are computed according to these rules.
flag(::Type{MOI.EqualTo{Float64}}) = 0x1
incompatible_mask(::Type{MOI.EqualTo{Float64}}) = 0x2f
flag(::Type{MOI.GreaterThan{Float64}}) = 0x2
incompatible_mask(::Type{MOI.GreaterThan{Float64}}) = 0x2b
flag(::Type{MOI.LessThan{Float64}}) = 0x4
incompatible_mask(::Type{MOI.LessThan{Float64}}) = 0x2d
flag(::Type{MOI.Interval{Float64}}) = 0x8
incompatible_mask(::Type{MOI.Interval{Float64}}) = 0x2f
flag(::Type{MOI.Integer}) = 0x10
incompatible_mask(::Type{MOI.Integer}) = 0x30
flag(::Type{MOI.ZeroOne}) = 0x20
incompatible_mask(::Type{MOI.ZeroOne}) = 0x3f
flag(::Type{<:VectorCone}) = 0x40
incompatible_mask(::Type{<:VectorCone}) = 0x40

function set_flag(model::MosekModel, vi::MOI.VariableIndex, S::Type)
    model.x_constraints[vi.value] |= flag(S)
end
function unset_flag(model::MosekModel, vi::MOI.VariableIndex, S::Type)
    model.x_constraints[vi.value] &= ~flag(S)
end
function has_flag(model::MosekModel, vi::MOI.VariableIndex, S::Type)
    return !iszero(model.x_constraints[vi.value] & flag(S))
end

###############################################################################
# MOI #########################################################################
###############################################################################

const ScalarLinearDomain = Union{MOI.LessThan{Float64},
                                 MOI.GreaterThan{Float64},
                                 MOI.EqualTo{Float64},
                                 MOI.Interval{Float64}}
const ScalarIntegerDomain = Union{MOI.ZeroOne, MOI.Integer}

## Add ########################################################################
###############################################################################

MOI.supports_constraint(m::MosekModel, ::Type{<:Union{MOI.SingleVariable, MOI.ScalarAffineFunction}}, ::Type{<:ScalarLinearDomain}) = true
MOI.supports_constraint(m::MosekModel, ::Type{MOI.VectorOfVariables}, ::Type{<:Union{VectorCone, MOI.PositiveSemidefiniteConeTriangle}}) = true
MOI.supports_constraint(m::MosekModel, ::Type{MOI.SingleVariable}, ::Type{<:ScalarIntegerDomain}) = true

## Affine Constraints #########################################################
##################### lc ≤ Ax ≤ uc ############################################

function MOI.add_constraint(m  ::MosekModel,
                            axb::MOI.ScalarAffineFunction{Float64},
                            dom::D) where {D <: MOI.AbstractScalarSet}

    if !iszero(axb.constant)
        throw(MOI.ScalarFunctionConstantNotZero{Float64, typeof(axb), D}(axb.constant))
    end

    # Duplicate indices not supported
    axb = MOIU.canonical(axb)

    N = 1
    conid = allocateconstraints(m, N)
    ci = MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, D}(conid)
    r = row(m, ci)
    set_row(m, r, axb.terms)

    add_bound(m, r, dom)

    return ci
end

## Variable Constraints #######################################################
####################### lx ≤ x ≤ u ############################################
#######################      x ∈ K ############################################

# We allow following. Each variable can have
# - at most most upper and one lower bound
# - belong to at most one non-semidefinite cone
# - any number of semidefinite cones, which are implemented as ordinary constraints
# This is when things get a bit funky; By default a variable has no
# bounds, i.e. "free". Adding a `GreaterThan`
# constraint causes it to have a defined lower bound but no upper
# bound, allowing a `LessThan` constraint to be
# added later. Adding a `Interval` constraint defines both upper and
# lower bounds.

cone_type(::Type{MOI.ExponentialCone})        = MSK_CT_PEXP
cone_type(::Type{MOI.DualExponentialCone})    = MSK_CT_DEXP
cone_type(::Type{MOI.PowerCone{Float64}})     = MSK_CT_PPOW
cone_type(::Type{MOI.DualPowerCone{Float64}}) = MSK_CT_DPOW
cone_type(::Type{MOI.SecondOrderCone})        = MSK_CT_QUAD
cone_type(::Type{MOI.RotatedSecondOrderCone}) = MSK_CT_RQUAD

function MOI.add_constraint(
    m   :: MosekModel,
    xs  :: MOI.SingleVariable,
    dom :: D) where {D <: MOI.AbstractScalarSet}

    msk_idx = mosek_index(m, xs.variable)
    if !(msk_idx isa ColumnIndex)
        error("Cannot add $D constraint on a matrix variable")
    end

    if !iszero(incompatible_mask(D) & m.x_constraints[xs.variable.value])
        error("Cannot put multiple bound sets of the same type on a variable")
    end

    set_flag(m, xs.variable, D)

    add_variable_constraint(m, msk_idx, dom)

    return MOI.ConstraintIndex{MOI.SingleVariable, D}(xs.variable.value)
end

function MOI.add_constraint(m::MosekModel, xs::MOI.VectorOfVariables,
                            dom::D) where {D<:VectorCone}
    if any(vi -> is_matrix(m, vi), xs.variables)
        error("Cannot add $D constraint on a matrix variable")
    end
    cols = ColumnIndices(reorder(columns(m, xs.variables).values, D))

    if !all(vi -> iszero(incompatible_mask(D) & m.x_constraints[vi.value]), xs.variables)
        error("Cannot multiple bound sets of the same type to a variable")
    end

    N = MOI.dimension(dom)
    id = add_cone(m, cols, dom)

    return MOI.ConstraintIndex{MOI.VectorOfVariables, D}(id)
end

################################################################################
################################################################################

function MOI.add_constraint(m  ::MosekModel,
                            xs ::MOI.VectorOfVariables,
                            dom::MOI.PositiveSemidefiniteConeTriangle)
    N = dom.side_dimension
    if N < 2
        error("Invalid dimension for semidefinite constraint, got $N which is ",
              "smaller than the minimum dimension 2.")
    end
    appendbarvars(m.task, [Int32(N)])
    push!(m.sd_dim, N)
    id = length(m.sd_dim)
    sd_row  = Dict{Int, Vector{Int32}}()
    sd_col  = Dict{Int, Vector{Int32}}()
    sd_coef = Dict{Int, Vector{Float64}}()
    obj_row = Int32[]
    obj_col = Int32[]
    obj_coef = Float64[]
    k = 0
    for i in 1:N
        for j in 1:i
            k += 1
            vi = xs.variables[k]
            m.x_sd[vi.value] = MatrixIndex(id, i, j)

            if variable_type(m, vi) == MatrixVariable
                error("Variable $vi cannot be part of two matrix variables.")
            elseif variable_type(m, vi) == ScalarVariable
                # The variable has already been used, we need replace the
                # coefficients in the `A` matrix by matrix coefficients
                nnz, rows, vals = getacol(m.task, column(m, vi).value)
                @assert nnz == length(rows) == length(vals)
                for ii in 1:nnz
                    if !haskey(sd_row, rows[ii])
                        sd_row[rows[ii]] = Int32[]
                        sd_col[rows[ii]] = Int32[]
                        sd_coef[rows[ii]] = Float64[]
                    end
                    push!(sd_row[rows[ii]], i)
                    push!(sd_col[rows[ii]], j)
                    push!(sd_coef[rows[ii]], i == j ? vals[ii] : vals[ii] / 2)
                end
                cj = getcj(m.task, column(m, vi).value)
                if !iszero(cj)
                    push!(obj_row, i)
                    push!(obj_col, j)
                    push!(obj_coef, i == j ? cj : cj / 2)
                end
                MOI.delete(m, vi)
            end
            m.x_type[vi.value] = MatrixVariable
        end
    end
    for row in eachindex(sd_row)
        sid = appendsparsesymmat(m.task, N, sd_row[row], sd_col[row],
                                sd_coef[row])
        putbaraij(m.task, row, id, [sid], [1.0])
    end
    if !isempty(obj_row)
        sid = appendsparsesymmat(m.task, N, obj_row, obj_col, obj_coef)
        putbarcj(m.task, id, [sid], [1.0])
    end
    @assert k == length(xs.variables)


    conref = MOI.ConstraintIndex{MOI.VectorOfVariables, MOI.PositiveSemidefiniteConeTriangle}(id)

    return conref
end

## Get ########################################################################
###############################################################################

_variable(ci::MOI.ConstraintIndex{MOI.SingleVariable}) = MOI.VariableIndex(ci.value)
function MOI.get(m::MosekModel, ::MOI.ConstraintFunction,
                 ci::MOI.ConstraintIndex{MOI.SingleVariable}) where S <: ScalarLinearDomain
    MOI.throw_if_not_valid(m, ci)
    return MOI.SingleVariable(_variable(ci))
end
function MOI.get(m::MosekModel, ::MOI.ConstraintSet,
                 ci::MOI.ConstraintIndex{MOI.SingleVariable, S}) where S <: ScalarIntegerDomain
    return S()
end
function MOI.get(m::MosekModel, ::MOI.ConstraintSet,
                 ci::MOI.ConstraintIndex{MOI.SingleVariable, S}) where S <: ScalarLinearDomain
    MOI.throw_if_not_valid(m, ci)
    sv = MOI.get(m, MOI.ConstraintFunction(), ci)
    return get_variable_constraint(m, sv.variable, ci)
end

function MOI.get(m::MosekModel, ::MOI.ConstraintFunction,
                 ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},
                                         <:ScalarLinearDomain})
    nnz, cols, vals = getarow(m.task, row(m, ci))
    @assert nnz == length(cols) == length(vals)
    terms = MOI.ScalarAffineTerm{Float64}[
        MOI.ScalarAffineTerm(vals[i], index_of_column(m, cols[i])) for i in 1:nnz]
    # TODO add matrix terms
    return MOI.ScalarAffineFunction(terms, 0.0)
end
function MOI.get(m::MosekModel, ::MOI.ConstraintSet,
                 ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},
                                         S}) where S
    return get_bound(m, ci)
end

## Modify #####################################################################
###############################################################################

### SET
chgbound(bl::Float64,bu::Float64,k::Float64,dom :: MOI.LessThan{Float64})    = bl,dom.upper-k
chgbound(bl::Float64,bu::Float64,k::Float64,dom :: MOI.GreaterThan{Float64}) = dom.lower-k,bu
chgbound(bl::Float64,bu::Float64,k::Float64,dom :: MOI.EqualTo{Float64})     = dom.value-k,dom.value-k
chgbound(bl::Float64,bu::Float64,k::Float64,dom :: MOI.Interval{Float64})    = dom.lower-k,dom.upper-k

function MOI.set(m::MosekModel,
                 ::MOI.ConstraintSet,
                 ci::MOI.ConstraintIndex{F,D},
                 dom::D) where {F<:MOI.SingleVariable,
                                D<:ScalarLinearDomain}
    col = column(m, _variable(ci))
    bk, bl, bu = getvarbound(m.task, col.value)
    bl, bu = chgbound(bl, bu, 0.0, dom)
    putvarbound(m.task, col.value, bk, bl, bu)
end


function MOI.set(m::MosekModel,
                 ::MOI.ConstraintSet,
                 cref::MOI.ConstraintIndex{F,D},
                 dom::D) where { F    <: MOI.ScalarAffineFunction,
                                 D    <: ScalarLinearDomain }
    cid = ref2id(cref)
    i = getindex(m.c_block,cid) # since we are in a scalar domain
    bk, bl, bu = getconbound(m.task,i)
    bl, bu = chgbound(bl,bu,0.0,dom)
    putconbound(m.task,i,bk,bl,bu)
end


### MODIFY
function MOI.modify(m   ::MosekModel,
                    c   ::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},D},
                    func::MOI.ScalarConstantChange{Float64}) where {D <: MOI.AbstractSet}
    if !iszero(func.new_constant)
        throw(MOI.ScalarFunctionConstantNotZero{Float64, MOI.ScalarAffineFunction{Float64}, D}(func.new_constant))
    end
end

function MOI.modify(m   ::MosekModel,
                    c   ::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}},
                    func::MOI.ScalarCoefficientChange{Float64})
    set_coefficient(m, row(m, c), func.variable, func.new_coefficient)
end

### TRANSFORM
function MOI.transform(m::MosekModel,
                       cref::MOI.ConstraintIndex{F,D},
                       newdom::D) where {F <: MOI.AbstractFunction,
                                         D <: MOI.AbstractSet}
    MOI.modify(m,cref,newdom)
    return cref
end

function MOI.transform(m::MosekModel,
                       cref::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},D1},
                       newdom::D2) where {D1 <: ScalarLinearDomain,
                                          D2 <: ScalarLinearDomain}
    F = MOI.ScalarAffineFunction{Float64}

    cid = ref2id(cref)

    r = row(m, cref)
    add_bound(m, r, newdom)

    newcref = MOI.ConstraintIndex{F, D2}(cid)
    return newcref
end

## Delete #####################################################################
###############################################################################

function MOI.is_valid(model::MosekModel,
                      ci::MOI.ConstraintIndex{<:MOI.ScalarAffineFunction{Float64},
                                              S}) where S<:ScalarLinearDomain
    return allocated(model.c_block, ci.value) && getconbound(model.task, row(model, ci))[1] == bound_key(S)
end
function MOI.delete(
    m::MosekModel,
    cref::MOI.ConstraintIndex{F,D}) where {F <: MOI.ScalarAffineFunction{Float64},
                                           D <: ScalarLinearDomain}
    MOI.throw_if_not_valid(m, cref)

    delete_name(m, cref)

    subi = getindexes(m.c_block, cref.value)

    n = length(subi)
    subi_i32 = convert(Vector{Int32}, subi)
    ptr = fill(Int64(0), n)
    putarowlist(m.task,subi_i32,ptr,ptr,Int32[],Float64[])
    b = fill(0.0,n)
    putconboundlist(m.task,subi_i32,fill(MSK_BK_FX,n),b,b)

    deleteblock(m.c_block, cref.value)
end

function MOI.is_valid(model::MosekModel,
                      ci::MOI.ConstraintIndex{MOI.SingleVariable,
                                              S}) where S<:Union{ScalarLinearDomain,
                                                                 ScalarIntegerDomain}
    return allocated(model.x_block, ci.value) && has_flag(model, _variable(ci), S)
end
function MOI.delete(
    m::MosekModel,
    ci::MOI.ConstraintIndex{MOI.SingleVariable, S}) where S<:Union{ScalarLinearDomain,
                                                                   ScalarIntegerDomain}
    MOI.throw_if_not_valid(m, ci)
    delete_name(m, ci)
    vi = _variable(ci)
    unset_flag(m, vi, S)
    delete_variable_constraint(m, column(m, vi), S)
end
function MOI.is_valid(model::MosekModel,
                      ci::MOI.ConstraintIndex{MOI.VectorOfVariables,
                                              S}) where S<:VectorCone
    return 1 ≤ ci.value ≤ getnumcone(model.task) &&
        getconeinfo(model.task, ci.value)[1] == cone_type(S)
end
function MOI.is_valid(model::MosekModel,
                      ci::MOI.ConstraintIndex{MOI.VectorOfVariables,
                                              MOI.PositiveSemidefiniteConeTriangle})
    # TODO add supports for deletion
    return 1 ≤ ci.value ≤ length(model.sd_dim)
end


## List #######################################################################
###############################################################################
function MOI.get(m::MosekModel, ::MOI.ListOfConstraintAttributesSet)
    set = MOI.AbstractConstraintAttribute[]
    if !isempty(m.constrnames)
        push!(set, MOI.ConstraintName())
    end
    # TODO add VariablePrimalStart when get is implemented on it
    return set
end

## Name #######################################################################
###############################################################################
function MOI.supports(::MosekModel, ::MOI.ConstraintName,
                      ::Type{<:MOI.ConstraintIndex})
    return true
end
function MOI.set(m::MosekModel, ::MOI.ConstraintName, ci::MOI.ConstraintIndex,
                 name ::AbstractString)
    delete_name(m, ci)
    if !haskey(m.constrnames, name)
        m.constrnames[name] = MOI.ConstraintIndex[]
    end
    push!(m.constrnames[name], ci)
    set_row_name(m, ci, name)
end
function MOI.get(m::MosekModel, ::MOI.ConstraintName,
                 ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}})
    # All rows should have same name so we take the first one
    return getconname(m.task, getindexes(m.c_block, ref2id(ci))[1])
end
function MOI.get(m::MosekModel, ::MOI.ConstraintName, ci::MOI.ConstraintIndex)
    return get(m.con_to_name, ci, "")
end
function MOI.get(m::MosekModel, CI::Type{<:MOI.ConstraintIndex}, name::String)
    #asgn, row = getconnameindex(m.task, name)
    #if iszero(asgn)
    #    return nothing
    #else
    #    # TODO how to recover function and constraint type ?
    #end
    if !haskey(m.constrnames, name)
        return nothing
    end
    cis = m.constrnames[name]
    if isempty(cis)
        return nothing
    end
    if !isone(length(cis))
        error("Multiple constraints have the name $name.")
    end
    ci = first(cis)
    if !(ci isa CI)
        return nothing
    end
    return ci
end
