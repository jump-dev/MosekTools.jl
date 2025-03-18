# Copyright (c) 2017: Ulf Worsøe, Mosek ApS
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

###############################################################################
# TASK ########################################################################
###### The `task` field should not be accessed outside this section. ##########

## Affine Constraints #########################################################
##################### lc ≤ Ax ≤ uc ############################################
const VectorCone = Union{
    MOI.SecondOrderCone,
    MOI.RotatedSecondOrderCone,
    MOI.PowerCone,
    MOI.DualPowerCone,
    MOI.ExponentialCone,
    MOI.DualExponentialCone,
}

# VectorConeDomain is a union of all domain types that can currently
# be mapped from Julia domain types to Mosek ACCs
const VectorConeDomain = Union{
    MOI.SecondOrderCone,
    MOI.RotatedSecondOrderCone,
    MOI.PowerCone,
    MOI.DualPowerCone,
    MOI.ExponentialCone,
    MOI.DualExponentialCone,
    MOI.Scaled{MOI.PositiveSemidefiniteConeTriangle},
}

function allocateconstraints(m::Optimizer, N::Int)
    numcon = Mosek.getnumcon(m.task)
    alloced = ensurefree(m.c_block, N)
    id = newblock(m.c_block, N)
    if alloced > 0
        Mosek.appendcons(m.task, alloced)
    end
    return id
end

function getconboundlist(t::Mosek.Task, subj::Vector{Int32})
    n = length(subj)
    bk = Vector{Mosek.Boundkey}(undef, n)
    bl = Vector{Float64}(undef, n)
    bu = Vector{Float64}(undef, n)
    for i in 1:n
        bki, bli, bui = Mosek.getconbound(t, subj[i])
        bk[i] = bki
        bl[i] = bli
        bu[i] = bui
    end
    return bk, bl, bu
end

# `putbaraij` and `putbarcj` need the whole matrix as a sum of sparse mat at once
function split_scalar_matrix(
    m::Optimizer,
    terms::Vector{MOI.ScalarAffineTerm{Float64}},
    set_sd::Function,
)
    cols = Int32[]
    values = Float64[]
    # `terms` is in canonical form so the variables belonging to the same
    # matrix appear adjacent to each other so we can reuse the vector for all
    # matrices. Allocating one vector for each matrix can cause performance
    # issues; see https://github.com/jump-dev/MosekTools.jl/issues/135
    current_matrix = -1
    sd_row = Int32[]
    sd_col = Int32[]
    sd_coef = Float64[]
    function add(col::ColumnIndex, coefficient::Float64)
        push!(cols, col.value)
        return push!(values, coefficient)
    end
    function add_sd()
        if current_matrix != -1
            @assert !isempty(sd_row)
            id = Mosek.appendsparsesymmat(
                m.task,
                m.sd_dim[current_matrix],
                sd_row,
                sd_col,
                sd_coef,
            )
            set_sd(current_matrix, [id], [1.0])
        end
    end
    function add(mat::MatrixIndex, coefficient::Float64)
        @assert mat.matrix != -1
        if current_matrix != mat.matrix
            add_sd()
            current_matrix = mat.matrix
            empty!(sd_row)
            empty!(sd_col)
            empty!(sd_coef)
        end
        coef = mat.row == mat.column ? coefficient : coefficient / 2
        push!(sd_row, mat.row)
        push!(sd_col, mat.column)
        return push!(sd_coef, coef)
    end
    for term in terms
        add(mosek_index(m, term.variable), term.coefficient)
    end
    add_sd()
    return cols, values
end

function set_row(
    task::Mosek.MSKtask,
    row::Int32,
    cols::ColumnIndices,
    values::Vector{Float64},
)
    return Mosek.putarow(task, row, cols.values, values)
end

function set_row(
    m::Optimizer,
    row::Int32,
    terms::Vector{MOI.ScalarAffineTerm{Float64}},
)
    cols, values = split_scalar_matrix(
        m,
        terms,
        (j, ids, coefs) -> Mosek.putbaraij(m.task, row, j, ids, coefs),
    )
    return set_row(m.task, row, ColumnIndices(cols), values)
end

function set_coefficients(
    task::Mosek.MSKtask,
    rows::Vector{Int32},
    cols::ColumnIndices,
    values::Vector{Float64},
)
    return Mosek.putaijlist(task, rows, cols.values, values)
end

function set_coefficients(
    task::Mosek.MSKtask,
    rows::Vector{Int32},
    col::ColumnIndex,
    values::Vector{Float64},
)
    n = length(rows)
    @assert n == length(values)
    return set_coefficient(
        task,
        rows,
        ColumnIndices(fill(col.value, n)),
        values,
    )
end

function set_coefficients(
    m::Optimizer,
    rows::Vector{Int32},
    vi::MOI.VariableIndex,
    values::Vector{Float64},
)
    return set_coefficient(m.task, rows, mosek_index(m, vi), values)
end

function set_coefficient(
    task::Mosek.MSKtask,
    row::Int32,
    col::ColumnIndex,
    value::Float64,
)
    return Mosek.putaij(task, row, col.value, value)
end

function set_coefficient(
    m::Optimizer,
    row::Int32,
    vi::MOI.VariableIndex,
    value::Float64,
)
    return set_coefficient(m.task, row, mosek_index(m, vi), value)
end

bound_key(::Type{MOI.GreaterThan{Float64}}) = Mosek.MSK_BK_LO
bound_key(::Type{MOI.LessThan{Float64}}) = Mosek.MSK_BK_UP
bound_key(::Type{MOI.EqualTo{Float64}}) = Mosek.MSK_BK_FX
bound_key(::Type{MOI.Interval{Float64}}) = Mosek.MSK_BK_RA

function _bounds(dom::MOI.Interval)
    bl = dom.lower
    bu = dom.upper
    bk = bound_key(typeof(dom))
    if bl < 0 && isinf(bl)
        if bu > 0 && isinf(bu)
            bk = Mosek.MSK_BK_FR
        else
            bk = Mosek.MSK_BK_UP
        end
    elseif bu > 0 && isinf(bu)
        bk = Mosek.MSK_BK_LO
    end
    return bk, bl, bu
end

function add_bound(m::Optimizer, row::Int32, dom::MOI.GreaterThan{Float64})
    return Mosek.putconbound(
        m.task,
        row,
        bound_key(typeof(dom)),
        dom.lower,
        dom.lower,
    )
end

function add_bound(m::Optimizer, row::Int32, dom::MOI.LessThan{Float64})
    return Mosek.putconbound(
        m.task,
        row,
        bound_key(typeof(dom)),
        dom.upper,
        dom.upper,
    )
end

function add_bound(m::Optimizer, row::Int32, dom::MOI.EqualTo{Float64})
    return Mosek.putconbound(
        m.task,
        row,
        bound_key(typeof(dom)),
        dom.value,
        dom.value,
    )
end

function add_bound(m::Optimizer, row::Int32, dom::MOI.Interval{Float64})
    return Mosek.putconbound(m.task, row, _bounds(dom)...)
end

function bounds_to_set(::Type{S}, bk, bl, bu) where {S}
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

function get_bound(
    m::Optimizer,
    ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S},
) where {S}
    return bounds_to_set(S, Mosek.getconbound(m.task, row(m, ci))...)
end

function get_bound(
    m::Optimizer,
    ci::MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64},S},
) where {S<:VectorConeDomain}
    domidx = Mosek.getaccdomain(m.task, ci.value)
    dt = Mosek.getdomaintype(m.task, domidx)
    if dt == Mosek.MSK_DOMAIN_QUADRATIC_CONE
        MOI.SecondOrderCone(solsize(m, ci))
    elseif dt == Mosek.MSK_DOMAIN_RQUADRATIC_CONE
        MOI.RotatedSecondOrderCone(solsize(m, ci))
    elseif dt == Mosek.MSK_DOMAIN_PRIMAL_EXP_CONE
        MOI.ExponentialCone()
    elseif dt == Mosek.MSK_DOMAIN_DUAL_EXP_CONE
        MOI.DualExponentialCone()
    elseif dt == Mosek.MSK_DOMAIN_PRIMAL_POWER_CONE
        (n, nleft) = Mosek.getpowerdomaininfo(m.task, domidx)
        if n != 3 || nleft != 2
            error("Incompatible power cone detected")
        end
        alpha = Mosek.getpowerdomainalpha(m.task, domidx)
        MOI.PowerCone(alpha[1])
    elseif dt == Mosek.MSK_DOMAIN_DUAL_POWER_CONE
        (n, nleft) = Mosek.getpowerdomaininfo(m.task, domidx)
        if n != 3 || nleft != 2
            error("Incompatible power cone detected")
        end
        alpha = Mosek.getpowerdomainalpha(m.task, domidx)
        MOI.DualPowerCone(alpha[1])
        # elseif dt == Mosek.MSK_DOMAIN_R
        # elseif dt == Mosek.MSK_DOMAIN_RZERO
        # elseif dt == Mosek.MSK_DOMAIN_RPLUS
        # elseif dt == Mosek.MSK_DOMAIN_RMINUS
        # elseif dt == Mosek.MSK_DOMAIN_PRIMAL_GEO_MEAN_CONE
        # elseif dt == Mosek.MSK_DOMAIN_DUAL_GEO_MEAN_CONE
    elseif dt == Mosek.MSK_DOMAIN_SVEC_PSD_CONE
        return MOI.Scaled{MOI.PositiveSemidefiniteConeTriangle}(
            MOI.Utilities.side_dimension_for_vectorized_dimension(
                solsize(m, ci),
            ),
        )
    else
        error("Incompatible cone detected")
    end
end

## Variable Constraints #######################################################
####################### lx ≤ x ≤ u ############################################
#######################      x ∈ K ############################################

function getvarboundlist(t::Mosek.Task, subj::Vector{Int32})
    n = length(subj)
    bk = Vector{Mosek.Boundkey}(undef, n)
    bl = Vector{Float64}(undef, n)
    bu = Vector{Float64}(undef, n)
    for i in 1:n
        bki, bli, bui = Mosek.getvarbound(t, subj[i])
        bk[i] = bki
        bl[i] = bli
        bu[i] = bui
    end
    return bk, bl, bu
end

function delete_variable_constraint(
    m::Optimizer,
    col::ColumnIndex,
    ::Type{<:Union{MOI.Interval,MOI.EqualTo}},
)
    return Mosek.putvarbound(m.task, col.value, Mosek.MSK_BK_FR, 0.0, 0.0)
end

function delete_variable_constraint(
    m::Optimizer,
    col::ColumnIndex,
    ::Type{MOI.Integer},
)
    return Mosek.putvartype(m.task, col.value, Mosek.MSK_VAR_TYPE_CONT)
end

function delete_variable_constraint(
    m::Optimizer,
    col::ColumnIndex,
    ::Type{MOI.LessThan{Float64}},
)
    bk, lo, up = Mosek.getvarbound(m.task, col.value)
    if bk == Mosek.MSK_BK_UP
        bk = Mosek.MSK_BK_FR
    else
        @assert bk == Mosek.MSK_BK_RA
        bk = Mosek.MSK_BK_LO
    end
    return Mosek.putvarbound(m.task, col.value, bk, lo, 0.0)
end

function delete_variable_constraint(
    m::Optimizer,
    col::ColumnIndex,
    ::Type{MOI.GreaterThan{Float64}},
)
    bk, lo, up = Mosek.getvarbound(m.task, col.value)
    if bk == Mosek.MSK_BK_LO
        bk = Mosek.MSK_BK_FR
    else
        @assert bk == Mosek.MSK_BK_RA
        bk = Mosek.MSK_BK_UP
    end
    return Mosek.putvarbound(m.task, col.value, bk, 0.0, up)
end

function add_variable_constraint(
    m::Optimizer,
    col::ColumnIndex,
    dom::MOI.Interval,
)
    return Mosek.putvarbound(m.task, col.value, _bounds(dom)...)
end

function add_variable_constraint(
    m::Optimizer,
    col::ColumnIndex,
    dom::MOI.EqualTo,
)
    return Mosek.putvarbound(
        m.task,
        col.value,
        Mosek.MSK_BK_FX,
        dom.value,
        dom.value,
    )
end

function add_variable_constraint(m::Optimizer, col::ColumnIndex, ::MOI.Integer)
    return Mosek.putvartype(m.task, col.value, Mosek.MSK_VAR_TYPE_INT)
end

function add_variable_constraint(
    m::Optimizer,
    col::ColumnIndex,
    dom::MOI.LessThan,
)
    bk, lo, up = Mosek.getvarbound(m.task, col.value)
    if bk == Mosek.MSK_BK_FR
        bk = Mosek.MSK_BK_UP
    else
        @assert bk == Mosek.MSK_BK_LO
        bk = Mosek.MSK_BK_RA
    end
    return Mosek.putvarbound(m.task, col.value, bk, lo, dom.upper)
end

function add_variable_constraint(
    m::Optimizer,
    col::ColumnIndex,
    dom::MOI.GreaterThan,
)
    bk, lo, up = Mosek.getvarbound(m.task, col.value)
    if bk == Mosek.MSK_BK_FR
        bk = Mosek.MSK_BK_LO
    else
        @assert bk == Mosek.MSK_BK_UP
        bk = Mosek.MSK_BK_RA
    end
    return Mosek.putvarbound(m.task, col.value, bk, dom.lower, up)
end

function get_variable_constraint(
    m::Optimizer,
    col::ColumnIndex,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,S},
) where {S}
    return bounds_to_set(S, Mosek.getvarbound(m.task, col.value)...)
end

function get_variable_constraint(
    m::Optimizer,
    vi::MOI.VariableIndex,
    ci::MOI.ConstraintIndex,
)
    return get_variable_constraint(m, mosek_index(m, vi), ci)
end

cone_parameter(dom::MOI.PowerCone{Float64}) = dom.exponent
cone_parameter(dom::MOI.DualPowerCone{Float64}) = dom.exponent
cone_parameter(dom::C) where {C<:MOI.AbstractSet} = 0.0

function add_cone(m::Optimizer, cols::ColumnIndices, set)
    Mosek.appendcone(
        m.task,
        cone_type(typeof(set)),
        cone_parameter(set),
        cols.values,
    )
    id = Mosek.getnumcone(m.task)
    if DEBUG
        Mosek.putconename(m.task, id, "$id")
    end
    return id
end

## Name #######################################################################
###############################################################################

function set_row_name(task::Mosek.MSKtask, row::Int32, name::String)
    return Mosek.putconname(task, row, name)
end

function set_row_name(
    m::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}},
    name::AbstractString,
)
    return set_row_name(m.task, row(m, c), name)
end

function set_row_name(
    m::Optimizer,
    c::MOI.ConstraintIndex,
    name::AbstractString,
)
    # Fallback for `SingleVariable` and `VectorOfVariables`.
    return m.con_to_name[c] = name
end

function delete_name(m::Optimizer, ci::MOI.ConstraintIndex)
    name = MOI.get(m, MOI.ConstraintName(), ci)
    if !isempty(name)
        cis = m.constrnames[name]
        deleteat!(cis, findfirst(isequal(ci), cis))
    end
    return
end

###############################################################################
# INDEXING ####################################################################
###############################################################################

function row(
    m::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}},
)::Int32
    return getindex(m.c_block, c.value)
end

function columns(m::Optimizer, ci::MOI.ConstraintIndex{MOI.VectorOfVariables})
    coneidx = cone_id(m, ci)
    if coneidx < 1 || coneidx > Mosek.getnumcone(m.task)
        throw(MOI.InvalidIndex(ci))
    end
    return ColumnIndices(Mosek.getcone(m.task, coneidx)[4])
end

function rows(
    m::Optimizer,
    ci::MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64}},
)
    if !(ci.value in keys(m.F_rows))
        throw(MOI.InvalidIndex(ci))
    end
    return m.F_rows[ci.value]
end

function appendconedomain(t::Mosek.Task, ::Int, dom::MOI.ExponentialCone)
    return Mosek.appendprimalexpconedomain(t)
end

function appendconedomain(t::Mosek.Task, ::Int, dom::MOI.DualExponentialCone)
    return Mosek.appenddualexpconedomain(t)
end

function appendconedomain(t::Mosek.Task, ::Int, dom::MOI.PowerCone{Float64})
    return Mosek.appendprimalpowerconedomain(
        t,
        3,
        Float64[dom.exponent, 1.0-dom.exponent],
    )
end

function appendconedomain(t::Mosek.Task, ::Int, dom::MOI.DualPowerCone{Float64})
    return Mosek.appenddualpowerconedomain(
        t,
        3,
        Float64[dom.exponent, 1.0-dom.exponent],
    )
end

function appendconedomain(t::Mosek.Task, n::Int, ::MOI.SecondOrderCone)
    return Mosek.appendquadraticconedomain(t, n)
end

function appendconedomain(t::Mosek.Task, n::Int, ::MOI.RotatedSecondOrderCone)
    return Mosek.appendrquadraticconedomain(t, n)
end

function appendconedomain(
    t::Mosek.Task,
    n::Int,
    ::MOI.Scaled{MOI.PositiveSemidefiniteConeTriangle},
)
    return Mosek.appendsvecpsdconedomain(t, n)
end

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
#flag(::Type{<:VectorCone}) = 0x40 # FIXME unused
incompatible_mask(::Type{<:VectorCone}) = 0x40

function set_flag(model::Optimizer, vi::MOI.VariableIndex, S::Type)
    return model.x_constraints[vi.value] |= flag(S)
end

function unset_flag(model::Optimizer, vi::MOI.VariableIndex, S::Type)
    return model.x_constraints[vi.value] &= ~flag(S)
end

function has_flag(model::Optimizer, vi::MOI.VariableIndex, S::Type)
    return !iszero(model.x_constraints[vi.value] & flag(S))
end

###############################################################################
# MOI #########################################################################
###############################################################################

const ScalarLinearDomain = Union{
    MOI.LessThan{Float64},
    MOI.GreaterThan{Float64},
    MOI.EqualTo{Float64},
    MOI.Interval{Float64},
}

## Add ########################################################################
###############################################################################

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{<:Union{MOI.VariableIndex,MOI.ScalarAffineFunction}},
    ::Type{<:ScalarLinearDomain},
)
    return true
end

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VectorOfVariables},
    ::Type{<:VectorCone},
)
    return true
end

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VariableIndex},
    ::Type{MOI.Integer},
)
    return true
end

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VectorAffineFunction{Float64}},
    ::Type{<:VectorConeDomain},
)
    return true
end

function MOI.supports_add_constrained_variables(
    ::Optimizer,
    ::Type{MOI.PositiveSemidefiniteConeTriangle},
)
    return true
end

## Affine Constraints #########################################################
##################### lc ≤ Ax ≤ uc ############################################

function MOI.add_constraint(
    m::Optimizer,
    axb::MOI.ScalarAffineFunction{Float64},
    dom::D,
) where {D<:MOI.AbstractScalarSet}
    if !iszero(axb.constant)
        throw(
            MOI.ScalarFunctionConstantNotZero{Float64,typeof(axb),D}(
                axb.constant,
            ),
        )
    end
    # Duplicate indices not supported
    axb = MOI.Utilities.canonical(axb)
    N = 1
    conid = allocateconstraints(m, N)
    ci = MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},D}(conid)
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

cone_type(::Type{MOI.ExponentialCone}) = Mosek.MSK_CT_PEXP
cone_type(::Type{MOI.DualExponentialCone}) = Mosek.MSK_CT_DEXP
cone_type(::Type{MOI.PowerCone{Float64}}) = Mosek.MSK_CT_PPOW
cone_type(::Type{MOI.DualPowerCone{Float64}}) = Mosek.MSK_CT_DPOW
cone_type(::Type{MOI.SecondOrderCone}) = Mosek.MSK_CT_QUAD
cone_type(::Type{MOI.RotatedSecondOrderCone}) = Mosek.MSK_CT_RQUAD

function MOI.add_constraint(
    m::Optimizer,
    xs::MOI.VariableIndex,
    dom::D,
) where {D<:Union{ScalarLinearDomain,MOI.Integer}}
    msk_idx = mosek_index(m, xs)
    if !(msk_idx isa ColumnIndex)
        error("Cannot add $D constraint on a matrix variable")
    end
    if !iszero(incompatible_mask(D) & m.x_constraints[xs.value])
        error("Cannot put multiple bound sets of the same type on a variable")
    end
    set_flag(m, xs, D)
    add_variable_constraint(m, msk_idx, dom)
    return MOI.ConstraintIndex{MOI.VariableIndex,D}(xs.value)
end

function MOI.add_constraint(
    m::Optimizer,
    xs::MOI.VectorOfVariables,
    dom::D,
) where {D<:VectorCone}
    if any(vi -> is_matrix(m, vi), xs.variables)
        error("Cannot add $D constraint on a matrix variable")
    end
    cols = ColumnIndices(reorder(columns(m, xs.variables).values, D, true))
    if !all(
        vi -> iszero(incompatible_mask(D) & m.x_constraints[vi.value]),
        xs.variables,
    )
        error("Cannot multiple bound sets of the same type to a variable")
    end
    id = add_cone(m, cols, dom)
    idx = first(xs.variables).value
    for vi in xs.variables
        m.variable_to_vector_constraint_id[vi.value] = -idx
    end
    m.variable_to_vector_constraint_id[idx] = id
    return MOI.ConstraintIndex{MOI.VectorOfVariables,D}(idx)
end

function MOI.add_constraint(
    m::Optimizer,
    func::MOI.VectorAffineFunction{Float64},
    dom::D,
) where {D<:VectorConeDomain}
    # if any(vi -> is_matrix(m, vi), xs.variables)
    #     error("Cannot add $D constraint on a matrix variable")
    # end
    axbs = MOI.Utilities.canonical(func)
    let acci = Mosek.getnumacc(m.task) + 1,
        afei = Mosek.getnumafe(m.task),
        b = -reorder(axbs.constants, D, true),
        num = length(axbs.constants),
        nnz = length(axbs.terms),
        domi = appendconedomain(m.task, num, dom)

        m.F_rows[acci] = afei .+ eachindex(b)
        Mosek.appendafes(m.task, num)
        Mosek.appendaccseq(m.task, domi, afei + 1, b)
        rsubi = Vector{Int64}()
        sizehint!(rsubi, nnz)
        rsubj = Vector{Int32}()
        sizehint!(rsubj, nnz)
        rcof = Vector{Float64}()
        sizehint!(rcof, nnz)
        rbarsubi = Vector{Int64}()
        sizehint!(rbarsubi, nnz)
        rbarsubj = Vector{Int32}()
        sizehint!(rbarsubj, nnz)
        rbarsubk = Vector{Int64}()
        sizehint!(rbarsubk, nnz)
        rbarsubl = Vector{Int64}()
        sizehint!(rbarsubl, nnz)
        rbarcof = Vector{Float64}()
        sizehint!(rbarcof, nnz)
        function add(row::Int, col::ColumnIndex, coefficient::Float64)
            push!(rsubi, row)
            push!(rsubj, col.value)
            return push!(rcof, coefficient)
        end
        function add(row::Int, mat::MatrixIndex, coefficient::Float64)
            push!(rbarsubi, row)
            push!(rbarsubj, mat.matrix)
            push!(rbarsubk, mat.row)
            push!(rbarsubl, mat.column)
            return push!(
                rbarcof,
                mat.row == mat.column ? coefficient : coefficient / 2,
            )
        end
        for term in axbs.terms
            add(
                reorder(term.output_index, dom, true) + afei,
                mosek_index(m, term.scalar_term.variable),
                term.scalar_term.coefficient,
            )
        end
        if !isempty(rsubi) # Mosek segfaults otherwise, see https://github.com/jump-dev/MosekTools.jl/actions/runs/3243196430/jobs/5317555832#step:7:132
            Mosek.putafefentrylist(m.task, rsubi, rsubj, rcof)
        end
        if !isempty(rbarsubi)
            Mosek.putafebarfblocktriplet(
                m.task,
                rbarsubi,
                rbarsubj,
                rbarsubk,
                rbarsubl,
                rbarcof,
            )
        end
        MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64},D}(acci)
    end
    # cols = ColumnIndices(reorder(columns(m, xs.variables).values, D))
    # id = add_cone(m, cols, dom)
    # idx = first(xs.variables).value
    # for vi in xs.variables
    #     m.variable_to_vector_constraint_id[vi.value] = -idx
    # end
    # m.variable_to_vector_constraint_id[idx] = id
    # ci = MOI.ConstraintIndex{MOI.VectorOfVariables, D}(idx)
    # return ci
end

function cone_id(
    model::Optimizer,
    ci::MOI.ConstraintIndex{MOI.VectorOfVariables},
)
    return model.variable_to_vector_constraint_id[ci.value]
end

################################################################################
################################################################################

function MOI.add_constrained_variables(
    m::Optimizer,
    dom::MOI.PositiveSemidefiniteConeTriangle,
)
    N = MOI.side_dimension(dom)
    if N < 1
        error(
            "Invalid dimension for semidefinite constraint, got $N which is ",
            "smaller than the minimum dimension 1.",
        )
    end
    Mosek.appendbarvars(m.task, [Int32(N)])
    push!(m.sd_dim, N)
    id = length(m.sd_dim)
    vis =
        [new_variable_index(m, MatrixIndex(id, i, j)) for i in 1:N for j in 1:i]
    con_idx = MOI.ConstraintIndex{
        MOI.VectorOfVariables,
        MOI.PositiveSemidefiniteConeTriangle,
    }(
        id,
    )
    return vis, con_idx
end

## Get ########################################################################
###############################################################################

function _variable(ci::MOI.ConstraintIndex{MOI.VariableIndex})
    return MOI.VariableIndex(ci.value)
end

function MOI.get(
    m::Optimizer,
    ::MOI.ConstraintFunction,
    ci::MOI.ConstraintIndex{MOI.VariableIndex},
)
    MOI.throw_if_not_valid(m, ci)
    return _variable(ci)
end

function MOI.get(
    m::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,S},
) where {S<:MOI.Integer}
    MOI.throw_if_not_valid(m, ci)
    return S()
end

function MOI.get(
    m::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,S},
) where {S<:ScalarLinearDomain}
    MOI.throw_if_not_valid(m, ci)
    vi = MOI.get(m, MOI.ConstraintFunction(), ci)
    return get_variable_constraint(m, vi, ci)
end

function MOI.set(
    ::Optimizer,
    ::MOI.ConstraintFunction,
    ci::MOI.ConstraintIndex{MOI.VariableIndex},
    ::MOI.VariableIndex,
)
    throw(MOI.SettingVariableIndexNotAllowed())
end

function MOI.get(
    m::Optimizer,
    attr::MOI.ConstraintFunction,
    ci::MOI.ConstraintIndex{
        MOI.ScalarAffineFunction{Float64},
        <:ScalarLinearDomain,
    },
)
    MOI.throw_if_not_valid(m, ci)
    if length(m.sd_dim) > 0
        throw(MOI.GetAttributeNotAllowed(attr))
    end
    nnz, cols, vals = Mosek.getarow(m.task, row(m, ci))
    @assert nnz == length(cols) == length(vals)
    terms = MOI.ScalarAffineTerm{Float64}[
        MOI.ScalarAffineTerm(vals[i], index_of_column(m, cols[i])) for
        i in 1:nnz
    ]
    return MOI.ScalarAffineFunction(terms, 0.0)
end

function MOI.get(
    m::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}},
)
    MOI.throw_if_not_valid(m, ci)
    return get_bound(m, ci)
end

function MOI.get(
    m::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64}},
)
    MOI.throw_if_not_valid(m, ci)
    return get_bound(m, ci)
end

function MOI.get(
    m::Optimizer,
    ::MOI.ConstraintFunction,
    ci::MOI.ConstraintIndex{MOI.VectorOfVariables,S},
) where {S<:VectorCone}
    return MOI.VectorOfVariables([
        index_of_column(m, col) for
        col in reorder(columns(m, ci).values, S, true)
    ])
end

function MOI.get(
    m::Optimizer,
    attr::MOI.ConstraintFunction,
    ci::MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64},S},
) where {S<:VectorConeDomain}
    if length(m.sd_dim) > 0
        # Cannot get function if there are matrix variables
        throw(MOI.GetAttributeNotAllowed(attr))
    end
    r = rows(m, ci)
    (frow, fcol, fval) = Mosek.getaccftrip(m.task)
    constants = Mosek.getaccb(m.task, ci.value)
    set = MOI.Utilities.set_with_dimension(S, length(r))
    terms = [
        MOI.VectorAffineTerm(
            reorder(frow[i] - first(r) + 1, set, false),
            MOI.ScalarAffineTerm(fval[i], index_of_column(m, fcol[i])),
        ) for i in eachindex(frow) if frow[i] in r
    ]
    return MOI.VectorAffineFunction(terms, -reorder(constants, S, false))
end

function type_cone(ct)
    if ct == Mosek.MSK_CT_PEXP
        return MOI.ExponentialCone
    elseif ct == Mosek.MSK_CT_DEXP
        return MOI.DualExponentialCone
    elseif ct == Mosek.MSK_CT_PPOW
        return MOI.PowerCone{Float64}
    elseif ct == Mosek.MSK_CT_DPOW
        return MOI.DualPowerCone{Float64}
    elseif ct == Mosek.MSK_CT_QUAD
        return MOI.SecondOrderCone
    elseif ct == Mosek.MSK_CT_RQUAD
        return MOI.RotatedSecondOrderCone
    else
        error("Unrecognized Mosek cone type `$ct`.")
    end
end

cone(::Type{MOI.ExponentialCone}, conepar, nummem) = MOI.ExponentialCone()

function cone(::Type{MOI.DualExponentialCone}, conepar, nummem)
    return MOI.DualExponentialCone()
end

cone(::Type{MOI.PowerCone{Float64}}, conepar, nummem) = MOI.PowerCone(conepar)

function cone(::Type{MOI.DualPowerCone{Float64}}, conepar, nummem)
    return MOI.DualPowerCone(conepar)
end

cone(::Type{MOI.SecondOrderCone}, conepar, nummem) = MOI.SecondOrderCone(nummem)

function cone(::Type{MOI.RotatedSecondOrderCone}, conepar, nummem)
    return MOI.RotatedSecondOrderCone(nummem)
end

function MOI.get(
    m::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.VectorOfVariables,<:VectorCone},
)
    MOI.throw_if_not_valid(m, ci)
    ct, conepar, nummem = Mosek.getconeinfo(m.task, cone_id(m, ci))
    return cone(type_cone(ct), conepar, nummem)
end

## Modify #####################################################################
###############################################################################

### SET
function chgbound(
    bl::Float64,
    bu::Float64,
    k::Float64,
    dom::MOI.LessThan{Float64},
)
    return bl, dom.upper - k
end

function chgbound(
    bl::Float64,
    bu::Float64,
    k::Float64,
    dom::MOI.GreaterThan{Float64},
)
    return dom.lower - k, bu
end

function chgbound(
    bl::Float64,
    bu::Float64,
    k::Float64,
    dom::MOI.EqualTo{Float64},
)
    return dom.value - k, dom.value - k
end

function chgbound(
    bl::Float64,
    bu::Float64,
    k::Float64,
    dom::MOI.Interval{Float64},
)
    return dom.lower - k, dom.upper - k
end

function MOI.set(
    m::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,D},
    dom::D,
) where {D<:ScalarLinearDomain}
    col = column(m, _variable(ci))
    bk, bl, bu = Mosek.getvarbound(m.task, col.value)
    bl, bu = chgbound(bl, bu, 0.0, dom)
    return Mosek.putvarbound(m.task, col.value, bk, bl, bu)
end

function MOI.set(
    m::Optimizer,
    ::MOI.ConstraintSet,
    cref::MOI.ConstraintIndex{<:MOI.ScalarAffineFunction,D},
    dom::D,
) where {D<:ScalarLinearDomain}
    cid = ref2id(cref)
    i = getindex(m.c_block, cid) # since we are in a scalar domain
    bk, bl, bu = Mosek.getconbound(m.task, i)
    bl, bu = chgbound(bl, bu, 0.0, dom)
    return Mosek.putconbound(m.task, i, bk, bl, bu)
end

### MODIFY
function MOI.modify(
    m::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},D},
    func::MOI.ScalarConstantChange{Float64},
) where {D<:MOI.AbstractSet}
    if !iszero(func.new_constant)
        throw(
            MOI.ScalarFunctionConstantNotZero{
                Float64,
                MOI.ScalarAffineFunction{Float64},
                D,
            }(
                func.new_constant,
            ),
        )
    end
end

# Is there a way to do this in Mosek API ? I haven't check so here is an error for now:
const _MODIFY_PSD_VAR_ERROR = "Modifying the coefficient of the variable correspond to an entry of a PSD matrix is not supported"

function MOI.modify(
    m::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}},
    func::MOI.ScalarCoefficientChange{Float64},
)
    if is_matrix(m, func.variable)
        throw(MOI.ModifyConstraintNotAllowed(c, func, _MODIFY_PSD_VAR_ERROR))
    end
    return set_coefficient(m, row(m, c), func.variable, func.new_coefficient)
end

### TRANSFORM
function MOI.transform(
    m::Optimizer,
    cref::MOI.ConstraintIndex{F,D},
    newdom::D,
) where {F<:MOI.AbstractFunction,D<:MOI.AbstractSet}
    MOI.modify(m, cref, newdom)
    return cref
end

function MOI.transform(
    m::Optimizer,
    cref::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},D1},
    newdom::D2,
) where {D1<:ScalarLinearDomain,D2<:ScalarLinearDomain}
    F = MOI.ScalarAffineFunction{Float64}
    cid = ref2id(cref)
    r = row(m, cref)
    add_bound(m, r, newdom)
    return MOI.ConstraintIndex{F,D2}(cid)
end

## Delete #####################################################################
###############################################################################

function MOI.is_valid(
    model::Optimizer,
    ci::MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64},S},
) where {S<:VectorConeDomain}
    numacc = Mosek.getnumacc(model.task)
    if ci.value < 1 || ci.value > numacc
        false
    else
        domidx = Mosek.getaccdomain(model.task, ci.value)
        if domidx == 1
            false
        else
            dt = Mosek.getdomaintype(model.task, domidx)
            (
                (
                    dt == Mosek.MSK_DOMAIN_QUADRATIC_CONE &&
                    S == MOI.SecondOrderCone
                ) ||
                (
                    dt == Mosek.MSK_DOMAIN_RQUADRATIC_CONE &&
                    S == MOI.RotatedSecondOrderCone
                ) ||
                (
                    dt == Mosek.MSK_DOMAIN_PRIMAL_EXP_CONE &&
                    S == MOI.ExponentialCone
                ) ||
                (
                    dt == Mosek.MSK_DOMAIN_DUAL_EXP_CONE &&
                    S == MOI.DualExponentialCone
                ) ||
                (
                    dt == Mosek.MSK_DOMAIN_PRIMAL_POWER_CONE &&
                    S == MOI.PowerCone{Float64}
                ) ||
                (
                    dt == Mosek.MSK_DOMAIN_DUAL_POWER_CONE &&
                    S == MOI.DualPowerCone{Float64}
                ) ||
                (
                    dt == Mosek.MSK_DOMAIN_SVEC_PSD_CONE &&
                    S == MOI.Scaled{MOI.PositiveSemidefiniteConeTriangle}
                )
            )
        end
    end
end

# Commenting out as it doesn't work (gives a MethodError)

# Deleting a constraint block means clearing non-zeros from the its
# AFE rows and resetting the underlying ACC to an empty domain. We do
# not reclaim the AFEs.
# function MOI.delete(m::Optimizer,
#                     cref::MOI.ConstraintIndex{F,D}) where {F <: MOI.VectorAffineFunction{Float64},
#                                                            D <: VectorConeDomain}
#     MOI.throw_if_not_valid(m, cref)
#     Mosek.putaccname(m.task,cref.value,"")
#     afeidxs = Mosek.getaccafeidxlist(m.task,cref.value)
#     # clear afe non-zeros, but don't delete or reuse afe idx
#     # FIXME gives a MethodError
#     Mosek.putafefrowlist(afeidxs,zeros(Int32,length(afeidxs)),zeros(Int64,length(afeidxs)),Int32[],Float64[])
#     putaccdom(m.task,
#               cref.value,
#               1, # the empty zero domain,
#               Int64[],
#               Float64[])
# end

function MOI.is_valid(
    model::Optimizer,
    ci::MOI.ConstraintIndex{<:MOI.ScalarAffineFunction{Float64},S},
) where {S<:ScalarLinearDomain}
    return allocated(model.c_block, ci.value) &&
           Mosek.getconbound(model.task, row(model, ci))[1] == bound_key(S)
end

function MOI.delete(
    m::Optimizer,
    cref::MOI.ConstraintIndex{F,D},
) where {F<:MOI.ScalarAffineFunction{Float64},D<:ScalarLinearDomain}
    MOI.throw_if_not_valid(m, cref)
    delete_name(m, cref)
    subi = getindexes(m.c_block, cref.value)
    n = length(subi)
    subi_i32 = convert(Vector{Int32}, subi)
    ptr = fill(Int64(0), n)
    Mosek.putarowlist(m.task, subi_i32, ptr, ptr, Int32[], Float64[])
    b = fill(0.0, n)
    Mosek.putconboundlist(m.task, subi_i32, fill(Mosek.MSK_BK_FX, n), b, b)
    return deleteblock(m.c_block, cref.value)
end

function MOI.is_valid(
    model::Optimizer,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,S},
) where {S<:Union{ScalarLinearDomain,MOI.Integer}}
    return allocated(model.x_block, ci.value) &&
           has_flag(model, _variable(ci), S)
end

function MOI.delete(
    m::Optimizer,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,S},
) where {S<:Union{ScalarLinearDomain,MOI.Integer}}
    MOI.throw_if_not_valid(m, ci)
    delete_name(m, ci)
    vi = _variable(ci)
    unset_flag(m, vi, S)
    return delete_variable_constraint(m, column(m, vi), S)
end

function MOI.is_valid(
    model::Optimizer,
    ci::MOI.ConstraintIndex{MOI.VectorOfVariables,S},
) where {S<:VectorCone}
    if !(ci.value in eachindex(model.variable_to_vector_constraint_id))
        return false
    end
    id = cone_id(model, ci)
    return 1 ≤ id ≤ Mosek.getnumcone(model.task) &&
           Mosek.getconeinfo(model.task, id)[1] == cone_type(S)
end

function MOI.delete(
    model::Optimizer,
    ci::MOI.ConstraintIndex{MOI.VectorOfVariables,<:VectorCone},
)
    id = cone_id(model, ci)
    for vi in MOI.get(model, MOI.ConstraintFunction(), ci).variables
        model.variable_to_vector_constraint_id[vi.value] = 0
    end
    Mosek.removecones(model.task, [id])
    # The conic constraints with id higher than `id` are renumbered.
    return map!(
        i -> i > id ? i - 1 : i,
        model.variable_to_vector_constraint_id,
        model.variable_to_vector_constraint_id,
    )
end

function MOI.is_valid(
    model::Optimizer,
    ci::MOI.ConstraintIndex{
        MOI.VectorOfVariables,
        MOI.PositiveSemidefiniteConeTriangle,
    },
)
    # TODO add supports for deletion
    return 1 ≤ ci.value ≤ length(model.sd_dim)
end

## List #######################################################################
###############################################################################
function MOI.get(m::Optimizer, ::MOI.ListOfConstraintAttributesSet)
    set = MOI.AbstractConstraintAttribute[]
    if !isempty(m.constrnames)
        push!(set, MOI.ConstraintName())
    end
    # TODO add VariablePrimalStart when get is implemented on it
    return set
end

## Name #######################################################################
###############################################################################
function MOI.supports(
    ::Optimizer,
    ::MOI.ConstraintName,
    ::Type{<:MOI.ConstraintIndex},
)
    return true
end

function MOI.set(
    m::Optimizer,
    ::MOI.ConstraintName,
    ci::MOI.ConstraintIndex,
    name::AbstractString,
)
    delete_name(m, ci)
    if !haskey(m.constrnames, name)
        m.constrnames[name] = MOI.ConstraintIndex[]
    end
    push!(m.constrnames[name], ci)
    return set_row_name(m, ci, name)
end

function MOI.get(
    m::Optimizer,
    ::MOI.ConstraintName,
    ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}},
)
    # All rows should have same name so we take the first one
    return Mosek.getconname(m.task, getindexes(m.c_block, ref2id(ci))[1])
end

function MOI.get(m::Optimizer, ::MOI.ConstraintName, ci::MOI.ConstraintIndex)
    return get(m.con_to_name, ci, "")
end

function MOI.get(m::Optimizer, CI::Type{<:MOI.ConstraintIndex}, name::String)
    #asgn, row = Mosek.getconnameindex(m.task, name)
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
