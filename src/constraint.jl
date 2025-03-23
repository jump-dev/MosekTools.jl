# Copyright (c) 2017: Ulf Worsøe, Mosek ApS
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

## Affine Constraints #########################################################
##################### lc <= Ax <= uc ############################################
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

_bound_key(::Type{MOI.GreaterThan{Float64}}) = Mosek.MSK_BK_LO
_bound_key(::Type{MOI.LessThan{Float64}}) = Mosek.MSK_BK_UP
_bound_key(::Type{MOI.EqualTo{Float64}}) = Mosek.MSK_BK_FX
_bound_key(::Type{MOI.Interval{Float64}}) = Mosek.MSK_BK_RA

function _bounds(dom::MOI.Interval)
    bl = dom.lower
    bu = dom.upper
    bk = _bound_key(typeof(dom))
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

function _putconbound(m::Optimizer, row::Int32, dom::MOI.GreaterThan{Float64})
    Mosek.putconbound(
        m.task,
        row,
        _bound_key(typeof(dom)),
        dom.lower,
        dom.lower,
    )
    return
end

function _putconbound(m::Optimizer, row::Int32, dom::MOI.LessThan{Float64})
    Mosek.putconbound(
        m.task,
        row,
        _bound_key(typeof(dom)),
        dom.upper,
        dom.upper,
    )
    return
end

function _putconbound(m::Optimizer, row::Int32, dom::MOI.EqualTo{Float64})
    Mosek.putconbound(
        m.task,
        row,
        _bound_key(typeof(dom)),
        dom.value,
        dom.value,
    )
    return
end

function _putconbound(m::Optimizer, row::Int32, dom::MOI.Interval{Float64})
    bk, bl, bu = _bounds(dom)
    Mosek.putconbound(m.task, row, bk, bl, bu)
    return
end

function _bounds_to_set(::Type{S}, bk, bl, bu) where {S}
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

function _get_bound(
    m::Optimizer,
    ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S},
) where {S}
    return _bounds_to_set(S, Mosek.getconbound(m.task, row(m, ci))...)
end

function _get_bound(
    m::Optimizer,
    ci::MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64},S},
) where {S<:VectorConeDomain}
    domidx = Mosek.getaccdomain(m.task, ci.value)
    dt = Mosek.getdomaintype(m.task, domidx)
    if dt == Mosek.MSK_DOMAIN_QUADRATIC_CONE
        return MOI.SecondOrderCone(solsize(m, ci))
    elseif dt == Mosek.MSK_DOMAIN_RQUADRATIC_CONE
        return MOI.RotatedSecondOrderCone(solsize(m, ci))
    elseif dt == Mosek.MSK_DOMAIN_PRIMAL_EXP_CONE
        return MOI.ExponentialCone()
    elseif dt == Mosek.MSK_DOMAIN_DUAL_EXP_CONE
        return MOI.DualExponentialCone()
    elseif dt == Mosek.MSK_DOMAIN_PRIMAL_POWER_CONE
        n, _ = Mosek.getpowerdomaininfo(m.task, domidx)
        @assert n == 3
        alpha = Mosek.getpowerdomainalpha(m.task, domidx)
        return MOI.PowerCone(alpha[1])
    elseif dt == Mosek.MSK_DOMAIN_DUAL_POWER_CONE
        n, _ = Mosek.getpowerdomaininfo(m.task, domidx)
        @assert n == 3
        alpha = Mosek.getpowerdomainalpha(m.task, domidx)
        return MOI.DualPowerCone(alpha[1])
    else
        @assert dt == Mosek.MSK_DOMAIN_SVEC_PSD_CONE
        return MOI.Scaled{MOI.PositiveSemidefiniteConeTriangle}(
            MOI.Utilities.side_dimension_for_vectorized_dimension(
                solsize(m, ci),
            ),
        )
    end
end

## Variable Constraints #######################################################
####################### lx <= x <= u ############################################
#######################      x ∈ K ############################################

function _delete_variable_constraint(
    m::Optimizer,
    col::ColumnIndex,
    ::Type{<:Union{MOI.Interval,MOI.EqualTo}},
)
    Mosek.putvarbound(m.task, col.value, Mosek.MSK_BK_FR, 0.0, 0.0)
    return
end

function _delete_variable_constraint(
    m::Optimizer,
    col::ColumnIndex,
    ::Type{MOI.Integer},
)
    Mosek.putvartype(m.task, col.value, Mosek.MSK_VAR_TYPE_CONT)
    return
end

function _delete_variable_constraint(
    m::Optimizer,
    col::ColumnIndex,
    ::Type{MOI.LessThan{Float64}},
)
    bk, lo, _ = Mosek.getvarbound(m.task, col.value)
    if bk == Mosek.MSK_BK_UP
        Mosek.putvarbound(m.task, col.value, Mosek.MSK_BK_FR, lo, 0.0)
    else
        @assert bk == Mosek.MSK_BK_RA
        Mosek.putvarbound(m.task, col.value, Mosek.MSK_BK_LO, lo, 0.0)
    end
    return
end

function _delete_variable_constraint(
    m::Optimizer,
    col::ColumnIndex,
    ::Type{MOI.GreaterThan{Float64}},
)
    bk, _, up = Mosek.getvarbound(m.task, col.value)
    if bk == Mosek.MSK_BK_LO
        Mosek.putvarbound(m.task, col.value, Mosek.MSK_BK_FR, 0.0, up)
    else
        @assert bk == Mosek.MSK_BK_RA
        Mosek.putvarbound(m.task, col.value, Mosek.MSK_BK_UP, 0.0, up)
    end
    return
end

function _add_variable_constraint(
    m::Optimizer,
    col::ColumnIndex,
    dom::MOI.Interval,
)
    bk, bl, bu = _bounds(dom)
    Mosek.putvarbound(m.task, col.value, bk, bl, bu)
    return
end

function _add_variable_constraint(
    m::Optimizer,
    col::ColumnIndex,
    dom::MOI.EqualTo,
)
    Mosek.putvarbound(m.task, col.value, Mosek.MSK_BK_FX, dom.value, dom.value)
    return
end

function _add_variable_constraint(m::Optimizer, col::ColumnIndex, ::MOI.Integer)
    Mosek.putvartype(m.task, col.value, Mosek.MSK_VAR_TYPE_INT)
    return
end

function _add_variable_constraint(
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
    Mosek.putvarbound(m.task, col.value, bk, lo, dom.upper)
    return
end

function _add_variable_constraint(
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
    Mosek.putvarbound(m.task, col.value, bk, dom.lower, up)
    return
end

function get_variable_constraint(
    m::Optimizer,
    col::ColumnIndex,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,S},
) where {S}
    return _bounds_to_set(S, Mosek.getvarbound(m.task, col.value)...)
end

function get_variable_constraint(
    m::Optimizer,
    vi::MOI.VariableIndex,
    ci::MOI.ConstraintIndex,
)
    return get_variable_constraint(m, mosek_index(m, vi), ci)
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

flag(::Type{MOI.EqualTo{Float64}}) = 0x1
flag(::Type{MOI.GreaterThan{Float64}}) = 0x2
flag(::Type{MOI.LessThan{Float64}}) = 0x4
flag(::Type{MOI.Interval{Float64}}) = 0x8
flag(::Type{MOI.Integer}) = 0x10

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
##################### lc <= Ax <= uc ############################################

function MOI.add_constraint(
    m::Optimizer,
    f::MOI.ScalarAffineFunction{Float64},
    set::S,
) where {S<:MOI.AbstractScalarSet}
    if !iszero(f.constant)
        throw(
            MOI.ScalarFunctionConstantNotZero{Float64,typeof(f),S}(f.constant),
        )
    end
    # Duplicate indices not supported
    f = MOI.Utilities.canonical(f)
    conid = allocateconstraints(m, 1)
    ci = MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S}(conid)
    r = row(m, ci)
    cols, values = split_scalar_matrix(
        m,
        f.terms,
        (j, ids, coefs) -> Mosek.putbaraij(m.task, r, j, ids, coefs),
    )
    Mosek.putarow(m.task, r, ColumnIndices(cols).values, values)
    _putconbound(m, r, set)
    return ci
end

## Variable Constraints #######################################################
####################### lx <= x <= u ############################################
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

function _check_bound_compat(
    m::Optimizer,
    x::MOI.VariableIndex,
    set::MOI.LessThan{Float64},
)
    if has_flag(m, x, MOI.EqualTo{Float64})
        throw(MOI.UpperBoundAlreadySet{MOI.EqualTo{Float64},typeof(set)}(x))
    elseif has_flag(m, x, MOI.Interval{Float64})
        throw(MOI.UpperBoundAlreadySet{MOI.Interval{Float64},typeof(set)}(x))
    elseif has_flag(m, x, MOI.LessThan{Float64})
        throw(MOI.UpperBoundAlreadySet{MOI.LessThan{Float64},typeof(set)}(x))
    end
    return
end

function _check_bound_compat(
    m::Optimizer,
    x::MOI.VariableIndex,
    set::MOI.GreaterThan{Float64},
)
    if has_flag(m, x, MOI.EqualTo{Float64})
        throw(MOI.LowerBoundAlreadySet{MOI.EqualTo{Float64},typeof(set)}(x))
    elseif has_flag(m, x, MOI.Interval{Float64})
        throw(MOI.LowerBoundAlreadySet{MOI.Interval{Float64},typeof(set)}(x))
    elseif has_flag(m, x, MOI.GreaterThan{Float64})
        throw(MOI.LowerBoundAlreadySet{MOI.GreaterThan{Float64},typeof(set)}(x))
    end
    return
end

function _check_bound_compat(
    m::Optimizer,
    x::MOI.VariableIndex,
    set::Union{MOI.EqualTo,MOI.Interval},
)
    if has_flag(m, x, MOI.EqualTo{Float64})
        throw(MOI.LowerBoundAlreadySet{MOI.EqualTo{Float64},typeof(set)}(x))
    elseif has_flag(m, x, MOI.Interval{Float64})
        throw(MOI.LowerBoundAlreadySet{MOI.Interval{Float64},typeof(set)}(x))
    elseif has_flag(m, x, MOI.LessThan{Float64})
        throw(MOI.UpperBoundAlreadySet{MOI.LessThan{Float64},typeof(set)}(x))
    elseif has_flag(m, x, MOI.GreaterThan{Float64})
        throw(MOI.LowerBoundAlreadySet{MOI.GreaterThan{Float64},typeof(set)}(x))
    end
    return
end

_check_bound_compat(::Optimizer, ::MOI.VariableIndex, ::MOI.Integer) = nothing

function MOI.add_constraint(
    m::Optimizer,
    x::MOI.VariableIndex,
    set::S,
) where {S<:Union{ScalarLinearDomain,MOI.Integer}}
    index = mosek_index(m, x)
    if index isa MatrixIndex
        msg = "Cannot add $S constraint on a matrix variable."
        throw(MOI.AddConstraintNotAllowed{MOI.VariableIndex,S}(msg))
    end
    _check_bound_compat(m, x, set)
    set_flag(m, x, S)
    _add_variable_constraint(m, index, set)
    return MOI.ConstraintIndex{MOI.VariableIndex,S}(x.value)
end

_cone_type(::Type{MOI.ExponentialCone}) = Mosek.MSK_CT_PEXP
_cone_type(::Type{MOI.DualExponentialCone}) = Mosek.MSK_CT_DEXP
_cone_type(::Type{MOI.PowerCone{Float64}}) = Mosek.MSK_CT_PPOW
_cone_type(::Type{MOI.DualPowerCone{Float64}}) = Mosek.MSK_CT_DPOW
_cone_type(::Type{MOI.SecondOrderCone}) = Mosek.MSK_CT_QUAD
_cone_type(::Type{MOI.RotatedSecondOrderCone}) = Mosek.MSK_CT_RQUAD

_cone_parameter(dom::MOI.PowerCone{Float64}) = dom.exponent
_cone_parameter(dom::MOI.DualPowerCone{Float64}) = dom.exponent
_cone_parameter(::MOI.AbstractSet) = 0.0

function MOI.add_constraint(
    m::Optimizer,
    xs::MOI.VectorOfVariables,
    dom::D,
) where {D<:VectorCone}
    if any(vi -> is_matrix(m, vi), xs.variables)
        msg = "Cannot add $D constraint on a matrix variable"
        throw(MOI.AddConstraintNotAllowed{MOI.VectorOfVariables,D}(msg))
    end
    cols = ColumnIndices(reorder(columns(m, xs.variables).values, D, true))
    Mosek.appendcone(m.task, _cone_type(D), _cone_parameter(dom), cols.values)
    id = Mosek.getnumcone(m.task)
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
        return MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64},D}(acci)
    end
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
    dom::S,
) where {S<:MOI.PositiveSemidefiniteConeTriangle}
    N = MOI.side_dimension(dom)
    @assert N >= 1
    Mosek.appendbarvars(m.task, Int32[N])
    push!(m.sd_dim, N)
    id = length(m.sd_dim)
    x = [new_variable_index(m, MatrixIndex(id, i, j)) for i in 1:N for j in 1:i]
    return x, MOI.ConstraintIndex{MOI.VectorOfVariables,S}(id)
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
    ci::MOI.ConstraintIndex{MOI.VariableIndex,MOI.Integer},
)
    MOI.throw_if_not_valid(m, ci)
    return MOI.Integer()
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
    return throw(MOI.SettingVariableIndexNotAllowed())
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
        MOI.ScalarAffineTerm(v, _col_to_index(m, c)) for
        (v, c) in zip(vals, cols)
    ]
    return MOI.ScalarAffineFunction(terms, 0.0)
end

function MOI.get(
    m::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}},
)
    MOI.throw_if_not_valid(m, ci)
    return _get_bound(m, ci)
end

function MOI.get(
    m::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64}},
)
    MOI.throw_if_not_valid(m, ci)
    return _get_bound(m, ci)
end

function MOI.get(
    m::Optimizer,
    ::MOI.ConstraintFunction,
    ci::MOI.ConstraintIndex{MOI.VectorOfVariables,S},
) where {S<:VectorCone}
    cols = reorder(columns(m, ci).values, S, true)
    return MOI.VectorOfVariables(_col_to_index.(m, cols))
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
    frow, fcol, fval = Mosek.getaccftrip(m.task)
    constants = Mosek.getaccb(m.task, ci.value)
    set = MOI.Utilities.set_with_dimension(S, length(r))
    terms = MOI.VectorAffineTerm{Float64}[]
    for (frowi, fcoli, fvali) in zip(frow, fcol, fval)
        if frowi in r
            row = reorder(frowi - first(r) + 1, set, false)
            term = MOI.ScalarAffineTerm(fvali, _col_to_index(m, fcoli))
            push!(terms, MOI.VectorAffineTerm(row, term))
        end
    end
    return MOI.VectorAffineFunction(terms, -reorder(constants, S, false))
end

function _type_cone(ct)
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
    else
        @assert ct == Mosek.MSK_CT_RQUAD
        return MOI.RotatedSecondOrderCone
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
    return cone(_type_cone(ct), conepar, nummem)
end

## Modify #####################################################################
###############################################################################

_change_bound(bl, bu, dom::MOI.LessThan{Float64}) = bl, dom.upper

_change_bound(bl, bu, dom::MOI.GreaterThan{Float64}) = dom.lower, bu

_change_bound(bl, bu, dom::MOI.EqualTo{Float64}) = dom.value, dom.value

_change_bound(bl, bu, dom::MOI.Interval{Float64}) = dom.lower, dom.upper

function MOI.set(
    m::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,D},
    dom::D,
) where {D<:ScalarLinearDomain}
    col = column(m, _variable(ci))
    bk, bl, bu = Mosek.getvarbound(m.task, col.value)
    bl, bu = _change_bound(bl, bu, dom)
    Mosek.putvarbound(m.task, col.value, bk, bl, bu)
    return
end

function MOI.set(
    m::Optimizer,
    ::MOI.ConstraintSet,
    cref::MOI.ConstraintIndex{<:MOI.ScalarAffineFunction,D},
    dom::D,
) where {D<:ScalarLinearDomain}
    i = getindex(m.c_block, cref.value) # since we are in a scalar domain
    bk, bl, bu = Mosek.getconbound(m.task, i)
    bl, bu = _change_bound(bl, bu, dom)
    Mosek.putconbound(m.task, i, bk, bl, bu)
    return
end

### MODIFY

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
    col = mosek_index(m, func.variable)::ColumnIndex
    Mosek.putaij(m.task, row(m, c), col.value, func.new_coefficient)
    return
end

### TRANSFORM

function MOI.transform(
    m::Optimizer,
    cref::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},D1},
    newdom::D2,
) where {D1<:ScalarLinearDomain,D2<:ScalarLinearDomain}
    F = MOI.ScalarAffineFunction{Float64}
    r = row(m, cref)
    _putconbound(m, r, newdom)
    return MOI.ConstraintIndex{F,D2}(cref.value)
end

## Delete #####################################################################
###############################################################################

_domain(::Type{MOI.SecondOrderCone}) = Mosek.MSK_DOMAIN_QUADRATIC_CONE

_domain(::Type{MOI.RotatedSecondOrderCone}) = Mosek.MSK_DOMAIN_RQUADRATIC_CONE

_domain(::Type{MOI.ExponentialCone}) = Mosek.MSK_DOMAIN_PRIMAL_EXP_CONE

_domain(::Type{MOI.DualExponentialCone}) = Mosek.MSK_DOMAIN_DUAL_EXP_CONE

_domain(::Type{MOI.PowerCone{Float64}}) = Mosek.MSK_DOMAIN_PRIMAL_POWER_CONE

_domain(::Type{MOI.DualPowerCone{Float64}}) = Mosek.MSK_DOMAIN_DUAL_POWER_CONE

function _domain(::Type{MOI.Scaled{MOI.PositiveSemidefiniteConeTriangle}})
    return Mosek.MSK_DOMAIN_SVEC_PSD_CONE
end

function MOI.is_valid(
    model::Optimizer,
    ci::MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64},S},
) where {S<:VectorConeDomain}
    numacc = Mosek.getnumacc(model.task)
    if !(1 <= ci.value <= numacc)
        return false
    end
    domidx = Mosek.getaccdomain(model.task, ci.value)
    return Mosek.getdomaintype(model.task, domidx) == _domain(S)
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
           Mosek.getconbound(model.task, row(model, ci))[1] == _bound_key(S)
end

function MOI.delete(
    m::Optimizer,
    cref::MOI.ConstraintIndex{F,D},
) where {F<:MOI.ScalarAffineFunction{Float64},D<:ScalarLinearDomain}
    MOI.throw_if_not_valid(m, cref)
    delete!(m.con_to_name, cref)
    m.name_to_con = nothing
    subi = getindexes(m.c_block, cref.value)
    n = length(subi)
    subi_i32 = convert(Vector{Int32}, subi)
    ptr = fill(Int64(0), n)
    Mosek.putarowlist(m.task, subi_i32, ptr, ptr, Int32[], Float64[])
    b = fill(0.0, n)
    Mosek.putconboundlist(m.task, subi_i32, fill(Mosek.MSK_BK_FX, n), b, b)
    deleteblock(m.c_block, cref.value)
    return
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
    vi = _variable(ci)
    unset_flag(m, vi, S)
    _delete_variable_constraint(m, column(m, vi), S)
    return
end

function MOI.is_valid(
    model::Optimizer,
    ci::MOI.ConstraintIndex{MOI.VectorOfVariables,S},
) where {S<:VectorCone}
    if !(ci.value in eachindex(model.variable_to_vector_constraint_id))
        return false
    end
    id = cone_id(model, ci)
    return 1 <= id <= Mosek.getnumcone(model.task) &&
           Mosek.getconeinfo(model.task, id)[1] == _cone_type(S)
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
    map!(
        i -> i > id ? i - 1 : i,
        model.variable_to_vector_constraint_id,
        model.variable_to_vector_constraint_id,
    )
    return
end

function MOI.is_valid(
    model::Optimizer,
    ci::MOI.ConstraintIndex{
        MOI.VectorOfVariables,
        MOI.PositiveSemidefiniteConeTriangle,
    },
)
    # TODO add supports for deletion
    return 1 <= ci.value <= length(model.sd_dim)
end

## List #######################################################################
###############################################################################
function MOI.get(
    m::Optimizer,
    ::MOI.ListOfConstraintAttributesSet{F,S},
) where {F,S}
    set = MOI.AbstractConstraintAttribute[]
    for ci in keys(m.con_to_name)
        if ci isa MOI.ConstraintIndex{F,S}
            push!(set, MOI.ConstraintName())
            break
        end
    end
    return set
end

## Name #######################################################################
###############################################################################

function _putconname(
    m::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}},
    name::AbstractString,
)
    Mosek.putconname(m.task, row(m, c), name)
    return
end

_putconname(::Optimizer, ::MOI.ConstraintIndex, ::AbstractString) = nothing

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
    m.con_to_name[ci] = name
    m.name_to_con = nothing
    _putconname(m, ci, name)
    return
end

function MOI.get(m::Optimizer, ::MOI.ConstraintName, ci::MOI.ConstraintIndex)
    return get(m.con_to_name, ci, "")
end

function _rebuild_name_to_constraint_index(m::Optimizer)
    m.name_to_con = Dict{String,Union{Nothing,MOI.ConstraintIndex}}()
    for (con, name) in m.con_to_name
        if isempty(name)
            continue
        end
        if haskey(m.name_to_con, name)
            m.name_to_con[name] = nothing
        else
            m.name_to_con[name] = con
        end
    end
    return
end

function MOI.get(
    m::Optimizer,
    ::Type{CI},
    name::String,
) where {CI<:MOI.ConstraintIndex}
    if m.name_to_con === nothing
        _rebuild_name_to_constraint_index(m)
    end
    index = get(m.name_to_con, name, missing)
    if ismissing(index)
        return nothing  # No name exists
    elseif isnothing(index)
        error("Multiple constraints have the name: $name")
    elseif !(index isa CI)
        return nothing  # A name exists, but not for this type
    end
    return index
end

function MOI.supports(
    ::Optimizer,
    ::MOI.ConstraintName,
    ::Type{<:MOI.ConstraintIndex{MOI.VariableIndex}},
)
    # Names are not defined for variable constraints
    return throw(MOI.VariableIndexConstraintNameError())
end

function MOI.set(
    ::Optimizer,
    ::MOI.ConstraintName,
    ::MOI.ConstraintIndex{MOI.VariableIndex},
    ::AbstractString,
)
    # Names are not defined for variable constraints
    return throw(MOI.VariableIndexConstraintNameError())
end
