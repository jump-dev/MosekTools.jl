###############################################################################
# TASK ########################################################################
###### The `task` field should not be accessed outside this section. ##########

function set_primal_start(task::Mosek.MSKtask, col::ColumnIndex, value::Float64)
    xx = Float64[value]
    for sol in [MSK_SOL_BAS, MSK_SOL_ITG]
        putxxslice(task, sol, col.value, col.value + Int32(1), xx)
    end
end
function set_primal_start(m::Optimizer, vi::MOI.VariableIndex, value::Float64)
    set_primal_start(m.task, mosek_index(m, vi), value)
end
function set_primal_start(task::Mosek.MSKtask, cols::ColumnIndices,
                          values::Vector{Float64})
    for sol in [MSK_SOL_BAS, MSK_SOL_ITG]
        if solutiondef(task, sol)
            xx = getxx(task, sol)
        else
            xx = zeros(Float64, getnumvar(task))
        end
        xx[cols.values] = values
        putxx(task, sol, xx)
    end
end
function set_primal_start(m::Optimizer, vis::Vector{MOI.VariableIndex},
                          values::Vector{Float64})
    if all(vi -> is_scalar(m, vi), vis)
        set_primal_start(m.task, columns(m, vis), values)
    else
        for (vi, value) in zip(vis, values)
            set_primal_start(m.task, mosek_index(m, vi), value)
        end
    end
end

###############################################################################
## INDEXING ###################################################################
###############################################################################

function variable_primal(m::Optimizer, N, col::ColumnIndex)
    return m.solutions[N].xx[col.value]
end
function variable_primal(m::Optimizer, N, mat::MatrixIndex)
    d = m.sd_dim[mat.matrix]
    r = d - mat.column + 1
    #   #entries full Δ       #entries right Δ      #entries above in lower Δ
    k = div((d + 1) * d, 2) - div((r + 1) * r, 2) + (mat.row - mat.column + 1)
    return m.solutions[N].barxj[mat.matrix][k]
end
function variable_primal(m::Optimizer, N, vi::MOI.VariableIndex)
    return variable_primal(m, N, mosek_index(m, vi))
end

###############################################################################
# MOI #########################################################################
###############################################################################

#### objective
function MOI.get(m::Optimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(m, attr)
    return getprimalobj(m.task, m.solutions[attr.result_index].whichsol)
end
function MOI.get(m::Optimizer, attr::MOI.DualObjectiveValue)
    MOI.check_result_index_bounds(m, attr)
    return getdualobj(m.task, m.solutions[attr.result_index].whichsol)
end

MOI.get(m::Optimizer, ::MOI.ObjectiveBound) =
    if solutiondef(m.task,MSK_SOL_ITG)
        getdouinf(m.task,MSK_DINF_MIO_OBJ_BOUND)
    elseif solutiondef(m.task,MSK_SOL_ITR)
        getprimalobj(m.task,MSK_SOL_ITR)
    elseif solutiondef(m.task,MSK_SOL_BAS)
        getprimalobj(m.task,MSK_SOL_BAS)
    else
        return 0.0
    end

MOI.get(m::Optimizer, ::MOI.RelativeGap) =
    if solutiondef(m.task,MSK_SOL_ITG)
        getdouinf(m.task,MSK_DINF_MIO_OBJ_REL_GAP)
    else
        0.0
    end

MOI.get(m::Optimizer, ::MOI.SolveTimeSec) = getdouinf(m.task,MSK_DINF_OPTIMIZER_TIME)


# NOTE: The MOSEK interface currently only supports Min and Max
function MOI.get(model::Optimizer, ::MOI.ObjectiveSense)
    if model.feasibility
        return MOI.FEASIBILITY_SENSE
    else
        sense = getobjsense(model.task)
        if sense == MSK_OBJECTIVE_SENSE_MINIMIZE
            MOI.MIN_SENSE
        else
            MOI.MAX_SENSE
        end
    end
end

function MOI.set(model::Optimizer,
                 attr::MOI.ObjectiveSense,
                 sense::MOI.OptimizationSense)
    if sense == MOI.MIN_SENSE
        model.feasibility = false
        putobjsense(model.task,MSK_OBJECTIVE_SENSE_MINIMIZE)
    elseif sense == MOI.MAX_SENSE
        model.feasibility = false
        putobjsense(model.task,MSK_OBJECTIVE_SENSE_MAXIMIZE)
    else
        @assert sense == MOI.FEASIBILITY_SENSE
        model.feasibility = true
        model.has_objective = false
        putobjsense(model.task,MSK_OBJECTIVE_SENSE_MINIMIZE)
        MOI.set(model,
                MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
                MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Float64}[], 0.0))
    end
end

#### Solver/Solution information

function MOI.get(m::Optimizer,attr::MOI.SimplexIterations)
    miosimiter = getlintinf(m.task,MSK_LIINF_MIO_SIMPLEX_ITER)
    if miosimiter > 0
        Int(miosimiter)
    else
        Int(getintinf(m.task,MSK_IINF_SIM_PRIMAL_ITER) + getintinf(m.task,MSK_IINF_SIM_DUAL_ITER))
    end
end

function MOI.get(m::Optimizer,attr::MOI.BarrierIterations)
    miosimiter = getlintinf(m.task,MSK_LIINF_MIO_INTPNT_ITER)
    if miosimiter > 0
        Int(miosimiter)
    else
        Int(getintinf(m.task,MSK_IINF_INTPNT_ITER))
    end
end

function MOI.get(m::Optimizer,attr::MOI.NodeCount)
        Int(getintinf(m.task,MSK_IINF_MIO_NUM_BRANCH))
end

MOI.get(m::Optimizer,attr::MOI.RawSolver) = m.task

MOI.get(m::Optimizer,attr::MOI.ResultCount) = length(m.solutions)

#### Problem information

function MOI.get(model::Optimizer,
                 ::MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64},
                                           S}) where S <: ScalarLinearDomain
    F = MOI.ScalarAffineFunction{Float64}
    return count(id -> MOI.is_valid(model, MOI.ConstraintIndex{F, S}(id)),
                 allocatedlist(model.c_block))
end
function MOI.get(model::Optimizer,
                 ::MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{Float64},
                                               S}) where S <: ScalarLinearDomain
    F = MOI.ScalarAffineFunction{Float64}
    ids = filter(id -> MOI.is_valid(model, MOI.ConstraintIndex{F, S}(id)),
                 allocatedlist(model.c_block))
    return [MOI.ConstraintIndex{F, S}(id) for id in ids]
end
function MOI.get(model::Optimizer,
                 ::MOI.NumberOfConstraints{MOI.VariableIndex, S}) where S<:Union{ScalarLinearDomain,
                                                                                  MOI.Integer}
    F = MOI.VariableIndex
    return count(id -> MOI.is_valid(model, MOI.ConstraintIndex{F, S}(id)),
                 allocatedlist(model.x_block))
end
function MOI.get(model::Optimizer,
                 ::MOI.ListOfConstraintIndices{MOI.VariableIndex, S}) where S<:Union{ScalarLinearDomain,
                                                                                      MOI.Integer}
    F = MOI.VariableIndex
    ids = filter(id -> MOI.is_valid(model, MOI.ConstraintIndex{F, S}(id)),
                 allocatedlist(model.x_block))
    return [MOI.ConstraintIndex{F, S}(id) for id in ids]
end
function MOI.get(model::Optimizer,
                 ::MOI.NumberOfConstraints{MOI.VectorOfVariables, S}) where S<:VectorCone
    F = MOI.VectorOfVariables
    return count(i -> MOI.is_valid(model, MOI.ConstraintIndex{F, S}(i)),
                 eachindex(model.variable_to_vector_constraint_id))
end
function MOI.get(model::Optimizer,
                 ::MOI.NumberOfConstraints{MOI.VectorAffineFunction{Float64}, S}) where S<:VectorConeDomain
    F = MOI.VectorAffineFunction{Float64}
    numacc = getnumacc(model.task)
    return count(i -> MOI.is_valid(model, MOI.ConstraintIndex{F, S}(i)), 1:numacc)
end
function MOI.get(model::Optimizer,
                 ::MOI.ListOfConstraintIndices{MOI.VectorOfVariables, S}) where S<:VectorCone
    F = MOI.VectorOfVariables
    ids = filter(i -> MOI.is_valid(model, MOI.ConstraintIndex{F, S}(i)),
                 eachindex(model.variable_to_vector_constraint_id))
    return [MOI.ConstraintIndex{F, S}(id) for id in ids]
end
function MOI.get(model::Optimizer,
                 ::MOI.NumberOfConstraints{MOI.VectorOfVariables,
                                           MOI.PositiveSemidefiniteConeTriangle})
    # TODO this only works because deletion of PSD constraints is not supported yet
    return length(model.sd_dim)
end
function MOI.get(model::Optimizer,
                 ::MOI.ListOfConstraintIndices{MOI.VectorOfVariables,
                                               MOI.PositiveSemidefiniteConeTriangle})
    # TODO this only works because deletion of PSD constraints is not supported yet
    return [MOI.ConstraintIndex{MOI.VectorOfVariables,
                                MOI.PositiveSemidefiniteConeTriangle}(id) for id in 1:length(model.sd_dim)]
end

function MOI.get(model::Optimizer,
                 ::MOI.ListOfConstraintTypesPresent)
    list = Tuple{DataType, DataType}[]
    F = MOI.VariableIndex
    for D in [MOI.LessThan{Float64}, MOI.GreaterThan{Float64},
              MOI.EqualTo{Float64}, MOI.Interval{Float64},
              MOI.Integer]
        if !iszero(MOI.get(model, MOI.NumberOfConstraints{F, D}()))
            push!(list, (F, D))
        end
    end
    F = MOI.ScalarAffineFunction{Float64}
    for D in [MOI.LessThan{Float64}, MOI.GreaterThan{Float64},
              MOI.EqualTo{Float64}, MOI.Interval{Float64}]
        if !iszero(MOI.get(model, MOI.NumberOfConstraints{F, D}()))
            push!(list, (F, D))
        end
    end
    F = MOI.VectorOfVariables
    for D in [MOI.SecondOrderCone, MOI.RotatedSecondOrderCone,
              MOI.PowerCone{Float64}, MOI.DualPowerCone{Float64},
              MOI.ExponentialCone, MOI.DualExponentialCone,
              MOI.PositiveSemidefiniteConeTriangle]
        if !iszero(MOI.get(model, MOI.NumberOfConstraints{F, D}()))
            push!(list, (F, D))
        end
    end
    return list
end

#### Warm start values

MOI.supports(::Optimizer, ::MOI.VariablePrimalStart, ::Type{MOI.VariableIndex}) = true

function MOI.set(m::Optimizer, ::MOI.VariablePrimalStart,
                 v::MOI.VariableIndex, ::Nothing)
    set_primal_start(m, v, 0.0)
end
function MOI.set(m::Optimizer, ::MOI.VariablePrimalStart,
                 v::MOI.VariableIndex, val::Float64)
    set_primal_start(m, v, val)
end

function MOI.set(m::Optimizer, ::MOI.VariablePrimalStart,
                 vis::Vector{MOI.VariableIndex}, values::Vector{Float64})
    set_primal_start(m, vis, values)
end

# function MOI.set(m::Optimizer,attr::MOI.ConstraintDualStart, vs::Vector{MOI.ConstraintIndex}, vals::Vector{Float64})
#     subj = columns(vs)

#     for sol in [ MSK_SOL_BAS, MSK_SOL_ITG ]
#         if solutiondef(m.task,sol)
#             xx = getxx(m.task,sol)
#             xx[subj] = vals
#             putxx(m.task,sol,xx)
#         else
#             xx = zeros(Float64,getnumvar(m.task))
#             xx[subj] = vals
#             putxx(m.task,sol,xx)
#         end
#     end
# end

#### Variable solution values

function MOI.get(m::Optimizer, attr::MOI.VariablePrimal, vi::MOI.VariableIndex)
    MOI.check_result_index_bounds(m, attr)
    return variable_primal(m, attr.result_index, vi)
end
function MOI.get!(output::Vector{Float64}, m::Optimizer,
                  attr::MOI.VariablePrimal, vs::Vector{MOI.VariableIndex})
    MOI.check_result_index_bounds(m, attr)
    @assert eachindex(output) == eachindex(vs)
    for i in eachindex(output)
        output[i] = MOI.get(m, attr, vs[i])
    end
end
function MOI.get(m::Optimizer, attr::MOI.VariablePrimal,
                 vs::Vector{MOI.VariableIndex})
    MOI.check_result_index_bounds(m, attr)
    output = Vector{Float64}(undef, length(vs))
    MOI.get!(output, m, attr, vs)
    return output
end

#### Variable basis status

function _basis_status(x::Stakey)
    if x == MSK_SK_BAS
        return BASIC
    elseif x == MSK_SK_SUPBAS
        return SUPER_BASIC
    elseif x == MSK_SK_LOW
        return NONBASIC_AT_LOWER
    elseif x == MSK_SK_UPR
        return NONBASIC_AT_UPPER
        # FIXME which one is `NONBASIC` ?
    elseif x == MSK_SK_UNK
        msg = "The status for the constraint or variable is unknown"
    elseif x == MSK_SK_FIX
        msg = "The constraint or variable is fixed"
    else
        @assert x == MSK_SK_INF
        msg = "The constraint or variable is infeasible in the bounds"
    end
    error("Mosek basic status `$msg` has no equivalent in the `MOI.BasisStatusCode` enum.")
end

function MOI.get(m::Optimizer, attr::MOI.VariableBasisStatus, col::ColumnIndex)
    return _basis_status(m.solutions[attr.result_index].xxstatus[col.value])
end

function MOI.get(::Optimizer, attr::MOI.VariableBasisStatus, mat::MatrixIndex)
    error("$attr not supported for PSD variable $mat")
end

function MOI.get(m::Optimizer, attr::MOI.VariableBasisStatus, vi::MOI.VariableIndex)
    MOI.check_result_index_bounds(m, attr)
    return MOI.get(m, attr, mosek_index(m, vi))
end

#### Constraint solution values

function MOI.get(
    m     ::Optimizer,
    attr  ::MOI.ConstraintPrimal,
    ci  ::MOI.ConstraintIndex{MOI.VariableIndex,D}) where D
    MOI.check_result_index_bounds(m, attr)
    col = column(m, _variable(ci))
    return m.solutions[attr.result_index].xx[col.value]
end

# Semidefinite domain for a variable
function MOI.get!(
    output::Vector{Float64},
    m     ::Optimizer,
    attr  ::MOI.ConstraintPrimal,
    ci  ::MOI.ConstraintIndex{MOI.VectorOfVariables,
                              MOI.PositiveSemidefiniteConeTriangle})
    MOI.check_result_index_bounds(m, attr)
    whichsol = getsolcode(m,attr.result_index)
    output[1:length(output)] = reorder(getbarxj(m.task, whichsol, ci.value),
                                       MOI.PositiveSemidefiniteConeTriangle,
                                       false)
end

# Any other domain for variable vector
function MOI.get!(
    output::Vector{Float64},
    m     ::Optimizer,
    attr  ::MOI.ConstraintPrimal,
    ci  ::MOI.ConstraintIndex{MOI.VectorOfVariables,D}) where D
    MOI.check_result_index_bounds(m, attr)
    cols = columns(m, ci)
    output[1:length(output)] = reorder(m.solutions[attr.result_index].xx[cols.values], D, false)
end

function MOI.get(m     ::Optimizer,
                 attr  ::MOI.ConstraintPrimal,
                 ci  ::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},D}) where D
    MOI.check_result_index_bounds(m, attr)
    cid = ref2id(ci)
    subi = getindex(m.c_block,cid)
    return m.solutions[attr.result_index].xc[subi]
end

function _variable_constraint_dual(sol::MosekSolution, col::ColumnIndex,
                                   ::Type{<:Union{MOI.Interval{Float64}, MOI.EqualTo{Float64}}})
    return sol.slx[col.value] - sol.sux[col.value]
end
function _variable_constraint_dual(sol::MosekSolution, col::ColumnIndex, ::Type{MOI.GreaterThan{Float64}})
    return sol.slx[col.value]
end
function _variable_constraint_dual(sol::MosekSolution, col::ColumnIndex, ::Type{MOI.LessThan{Float64}})
    return -sol.sux[col.value]
end
function MOI.get(m::Optimizer, attr::MOI.ConstraintDual,
                 ci::MOI.ConstraintIndex{MOI.VariableIndex, S}) where S <: ScalarLinearDomain
    MOI.check_result_index_bounds(m, attr)
    col = column(m, _variable(ci))
    dual = _variable_constraint_dual(m.solutions[attr.result_index], col, S)
    if getobjsense(m.task) == MSK_OBJECTIVE_SENSE_MINIMIZE
        return dual
    else
        return -dual
    end
end

getsolcode(m::Optimizer, N) = m.solutions[N].whichsol

# The dual or primal of an SDP variable block is returned in lower triangular
# form but the constraint is in upper triangular form.
function reorder(k::Integer, set::MOI.Scaled{MOI.PositiveSemidefiniteConeTriangle}, moi_to_mosek::Bool)
    # `i` is the row in columnwise upper triangular form
    # the returned value is in columnwise lower triangular form
    if !moi_to_mosek
        # If we reverse things, the columnwise lower triangular becomes a
        # columnwise upper triangular so we can use `trimap`
        k = MOI.dimension(set) - k + 1
    end
    j = div(1 + isqrt(8k - 7), 2)
    i = k - div((j - 1) * j, 2)
    d = MOI.side_dimension(set)
    @assert 0 < j <= d
    @assert 0 < i <= j
    k = MOI.Utilities.trimap(d - j + 1, d - i + 1)
    if moi_to_mosek
        k = MOI.dimension(set) - k + 1
    end
    return k
end
function reorder(x::AbstractVector, ::Type{<:Union{MOI.Scaled{MOI.PositiveSemidefiniteConeTriangle},MOI.PositiveSemidefiniteConeTriangle}}, moi_to_mosek::Bool)
    n = MOI.Utilities.side_dimension_for_vectorized_dimension(length(x))
    @assert length(x) == MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(n))
    y = similar(x)
    k = 0
    for j in 1:n, i in j:n
        k += 1
        o = div((i - 1) * i, 2) + j
        if moi_to_mosek
            y[k] = x[o]
        else
            y[o] = x[k]
        end
    end
    @assert k == length(x)
    return y
end


const ExpCones = Union{
    MOI.ExponentialCone,
    MOI.DualExponentialCone,
}

function reorder(i::Integer, ::ExpCones, ::Bool)
    return (3:-1:1)[i]
end
function reorder(x::AbstractVector, ::Union{ExpCones,Type{<:ExpCones}}, ::Bool)
    return [x[3], x[2], x[1]]
end

const NoReorder = Union{
    MOI.SecondOrderCone,
    MOI.RotatedSecondOrderCone,
    MOI.PowerCone,
    MOI.DualPowerCone,
}

reorder(x, ::Union{NoReorder,Type{<:NoReorder}}, ::Bool) = x


# Semidefinite domain for a variable
function MOI.get!(
    output::Vector{Float64},
    m     ::Optimizer,
    attr  ::MOI.ConstraintDual,
    ci  ::MOI.ConstraintIndex{MOI.VectorOfVariables, MOI.PositiveSemidefiniteConeTriangle})
    MOI.check_result_index_bounds(m, attr)
    whichsol = getsolcode(m,attr.result_index)
    # It is in fact a real constraint and cid is the id of an ordinary constraint
    dual = reorder(getbarsj(m.task, whichsol, ci.value),
                   MOI.PositiveSemidefiniteConeTriangle, false)
    if (getobjsense(m.task) == MSK_OBJECTIVE_SENSE_MINIMIZE)
        output[1:length(output)] .= dual
    else
        output[1:length(output)] .= -dual
    end
end

# Any other domain for variable vector
function MOI.get!(
    output::Vector{Float64},
    m     ::Optimizer,
    attr  ::MOI.ConstraintDual,
    ci  ::MOI.ConstraintIndex{MOI.VectorOfVariables,D}) where D
    MOI.check_result_index_bounds(m, attr)

    xcid = ref2id(ci)
    @assert(xcid > 0)

    cols = columns(m, ci)

    idx = reorder(1:length(output), D, true)

    if (getobjsense(m.task) == MSK_OBJECTIVE_SENSE_MINIMIZE)
        output[idx] = m.solutions[attr.result_index].snx[cols.values]
    else
        output[idx] = -m.solutions[attr.result_index].snx[cols.values]
    end
end

function MOI.get!(
    output::Vector{Float64},
    m     ::Optimizer,
    attr  ::MOI.ConstraintDual,
    ci  ::MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64},D}) where D
    MOI.check_result_index_bounds(m, attr)

    r = rows(m, ci)

    idx = reorder(1:length(output), D, true)

    if (getobjsense(m.task) == MSK_OBJECTIVE_SENSE_MINIMIZE)
        output[idx] = m.solutions[attr.result_index].doty[r]
    else
        output[idx] = -m.solutions[attr.result_index].doty[r]
    end
end

function MOI.get(m     ::Optimizer,
                 attr  ::MOI.ConstraintDual,
                 ci  ::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},D}) where D
    MOI.check_result_index_bounds(m, attr)

    cid = ref2id(ci)
    subi = getindex(m.c_block, cid)

    if getobjsense(m.task) == MSK_OBJECTIVE_SENSE_MINIMIZE
        m.solutions[attr.result_index].y[subi]
    else
        - m.solutions[attr.result_index].y[subi]
    end
end

solsize(m::Optimizer, ::MOI.ConstraintIndex{<:MOI.AbstractScalarFunction}) = 1
function solsize(m::Optimizer, ci::MOI.ConstraintIndex{MOI.VectorOfVariables})
    return getconeinfo(m.task, cone_id(m, ci))[3]
end
function solsize(m::Optimizer, ci::MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64}})
    return length(rows(m, ci))
end
function solsize(m::Optimizer,
                 ci::MOI.ConstraintIndex{MOI.VectorOfVariables,
                                         MOI.PositiveSemidefiniteConeTriangle})
    d = m.sd_dim[ci.value]
    return MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(d))
end

function MOI.get(m::Optimizer,
                 attr::Union{MOI.ConstraintPrimal, MOI.ConstraintDual},
                 ci::MOI.ConstraintIndex{<:MOI.AbstractVectorFunction})
    MOI.check_result_index_bounds(m, attr)
    cid = ref2id(ci)
    output = Vector{Float64}(undef, solsize(m, ci))
    MOI.get!(output, m, attr, ci)
    return output
end
function MOI.get(
    m::Optimizer, # FIXME does Mosek provide this ?
    attr::MOI.ConstraintPrimal,
    ci::MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64}},
)
    return MOI.Utilities.get_fallback(m, attr, ci)
end

#### Status codes
function MOI.get(m::Optimizer, attr::MOI.RawStatusString)
    if     m.trm === nothing
        return "MOI.OPTIMIZE_NOT_CALLED"
    elseif m.trm == MSK_RES_OK
        return join([Mosek.tostr(sol.solsta) for sol in m.solutions], ", ")
    else
        return Mosek.tostr(m.trm)
    end

end
function MOI.get(m::Optimizer, attr::MOI.TerminationStatus)
    if     m.trm === nothing
        MOI.OPTIMIZE_NOT_CALLED
    elseif m.trm == MSK_RES_OK
        # checking `any(sol -> sol.solsta == MSK_SOL_STA_PRIM_INFEAS_CER, m.solutions)`
        # doesn't work for MIP as there is not certificate, i.e. the solutions status is
        # `UNKNOWN`, only the problem status is `INFEAS`.
        if any(sol -> sol.prosta in [ MSK_PRO_STA_PRIM_INFEAS,
                                      MSK_PRO_STA_ILL_POSED], m.solutions)
            MOI.INFEASIBLE
        elseif any(sol -> sol.prosta == MSK_PRO_STA_DUAL_INFEAS, m.solutions)
            MOI.DUAL_INFEASIBLE
        elseif any(sol -> sol.prosta == MSK_PRO_STA_PRIM_INFEAS_OR_UNBOUNDED,m.solutions)
            MOI.INFEASIBLE_OR_UNBOUNDED
        elseif any(sol -> sol.solsta in [ MSK_SOL_STA_OPTIMAL,MSK_SOL_STA_INTEGER_OPTIMAL ], m.solutions)
            MOI.OPTIMAL
        else
            MOI.OTHER_ERROR # ??
        end
    elseif m.trm == MSK_RES_TRM_MAX_ITERATIONS
        MOI.ITERATION_LIMIT
    elseif m.trm == MSK_RES_TRM_MAX_TIME
        MOI.TIME_LIMIT
    elseif m.trm == MSK_RES_TRM_OBJECTIVE_RANGE
        MOI.OBJECTIVE_LIMIT
    elseif m.trm == MSK_RES_TRM_MIO_NUM_RELAXS
        MOI.OTHER_LIMIT
    elseif m.trm == MSK_RES_TRM_MIO_NUM_BRANCHES
        MOI.NODE_LIMIT
    elseif m.trm == MSK_RES_TRM_NUM_MAX_NUM_INT_SOLUTIONS
        MOI.SOLUTION_LIMIT
    elseif m.trm == MSK_RES_TRM_STALL
        MOI.SLOW_PROGRESS
    elseif m.trm == MSK_RES_TRM_USER_CALLBACK
        MOI.INTERRUPTED
    elseif m.trm == MSK_RES_TRM_MAX_NUM_SETBACKS
        MOI.OTHER_LIMIT
    elseif m.trm == MSK_RES_TRM_NUMERICAL_PROBLEM
        MOI.SLOW_PROGRESS
    elseif m.trm == MSK_RES_TRM_INTERNAL
        MOI.OTHER_ERROR
    elseif m.trm == MSK_RES_TRM_INTERNAL_STOP
        MOI.OTHER_ERROR
    else
        MOI.OTHER_ERROR
    end
end

function MOI.get(m::Optimizer, attr::MOI.PrimalStatus)
    if attr.result_index > MOI.get(m, MOI.ResultCount())
        return MOI.NO_SOLUTION
    end
    solsta = m.solutions[attr.result_index].solsta
    if     solsta == MSK_SOL_STA_UNKNOWN
        MOI.UNKNOWN_RESULT_STATUS
    elseif solsta == MSK_SOL_STA_OPTIMAL
        MOI.FEASIBLE_POINT
    elseif solsta == MSK_SOL_STA_PRIM_FEAS
        MOI.FEASIBLE_POINT
    elseif solsta == MSK_SOL_STA_DUAL_FEAS
        MOI.UNKNOWN_RESULT_STATUS
    elseif solsta == MSK_SOL_STA_PRIM_AND_DUAL_FEAS
        MOI.FEASIBLE_POINT
    elseif solsta == MSK_SOL_STA_PRIM_INFEAS_CER
        MOI.NO_SOLUTION
    elseif solsta == MSK_SOL_STA_DUAL_INFEAS_CER
        MOI.INFEASIBILITY_CERTIFICATE
    elseif solsta == MSK_SOL_STA_PRIM_ILLPOSED_CER
        MOI.NO_SOLUTION
    elseif solsta == MSK_SOL_STA_DUAL_ILLPOSED_CER
        MOI.REDUCTION_CERTIFICATE
    elseif solsta == MSK_SOL_STA_INTEGER_OPTIMAL
        MOI.FEASIBLE_POINT
    else
        MOI.UNKNOWN_RESULT_STATUS
    end
end

function MOI.get(m::Optimizer,attr::MOI.DualStatus)
    if attr.result_index > MOI.get(m, MOI.ResultCount())
        return MOI.NO_SOLUTION
    end
    solsta = m.solutions[attr.result_index].solsta
    if     solsta == MSK_SOL_STA_UNKNOWN
        MOI.UNKNOWN_RESULT_STATUS
    elseif solsta == MSK_SOL_STA_OPTIMAL
        MOI.FEASIBLE_POINT
    elseif solsta == MSK_SOL_STA_PRIM_FEAS
        MOI.UNKNOWN_RESULT_STATUS
    elseif solsta == MSK_SOL_STA_DUAL_FEAS
        MOI.FEASIBLE_POINT
    elseif solsta == MSK_SOL_STA_PRIM_AND_DUAL_FEAS
        MOI.FEASIBLE_POINT
    elseif solsta == MSK_SOL_STA_PRIM_INFEAS_CER
        MOI.INFEASIBILITY_CERTIFICATE
    elseif solsta == MSK_SOL_STA_DUAL_INFEAS_CER
        MOI.NO_SOLUTION
    elseif solsta == MSK_SOL_STA_PRIM_ILLPOSED_CER
        MOI.REDUCTION_CERTIFICATE
    elseif solsta == MSK_SOL_STA_DUAL_ILLPOSED_CER
        MOI.NO_SOLUTION
    elseif solsta == MSK_SOL_STA_INTEGER_OPTIMAL
        MOI.NO_SOLUTION
    else
        MOI.UNKNOWN_RESULT_STATUS
    end
end

function MOI.Utilities.substitute_variables(
    ::F,
    x::Mosek.MosekEnum,
) where {F<:Function}
    return x
end
