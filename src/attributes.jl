
#### solver attributes
#MOI.get(m::Union{MosekSolver,MosekModel},::MOI.SupportsDuals) = true
#MOI.get(m::Union{MosekSolver,MosekModel},::MOI.SupportsAddConstraintAfterSolve) = true
#MOI.get(m::Union{MosekSolver,MosekModel},::MOI.SupportsAddVariableAfterSolve) = true
#MOI.get(m::Union{MosekSolver,MosekModel},::MOI.SupportsDeleteConstraint) = true
#MOI.get(m::Union{MosekSolver,MosekModel},::MOI.SupportsDeleteVariable) = true
#MOI.get(m::Union{MosekSolver,MosekModel},::MOI.SupportsConicThroughQuadratic) = false # though actually the solver does

#### objective
MOI.get(m::MosekModel,attr::MOI.ObjectiveValue) = getprimalobj(m.task,m.solutions[attr.resultindex].whichsol)

MOI.get(m::MosekModel,attr::MOI.ObjectiveBound) = getdouinf(m.task,MSK_DINF_MIO_OBJ_BOUND)

MOI.get(m::MosekModel,attr::MOI.RelativeGap) = getdouinf(m.task,MSK_DINF_MIO_OBJ_REL_GAP)

MOI.get(m::MosekModel,attr::MOI.SolveTime) = getdouinf(m.task,MSK_DINF_OPTIMIZER_TIME)


# NOTE: The MOSEK interface currently only supports Min and Max
function MOI.get(model::MosekModel, ::MOI.ObjectiveSense)
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

function MOI.set(model::MosekModel,
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
        putobjsense(model.task,MSK_OBJECTIVE_SENSE_MINIMIZE)
        MOI.set(model,
                MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
                MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Float64}[], 0.0))
    end
end

#### Solver/Solution information

function MOI.get(m::MosekModel,attr::MOI.SimplexIterations)
    miosimiter = getlintinf(m.task,MSK_LIINF_MIO_SIMPLEX_ITER)
    if miosimiter > 0
        Int(miosimiter)
    else
        Int(getintinf(m.task,MSK_IINF_SIM_PRIMAL_ITER) + getintinf(m.task,MSK_IINF_SIM_DUAL_ITER))
    end
end

function MOI.get(m::MosekModel,attr::MOI.BarrierIterations)
    miosimiter = getlintinf(m.task,MSK_LIINF_MIO_INTPNT_ITER)
    if miosimiter > 0
        Int(miosimiter)
    else
        Int(getintinf(m.task,MSK_IINF_INTPNT_ITER))
    end
end

function MOI.get(m::MosekModel,attr::MOI.NodeCount)
        Int(getintinf(m.task,MSK_IINF_MIO_NUM_BRANCH))
end

MOI.get(m::MosekModel,attr::MOI.RawSolver) = m.task

MOI.get(m::MosekModel,attr::MOI.ResultCount) = length(m.solutions)

#### Problem information

function MOI.get(m::MosekModel, ::MOI.NumberOfConstraints{F,D}) where {F,D}
    return length(select(m.constrmap, F, D))
end

function MOI.get(model::MosekModel,
                 ::MOI.ListOfConstraints)
    list = Tuple{DataType, DataType}[]
    for F in [MOI.SingleVariable, MOI.ScalarAffineFunction{Float64}]
        for D in [MOI.LessThan{Float64}, MOI.GreaterThan{Float64},
                  MOI.EqualTo{Float64}, MOI.Interval{Float64}]
            if !isempty(select(model.constrmap, F, D))
                push!(list, (F, D))
            end
        end
    end
    F = MOI.VectorOfVariables
    for D in [MOI.Nonpositives, MOI.Nonnegatives, MOI.Reals, MOI.Zeros,
              MOI.SecondOrderCone, MOI.RotatedSecondOrderCone,
              MOI.PowerCone{Float64}, MOI.DualPowerCone{Float64},
              MOI.ExponentialCone, MOI.DualExponentialCone,
              MOI.PositiveSemidefiniteConeTriangle]
        if !isempty(select(model.constrmap, F, D))
            push!(list, (F, D))
        end
    end
    F = MOI.VectorAffineFunction{Float64}
    for D in [MOI.Nonpositives, MOI.Nonnegatives, MOI.Reals, MOI.Zeros,
              MOI.PositiveSemidefiniteConeTriangle]
        if !isempty(select(model.constrmap, F, D))
            push!(list, (F, D))
        end
    end
    return list
end

function MOI.get(model::MosekModel,
                 ::MOI.ListOfConstraintIndices{F, D}) where {F, D}
    return MOI.ConstraintIndex{F, D}.(keys(select(model.constrmap, F, D)))
end

#### Warm start values

function MOI.set(m::MosekModel, attr::MOI.VariablePrimalStart,
                 v::MOI.VariableIndex, val::Float64)
    col = column(m, v)

    xx = Float64[val]
    for sol in [ MSK_SOL_BAS, MSK_SOL_ITG ]
        putxxslice(m.task, sol, Int(col), Int(col + 1), xx)
    end
end

function MOI.set(m::MosekModel, ::MOI.VariablePrimalStart,
                 vs::Vector{MOI.VariableIndex}, vals::Vector{Float64})
    cols = columns(m, vs)

    for sol in [ MSK_SOL_BAS, MSK_SOL_ITG ]
        if solutiondef(m.task,sol)
            xx = getxx(m.task,sol)
            xx[cols] = vals
            putxx(m.task,sol,xx)
        else
            xx = zeros(Float64, getnumvar(m.task))
            xx[cols] = vals
            putxx(m.task,sol,xx)
        end
    end
end

# function MOI.set(m::MosekModel,attr::MOI.ConstraintDualStart, vs::Vector{MOI.ConstraintIndex}, vals::Vector{Float64})
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

function MOI.get!(output::Vector{Float64}, m::MosekModel,
                  attr::MOI.VariablePrimal, vs::Vector{MOI.VariableIndex})
    output[1:length(output)] = m.solutions[attr.N].xx[columns(m, vs)]
end

function MOI.get(m::MosekModel, attr::MOI.VariablePrimal,
                 vs::Vector{MOI.VariableIndex})
    output = Vector{Float64}(undef, length(vs))
    MOI.get!(output, m, attr, vs)
    return output
end

function MOI.get(m::MosekModel, attr::MOI.VariablePrimal, vi::MOI.VariableIndex)
    return m.solutions[attr.N].xx[column(m, vi)]
end


#### Constraint solution values

function MOI.get(
    m     ::MosekModel,
    attr  ::MOI.ConstraintPrimal,
    ci  ::MOI.ConstraintIndex{MOI.SingleVariable,D}) where D

    conid = ref2id(ci)
    idx  = getindex(m.xc_block, conid)
    subj  = m.xc_idxs[idx]

    m.solutions[attr.N].xx[subj]
end

# Semidefinite domain for a variable
function MOI.get!(
    output::Vector{Float64},
    m     ::MosekModel,
    attr  ::MOI.ConstraintPrimal,
    ci  ::MOI.ConstraintIndex{MOI.VectorOfVariables,MOI.PositiveSemidefiniteConeTriangle})

    whichsol = getsolcode(m,attr.N)
    cid = ref2id(ci)
    @assert(cid < 0)

    # It is in fact a real constraint and cid is the id of an ordinary constraint
    barvaridx = - m.c_block_slack[-cid]
    output[1:length(output)] = reorder(getbarxj(m.task,whichsol,barvaridx), PositiveSemidefiniteCone)
end

# Any other domain for variable vector
function MOI.get!(
    output::Vector{Float64},
    m     ::MosekModel,
    attr  ::MOI.ConstraintPrimal,
    ci  ::MOI.ConstraintIndex{MOI.VectorOfVariables,D}) where D

    xcid = ref2id(ci)
    @assert(xcid > 0)

    mask = m.xc_bounds[xcid]
    idxs = getindexes(m.xc_block, xcid)
    subj = m.xc_idxs[idxs]

    output[1:length(output)] = reorder(m.solutions[attr.N].xx[subj], D)
end

function MOI.get(m     ::MosekModel,
                 attr  ::MOI.ConstraintPrimal,
                 ci  ::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},D}) where D
    cid = ref2id(ci)
    subi = getindex(m.c_block,cid)
    m.solutions[attr.N].xc[subi] + m.c_constant[subi]
end



function MOI.get!(
    output::Vector{Float64},
    m     ::MosekModel,
    attr  ::MOI.ConstraintPrimal,
    ci  ::MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64},D}) where D

    cid = ref2id(ci)
    subi = getindexes(m.c_block, cid)

    if     m.c_block_slack[cid] == 0 # no slack
        output[1:length(output)] = m.solutions[attr.N].xc[subi] + m.c_constant[subi]
    else # psd slack
        @assert m.c_block_slack[cid] < 0
        xid = - m.c_block_slack[cid]
        output[1:length(output)] = reorder(getbarxj(m.task,m.solutions[attr.N].whichsol,Int32(xid)), PositiveSemidefiniteCone)
    end
end













function MOI.get(m::MosekModel, attr::MOI.ConstraintDual,
                 ci::MOI.ConstraintIndex{MOI.SingleVariable})
    xcid = ref2id(ci)
    idx = getindex(m.xc_block, xcid) # variable id

    @assert(blocksize(m.xc_block,xcid) > 0)

    subj  = getindex(m.x_block, m.xc_idxs[idx])
    if (getobjsense(m.task) == MSK_OBJECTIVE_SENSE_MINIMIZE)
        if     m.xc_bounds[xcid] & boundflag_lower != 0 && m.xc_bounds[xcid] & boundflag_upper != 0
            m.solutions[attr.N].slx[subj] - m.solutions[attr.N].sux[subj]
        elseif (m.xc_bounds[xcid] & boundflag_lower) != 0
            m.solutions[attr.N].slx[subj]
        elseif (m.xc_bounds[xcid] & boundflag_upper) != 0
            - m.solutions[attr.N].sux[subj]
        elseif (m.xc_bounds[xcid] & boundflag_cone) != 0
            m.solutions[attr.N].snx[subj]
        else
            error("Dual value available for this constraint")
        end
    else
        if     m.xc_bounds[xcid] & boundflag_lower != 0 && m.xc_bounds[xcid] & boundflag_upper != 0
            m.solutions[attr.N].sux[subj] - m.solutions[attr.N].slx[subj]
        elseif (m.xc_bounds[xcid] & boundflag_lower) != 0
            - m.solutions[attr.N].slx[subj]
        elseif (m.xc_bounds[xcid] & boundflag_upper) != 0
            m.solutions[attr.N].sux[subj]
        elseif (m.xc_bounds[xcid] & boundflag_cone) != 0
            - m.solutions[attr.N].snx[subj]
        else
            error("Dual value available for this constraint")
        end
    end
end

getsolcode(m::MosekModel, N) = m.solutions[N].whichsol

# The dual or primal of an SDP variable block is returned in lower triangular
# form but the constraint is in upper triangular form.
function reorder(x, ::Type{PositiveSemidefiniteCone})
    n = div(isqrt(1 + 8length(x)) - 1, 2)
    @assert length(x) == MOI.dimension(PositiveSemidefiniteCone(n))
    y = similar(x)
    k = 0
    for j in 1:n, i in j:n
        k += 1
        y[div((i - 1) * i, 2) + j] = x[k]

    end
    @assert k == length(x)
    return y
end

function reorder(x, ::Type{<:Union{MOI.ExponentialCone,
                                   MOI.DualExponentialCone}})
    return [x[3], x[2], x[1]]
end
reorder(x, ::Type{<:MOI.AbstractVectorSet}) = x



# Semidefinite domain for a variable
function MOI.get!(
    output::Vector{Float64},
    m     ::MosekModel,
    attr  ::MOI.ConstraintDual,
    ci  ::MOI.ConstraintIndex{MOI.VectorOfVariables,MOI.PositiveSemidefiniteConeTriangle})

    whichsol = getsolcode(m,attr.N)
    cid = ref2id(ci)
    @assert(cid < 0)

    # It is in fact a real constraint and cid is the id of an ordinary constraint
    barvaridx = - m.c_block_slack[-cid]
    dual = reorder(getbarsj(m.task,whichsol,barvaridx), MOI.PositiveSemidefiniteConeTriangle)
    if (getobjsense(m.task) == MSK_OBJECTIVE_SENSE_MINIMIZE)
        output[1:length(output)] = dual
    else
        output[1:length(output)] = -dual
    end
end

# Any other domain for variable vector
function MOI.get!(
    output::Vector{Float64},
    m     ::MosekModel,
    attr  ::MOI.ConstraintDual,
    ci  ::MOI.ConstraintIndex{MOI.VectorOfVariables,D}) where D

    xcid = ref2id(ci)
    @assert(xcid > 0)

    mask = m.xc_bounds[xcid]
    idxs = getindexes(m.xc_block,xcid)
    subj = m.xc_idxs[idxs]

    idx = reorder(1:length(output), D)

    if (getobjsense(m.task) == MSK_OBJECTIVE_SENSE_MINIMIZE)
        if     mask & boundflag_lower != 0 && mask & boundflag_upper != 0
            output[idx] = m.solutions[attr.N].slx[subj] - m.solutions[attr.N].sux[subj]
        elseif (mask & boundflag_lower) != 0
            output[idx] = m.solutions[attr.N].slx[subj]
        elseif (mask & boundflag_upper) != 0
            output[idx] = - m.solutions[attr.N].sux[subj]
        elseif (mask & boundflag_cone) != 0
            output[idx] = m.solutions[attr.N].snx[subj]
        else
            error("Dual value available for this constraint")
        end
    else
        if     mask & boundflag_lower != 0 && mask & boundflag_upper != 0
            output[idx] = m.solutions[attr.N].sux[subj] - m.solutions[attr.N].slx[subj]
        elseif (mask & boundflag_lower) != 0
            output[idx] = - m.solutions[attr.N].slx[subj]
        elseif (mask & boundflag_upper) != 0
            output[idx] = m.solutions[attr.N].sux[subj]
        elseif (mask & boundflag_cone) != 0
            output[idx] = - m.solutions[attr.N].snx[subj]
        else
            error("Dual value available for this constraint")
        end
    end
end

function MOI.get(m     ::MosekModel,
                 attr  ::MOI.ConstraintDual,
                 ci  ::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},D}) where D

    cid = ref2id(ci)
    subi = getindex(m.c_block, cid)

    if getobjsense(m.task) == MSK_OBJECTIVE_SENSE_MINIMIZE
        m.solutions[attr.N].y[subi]
    else
        - m.solutions[attr.N].y[subi]
    end
end


function MOI.get!(
    output::Vector{Float64},
    m     ::MosekModel,
    attr  ::MOI.ConstraintDual,
    ci  ::MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64},D}) where { D <: MOI.AbstractSet }

    cid = ref2id(ci)
    subi = getindexes(m.c_block, cid)

    if getobjsense(m.task) == MSK_OBJECTIVE_SENSE_MINIMIZE
        if     m.c_block_slack[cid] == 0 # no slack
            output[1:length(output)] = m.solutions[attr.N].y[subi]
        else # psd slack
            @assert m.c_block_slack[cid] < 0
            whichsol = getsolcode(m,attr.N)
            xid = - m.c_block_slack[cid]
            output[1:length(output)] = reorder(getbarsj(m.task,whichsol,Int32(xid)), PositiveSemidefiniteCone)
        end
    else
        if     m.c_block_slack[cid] == 0 # no slack
            output[1:length(output)] = - m.solutions[attr.N].y[subi]
        else # psd slack
            @assert m.c_block_slack[cid] < 0
            whichsol = getsolcode(m,attr.N)
            xid = - m.c_block_slack[cid]
            output[1:length(output)] = reorder(-getbarsj(m.task,whichsol,Int32(xid)), PositiveSemidefiniteCone)
        end
    end
end














solsize(m::MosekModel, ::MOI.ConstraintIndex{<:MOI.AbstractScalarFunction}) = 1
function solsize(m::MosekModel, ci::MOI.ConstraintIndex{MOI.VectorOfVariables})
    cid = ref2id(ci)
    if cid < 0
        blocksize(m.c_block,-cid)
    else
        blocksize(m.xc_block,cid)
    end
end

function solsize(m::MosekModel,
                 ci::MOI.ConstraintIndex{<:MOI.VectorAffineFunction})
    return blocksize(m.c_block, ref2id(ci))
end

function MOI.get(m::MosekModel,
                 attr::Union{MOI.ConstraintPrimal, MOI.ConstraintDual},
                 ci::MOI.ConstraintIndex{<:MOI.AbstractVectorFunction})
    cid = ref2id(ci)
    output = Vector{Float64}(undef, solsize(m, ci))
    MOI.get!(output, m, attr, ci)
    return output
end






#### Status codes
function MOI.get(m::MosekModel,attr::MOI.TerminationStatus)
    if     m.trm === nothing
        MOI.OPTIMIZE_NOT_CALLED
    elseif m.trm == MSK_RES_OK
        if any(sol -> sol.solsta == MSK_SOL_STA_PRIM_INFEAS_CER, m.solutions)
            MOI.INFEASIBLE
        elseif any(sol -> sol.solsta == MSK_SOL_STA_DUAL_INFEAS_CER,
                   m.solutions)
            MOI.DUAL_INFEASIBLE
        else
            MOI.OPTIMAL
        end
    elseif m.trm == MSK_RES_TRM_MAX_ITERATIONS
        MOI.ITERATION_LIMIT
    elseif m.trm == MSK_RES_TRM_MAX_TIME
        MOI.TIME_LIMIT
    elseif m.trm == MSK_RES_TRM_OBJECTIVE_RANGE
        MOI.OBJECTIVE_LIMIT
    elseif m.trm == MSK_RES_TRM_MIO_NEAR_REL_GAP
        MOI.ALMOST_OPTIMAL
    elseif m.trm == MSK_RES_TRM_MIO_NEAR_ABS_GAP
        MOI.ALMOST_OPTIMAL
    elseif m.trm == MSK_RES_TRM_MIO_NUM_RELAXS
        MOI.OTHER_LIMIT
    elseif m.trm == MSK_RES_TRM_MIO_NUM_BRANCHES
        MOI.NODE_LIMIT
    elseif m.trm == MSK_RES_TRM_NUM_MAX_NUM_INT_SOLUTIONS
        MOI.SOLUTION_LIMIT
    elseif m.trm == MSK_RES_TRM_STALL
        println("STALL")
        MOI.SLOW_PROGRESS
    elseif m.trm == MSK_RES_TRM_USER_CALLBACK
        MOI.INTERRUPTED
    elseif m.trm == MSK_RES_TRM_MAX_NUM_SETBACKS
        MOI.OTHER_LIMIT
    elseif m.trm == MSK_RES_TRM_NUMERICAL_PROBLEM
        println("NUMERICAL_PROBLEM")
        MOI.SLOW_PROGRESS
    elseif m.trm == MSK_RES_TRM_INTERNAL
        MOI.OTHER_ERROR
    elseif m.trm == MSK_RES_TRM_INTERNAL_STOP
        MOI.OTHER_ERROR
    else
        MOI.OTHER_ERROR
    end
end

function MOI.get(m::MosekModel, attr::MOI.PrimalStatus)
    solsta = m.solutions[attr.N].solsta
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

function MOI.get(m::MosekModel,attr::MOI.DualStatus)
    solsta = m.solutions[attr.N].solsta
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
