module MosekTools

import MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities

using Mosek
export Mosek
# Allows the user to use `Mosek.Optimizer` instead of `MosekTools.Optimizer`
# for convenience and for consitency with other solvers where the syntax is
# `SolverName.Optimizer`.
function Mosek.Optimizer(; kws...)
    return MosekTools.Optimizer(; kws...)
end

include("LinkedInts.jl")

const DEBUG = false

struct MosekSolution
    whichsol :: Soltype
    solsta   :: Solsta
    prosta   :: Prosta

    xxstatus :: Vector{Stakey}
    xx       :: Vector{Float64}
    barxj    :: Vector{Vector{Float64}}
    slx      :: Vector{Float64}
    sux      :: Vector{Float64}
    snx      :: Vector{Float64}

    cstatus  :: Vector{Stakey}
    xc       :: Vector{Float64}
    slc      :: Vector{Float64}
    suc      :: Vector{Float64}
    y        :: Vector{Float64}
end

struct ColumnIndex
    value::Int32
end

struct ColumnIndices
    values::Vector{Int32}
end

struct MatrixIndex
    matrix::Int32
    row::Int32
    column::Int32
    function MatrixIndex(matrix::Integer, row::Integer, column::Integer)
        @assert column ≤ row
        new(matrix, row, column)
    end
end

"""
    Optimizer <: MathOptInterface.AbstractModel

Linear variables and constraint can be deleted. For some reason MOSEK
does not support deleting PSD variables.

Note also that adding variables and constraints will permanently add
some (currently between 1 and 3) Int64s that a `delete!` will not
remove. This ensures that Indices (Variable and constraint) that
are deleted are thereafter invalid.
"""
mutable struct Optimizer  <: MOI.AbstractOptimizer
    task :: Mosek.MSKtask
    ## Options passed in `Mosek.Optimizer` that are used to create a new task
    ## in `MOI.empty!`:
    # Should Mosek output be ignored or printed ?
    be_quiet :: Bool
    # Integer parameters, i.e. parameters starting with `MSK_IPAR_`
    ipars :: Dict{String, Int32}
    # Floating point parameters, i.e. parameters starting with `MSK_DPAR_`
    dpars :: Dict{String, Float64}
    # String parameters, i.e. parameters starting with `MSK_SPAR_`
    spars :: Dict{String, AbstractString}

    has_variable_names::Bool
    constrnames :: Dict{String, Vector{MOI.ConstraintIndex}}
    # Mosek only support names for `MOI.ScalarAffineFunction` so we need a
    # fallback for `SingleVariable` and `VectorOfVariables`.
    con_to_name :: Dict{MOI.ConstraintIndex, String}

    # For each MOI index of variables, gives the flags of constraints present
    # The SingleVariable constraints added cannot just be inferred from getvartype
    # and getvarbound so we need to keep them here so implement `MOI.is_valid`
    x_constraints::Vector{UInt8}

    """
        The total length of `x_block` matches the number of variables in
    the underlying task, and the number of blocks corresponds to the
    number variables allocated in the Model.
    """
    x_block::LinkedInts

    """
    One entry per scalar variable in the task indicating in which semidefinite
    block it is and at which index.
    MOI index -> MatrixIndex
    """
    x_sd::Vector{MatrixIndex}

    sd_dim::Vector{Int}

    ###########################
    """
    One scalar entry per constraint in the underlying task. One block
    per constraint allocated in the Model.
    """
    c_block :: LinkedInts

    # i -> 0: Not in a VectorOfVariables constraint
    # i -> +j: In `MOI.ConstraintIndex{MOI.VectorOfVariables, ?}(j)`
    # i -> -j: In `MOI.VectorOfVariables` constraint with `MOI.VariableIndex(j)` as first variable
    variable_to_vector_constraint_id::Vector{Int32}

    ###########################
    trm :: Union{Nothing, Rescode}
    solutions :: Vector{MosekSolution}

    ###########################
    """
    Indicating whether the objective sense is MOI.FEASIBILITY_SENSE. It is
    encoded as a MOI.MIN_SENSE with a zero objective internally but this allows
    MOI.get(::Optimizer, ::ObjectiveSense) to still return the right value
    """
    feasibility :: Bool

    fallback :: Union{String, Nothing}

    function Optimizer(; kws...)
        optimizer = new(
            maketask(), # task
            false, # be_quiet
            Dict{String, Int32}(), # ipars
            Dict{String, Float64}(), # dpars
            Dict{String, AbstractString}(), # spars
            false, # has_variable_names
            Dict{String, Vector{MOI.ConstraintIndex}}(), # constrnames
            Dict{MOI.ConstraintIndex, String}(), # con_to_name
            UInt8[], # x_constraints
            LinkedInts(),# x_block
            MatrixIndex[], # x_sd
            Int[], # sd_dim
            LinkedInts(), # c_block
            Int32[], # variable_to_vector_constraint_id
            nothing,# trm
            MosekSolution[],
            true, # feasibility_sense
            nothing,
        )
        Mosek.putstreamfunc(optimizer.task, Mosek.MSK_STREAM_LOG, m -> print(m))
        if length(kws) > 0
            @warn("""Passing optimizer attributes as keyword arguments to
            Mosek.Optimizer is deprecated. Use
                MOI.set(model, MOI.RawParameter("key"), value)
            or
                JuMP.set_optimizer_attribute(model, "key", value)
            instead.
            """)
        end
        for (option, value) in kws
            MOI.set(optimizer, MOI.RawParameter(string(option)), value)
        end
        return optimizer
    end
end

struct IntegerParameter <: MOI.AbstractOptimizerAttribute
    name::String
end
function MOI.set(m::Optimizer, p::IntegerParameter, value)
    m.ipars[p.name] = value
    Mosek.putnaintparam(m.task, p.name, value)
end
function MOI.get(m::Optimizer, p::IntegerParameter)
    Mosek.getnaintparam(m.task, p.name)
end

struct DoubleParameter <: MOI.AbstractOptimizerAttribute
    name::String
end
function MOI.set(m::Optimizer, p::DoubleParameter, value)
    m.dpars[p.name] = value
    Mosek.putnadouparam(m.task, p.name, value)
end
function MOI.get(m::Optimizer, p::DoubleParameter)
    Mosek.getnadouparam(m.task, p.name)
end

struct StringParameter <: MOI.AbstractOptimizerAttribute
    name::String
end
function MOI.set(m::Optimizer, p::StringParameter, value::AbstractString)
    m.spars[p.name] = value
    Mosek.putnastrparam(m.task, p.name, value)
end
function MOI.get(m::Optimizer, p::StringParameter)
    # We need to give the maximum length of the value of the parameter.
    # 255 should be ok in most cases.
    len, str = Mosek.getnastrparam(m.task, p.name, 255)
    return str
end

function MOI.set(m::Optimizer, p::MOI.RawParameter, value)
    if p.name == "QUIET"
        if m.be_quiet != convert(Bool, value)
            m.be_quiet = !m.be_quiet
            if m.be_quiet
                Mosek.putstreamfunc(m.task, Mosek.MSK_STREAM_LOG,
                                    m -> begin end)
            else
                Mosek.putstreamfunc(m.task, Mosek.MSK_STREAM_LOG,
                                    m -> print(m))
            end
        end
    elseif p.name == "fallback"
        m.fallback = value
    else
        if startswith(p.name, "MSK_IPAR_")
            par = IntegerParameter(p.name)
        elseif startswith(p.name, "MSK_DPAR_")
            par = DoubleParameter(p.name)
        elseif startswith(p.name, "MSK_SPAR_")
            par = StringParameter(p.name)
        elseif isa(value, Integer)
            par = IntegerParameter("MSK_IPAR_" * p.name)
        elseif isa(value, AbstractFloat)
            par = DoubleParameter("MSK_DPAR_" * p.name)
        elseif isa(value, AbstractString)
            par = StringParameter("MSK_SPAR_" * p.name)
        else
            error("Value $value for parameter $(p.name) has unrecognized type")
        end
        MOI.set(m, par, value)
    end
end

function MOI.get(m::Optimizer, p::MOI.RawParameter)
    if p.name == "QUIET"
        return m.be_quiet
    elseif p.name == "fallback"
        return m.fallback
    else
        if startswith(p.name, "MSK_IPAR_")
            par = IntegerParameter(p.name)
        elseif startswith(p.name, "MSK_DPAR_")
            par = DoubleParameter(p.name)
        elseif startswith(p.name, "MSK_SPAR_")
            par = StringParameter(p.name)
        else
            error("The parameter $(p.name) should start by `MSK_IPAR_`, `MSK_DPAR_` or `MSK_SPAR_`.")
        end
        MOI.get(m, par)
    end
end

MOI.supports(::Optimizer, ::MOI.Silent) = true
function MOI.set(model::Optimizer, ::MOI.Silent, value::Bool)
    MOI.set(model, MOI.RawParameter("QUIET"), value)
end
function MOI.get(model::Optimizer, ::MOI.Silent)
    MOI.get(model, MOI.RawParameter("QUIET"))
end

MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true
function MOI.set(model::Optimizer, ::MOI.TimeLimitSec, value::Real)
    MOI.set(model, MOI.RawParameter("MSK_DPAR_OPTIMIZER_MAX_TIME"), value)
end
function MOI.set(model::Optimizer, ::MOI.TimeLimitSec, ::Nothing)
    MOI.set(model, MOI.RawParameter("MSK_DPAR_OPTIMIZER_MAX_TIME"), -1.0)
end
function MOI.get(model::Optimizer, ::MOI.TimeLimitSec)
    value = MOI.get(model, MOI.RawParameter("MSK_DPAR_OPTIMIZER_MAX_TIME"))
    if value < 0.0
        return nothing
    else
        return value
    end
end

function matrix_solution(m::Optimizer, sol)
    return Vector{Float64}[getbarxj(m.task, sol, j) for j in 1:length(m.sd_dim)]
end

function MOI.optimize!(m::Optimizer)
    m.trm = if m.fallback == nothing; optimize(m.task) else optimize(m.task,m.fallback) end
    m.solutions = MosekSolution[]
    # If the problem is conic but a starting value is set,
    # `MSK_SOL_ITG` simply contains the starting value and
    # `MSK_SOL_ITR` contains the actual solution found.
    # We want `VariablePrimal(1)`, ... to be `MSK_SOL_ITR` in that case
    # so puting it first should solve the issue in that case, see
    # https://github.com/jump-dev/MosekTools.jl/issues/66
    if solutiondef(m.task,MSK_SOL_ITR)
        push!(m.solutions,
              MosekSolution(MSK_SOL_ITR,
                            getsolsta(m.task,MSK_SOL_ITR),
                            getprosta(m.task,MSK_SOL_ITR),
                            getskx(m.task,MSK_SOL_ITR),
                            getxx(m.task,MSK_SOL_ITR),
                            matrix_solution(m, MSK_SOL_ITR),
                            getslx(m.task,MSK_SOL_ITR),
                            getsux(m.task,MSK_SOL_ITR),
                            getsnx(m.task,MSK_SOL_ITR),
                            getskc(m.task,MSK_SOL_ITR),
                            getxc(m.task,MSK_SOL_ITR),
                            getslc(m.task,MSK_SOL_ITR),
                            getsuc(m.task,MSK_SOL_ITR),
                            gety(m.task,MSK_SOL_ITR)))
    end
    if solutiondef(m.task,MSK_SOL_ITG)
        push!(m.solutions,
              MosekSolution(MSK_SOL_ITG,
                            getsolsta(m.task, MSK_SOL_ITG),
                            getprosta(m.task, MSK_SOL_ITG),
                            getskx(m.task, MSK_SOL_ITG),
                            getxx(m.task, MSK_SOL_ITG),
                            matrix_solution(m, MSK_SOL_ITG),
                            Float64[],
                            Float64[],
                            Float64[],
                            getskc(m.task, MSK_SOL_ITG),
                            getxc(m.task, MSK_SOL_ITG),
                            Float64[],
                            Float64[],
                            Float64[]))
    end
    if solutiondef(m.task,MSK_SOL_BAS)
        push!(m.solutions,
              MosekSolution(MSK_SOL_BAS,
                            getsolsta(m.task,MSK_SOL_BAS),
                            getprosta(m.task,MSK_SOL_BAS),
                            getskx(m.task,MSK_SOL_BAS),
                            getxx(m.task,MSK_SOL_BAS),
                            matrix_solution(m, MSK_SOL_BAS),
                            getslx(m.task,MSK_SOL_BAS),
                            getsux(m.task,MSK_SOL_BAS),
                            Float64[],
                            getskc(m.task,MSK_SOL_BAS),
                            getxc(m.task,MSK_SOL_BAS),
                            getslc(m.task,MSK_SOL_BAS),
                            getsuc(m.task,MSK_SOL_BAS),
                            gety(m.task,MSK_SOL_BAS)))
    end
end

MOI.supports(::Optimizer, ::MOI.Name) = true
function MOI.set(m::Optimizer, ::MOI.Name, name::String)
    puttaskname(m.task, name)
end
function MOI.get(m::Optimizer, ::MOI.Name)
    return gettaskname(m.task)
end
function MOI.get(m::Optimizer, ::MOI.ListOfModelAttributesSet)
    set = MOI.AbstractModelAttribute[]
    if !m.feasibility
        push!(set, MOI.ObjectiveSense())
        push!(set, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
    end
    set = [MOI.ObjectiveSense(), MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}()]
    if !isempty(MOI.get(m, MOI.Name()))
        push!(set, MOI.Name())
    end
    return set
end

function MOI.is_empty(m::Optimizer)
    getnumvar(m.task) == 0 && getnumcon(m.task) == 0 && getnumcone(m.task) == 0 && getnumbarvar(m.task) == 0
end

function MOI.empty!(model::Optimizer)
    model.task               = maketask()
    for (name, value) in model.ipars
        Mosek.putnaintparam(model.task, name, value)
    end
    for (name, value) in model.dpars
        Mosek.putnadouparam(model.task, name, value)
    end
    for (name, value) in model.spars
        Mosek.putnastrparam(model.task, name, value)
    end
    if !model.be_quiet
        Mosek.putstreamfunc(model.task, Mosek.MSK_STREAM_LOG, m -> print(m))
    end
    model.has_variable_names = false
    empty!(model.constrnames)
    empty!(model.con_to_name)
    empty!(model.x_constraints)
    model.x_block            = LinkedInts()
    empty!(model.x_sd)
    empty!(model.sd_dim)
    model.c_block            = LinkedInts()
    empty!(model.variable_to_vector_constraint_id)
    model.trm                = nothing
    empty!(model.solutions)
    model.feasibility        = true
end

MOI.get(::Optimizer, ::MOI.SolverName) = "Mosek"

MOIU.supports_default_copy_to(::Optimizer, copy_names::Bool) = true
function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike; kws...)
    return MOIU.automatic_copy_to(dest, src; kws...)
end

function MOI.write_to_file(m::Optimizer, filename :: String)
    putintparam(m.task,MSK_IPAR_OPF_WRITE_SOLUTIONS, MSK_ON)
    writedata(m.task,filename)
end

# For linear objectives we accept:
# EITER affine left-hand side and ranged, unbounded, half-open, fixed (equality), PSD or SOC domains
# OR affine and quadratic left-hand side, and ranged, unbounded, half-open, fixed (equality) domains (quadratic constraints must be unbounded or half-open)
#
# For non-quadratic problems we allow binary and integer variables (but not constraints)
#function supportsconstraints(m::MosekSolver, constraint_types) :: Bool
#    for (fun,dom) in constraint_types
#        if  fun in [MOI.ScalarAffineFunction{Float64},
#                    MOI.SingleVariable,
#                    MOI.VectorOfVariables] &&
#            dom in [MOI.GreaterThan{Float64},
#                    MOI.LessThan{Float64},
#                    MOI.EqualTo{Float64},
#                    MOI.Interval{Float64},
#                    MOI.SecondOrderCone,
#                    MOI.RotatedSecondOrderCone,
#                    MOI.PositiveSemidefiniteConeTriangle,
#                    MOI.PositiveSemidefiniteConeScaled ]
#            # ok
#        elseif dom == MOI.Integer && fun in [MOI.SingleVariable, MOI.VectorOfVariables]
#            # ok
#        else
#            return false
#        end
#    end
#    true
#end

ref2id(vi::MOI.VariableIndex)::Int = vi.value
ref2id(ci::MOI.ConstraintIndex)::Int = ci.value

include("objective.jl")
include("variable.jl")
include("constraint.jl")
include("attributes.jl")

end # module
