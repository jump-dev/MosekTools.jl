# Copyright (c) 2017: Ulf Worsøe, Mosek ApS
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module MosekTools

import MathOptInterface as MOI
import Mosek

export Mosek

# Allows the user to use `Mosek.Optimizer` instead of `MosekTools.Optimizer`
# for convenience and for consitency with other solvers where the syntax is
# `SolverName.Optimizer`.
Mosek.Optimizer(; kwargs...) = MosekTools.Optimizer(; kwargs...)

include("LinkedInts.jl")

const DEBUG = false

struct MosekSolution
    whichsol::Mosek.Soltype
    solsta::Mosek.Solsta
    prosta::Mosek.Prosta

    xxstatus::Vector{Mosek.Stakey}
    xx::Vector{Float64}
    barxj::Vector{Vector{Float64}}
    slx::Vector{Float64}
    sux::Vector{Float64}
    snx::Vector{Float64}
    doty::Vector{Float64}

    cstatus::Vector{Mosek.Stakey}
    xc::Vector{Float64}
    slc::Vector{Float64}
    suc::Vector{Float64}
    y::Vector{Float64}
end

struct ColumnIndex
    value::Int32
end

struct ColumnIndices
    values::Vector{Int32}
end

struct MatrixIndex
    # `-1` means it has been deleted (hence it was a scalar variable since deleting a matrix variable is not supported)
    # `0` means it is a scalar variable
    # `> 0` means it is a matrix variable part of the `matrix`th block
    matrix::Int32
    # `row` in the lower-triangular part of the `matrix`th block
    row::Int32
    # `column` in the lower-triangular part of the `matrix`th block
    column::Int32
    function MatrixIndex(matrix::Integer, row::Integer, column::Integer)
        # Since it is in the lower-triangular part:
        @assert column ≤ row
        return new(matrix, row, column)
    end
end

"""
    Optimizer <: MOI.AbstractOptimizer

Linear variables and constraint can be deleted. For some reason MOSEK
does not support deleting PSD variables.

Note also that adding variables and constraints will permanently add
some (currently between 1 and 3) Int64s that a `delete!` will not
remove. This ensures that Indices (Variable and constraint) that
are deleted are thereafter invalid.
"""
mutable struct Optimizer <: MOI.AbstractOptimizer
    task::Mosek.MSKtask
    ## Options passed in `Mosek.Optimizer` that are used to create a new task
    ## in `MOI.empty!`:
    # Should Mosek output be ignored or printed ?
    be_quiet::Bool
    # Integer parameters, i.e. parameters starting with `MSK_IPAR_`
    ipars::Dict{String,Int32}
    # Floating point parameters, i.e. parameters starting with `MSK_DPAR_`
    dpars::Dict{String,Float64}
    # String parameters, i.e. parameters starting with `MSK_SPAR_`
    spars::Dict{String,AbstractString}
    has_variable_names::Bool

    # Mappings for MOI.ConstraintName
    con_to_name::Dict{MOI.ConstraintIndex,String}
    # name_to_con is built lazily
    name_to_con::Union{Nothing,Dict{String,Union{MOI.ConstraintIndex,Nothing}}}

    # For each MOI index of variables, gives the flags of constraints present
    # The SingleVariable constraints added cannot just be inferred from Mosek.getvartype
    # and Mosek.getvarbound so we need to keep them here so implement `MOI.is_valid`
    x_constraints::Vector{UInt8}

    F_rows::Dict{Int,UnitRange{Int}} # TODO can it be obtained from Mosek ?

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
    c_block::LinkedInts

    # i -> 0: Not in a VectorOfVariables constraint
    # i -> +j: In `MOI.ConstraintIndex{MOI.VectorOfVariables, ?}(j)`
    # i -> -j: In `MOI.VectorOfVariables` constraint with `MOI.VariableIndex(j)` as first variable
    variable_to_vector_constraint_id::Vector{Int32}

    ###########################
    trm::Union{Nothing,Mosek.Rescode}
    solutions::Vector{MosekSolution}

    ###########################
    """
    Indicating whether the objective sense is MOI.FEASIBILITY_SENSE. It is
    encoded as a MOI.MIN_SENSE with a zero objective internally but this allows
    MOI.get(::Optimizer, ::ObjectiveSense) to still return the right value
    """
    feasibility::Bool
    """
    Indicates whether there is an objective set.
    If `has_objective` is `false` then Mosek has a zero objective internally.
    This affects `MOI.ListOfModelAttributesSet`.
    """
    has_objective::Bool
    """
    Indicates whether there was any PSD variables used when setting the objective.
    When resetting it, I don't know how to remove the contributions from these
    variables so we need to keep this boolean to throw in that case.
    """
    has_psd_in_objective::Bool

    fallback::Union{String,Nothing}

    function Optimizer(; kwargs...)
        optimizer = new(
            Mosek.maketask(), # task
            false, # be_quiet
            Dict{String,Int32}(), # ipars
            Dict{String,Float64}(), # dpars
            Dict{String,AbstractString}(), # spars
            false, # has_variable_names
            Dict{MOI.ConstraintIndex,String}(),  # con_to_name
            nothing,                             # name_to_con
            UInt8[], # x_constraints
            Dict{Int,UnitRange{Int}}(),
            LinkedInts(),# x_block
            MatrixIndex[], # x_sd
            Int[], # sd_dim
            LinkedInts(), # c_block
            Int32[], # variable_to_vector_constraint_id
            nothing,# trm
            MosekSolution[],
            true, # feasibility_sense
            false, # has_objective
            false, # has_psd_in_objective
            nothing,
        )
        Mosek.appendrzerodomain(optimizer.task, 0)
        Mosek.putstreamfunc(optimizer.task, Mosek.MSK_STREAM_LOG, print)
        if length(kwargs) > 0
            @warn("""Passing optimizer attributes as keyword arguments to
            Mosek.Optimizer is deprecated. Use
                MOI.set(model, MOI.RawOptimizerAttribute("key"), value)
            or
                JuMP.set_optimizer_attribute(model, "key", value)
            instead.
            """)
        end
        for (option, value) in kwargs
            MOI.set(optimizer, MOI.RawOptimizerAttribute(string(option)), value)
        end
        return optimizer
    end
end

# MosekTools.IntegerParameter

struct IntegerParameter <: MOI.AbstractOptimizerAttribute
    name::String
end

function MOI.set(m::Optimizer, p::IntegerParameter, value)
    m.ipars[p.name] = value
    Mosek.putnaintparam(m.task, p.name, value)
    return
end

function MOI.get(m::Optimizer, p::IntegerParameter)
    return Mosek.getnaintparam(m.task, p.name)
end

# MosekTools.DoubleParameter

struct DoubleParameter <: MOI.AbstractOptimizerAttribute
    name::String
end

function MOI.set(m::Optimizer, p::DoubleParameter, value)
    m.dpars[p.name] = value
    Mosek.putnadouparam(m.task, p.name, value)
    return
end

function MOI.get(m::Optimizer, p::DoubleParameter)
    return Mosek.getnadouparam(m.task, p.name)
end

# MosekTools.StringParameter

struct StringParameter <: MOI.AbstractOptimizerAttribute
    name::String
end

function MOI.set(m::Optimizer, p::StringParameter, value::AbstractString)
    m.spars[p.name] = value
    Mosek.putnastrparam(m.task, p.name, value)
    return
end

function MOI.get(m::Optimizer, p::StringParameter)
    # We need to give the maximum length of the value of the parameter.
    # 255 should be ok in most cases.
    _, str = Mosek.getnastrparam(m.task, p.name, 255)
    return str
end

# MOI.RawOptimizerAttribute

"""
Set optimizer parameters. Set MOSEK solver parameters, or one of the
additional parametes:

- "QUIET" (true|false), to enable or disable solver log output
- "fallback" (string), to set a solver server to use if no local license file was found,
"""
function MOI.set(m::Optimizer, p::MOI.RawOptimizerAttribute, value)
    if p.name == "QUIET"
        if m.be_quiet != convert(Bool, value)
            m.be_quiet = !m.be_quiet
            if m.be_quiet
                Mosek.putstreamfunc(m.task, Mosek.MSK_STREAM_LOG, m -> nothing)
            else
                Mosek.putstreamfunc(m.task, Mosek.MSK_STREAM_LOG, print)
            end
        end
    elseif p.name == "fallback"
        m.fallback = value
    elseif startswith(p.name, "MSK_IPAR_")
        MOI.set(m, IntegerParameter(p.name), value)
    elseif startswith(p.name, "MSK_DPAR_")
        MOI.set(m, DoubleParameter(p.name), value)
    elseif startswith(p.name, "MSK_SPAR_")
        MOI.set(m, StringParameter(p.name), value)
    elseif value isa Integer
        MOI.set(m, IntegerParameter("MSK_IPAR_$(p.name)"), value)
    elseif value isa AbstractFloat
        MOI.set(m, DoubleParameter("MSK_DPAR_$(p.name)"), value)
    elseif value isa AbstractString
        MOI.set(m, StringParameter("MSK_SPAR_$(p.name)"), value)
    else
        error("Value $value for parameter $(p.name) has unrecognized type")
    end
    return
end

function MOI.get(m::Optimizer, p::MOI.RawOptimizerAttribute)
    if p.name == "QUIET"
        return m.be_quiet
    elseif p.name == "fallback"
        return m.fallback
    elseif startswith(p.name, "MSK_IPAR_")
        return MOI.get(m, IntegerParameter(p.name))
    elseif startswith(p.name, "MSK_DPAR_")
        return MOI.get(m, DoubleParameter(p.name))
    elseif startswith(p.name, "MSK_SPAR_")
        return MOI.get(m, StringParameter(p.name))
    else
        error(
            "The parameter $(p.name) should start by `MSK_IPAR_`, `MSK_DPAR_` or `MSK_SPAR_`.",
        )
    end
end

# MOI.Silent

MOI.supports(::Optimizer, ::MOI.Silent) = true

function MOI.set(model::Optimizer, ::MOI.Silent, value::Bool)
    MOI.set(model, MOI.RawOptimizerAttribute("QUIET"), value)
    return
end

function MOI.get(model::Optimizer, ::MOI.Silent)
    return MOI.get(model, MOI.RawOptimizerAttribute("QUIET"))
end

# MOI.Silent

MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true

function MOI.set(model::Optimizer, ::MOI.TimeLimitSec, value::Real)
    MOI.set(
        model,
        MOI.RawOptimizerAttribute("MSK_DPAR_OPTIMIZER_MAX_TIME"),
        value,
    )
    return
end

function MOI.set(model::Optimizer, ::MOI.TimeLimitSec, ::Nothing)
    MOI.set(
        model,
        MOI.RawOptimizerAttribute("MSK_DPAR_OPTIMIZER_MAX_TIME"),
        -1.0,
    )
    return
end

function MOI.get(model::Optimizer, ::MOI.TimeLimitSec)
    value =
        MOI.get(model, MOI.RawOptimizerAttribute("MSK_DPAR_OPTIMIZER_MAX_TIME"))
    if value < 0.0
        return nothing
    end
    return value
end

function matrix_solution(m::Optimizer, sol)
    return Vector{Float64}[
        Mosek.getbarxj(m.task, sol, j) for j in 1:length(m.sd_dim)
    ]
end

function MOI.optimize!(m::Optimizer)
    # See https://github.com/jump-dev/MosekTools.jl/issues/70
    Mosek.putintparam(
        m.task,
        Mosek.MSK_IPAR_REMOVE_UNUSED_SOLUTIONS,
        Mosek.MSK_ON,
    )
    if m.fallback == nothing
        m.trm = Mosek.optimize(m.task)
    else
        m.trm = Mosek.optimize(m.task, m.fallback)
    end
    m.solutions = MosekSolution[]
    if Mosek.solutiondef(m.task, Mosek.MSK_SOL_ITR)
        push!(
            m.solutions,
            MosekSolution(
                Mosek.MSK_SOL_ITR,
                Mosek.getsolsta(m.task, Mosek.MSK_SOL_ITR),
                Mosek.getprosta(m.task, Mosek.MSK_SOL_ITR),
                Mosek.getskx(m.task, Mosek.MSK_SOL_ITR),
                Mosek.getxx(m.task, Mosek.MSK_SOL_ITR),
                matrix_solution(m, Mosek.MSK_SOL_ITR),
                Mosek.getslx(m.task, Mosek.MSK_SOL_ITR),
                Mosek.getsux(m.task, Mosek.MSK_SOL_ITR),
                Mosek.getsnx(m.task, Mosek.MSK_SOL_ITR),
                Mosek.getaccdotys(m.task, Mosek.MSK_SOL_ITR),
                Mosek.getskc(m.task, Mosek.MSK_SOL_ITR),
                Mosek.getxc(m.task, Mosek.MSK_SOL_ITR),
                Mosek.getslc(m.task, Mosek.MSK_SOL_ITR),
                Mosek.getsuc(m.task, Mosek.MSK_SOL_ITR),
                Mosek.gety(m.task, Mosek.MSK_SOL_ITR),
            ),
        )
    end
    if Mosek.solutiondef(m.task, Mosek.MSK_SOL_ITG)
        push!(
            m.solutions,
            MosekSolution(
                Mosek.MSK_SOL_ITG,
                Mosek.getsolsta(m.task, Mosek.MSK_SOL_ITG),
                Mosek.getprosta(m.task, Mosek.MSK_SOL_ITG),
                Mosek.getskx(m.task, Mosek.MSK_SOL_ITG),
                Mosek.getxx(m.task, Mosek.MSK_SOL_ITG),
                # See https://github.com/jump-dev/MosekTools.jl/issues/71
                Float64[], #matrix_solution(m, Mosek.MSK_SOL_ITG),
                Float64[],
                Float64[],
                Float64[],
                Float64[],
                Mosek.getskc(m.task, Mosek.MSK_SOL_ITG),
                Mosek.getxc(m.task, Mosek.MSK_SOL_ITG),
                Float64[],
                Float64[],
                Float64[],
            ),
        )
    end
    if Mosek.solutiondef(m.task, Mosek.MSK_SOL_BAS)
        push!(
            m.solutions,
            MosekSolution(
                Mosek.MSK_SOL_BAS,
                Mosek.getsolsta(m.task, Mosek.MSK_SOL_BAS),
                Mosek.getprosta(m.task, Mosek.MSK_SOL_BAS),
                Mosek.getskx(m.task, Mosek.MSK_SOL_BAS),
                Mosek.getxx(m.task, Mosek.MSK_SOL_BAS),
                # See https://github.com/jump-dev/MosekTools.jl/issues/71
                Float64[], #matrix_solution(m, Mosek.MSK_SOL_BAS),
                Mosek.getslx(m.task, Mosek.MSK_SOL_BAS),
                Mosek.getsux(m.task, Mosek.MSK_SOL_BAS),
                Float64[],
                Float64[],
                Mosek.getskc(m.task, Mosek.MSK_SOL_BAS),
                Mosek.getxc(m.task, Mosek.MSK_SOL_BAS),
                Mosek.getslc(m.task, Mosek.MSK_SOL_BAS),
                Mosek.getsuc(m.task, Mosek.MSK_SOL_BAS),
                Mosek.gety(m.task, Mosek.MSK_SOL_BAS),
            ),
        )
    end
    # We need to sort the solutions, so that an optimal one is first (if it exists).
    sort!(
        m.solutions;
        by = x -> x.solsta in
        [Mosek.MSK_SOL_STA_OPTIMAL, Mosek.MSK_SOL_STA_INTEGER_OPTIMAL],
        rev = true,
    )
    return
end

# MOI.Name

MOI.supports(::Optimizer, ::MOI.Name) = true

function MOI.set(m::Optimizer, ::MOI.Name, name::String)
    Mosek.puttaskname(m.task, name)
    return
end

function MOI.get(m::Optimizer, ::MOI.Name)
    return Mosek.gettaskname(m.task)
end

# MOI.ListOfModelAttributesSet

function MOI.get(m::Optimizer, ::MOI.ListOfModelAttributesSet)
    set = MOI.AbstractModelAttribute[]
    if !m.feasibility
        push!(set, MOI.ObjectiveSense())
    end
    if m.has_objective
        push!(set, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
    end
    if !isempty(MOI.get(m, MOI.Name()))
        push!(set, MOI.Name())
    end
    return set
end

# MOI.is_empty

function MOI.is_empty(m::Optimizer)
    return Mosek.getnumvar(m.task) == 0 &&
           Mosek.getnumcon(m.task) == 0 &&
           Mosek.getnumcone(m.task) == 0 &&
           Mosek.getnumbarvar(m.task) == 0 &&
           isempty(m.F_rows)
end

# MOI.empty!

function MOI.empty!(model::Optimizer)
    model.task = Mosek.maketask()
    Mosek.appendrzerodomain(model.task, 0)
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
    empty!(model.con_to_name)
    model.name_to_con = nothing
    empty!(model.F_rows)
    empty!(model.x_constraints)
    model.x_block = LinkedInts()
    empty!(model.x_sd)
    empty!(model.sd_dim)
    model.c_block = LinkedInts()
    empty!(model.variable_to_vector_constraint_id)
    model.trm = nothing
    empty!(model.solutions)
    model.feasibility = true
    model.has_objective = false
    model.has_psd_in_objective = false
    return
end

# MOI.SolverName

MOI.get(::Optimizer, ::MOI.SolverName) = "Mosek"

# MOI.SolverVersion

function MOI.get(::Optimizer, ::MOI.SolverVersion)
    major, minor, revision = Mosek.getversion()
    return string(VersionNumber(major, minor, revision))
end

# MOI.supports_incremental_interface

MOI.supports_incremental_interface(::Optimizer) = true

# MOI.copy_to

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike)
    return MOI.Utilities.default_copy_to(dest, src)
end

# MOI.write_to_file

function MOI.write_to_file(m::Optimizer, filename::String)
    Mosek.putintparam(m.task, Mosek.MSK_IPAR_OPF_WRITE_SOLUTIONS, Mosek.MSK_ON)
    Mosek.writedata(m.task, filename)
    return
end

ref2id(vi::MOI.VariableIndex)::Int = vi.value

ref2id(ci::MOI.ConstraintIndex)::Int = ci.value

include("objective.jl")
include("variable.jl")
include("constraint.jl")
include("attributes.jl")

end # module
