module MathOptInterfaceMosek

import MathOptInterface
const MOI = MathOptInterface
using Mosek
#using Mosek.Ext

export MosekOptimizer


include("LinkedInts.jl")

const DEBUG = false

mosek_block_type_unallocated = 0
mosek_block_type_zero   = 1
mosek_block_type_nonneg = 2
mosek_block_type_nonpos = 3
mosek_block_type_range  = 4
mosek_block_type_qcone  = 5
mosek_block_type_rqcone = 6
mosek_block_type_psd    = 7
mosek_block_type_integer = 8

problemtype_linear    = 0
problemtype_conic     = 1
problemtype_quadratic = 2

import MathOptInterface



boundflag_lower = 0x1
boundflag_upper = 0x2
boundflag_cone  = 0x4
boundflag_int   = 0x8
boundflag_all   = 0x0f


# Mapping of all constraint types to its index
struct ConstraintMap
    x_lessthan           :: Dict{Int64,Int} # (SingleVariable,LessThan) -> constraint number
    x_greaterthan        :: Dict{Int64,Int} # (SingleVariable,GreaterThan) -> constraint number
    x_equalto            :: Dict{Int64,Int} # (SingleVariable,EqualTo) -> constraint number
    x_interval           :: Dict{Int64,Int} # (SingleVariable,Interval) -> constraint number
    x_nonpositives       :: Dict{Int64,Int} # (SingleVariable,Nonpositives) -> constraint number
    x_nonnegatives       :: Dict{Int64,Int} # (SingleVariable,Nonnegatives) -> constraint number
    x_binary             :: Dict{Int64,Int} # (SingleVariable,ZeroOne) -> constraint number
    x_integer            :: Dict{Int64,Int} # (SingleVariable,Integer) -> constraint number

    xs_nonpositives      :: Dict{Int64,Int} # (VectorOfVariables,Nonpositives) -> constraint number
    xs_nonnegatives      :: Dict{Int64,Int} # (VectorOfVariables,Nonnegatives) -> constraint number
    xs_zeros             :: Dict{Int64,Int} # (VectorOfVariables,Zeros) -> constraint number
    xs_reals             :: Dict{Int64,Int} # (VectorOfVariables,Reals) -> constraint number
    xs_qcone             :: Dict{Int64,Int} # (VectorOfVariables,SecondOrderCone) -> constraint number
    xs_rqcone            :: Dict{Int64,Int} # (VectorOfVariables,RotatedSecondOrderCone) -> constraint number
    xs_psdconetriangle   :: Dict{Int64,Int} # (VectorOfVariables,PositiveSemidefiniteConeTriangle) -> constraint number

    axb_lessthan         :: Dict{Int64,Int} # (ScalarAffineFunction,LessThan) -> constraint number
    axb_greaterthan      :: Dict{Int64,Int} # (ScalarAffineFunction,GreaterThan) -> constraint number
    axb_equalto          :: Dict{Int64,Int} # (ScalarAffineFunction,EqualTo) -> constraint number
    axb_interval         :: Dict{Int64,Int} # (ScalarAffineFunction,Interval) -> constraint number
    axb_nonpositives     :: Dict{Int64,Int} # (ScalarAffineFunction,Nonpositives) -> constraint number
    axb_nonnegatives     :: Dict{Int64,Int} # (ScalarAffineFunction,Nonnegatives) -> constraint number
    axb_binary           :: Dict{Int64,Int} # (ScalarAffineFunction,ZeroOne) -> constraint number
    axb_integer          :: Dict{Int64,Int} # (ScalarAffineFunction,Integer) -> constraint number

    axbs_nonpositives    :: Dict{Int64,Int} # (VectorAffineFunction,Nonpositives) -> constraint number
    axbs_nonnegatives    :: Dict{Int64,Int} # (VectorAffineFunction,Nonnegatives) -> constraint number
    axbs_zeros           :: Dict{Int64,Int} # (VectorAffineFunction,Zeros) -> constraint number
    axbs_reals           :: Dict{Int64,Int} # (VectorAffineFunction,Reals) -> constraint number
    axbs_qcone           :: Dict{Int64,Int} # (VectorAffineFunction,SecondOrderCone) -> constraint number
    axbs_rqcone          :: Dict{Int64,Int} # (VectorAffineFunction,RotatedSecondOrderCone) -> constraint number
    axbs_psdconetriangle :: Dict{Int64,Int} # (VectorAffineFunction,PositiveSemidefiniteConeTriangle) -> constraint number
end

ConstraintMap() = ConstraintMap([Dict{Int64,Int}() for i in 1:30]...)
select(cm::ConstraintMap,::Type{MOI.SingleVariable},               ::Type{MOI.LessThan{Float64}}) =                cm.x_lessthan
select(cm::ConstraintMap,::Type{MOI.SingleVariable},               ::Type{MOI.GreaterThan{Float64}}) =             cm.x_greaterthan
select(cm::ConstraintMap,::Type{MOI.SingleVariable},               ::Type{MOI.EqualTo{Float64}}) =                 cm.x_equalto
select(cm::ConstraintMap,::Type{MOI.SingleVariable},               ::Type{MOI.Interval{Float64}}) =                cm.x_interval
select(cm::ConstraintMap,::Type{MOI.SingleVariable},               ::Type{MOI.Nonpositives}) =                     cm.x_nonpositives
select(cm::ConstraintMap,::Type{MOI.SingleVariable},               ::Type{MOI.Nonnegatives}) =                     cm.x_nonnegatives
select(cm::ConstraintMap,::Type{MOI.SingleVariable},               ::Type{MOI.ZeroOne}) =                          cm.x_binary
select(cm::ConstraintMap,::Type{MOI.SingleVariable},               ::Type{MOI.Integer}) =                          cm.x_integer
select(cm::ConstraintMap,::Type{MOI.VectorOfVariables},            ::Type{MOI.Nonpositives}) =                     cm.xs_nonpositives
select(cm::ConstraintMap,::Type{MOI.VectorOfVariables},            ::Type{MOI.Nonnegatives}) =                     cm.xs_nonnegatives
select(cm::ConstraintMap,::Type{MOI.VectorOfVariables},            ::Type{MOI.Zeros}) =                            cm.xs_zeros
select(cm::ConstraintMap,::Type{MOI.VectorOfVariables},            ::Type{MOI.Reals}) =                            cm.xs_reals
select(cm::ConstraintMap,::Type{MOI.VectorOfVariables},            ::Type{MOI.SecondOrderCone}) =                  cm.xs_qcone
select(cm::ConstraintMap,::Type{MOI.VectorOfVariables},            ::Type{MOI.RotatedSecondOrderCone}) =           cm.xs_rqcone
select(cm::ConstraintMap,::Type{MOI.VectorOfVariables},            ::Type{MOI.PositiveSemidefiniteConeTriangle}) = cm.xs_psdconetriangle
select(cm::ConstraintMap,::Type{MOI.ScalarAffineFunction{Float64}},::Type{MOI.LessThan{Float64}}) =                cm.axb_lessthan
select(cm::ConstraintMap,::Type{MOI.ScalarAffineFunction{Float64}},::Type{MOI.GreaterThan{Float64}}) =             cm.axb_greaterthan
select(cm::ConstraintMap,::Type{MOI.ScalarAffineFunction{Float64}},::Type{MOI.EqualTo{Float64}}) =                 cm.axb_equalto
select(cm::ConstraintMap,::Type{MOI.ScalarAffineFunction{Float64}},::Type{MOI.Interval{Float64}}) =                cm.axb_interval
select(cm::ConstraintMap,::Type{MOI.ScalarAffineFunction{Float64}},::Type{MOI.ZeroOne}) =                          cm.axb_binary
select(cm::ConstraintMap,::Type{MOI.ScalarAffineFunction{Float64}},::Type{MOI.Integer}) =                          cm.axb_integer
select(cm::ConstraintMap,::Type{MOI.ScalarAffineFunction{Float64}},::Type{MOI.Nonpositives}) =                     cm.axb_nonpositives
select(cm::ConstraintMap,::Type{MOI.ScalarAffineFunction{Float64}},::Type{MOI.Nonnegatives}) =                     cm.axb_nonnegatives
select(cm::ConstraintMap,::Type{MOI.VectorAffineFunction{Float64}},::Type{MOI.Nonpositives}) =                     cm.axbs_nonpositives
select(cm::ConstraintMap,::Type{MOI.VectorAffineFunction{Float64}},::Type{MOI.Nonnegatives}) =                     cm.axbs_nonnegatives
select(cm::ConstraintMap,::Type{MOI.VectorAffineFunction{Float64}},::Type{MOI.Zeros}) =                            cm.axbs_zeros
select(cm::ConstraintMap,::Type{MOI.VectorAffineFunction{Float64}},::Type{MOI.Reals}) =                            cm.axbs_reals
select(cm::ConstraintMap,::Type{MOI.VectorAffineFunction{Float64}},::Type{MOI.SecondOrderCone}) =                  cm.axbs_qcone
select(cm::ConstraintMap,::Type{MOI.VectorAffineFunction{Float64}},::Type{MOI.RotatedSecondOrderCone}) =           cm.axbs_rqcone
select(cm::ConstraintMap,::Type{MOI.VectorAffineFunction{Float64}},::Type{MOI.PositiveSemidefiniteConeTriangle}) = cm.axbs_psdconetriangle

Base.getindex{F,D}(cm::ConstraintMap,r :: MOI.ConstraintIndex{F,D}) = select(cm,F,D)[r.value]


struct MosekSolution
    whichsol :: Soltype
    solsta   :: Solsta
    prosta   :: Prosta

    xxstatus :: Vector{Stakey}
    xx       :: Vector{Float64}
    slx      :: Vector{Float64}
    sux      :: Vector{Float64}
    snx      :: Vector{Float64}

    cstatus  :: Vector{Stakey}
    xc       :: Vector{Float64}
    slc      :: Vector{Float64}
    suc      :: Vector{Float64}
    y        :: Vector{Float64}
end

"""
    MosekModel <: MathOptInterface.AbstractModel

Linear variables and constraint can be deleted. For some reason MOSEK
does not support deleting PSD variables.

Note also that adding variables and constraints will permanently add
some (currently between 1 and 3) Int64s that a `delete!` will not
remove. This ensures that Indices (Variable and constraint) that
are deleted are thereafter invalid.
"""
mutable struct MosekModel  <: MOI.AbstractOptimizer
    task :: Mosek.MSKtask

    """
    Number of variables explicitly created by user
    """
    publicnumvar :: Int

    """
    """
    constrmap :: ConstraintMap

    """
        The total length of `x_block` matches the number of variables in
    the underlying task, and the number of blocks corresponds to the
    number variables allocated in the Model.
    """
    x_block      :: LinkedInts

    """
    One entry per scalar variable in the task indicating which bound
    types are defined
    """
    x_boundflags :: Vector{Int}

    """
    One entry per scalar variable in the task defining the number of
    variable constraints imposed on that variable. It is not allowed
    to delete a variable without deleting the constraints on it first
    (to get around the problem of deleting roots in conic constraints).
    """
    x_numxc :: Vector{Int}

    """
    One entry per variable-constraint
    """
    xc_block     :: LinkedInts

    """
    One entry per variable-constraint block indicating which bound
    types it defines. The values are binary ORed `boundflag_...` values.
    """
    xc_bounds    :: Vector{UInt8}
    """
    One entry per variable-constraint block indicating which cone it belongs to, 0 if none.
    """
    xc_coneid    :: Vector{Int}

    """
    One entry per scalar variable-constraint, defining which native
    variables the bound block covers.
    """
    xc_idxs      :: Vector{Int}

    ###########################
    """
    One scalar entry per constraint in the underlying task. One block
    per constraint allocated in the Model.
    """
    c_block :: LinkedInts

    """
    One entry per allocated scalar constraint. Defines the fixed term
    on the left-hand for each scalar constraint.
    """
    c_constant :: Vector{Float64}

    """
    One entry per allocated scalar constraint.
    Each element is either
    - 0: meaning: no slack, when domain is defined directly as a bound,
    - positive: a `x_block` reference, e.g. for qcones, or
    - negative: Negated index of a PSD variable in the underlying task
    """
    c_block_slack   :: Vector{Int}
    """
    One entry per variable-constraint block indicating which cone it belongs to, 0 if none.
    """
    c_coneid   :: Vector{Int}

    ###########################
    conecounter :: Int

    ###########################
    trm :: Rescode
    solutions :: Vector{MosekSolution}
end

function MosekOptimizer(; kws...)
    t = maketask()
    try
        be_quiet = false
        for (option,val) in kws
            parname = string(option)
            if parname == "QUIET"
                be_quiet = be_quiet || convert(Bool,val)
            elseif startswith(parname, "MSK_IPAR_")
                Mosek.putnaintparam(t, parname, convert(Integer, val))
            elseif startswith(parname, "MSK_DPAR_")
                Mosek.putnadouparam(t, parname, convert(AbstractFloat, val))
            elseif startswith(parname, "MSK_SPAR_")
                Mosek.putnastrparam(t, parname, convert(AbstractString, val))
            elseif isa(val, Integer)
                parname = "MSK_IPAR_$parname"
                Mosek.putnaintparam(t, parname, val)
            elseif isa(val, AbstractFloat)
                parname = "MSK_DPAR_$parname"
                Mosek.putnadouparam(t, parname, val)
            elseif isa(val, AbstractString)
                parname = "MSK_SPAR_$parname"
                Mosek.putnastrparam(t, parname, val)
            else
                error("Value $val for parameter $option has unrecognized type")
            end
        end
        if ! be_quiet
            Mosek.putstreamfunc(t,Mosek.MSK_STREAM_LOG,m -> print(m))
        end
        MosekModel(t,# task
                   0, # public numvar
                   ConstraintMap(), # public constraints
                   LinkedInts(),# c_block
                   Int[], # x_boundflags
                   Int[], # x_numxc
                   LinkedInts(), # xc_block
                   UInt8[], # xc_bounds
                   Int[], # xc_coneid
                   Int[], # xc_idxs
                   LinkedInts(), # c_block
                   Float64[], # c_constant
                   Int[], # c_block_slack
                   Int[], # c_coneid
                   0, # cone counter
                   Mosek.MSK_RES_OK,
                   MosekSolution[]) # trm
    catch
        Mosek.deletetask(t)
        rethrow()
    end
end

function MOI.free!(m::MosekModel)
    Mosek.deletetask(m.task)
end

function MOI.isempty(m::MosekModel)
    getnumvar(m.task) == 0 && getnumcon(m.task) == 0 && getnumcone(m.task) == 0 && getnumbarvar(m.task) == 0
end

function MOI.optimize!(m::MosekModel)
    m.trm = optimize(m.task)
    m.solutions = MosekSolution[]
    if solutiondef(m.task,MSK_SOL_ITG)

        push!(m.solutions,
              MosekSolution(MSK_SOL_ITG,
                            getsolsta(m.task,MSK_SOL_ITG),
                            getprosta(m.task,MSK_SOL_ITG),
                            getskx(m.task,MSK_SOL_ITG),
                            getxx(m.task,MSK_SOL_ITG),
                            Float64[],
                            Float64[],
                            Float64[],
                            getskc(m.task,MSK_SOL_ITG),
                            getxc(m.task,MSK_SOL_ITG),
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
                            getslx(m.task,MSK_SOL_BAS),
                            getsux(m.task,MSK_SOL_BAS),
                            Float64[],
                            getskc(m.task,MSK_SOL_BAS),
                            getxc(m.task,MSK_SOL_BAS),
                            getslc(m.task,MSK_SOL_BAS),
                            getsuc(m.task,MSK_SOL_BAS),
                            gety(m.task,MSK_SOL_BAS)))
    end
    if solutiondef(m.task,MSK_SOL_ITR)
        push!(m.solutions,
              MosekSolution(MSK_SOL_ITR,
                            getsolsta(m.task,MSK_SOL_ITR),
                            getprosta(m.task,MSK_SOL_ITR),
                            getskx(m.task,MSK_SOL_ITR),
                            getxx(m.task,MSK_SOL_ITR),
                            getslx(m.task,MSK_SOL_ITR),
                            getsux(m.task,MSK_SOL_ITR),
                            getsnx(m.task,MSK_SOL_ITR),
                            getskc(m.task,MSK_SOL_ITR),
                            getxc(m.task,MSK_SOL_ITR),
                            getslc(m.task,MSK_SOL_ITR),
                            getsuc(m.task,MSK_SOL_ITR),
                            gety(m.task,MSK_SOL_ITR)))
    end
end


function MOI.empty!(m::MosekModel)
    m.task          = maketask()
    m.publicnumvar  = 0
    m.constrmap     = ConstraintMap()
    m.x_block       = LinkedInts()
    m.x_boundflags  = Int[]
    m.x_numxc       = Int[]
    m.xc_block      = LinkedInts()
    m.xc_bounds     = UInt8[]
    m.xc_coneid     = Int[]
    m.xc_idxs       = Int[]
    m.c_block       = LinkedInts()
    m.c_constant    = Float64[]
    m.c_block_slack = Int[]
    m.c_coneid      = Int[]
    m.conecounter   = 0
    m.trm           = Mosek.MSK_RES_OK
    m.solutions     = MosekSolution[]
end

function MOI.Utilities.copyconstraints!(dest::MosekModel, src::MOI.ModelLike, idxmap::MOI.Utilities.IndexMap, ::Type{F}, ::Type{S}) where {F<:MOI.AbstractFunction, S<:MOI.AbstractSet}
    # Copy constraints
    cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
    for ci_src in cis_src
        f_src  = MOI.get(src, MOI.ConstraintFunction(), ci_src)
        f_dest = MOI.Utilities.mapvariables(idxmap, f_src) # appears to be expensive
        s = MOI.get(src, MOI.ConstraintSet(), ci_src)
        ci_dest = MOI.addconstraint!(dest, f_dest, s)
        idxmap.conmap[ci_src] = ci_dest
    end

    nothing
    #return passattributes!(dest, src, idxmap, cis_src)
end


function MOI.copy!(dest::MosekModel, src::MOI.ModelLike; copynames=true)
    if ! MOI.isempty(dest)
        MOI.empty!(dest)
    end

    idxmap = MOI.Utilities.IndexMap()

    vis_src = MOI.get(src, MOI.ListOfVariableIndices())
    for vi in vis_src
        idxmap.varmap[vi] = MOI.addvariable!(dest)
    end

    MOI.set!(dest, MOI.ObjectiveSense(), MOI.get(src,MOI.ObjectiveSense()))
    if MOI.canget(src,MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
        MOI.set!(dest,
                 MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
                 MOI.get(src,MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}()))
    elseif MOI.canget(src,MOI.ObjectiveFunction{MOI.SingleVariable}())
        MOI.set!(dest,
                 MOI.ObjectiveFunction{MOI.SingleVariable}(),
                 MOI.get(src,MOI.ObjectiveFunction{MOI.SingleVariable}()))
    else
        # no objective functrion
    end

    # Copy constraints
    for (F, S) in MOI.get(src, MOI.ListOfConstraints()) # expensive?!
        # do the rest in copyconstraints! which is type stable
        res = MOI.Utilities.copyconstraints!(dest, src, idxmap, F, S)
        #res.status == MOI.CopySuccess || return res
    end

    return MOI.CopyResult(MOI.CopySuccess, "", idxmap)
end



function MOI.write(m::MosekModel, filename :: String)
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
#                    MOI.VectorAffineFunction{Float64},
#                    MOI.VectorOfVariables] &&
#            dom in [MOI.Zeros,
#                    MOI.Reals,
#                    MOI.Nonnegatives,
#                    MOI.Nonpositives,
#                    MOI.GreaterThan{Float64},
#                    MOI.LessThan{Float64},
#                    MOI.EqualTo{Float64},
#                    MOI.Interval{Float64},
#                    MOI.SecondOrderCone,
#                    MOI.RotatedSecondOrderCone,
#                    MOI.PositiveSemidefiniteConeTriangle,
#                    MOI.PositiveSemidefiniteConeScaled ]
#            # ok
#        elseif dom in [MOI.ZeroOne,
#                       MOI.Integer] &&
#                           fun in [MOI.SingleVariable,
#                                   MOI.VectorOfVariables]
#            # ok
#        else
#            return false
#        end
#    end
#    true
#end


#MOI.supportsproblem(m::MosekSolver, ::Type{MOI.SingleVariable},                constraint_types) :: Bool = supportsconstraints(m,constraint_types)
#MOI.supportsproblem(m::MosekSolver, ::Type{MOI.ScalarAffineFunction{Float64}}, constraint_types) :: Bool = supportsconstraints(m,constraint_types)
#MOI.supportsproblem{F}(m::MosekSolver, ::Type{F}, constraint_types) :: Bool = false

ref2id(ref :: MOI.VariableIndex) :: Int = Int(ref.value)

ref2id(ref :: MOI.ConstraintIndex) :: Int =
    if ref.value & 1 == 0
        Int(ref.value >> 1)
    else
        - Int(ref.value >> 1)
    end

function id2vref(id :: Int) :: MOI.VariableIndex
    assert(id > 0)
    MOI.VariableIndex(id)
end

include("objective.jl")
include("variable.jl")
include("constraint.jl")
include("attributes.jl")


export MosekModel


end # module
