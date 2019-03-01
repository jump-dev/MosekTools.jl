###############################################################################
# TASK ########################################################################
###### The `task` field should not be accessed outside this section. ##########

## Affine Constraints #########################################################
##################### lc ≤ Ax ≤ uc ############################################

function rows(m::MosekModel,
              c::MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64}})::Vector{Int32} # TODO shouldn't need to convert
    return getindexes(m.c_block, ref2id(c))
end
function row(m::MosekModel,
             c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}})::Int32
    return getindex(m.c_block, ref2id(c))
end

function allocateconstraints(m::MosekModel, N::Int)
    numcon = getnumcon(m.task)
    alloced = ensurefree(m.c_block,N)
    id = newblock(m.c_block, N)

    M = numblocks(m.c_block) - length(m.c_block_slack)
    if alloced > 0
        appendcons(m.task, alloced)
        append!(m.c_constant, zeros(Float64,alloced))
    end
    if M > 0
        append!(m.c_block_slack, zeros(Float64,M))
        append!(m.c_coneid, zeros(Float64,M))
    end
    id
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

# Put the linear left-hand side
function addlhsblock!(m        :: MosekModel,
                      conid    :: Int,
                      conidxs  :: Vector{Int},
                      terms    :: Vector{MOI.ScalarAffineTerm{Float64}})
    consubi = getindexes(m.c_block,conid)
    cols = Int32[column(m, term.variable_index).value for term in terms]

    N = length(consubi)
    nnz = length(terms)

    msk_subi = convert(Vector{Int32}, consubi)

    msk_rowptr = zeros(Int64, N+1)
    for i in conidxs
        msk_rowptr[i+1] += 1
    end
    msk_rowptr[1] = 1
    for i in 2:N+1
        msk_rowptr[i] += msk_rowptr[i-1]
    end

    msk_subj = Array{Int32}(undef,nnz)
    msk_cof  = Array{Float64}(undef,nnz)

    # sort by row
    for i in 1:nnz
        msk_subj[msk_rowptr[conidxs[i]]] = cols[i]
        msk_cof[msk_rowptr[conidxs[i]]]  = terms[i].coefficient
        msk_rowptr[conidxs[i]] += 1
    end

    for i in N:-1:1
        msk_rowptr[i+1] = msk_rowptr[i]
    end
    msk_rowptr[1] = 1

    putarowlist(m.task, msk_subi, msk_rowptr[1:N], msk_rowptr[2:N+1], msk_subj, msk_cof)
end

addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MOI.Reals)                = putconboundlist(m.task,convert(Vector{Int32},conidxs),fill(MSK_BK_FR,MOI.dimension(dom)),-constant,-constant)
addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MOI.Zeros)                = putconboundlist(m.task,convert(Vector{Int32},conidxs),fill(MSK_BK_FX,MOI.dimension(dom)),-constant,-constant)
addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MOI.Nonnegatives)         = putconboundlist(m.task,convert(Vector{Int32},conidxs),fill(MSK_BK_LO,MOI.dimension(dom)),-constant,-constant)
addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MOI.Nonpositives)         = putconboundlist(m.task,convert(Vector{Int32},conidxs),fill(MSK_BK_UP,MOI.dimension(dom)),-constant,-constant)
addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MOI.GreaterThan{Float64}) = putconbound(m.task,Int32(conidxs[1]),MSK_BK_LO,dom.lower-constant[1],dom.lower-constant[1])
addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MOI.LessThan{Float64})    = putconbound(m.task,Int32(conidxs[1]),MSK_BK_UP,dom.upper-constant[1],dom.upper-constant[1])
addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MOI.EqualTo{Float64})     = putconbound(m.task,Int32(conidxs[1]),MSK_BK_FX,dom.value-constant[1],dom.value-constant[1])
addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MOI.Interval{Float64})    = putconbound(m.task,Int32(conidxs[1]),MSK_BK_RA,dom.lower-constant[1],dom.upper-constant[1])

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


function add_variable_constraint(m::MosekModel, col::Int32, dom::MOI.Interval)
    putvarbound(m.task, col, MSK_BK_RA, dom.lower, dom.upper)
end
function add_variable_constraint(m::MosekModel, col::Int32, dom::MOI.EqualTo)
    putvarbound(m.task, col, MSK_BK_FX, dom.value, dom.value)
end
function add_variable_constraint(m::MosekModel, col::Int32, ::MOI.Integer)
    putvartype(m.task, col, MSK_VAR_TYPE_INT)
end
function add_variable_constraint(m::MosekModel, col::Int32, ::MOI.ZeroOne)
    putvartype(m.task, col, MSK_VAR_TYPE_INT)
    putvarbound(m.task, col, MSK_BK_RA, 0.0, 1.0)
end
function add_variable_constraint(m::MosekModel, col::Int32, dom::MOI.LessThan)
    bk, lo, up = getvarbound(m.task, col)
    bk = if (bk == MSK_BK_FR) MSK_BK_UP else MSK_BK_RA end
    putvarbound(m.task, Int32(col), bk, lo, dom.upper)
end
function add_variable_constraint(m::MosekModel, col::Int32,
                                 dom::MOI.GreaterThan)
    bk, lo, up = getvarbound(m.task, col)
    bk = if (bk == MSK_BK_FR) MSK_BK_LO else MSK_BK_RA end
    putvarbound(m.task, Int32(col), bk, dom.lower, up)
end

add_variable_constraint(::MosekModel, ::Vector{Int32}, ::MOI.Reals) = nothing
function add_variable_constraint(m::MosekModel, cols::Vector{Int32},
                                 dom::MOI.Zeros)
    n = length(cols)
    @assert n == MOI.dimension(dom)
    bnd = zeros(Float64, n)
    putvarboundlist(m.task, cols, fill(MSK_BK_FX, n), bnd, bnd)
end
function add_variable_constraint(m::MosekModel, cols::Vector{Int32},
                                 dom::MOI.Nonnegatives)
    n = length(cols)
    @assert n == MOI.dimension(dom)
    bkx = Vector{Boundkey}(undef, n)
    blx = zeros(Float64, n)
    bux = zeros(Float64, n)
    for (i, col) in enumerate(cols)
        bk, lo, up = getvarbound(m.task, col)
        bkx[i] = bk == MSK_BK_FR ? MSK_BK_LO : MSK_BK_RA
        bux[i] = up
    end
    putvarboundlist(m.task, cols, bkx, blx, bux)
end

function add_variable_constraint(m::MosekModel, cols::Vector{Int32},
                                 dom::MOI.Nonpositives)
    n = length(cols)
    @assert n == MOI.dimension(dom)
    bkx = Vector{Boundkey}(undef, n)
    blx = zeros(Float64, n)
    bux = zeros(Float64, n)
    for (i, col) in enumerate(cols)
        bk, lo, up = getvarbound(m.task, col)
        bkx[i] = bk == MSK_BK_FR ? MSK_BK_UP : MSK_BK_RA
        blx[i] = lo
    end
    putvarboundlist(m.task, cols, bkx, blx, bux)
end

cone_parameter(dom :: MOI.PowerCone{Float64})     = dom.exponent
cone_parameter(dom :: MOI.DualPowerCone{Float64}) = dom.exponent
cone_parameter(dom :: C) where C <: MOI.AbstractSet = 0.0

const VectorCone = Union{MOI.SecondOrderCone,
                         MOI.RotatedSecondOrderCone,
                         MOI.PowerCone,
                         MOI.DualPowerCone,
                         MOI.ExponentialCone,
                         MOI.DualExponentialCone}

function set_variable_domain(m::MosekModel, xcid::Int, cols::Vector{Int32},
                             dom::VectorCone)
    appendcone(m.task, cone_type(dom), cone_parameter(dom), cols)
    coneidx = getnumcone(m.task)
    m.conecounter += 1
    putconename(m.task,coneidx,"$(m.conecounter)")
    m.xc_coneid[xcid] = m.conecounter
end
function set_variable_domain(m::MosekModel, ::Int, cols::Vector{Int32},
                             dom::MOI.AbstractSet)
    add_variable_constraint(m, cols, dom)
end

## Name #######################################################################
###############################################################################

function set_row_name(task::Mosek.MSKtask, row::Int32, name::String)
    putconname(task, row, name)
end

function set_row_name(m::MosekModel,
                      c::MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64}},
                      name::AbstractString) where {D}
    for row in rows(m, c)
        set_row_name(m.task, row, name)
    end
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

function allocatevarconstraints(m :: MosekModel,
                                N :: Int)
    nalloc = ensurefree(m.xc_block, N)
    id = newblock(m.xc_block, N)

    M = numblocks(m.xc_block) - length(m.xc_bounds)
    if M > 0
        append!(m.xc_bounds, zeros(Float64, M))
        append!(m.xc_coneid, zeros(Float64, M))
    end
    if nalloc > 0
        append!(m.xc_idxs, zeros(Float64, nalloc))
    end

    return id
end

###############################################################################
# MOI #########################################################################
###############################################################################

const PositiveSemidefiniteCone = MOI.PositiveSemidefiniteConeTriangle
const LinearFunction = Union{MOI.SingleVariable,
                             MOI.VectorOfVariables,
                             MOI.ScalarAffineFunction,
                             MOI.VectorAffineFunction}
const AffineFunction = Union{MOI.ScalarAffineFunction,
                             MOI.VectorAffineFunction}

const ScalarLinearDomain = Union{MOI.LessThan{Float64},
                                 MOI.GreaterThan{Float64},
                                 MOI.EqualTo{Float64},
                                 MOI.Interval{Float64}}
const VectorLinearDomain = Union{MOI.Nonpositives,
                                 MOI.Nonnegatives,
                                 MOI.Reals,
                                 MOI.Zeros}
const LinearDomain = Union{ScalarLinearDomain,VectorLinearDomain}
const ScalarIntegerDomain = Union{MOI.ZeroOne, MOI.Integer}

## Add ########################################################################
###############################################################################

MOI.supports_constraint(m::MosekModel, ::Type{<:Union{MOI.SingleVariable, MOI.ScalarAffineFunction}}, ::Type{<:ScalarLinearDomain}) = true
MOI.supports_constraint(m::MosekModel, ::Type{<:Union{MOI.VectorOfVariables, MOI.VectorAffineFunction}}, ::Type{<:Union{PositiveSemidefiniteCone, VectorLinearDomain}}) = true
MOI.supports_constraint(m::MosekModel, ::Type{MOI.VectorOfVariables}, ::Type{<:VectorCone}) = true
MOI.supports_constraint(m::MosekModel, ::Type{MOI.SingleVariable}, ::Type{<:ScalarIntegerDomain}) = true

## Affine Constraints #########################################################
##################### lc ≤ Ax ≤ uc ############################################

function MOI.add_constraint(m   :: MosekModel,
                            axb :: MOI.ScalarAffineFunction{Float64},
                            dom :: D) where {D <: MOI.AbstractScalarSet}

    # Duplicate indices not supported
    axb = MOIU.canonical(axb)

    N = 1
    conid = allocateconstraints(m,N)
    addlhsblock!(m,
                 conid,
                 fill(1, length(axb.terms)),
                 axb.terms)

    if length(m.c_constant) < length(m.c_block)
        append!(m.c_constant,
                zeros(Float64, length(m.c_block) - length(m.c_constant)))
    end
    conidxs = getindexes(m.c_block, conid)
    m.c_constant[conidxs] .= axb.constant

    addbound!(m, conid, conidxs, Float64[axb.constant], dom)
    conref = MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, D}(UInt64(conid) << 1)
    select(m.constrmap, MOI.ScalarAffineFunction{Float64}, D)[conref.value] = conid

    return conref
end

function MOI.add_constraint(m   :: MosekModel,
                            axb :: MOI.VectorAffineFunction{Float64},
                            dom :: D) where { D <: VectorLinearDomain }

    # Duplicate indices not supported
    axb = MOIU.canonical(axb)

    N = MOI.dimension(dom)
    conid = allocateconstraints(m,N)
    addlhsblock!(m,
                 conid,
                 map(t -> t.output_index, axb.terms),
                 map(t -> t.scalar_term, axb.terms))

    conidxs = getindexes(m.c_block, conid)
    m.c_constant[conidxs] .= axb.constants

    addbound!(m, conid, conidxs, axb.constants, dom)
    conref = MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64}, D}(UInt64(conid) << 1)
    select(m.constrmap, MOI.VectorAffineFunction{Float64}, D)[conref.value] = conid

    return conref
end

function MOI.add_constraint(m   :: MosekModel,
                            axb :: MOI.VectorAffineFunction{Float64},
                            dom :: PSDCone) where { PSDCone <: PositiveSemidefiniteCone }

    # Duplicate indices not supported
    axb = MOIU.canonical(axb)

    N = dom.side_dimension
    NN = MOI.dimension(dom)

    conid = allocateconstraints(m,NN)
    addlhsblock!(m,
                 conid,
                 map(t -> t.output_index, axb.terms),
                 map(t -> t.scalar_term, axb.terms))
    conidxs = getindexes(m.c_block,conid)
    constant = axb.constants
    m.c_constant[conidxs] = constant

    addbound!(m, conid, conidxs, constant, dom)
    conref = MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64},PSDCone}(UInt64(conid) << 1)
    select(m.constrmap,MOI.VectorAffineFunction{Float64},PSDCone)[conref.value] = conid
    conref
end

## Variable Constraints #######################################################
####################### lx ≤ x ≤ u ############################################
#######################      x ∈ K ############################################

# We allow following. Each variable can have
# - at most most upper and one lower bound
# - belong to at most one non-semidefinite cone
# - any number of semidefinite cones, which are implemented as ordinary constraints
# This is when things get a bit funky; By default a variable has no
# bounds, i.e. "free". Adding a `GreaterThan` or `Nonnegatives`
# constraint causes it to have a defined lower bound but no upper
# bound, allowing a `LessThan` or ``Nonpositives` constraint to be
# added later. Adding a `Interval` constraint defines both upper and
# lower bounds. Adding a `Reals` constraint will effectively be the
# same as an interval ]-infty;infty[, in that it will define both
# upper and lower bounds, and not allow those to be set afterwards.

domain_type_mask(dom :: MOI.Reals)        = (boundflag_lower | boundflag_upper)
domain_type_mask(dom :: MOI.Interval)     = (boundflag_lower | boundflag_upper)
domain_type_mask(dom :: MOI.EqualTo)      = (boundflag_lower | boundflag_upper)
domain_type_mask(dom :: MOI.Zeros)        = (boundflag_lower | boundflag_upper)
domain_type_mask(dom :: MOI.GreaterThan)  = boundflag_lower
domain_type_mask(dom :: MOI.Nonnegatives) = boundflag_lower
domain_type_mask(dom :: MOI.LessThan)     = boundflag_upper
domain_type_mask(dom :: MOI.Nonpositives) = boundflag_upper

domain_type_mask(dom :: MOI.SecondOrderCone)        = boundflag_cone
domain_type_mask(dom :: MOI.RotatedSecondOrderCone) = boundflag_cone
domain_type_mask(dom :: MOI.ExponentialCone)        = boundflag_cone
domain_type_mask(dom :: MOI.PowerCone)              = boundflag_cone
domain_type_mask(dom :: MOI.DualExponentialCone)    = boundflag_cone
domain_type_mask(dom :: MOI.DualPowerCone)          = boundflag_cone

domain_type_mask(dom :: MOI.PositiveSemidefiniteConeTriangle) = 0

domain_type_mask(dom :: MOI.Integer) = boundflag_int
domain_type_mask(dom :: MOI.ZeroOne) = (boundflag_int | boundflag_upper | boundflag_lower)

is_positivesemidefinite_set(dom :: PositiveSemidefiniteCone) = true
is_positivesemidefinite_set(dom) = false
is_conic_set(dom :: VectorCone) = true
is_conic_set(dom) = false

cone_type(::MOI.ExponentialCone)        = MSK_CT_PEXP
cone_type(::MOI.DualExponentialCone)    = MSK_CT_DEXP
cone_type(::MOI.PowerCone)              = MSK_CT_PPOW
cone_type(::MOI.DualPowerCone)          = MSK_CT_DPOW
cone_type(::MOI.SecondOrderCone)        = MSK_CT_QUAD
cone_type(::MOI.RotatedSecondOrderCone) = MSK_CT_RQUAD

function MOI.add_constraint(
    m   :: MosekModel,
    xs  :: MOI.SingleVariable,
    dom :: D) where {D <: MOI.AbstractScalarSet}

    col = column(m, xs.variable).value

    mask = domain_type_mask(dom)
    if mask & m.x_boundflags[col] != 0
        error("Cannot put multiple bound sets of the same type on a variable")
    end

    xcid = allocatevarconstraints(m, 1)

    xc_sub = getindex(m.xc_block, xcid)

    m.xc_bounds[xcid]  = mask
    m.xc_idxs[xc_sub] = col

    add_variable_constraint(m, col, dom)

    m.x_boundflags[col] |= mask

    conref = MOI.ConstraintIndex{MOI.SingleVariable,D}(UInt64(xcid) << 1)

    select(m.constrmap,MOI.SingleVariable,D)[conref.value] = xcid
    conref
end

function MOI.add_constraint(m::MosekModel, xs::MOI.VectorOfVariables,
                            dom::D) where {D <: MOI.AbstractVectorSet}
    cols = reorder(columns(m, xs.variables).values, D)

    mask = domain_type_mask(dom)
    if any(mask .& m.x_boundflags[cols] .> 0)
        error("Cannot multiple bound sets of the same type to a variable")
    end

    N = MOI.dimension(dom)
    xcid = allocatevarconstraints(m, N)
    xc_sub = getindexes(m.xc_block, xcid)

    m.xc_bounds[xcid] = mask
    m.xc_idxs[xc_sub] = cols

    set_variable_domain(m, xcid, cols, dom)

    m.x_boundflags[cols] .|= mask

    conref = MOI.ConstraintIndex{MOI.VectorOfVariables, D}(UInt64(xcid) << 1)
    select(m.constrmap,MOI.VectorOfVariables, D)[conref.value] = xcid
    return conref
end

################################################################################
################################################################################

function MOI.add_constraint(m   :: MosekModel,
                            xs  :: MOI.VectorOfVariables,
                            dom :: D) where { D <: PositiveSemidefiniteCone }
    N = dom.side_dimension
    if N < 2
        error("Invalid dimension for semidefinite constraint, got $N which is ",
              "smaller than the minimum dimension 2.")
    end

    vars = xs.variables
    cols = columns(m, xs.variables).values

    mask = domain_type_mask(dom)
    if any(mask .& m.x_boundflags[cols[1]] .> 0)
        error("Cannot put multiple bound sets of the same type to a variable")
    end

    NN = MOI.dimension(dom)

    # The scalar variables `xs.variables` should not have been created, we
    # should have created a matrix variable. To fix this we delete the scalar
    # variables and we create a matrix variable and make the scalar variables
    # point to the matrix variable

    # Delete the matrix variables

    # Create matrix variable
    appendbarvars(m.task, Int32[N])
    barvaridx = getnumbarvar(m.task)

    # Redirect scalar variables to the matrix ones

    if length(cols) != NN
        error("Mismatching variable length for semidefinite constraint")
    end

    id = allocateconstraints(m,NN)

    subi = getindexes(m.c_block,id)

    subii32 = convert(Vector{Int32},subi)
    putaijlist(m.task,
               subii32,
               cols,
               ones(Float64, NN))

    addbound!(m, id, subi, zeros(Float64, NN), dom)

    # HACK: We need to return a negative to indicate that this is
    # not, in fact, a real variable constraint, but rather a real
    # constraint, but to resolve return value at compile time we
    # need to disguise it as a variable constraint.
    #id2cref{MOI.VectorOfVariables,D}(-id)

    conref = MOI.ConstraintIndex{MOI.VectorOfVariables,MOI.PositiveSemidefiniteConeTriangle}(UInt64((id << 1) | 1))
    select(m.constrmap,MOI.VectorOfVariables,MOI.PositiveSemidefiniteConeTriangle)[conref.value] = id

    return conref
end



function addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MOI.PositiveSemidefiniteConeTriangle)
    dim = dom.side_dimension
    appendbarvars(m.task,Int32[dim])
    barvaridx = getnumbarvar(m.task)

    idx = 1
    for i in 1:dim
        for j in 1:i
            matrixid = appendsparsesymmat(m.task, Int32(dim), Int32[i], Int32[j], Float64[1.0])
            if i == j
                putbaraij(m.task, conidxs[idx], barvaridx, Int[matrixid], Float64[-1.0])
            else
                putbaraij(m.task, conidxs[idx], barvaridx, Int[matrixid], Float64[-0.5])
            end
            #putconname(m.task,Int32(conidxs[idx]),"bar_slack[$i,$j]")
            idx += 1
        end
    end

    putconboundlist(m.task, convert(Vector{Int32}, conidxs),
                    fill(MSK_BK_FX, length(constant)), -constant, -constant)


    m.c_block_slack[conid] = -barvaridx
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
                 xcref::MOI.ConstraintIndex{F,D},
                 dom::D) where { F    <: MOI.SingleVariable,
                                 D    <: ScalarLinearDomain }
    xcid = ref2id(xcref)
    col = m.xc_idxs[getindex(m.xc_block, xcid)]
    bk, bl, bu = getvarbound(m.task, col)
    bl, bu = chgbound(bl, bu, 0.0, dom)
    putvarbound(m.task, col, bk, bl, bu)
end


function MOI.set(m::MosekModel,
                 ::MOI.ConstraintSet,
                 cref::MOI.ConstraintIndex{F,D},
                 dom::D) where { F    <: MOI.ScalarAffineFunction,
                                 D    <: ScalarLinearDomain }
    cid = ref2id(cref)
    i = getindex(m.c_block,cid) # since we are in a scalar domain
    bk,bl,bu = getconbound(m.task,i)
    bl,bu = chgbound(bl,bu,0.0,dom)
    putconbound(m.task,i,bk,bl,bu)
end


### MODIFY
function MOI.modify(m   ::MosekModel,
                    c   ::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},D},
                    func::MOI.ScalarConstantChange{Float64}) where {D <: MOI.AbstractSet}

    cid = ref2id(c)

    i = getindex(m.c_block, cid)
    bk,bl,bu = getconbound(m.task,i)
    bl += m.c_constant[i] - func.new_constant
    bu += m.c_constant[i] - func.new_constant
    m.c_constant[i] = func.new_constant
    putconbound(m.task,i,bk,bl,bu)
end

function MOI.modify(m   ::MosekModel,
                    c   ::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},D},
                    func::MOI.ScalarCoefficientChange{Float64}) where {D <: MOI.AbstractSet}
    cid = ref2id(c)

    i = getindex(m.c_block, cid)
    j = column(m, func.variable).value

    putaij(m.task, i, j, func.new_coefficient)
end

function MOI.modify(m::MosekModel,
                    c::MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64},D},
                    func::MOI.VectorConstantChange{Float64}) where {D <: MOI.AbstractSet}
    cid = ref2id(c)
    @assert(cid > 0)

    subi = getindexes(m.c_block, cid)
    bk = Vector{Int32}(undef,length(subi))
    bl = Vector{Float64}(undef,length(subi))
    bu = Vector{Float64}(undef,length(subi))

    bk,bl,bu = getconboundlist(m.task,convert(Vector{Int32},subi))
    bl += m.c_constant[subi] - func.new_constant
    bu += m.c_constant[subi] - func.new_constant
    m.c_constant[subi] = func.new_constant

    putconboundlist(m.task, convert(Vector{Int32}, subi), bk, bl, bu)
end

function MOI.modify(m::MosekModel,
                    c::MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64},D},
                    func::MOI.MultirowChange{Float64}) where {D <: MOI.AbstractSet}
    cid = ref2id(c)
    @assert(cid > 0)

    subi = getindexes(m.c_block, cid)[getindex.(func.new_coefficients, 1)]
    j = column(m, func.variable).value

    putaijlist(m.task,convert(Vector{Int32},subi),fill(j,length(subi)),getindex.(func.new_coefficients,2))
end

### TRANSFORM
function MOI.transform(m::MosekModel,
                       cref::MOI.ConstraintIndex{F,D},
                       newdom::D) where {F <: MOI.AbstractFunction,
                                         D <: MOI.AbstractSet}
    MOI.modify(m,cref,newdom)
    cref
end

function MOI.transform(m::MosekModel,
                       cref::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},D1},
                       newdom::D2) where {D1 <: ScalarLinearDomain,
                                          D2 <: ScalarLinearDomain}
    F = MOI.ScalarAffineFunction{Float64}

    cid = ref2id(cref)

    subi = getindexes(m.c_block,cid)

    addbound!(m,cid,subi,m.c_constant[subi], newdom)

    newcref = MOI.ConstraintIndex{F,D2}(UInt64(cid) << 1)
    delete!(select(m.constrmap,F,D1), cref.value)
    select(m.constrmap, F, D2)[newcref.value] = cid
    newcref
end

function MOI.transform(m::MosekModel,
                       cref::MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64},D1},
                       newdom::D2) where {D1 <: VectorLinearDomain,
                                          D2 <: VectorLinearDomain}
    F = MOI.VectorAffineFunction{Float64}

    cid = ref2id(cref)

    subi = getindexes(m.c_block,cid)

    addbound!(m,cid,subi,m.c_constant[subi], newdom)

    newcref = MOI.ConstraintIndex{F,D2}(UInt64(cid) << 1)
    delete!(select(m.constrmap,F,D1), cref.value)
    select(m.constrmap,F,D2)[newcref.value] = cid
    newcref
end

## Delete #####################################################################
###############################################################################

function MOI.is_valid(model::MosekModel,
                      ref::MOI.ConstraintIndex{F, D}) where {F, D}
    return haskey(select(model.constrmap, F, D), ref.value)
end

function MOI.delete(
    m::MosekModel,
    cref::MOI.ConstraintIndex{F,D}) where {F <: AffineFunction,
                                           D <: Union{ScalarLinearDomain,
                                                      VectorLinearDomain}}
    if !MOI.is_valid(m, cref)
        throw(MOI.InvalidIndex(cref))
    end

    delete_name(m, cref)

    delete!(select(m.constrmap, F, D), cref.value)

    cid = ref2id(cref)
    subi = getindexes(m.c_block,cid)

    n = length(subi)
    subi_i32 = convert(Vector{Int32},subi)
    ptr = fill(Int64(0),n)
    putarowlist(m.task,subi_i32,ptr,ptr,Int32[],Float64[])
    b = fill(0.0,n)
    putconboundlist(m.task,subi_i32,fill(MSK_BK_FX,n),b,b)

    m.c_constant[subi] .= 0.0
    deleteblock(m.c_block,cid)
end

function MOI.delete(
    m::MosekModel,
    cref::MOI.ConstraintIndex{F,D}) where {F <: Union{MOI.SingleVariable,
                                                      MOI.VectorOfVariables},
                                           D <: Union{ScalarLinearDomain,
                                                      VectorLinearDomain,
                                                      ScalarIntegerDomain}}
    if !MOI.is_valid(m, cref)
        throw(MOI.InvalidIndex(cref))
    end

    delete_name(m, cref)

    delete!(select(m.constrmap, F, D), cref.value)

    xcid = ref2id(cref)
    sub = getindexes(m.xc_block,xcid)

    subj = [ getindex(m.x_block,i) for i in sub ]
    N = length(subj)

    m.x_boundflags[subj] .&= ~m.xc_bounds[xcid]
    if m.xc_bounds[xcid] & boundflag_int != 0
        for i in 1:length(subj)
            putvartype(m.task,subj[i],MSK_VAR_TYPE_CONT)
        end
    end

    if m.xc_bounds[xcid] & boundflag_lower != 0 && m.xc_bounds[xcid] & boundflag_upper != 0
        bnd = fill(0.0, length(N))
        putvarboundlist(m.task,convert(Vector{Int32},subj),fill(MSK_BK_FR,N),bnd,bnd)
    elseif m.xc_bounds[xcid] & boundflag_lower != 0
        bnd = fill(0.0, length(N))
        bk,bl,bu = getvarboundlist(m.task, convert(Vector{Int32},subj))

        for i in 1:N
            if MSK_BK_RA == bk[i] || MSK_BK_FX == bk[i]
                bk[i] = MSK_BK_UP
            else
                bk[i] = MSK_BK_FR
            end
        end

        putvarboundlist(m.task,convert(Vector{Int32},subj),bk,bl,bu)
    elseif m.xc_bounds[xcid] & boundflag_upper != 0
        bnd = fill(0.0, length(N))
        bk,bl,bu = getvarboundlist(m.task, subj)

        for i in 1:N
            if MSK_BK_RA == bk[i] || MSK_BK_FX == bk[i]
                bk[i] = MSK_BK_LO
            else
                bk[i] = MSK_BK_FR
            end
        end

        putvarboundlist(m.task,convert(Vector{Int32},subj),bk,bl,bu)
    else
        @assert(false)
        # should not happen
    end

    m.x_numxc[subj] .-= 1
    m.xc_idxs[sub] .= 0
    m.xc_bounds[xcid] = 0

    deleteblock(m.xc_block, xcid)
end

## List #######################################################################
###############################################################################

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
                 ci::MOI.ConstraintIndex{<:Union{MOI.ScalarAffineFunction{Float64},
                                                 MOI.VectorAffineFunction{Float64}}})
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
