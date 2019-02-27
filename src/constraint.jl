const VectorCone = Union{MOI.SecondOrderCone,
                         MOI.RotatedSecondOrderCone,
                         MOI.PowerCone,
                         MOI.DualPowerCone,
                         MOI.ExponentialCone,
                         MOI.DualExponentialCone}
const PositiveSemidefiniteCone = MOI.PositiveSemidefiniteConeTriangle
const LinearFunction = Union{MOI.SingleVariable,
                             MOI.VectorOfVariables,
                             MOI.ScalarAffineFunction,
                             MOI.VectorAffineFunction}
const AffineFunction = Union{MOI.ScalarAffineFunction,
                             MOI.VectorAffineFunction}
const VariableFunction = Union{MOI.ScalarAffineFunction,
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
################################################################################
# ADD CONSTRAINT ###############################################################
################################################################################

MOI.supports_constraint(m::MosekModel, ::Type{<:Union{MOI.SingleVariable, MOI.ScalarAffineFunction}}, ::Type{<:ScalarLinearDomain}) = true
MOI.supports_constraint(m::MosekModel, ::Type{<:Union{MOI.VectorOfVariables, MOI.VectorAffineFunction}}, ::Type{<:Union{PositiveSemidefiniteCone, VectorLinearDomain}}) = true
MOI.supports_constraint(m::MosekModel, ::Type{MOI.VectorOfVariables}, ::Type{<:VectorCone}) = true
MOI.supports_constraint(m::MosekModel, ::Type{MOI.SingleVariable}, ::Type{<:ScalarIntegerDomain}) = true
#MOI.canaddconstraint(m::MosekModel, ::Type{<:Union{MOI.SingleVariable, MOI.ScalarAffineFunction}}, ::Type{<:ScalarLinearDomain}) = true
#MOI.canaddconstraint(m::MosekModel, ::Type{<:Union{MOI.VectorOfVariables, MOI.VectorAffineFunction}}, ::Type{<:Union{VectorCone, PositiveSemidefiniteCone, VectorLinearDomain}}) = true
#MOI.canaddconstraint(m::MosekModel, ::Type{MOI.SingleVariable}, ::Type{<:ScalarIntegerDomain}) = true

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

function trimapL(i, j, n)
    if i < j
        trimapL(j, i, n)
    else
        i + div((2n-j) * (j-1), 2)
    end
end

function trimapU(i::Integer, j::Integer)
    if i < j
        trimapU(j, i)
    else
        div((i-1)*i, 2) + j
    end
end

# Vectorized length for matrix dimension n
sympackedlen(n) = (n*(n+1)) >> 1
# Matrix dimension for vectorized length n
sympackeddim(n) = div(isqrt(1+8n) - 1, 2)

function _sympackedto(x, n, mapfrom, mapto)
    @assert length(x) == sympackedlen(n)
    y = similar(x)
    for i in 1:n, j in 1:i
        y[mapto(i, j)] = x[mapfrom(i, j)]
    end
    y
end
sympackedLtoU(x, n=sympackeddim(length(x))) = _sympackedto(x, n, (i, j) -> trimapL(i, j, n), trimapU)
sympackedUtoL(x, n=sympackeddim(length(x))) = _sympackedto(x, n, trimapU, (i, j) -> trimapL(i, j, n))

function sympackedUtoLidx(x::AbstractVector{<:Integer}, n)
    y = similar(x)
    map = sympackedLtoU(1:sympackedlen(n), n)
    for i in eachindex(y)
        y[i] = map[x[i]]
    end
    y
end

function MOI.add_constraint(m   :: MosekModel,
                            axb :: MOI.VectorAffineFunction{Float64},
                            dom :: PSDCone) where { PSDCone <: PositiveSemidefiniteCone }

    # Duplicate indices not supported
    axb = MOIU.canonical(axb)

    N = dom.side_dimension
    NN = sympackedlen(N)

    conid = allocateconstraints(m,NN)
    addlhsblock!(m,
                 conid,
                 sympackedUtoLidx(map(t -> t.output_index, axb.terms), N),
                 map(t -> t.scalar_term, axb.terms))
    conidxs = getindexes(m.c_block,conid)
    constant = sympackedUtoL(axb.constants, N)
    m.c_constant[conidxs] = constant

    addbound!(m,conid,conidxs,constant,dom)
    conref = MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64},PSDCone}(UInt64(conid) << 1)
    select(m.constrmap,MOI.VectorAffineFunction{Float64},PSDCone)[conref.value] = conid
    conref
end

################################################################################
# Variable constraints

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

function addvarconstr(m::MosekModel, subj::Int, dom::MOI.Interval)
    putvarbound(m.task, subj, MSK_BK_RA, dom.lower, dom.upper)
end
function addvarconstr(m::MosekModel, subj::Int, dom::MOI.EqualTo)
    putvarbound(m.task, subj, MSK_BK_FX, dom.value, dom.value)
end
function addvarconstr(m::MosekModel, subj::Int, dom::MOI.Integer)
    putvartype(m.task, subj, MSK_VAR_TYPE_INT)
end
function addvarconstr(m::MosekModel, subj::Int, dom::MOI.ZeroOne)
    putvartype(m.task, subj, MSK_VAR_TYPE_INT)
    putvarbound(m.task, subj, MSK_BK_RA, 0.0, 1.0)
end

function addvarconstr(m::MosekModel, subj::Int, dom::MOI.LessThan)
    bk, lo, up = getvarbound(m.task, subj)
    bk = if (bk == MSK_BK_FR) MSK_BK_UP else MSK_BK_RA end
    putvarbound(m.task, Int32(subj), bk, lo, dom.upper)
end

function addvarconstr(m::MosekModel, subj::Int, dom::MOI.GreaterThan)
    bk, lo, up = getvarbound(m.task, subj)
    bk = if (bk == MSK_BK_FR) MSK_BK_LO else MSK_BK_RA end
    putvarbound(m.task, Int32(subj), bk, dom.lower, up)
end

addvarconstr(m :: MosekModel, subj::Vector{Int}, dom :: MOI.Reals) = nothing
function addvarconstr(m :: MosekModel, subj::Vector{Int}, dom :: MOI.Zeros)
    bnd = zeros(Float64, length(subj))
    putvarboundlist(m.task, subj, fill(MSK_BK_FX, length(subj)), bnd, bnd)
end
function addvarconstr(m :: MosekModel, subj::Vector{Int}, dom :: MOI.Nonnegatives)
    bkx = Vector{Boundkey}(undef, length(subj))
    blx = zeros(Float64, length(subj))
    bux = zeros(Float64, length(subj))
    for (i,j) in enumerate(subj)
        bk,lo,up = getvarbound(m.task, j)
        bkx[i] = if (bk == MSK_BK_FR) MSK_BK_RA else MSK_BK_LO end
        bux[i] = up
    end
    putvarboundlist(m.task, subj, bkx, blx, bux)
end

function addvarconstr(m :: MosekModel, subj::Vector{Int}, dom :: MOI.Nonpositives)
    bkx = Vector{Boundkey}(undef, length(subj))
    blx = zeros(Float64, length(subj))
    bux = zeros(Float64, length(subj))
    for (i,j) in enumerate(subj)
        bk,lo,up = getvarbound(m.task, j)
        bkx[i] = if (bk == MSK_BK_FR) MSK_BK_RA else MSK_BK_UP end
        blx[i] = lo
    end
    putvarboundlist(m.task, subj, bkx, blx, bux)
end

abstractset2ct(dom::MOI.ExponentialCone)        = MSK_CT_PEXP
abstractset2ct(dom::MOI.DualExponentialCone)    = MSK_CT_DEXP
abstractset2ct(dom::MOI.PowerCone)              = MSK_CT_PPOW
abstractset2ct(dom::MOI.DualPowerCone)          = MSK_CT_DPOW
abstractset2ct(dom::MOI.SecondOrderCone)        = MSK_CT_QUAD
abstractset2ct(dom::MOI.RotatedSecondOrderCone) = MSK_CT_RQUAD

function MOI.add_constraint(
    m   :: MosekModel,
    xs  :: MOI.SingleVariable,
    dom :: D) where {D <: MOI.AbstractScalarSet}

    subj = getindex(m.x_block, ref2id(xs.variable))

    mask = domain_type_mask(dom)
    if mask & m.x_boundflags[subj] != 0
        error("Cannot put multiple bound sets of the same type on a variable")
    end

    xcid = allocatevarconstraints(m, 1)

    xc_sub = getindex(m.xc_block,xcid)

    m.xc_bounds[xcid]  = mask
    m.xc_idxs[xc_sub] = subj

    addvarconstr(m, subj, dom)

    m.x_boundflags[subj] |= mask

    conref = MOI.ConstraintIndex{MOI.SingleVariable,D}(UInt64(xcid) << 1)

    select(m.constrmap,MOI.SingleVariable,D)[conref.value] = xcid
    conref
end

function MOI.add_constraint(m   :: MosekModel,
                            xs  :: MOI.VectorOfVariables,
                            dom :: D) where { D <: PositiveSemidefiniteCone }
    N = dom.side_dimension
    vars = sympackedUtoL(xs.variables, N)
    subj = Vector{Int}(undef, length(vars))
    for i in 1:length(subj)
        getindexes(m.x_block, ref2id(vars[i]), subj, i)
    end

    mask = domain_type_mask(dom)
    if any(mask .& m.x_boundflags[subj[1]] .> 0)
        error("Cannot put multiple bound sets of the same type to a variable")
    end

    if N < 2
        error("Invalid dimension for semidefinite constraint")
    end

    NN = sympackedlen(N)

    if length(subj) != NN
        error("Mismatching variable length for semidefinite constraint")
    end

    id = allocateconstraints(m,NN)

    subi = getindexes(m.c_block,id)

    subii32 = convert(Vector{Int32},subi)
    putaijlist(m.task,
               subii32,
               convert(Vector{Int32},subj),
               ones(Float64,NN))

    addbound!(m,id,subi,zeros(Float64,NN),dom)

    #putconboundlist(m.task,subii32,fill(MSK_BK_FX,NN),zeros(Float64,NN),zeros(Float64,NN))
    #idx = 1
    #for j in 1:N
    #    for i in 1:N
    #        symmatidx = appendsparsesymmat(m.task,N,Int32[i],Int32[j],Float64[-1.0])
    #        putbaraij(m.task,subii32[idx],barvaridx,[symmatidx],Float64[1.0])
    #        idx += 1
    #    end
    #end

    # HACK: We need to return a negative to indicate that this is
    # not, in fact, a real variable constraint, but rather a real
    # constraint, but to resolve return value at compile time we
    # need to disguise it as a variable constraint.
    #id2cref{MOI.VectorOfVariables,D}(-id)

    conref = MOI.ConstraintIndex{MOI.VectorOfVariables,MOI.PositiveSemidefiniteConeTriangle}(UInt64((id << 1) | 1))
    select(m.constrmap,MOI.VectorOfVariables,MOI.PositiveSemidefiniteConeTriangle)[conref.value] = id

    conref
end



coneparfromset(dom :: MOI.PowerCone{Float64})     = dom.exponent
coneparfromset(dom :: MOI.DualPowerCone{Float64}) = dom.exponent
coneparfromset(dom :: C) where C <: MOI.AbstractSet = 0.0

function aux_setvardom(m::MosekModel,
                       xcid::Int,
                       subj::Vector{Int},
                       dom::D) where { D <: VectorCone }
    appendcone(m.task,abstractset2ct(dom),  coneparfromset(dom), subj)
    coneidx = getnumcone(m.task)
    m.conecounter += 1
    putconename(m.task,coneidx,"$(m.conecounter)")
    m.xc_coneid[xcid] = m.conecounter
end
function aux_setvardom(m::MosekModel, xcid::Int, subj::Vector{Int},
                       dom :: D) where {D <: MOI.AbstractSet}
    addvarconstr(m, subj, dom)
end

function MOI.add_constraint(m :: MosekModel, xs :: MOI.VectorOfVariables, dom :: D) where {D <: MOI.AbstractSet}
    subj = Vector{Int}(undef,length(xs.variables))
    for i in 1:length(subj)
        getindexes(m.x_block, ref2id(xs.variables[i]),subj,i)
    end

    mask = domain_type_mask(dom)
    if any(mask .& m.x_boundflags[subj] .> 0)
        error("Cannot multiple bound sets of the same type to a variable")
    end

    N = MOI.dimension(dom)
    xcid = allocatevarconstraints(m,N)
    xc_sub = getindexes(m.xc_block,xcid)

    m.xc_bounds[xcid] = mask
    m.xc_idxs[xc_sub] = subj

    aux_setvardom(m, xcid, subj, dom)

    m.x_boundflags[subj] .|= mask

    conref = MOI.ConstraintIndex{MOI.VectorOfVariables,D}(UInt64(xcid) << 1)
    select(m.constrmap,MOI.VectorOfVariables,D)[conref.value] = xcid
    conref
end

################################################################################
################################################################################


# Put the linear left-hand side
function addlhsblock!(m        :: MosekModel,
                      conid    :: Int,
                      conidxs  :: Vector{Int},
                      terms    :: Vector{MOI.ScalarAffineTerm{Float64}})
    consubi = getindexes(m.c_block,conid)
    subj = Vector{Int}(undef,length(terms))
    for i in 1:length(subj)
        getindexes(m.x_block,ref2id(terms[i].variable_index),subj,i)
    end

    N = length(consubi)
    nnz = length(terms)

    msk_subi = convert(Vector{Int32},consubi)

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
        msk_subj[msk_rowptr[conidxs[i]]] = subj[i]
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

# need a special case for this since MOI's variable order in PEXP cone is the reverse if Mosek's
function addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: D) where { D <: Union{MOI.ExponentialCone, MOI.DualExponentialCone} }
    N = MOI.dimension(dom)
    nalloc = ensurefree(m.x_block,N)

    varid = newblock(m.x_block,N)
    numvar = getnumvar(m.task)
    if nalloc > 0
        appendvars(m.task, length(m.x_block) - numvar)
        append!(m.x_boundflags, zeros(Int,length(m.x_block) - numvar))
        append!(m.x_numxc, zeros(Int,length(m.x_block) - numvar))
    end
    subj = getindexes(m.x_block,varid)
    subj = [ subj[3], subj[2], subj[1] ]

    putaijlist(m.task,conidxs,subj,-ones(Float64,N))
    putvarboundlist(m.task,subj,fill(MSK_BK_FR,N),zeros(Float64,N),zeros(Float64,N))
    putconboundlist(m.task,convert(Vector{Int32},conidxs),fill(MSK_BK_FX,N),-constant,-constant)

    m.c_block_slack[conid] = varid

    appendcone(m.task,abstractset2ct(dom),coneparfromset(dom),subj)
    coneidx = getnumcone(m.task)
    m.conecounter += 1
    #putconename(m.task,coneidx,"$(m.conecounter)")
    m.c_coneid[conid] = m.conecounter
end



function addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MOI.PositiveSemidefiniteConeTriangle)
    dim = dom.side_dimension
    appendbarvars(m.task,Int32[dim])
    barvaridx = getnumbarvar(m.task)

    idx = 1
    for j in 1:dim
        for i in j:dim
            matrixid = appendsparsesymmat(m.task,Int32(dim), Int32[i], Int32[j], Float64[1.0])
            if i == j
                putbaraij(m.task, conidxs[idx], barvaridx, Int[matrixid], Float64[-1.0])
            else
                putbaraij(m.task, conidxs[idx], barvaridx, Int[matrixid], Float64[-0.5])
            end
            #putconname(m.task,Int32(conidxs[idx]),"bar_slack[$i,$j]")
            idx += 1
        end
    end

    putconboundlist(m.task,convert(Vector{Int32},conidxs),fill(MSK_BK_FX,length(constant)),-constant,-constant)


    m.c_block_slack[conid] = -barvaridx
end

#function addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MOI.PositiveSemidefiniteConeScaled)
#    dim = MOI.dimension(dom)
#    #dim = floor(Int,-.5 + sqrt(.25+2*MOI.dimension(dom)))
#    n = MOI.dimension(dom)
#
#    appendbarvars(m.task,Int32[dim])
#    barvaridx = getnumbarvar(m.task)
#
#    idx = 1
#    for j in 1:dim
#        for i in j:dim
#            matrixid = appendsparsesymmat(m.task,Int32(dim), Int32[i], Int32[j], Float64[1.0])
#            if i == j
#                putbaraij(m.task, conidxs[idx], barvaridx, Int[matrixid], Float64[-1.0])
#            else
#                putbaraij(m.task, conidxs[idx], barvaridx, Int[matrixid], Float64[-sqrt(2.0)/2.0])
#            end
#            idx += 1
#        end
#    end
#
#    putconboundlist(m.task,convert(Vector{Int32},conidxs),fill(MSK_BK_FX,length(constant)),-constant,-constant)
#
#    m.c_block_slack[conid] = -barvaridx
#end



################################################################################
##  MODIFY #####################################################################
################################################################################

#MOI.canset(m::MosekModel, ::MOI.ConstraintSet, ::Type{MOI.ConstraintIndex{F,D}}) where { F <: LinearFunction, D <: LinearDomain } = true
#MOI.canset(m::MosekModel, ::MOI.ConstraintFunction, ::Type{MOI.ConstraintIndex{F,D}}) where { F <: Union{MOI.SingleVariable,MOI.ScalarAffineFunction}, D <: ScalarLinearDomain } = true
#MOI.canmodify(m::MosekModel, ::Type{MOI.ConstraintIndex{F,D}}, ::Type{<:MOI.AbstractFunctionModification}) where { F <: AffineFunction, D <: MOI.AbstractSet } = true

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
    j = m.xc_idxs[getindex(m.xc_block,xcid)]
    bk,bl,bu = getvarbound(m.task,j)
    bl,bu = chgbound(bl,bu,0.0,dom)
    putvarbound(m.task,j,bk,bl,bu)
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


function set_internal_name(m::MosekModel,c::MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64},D},name::AbstractString) where {D}
    cid = ref2id(c)
    for i in getindexes(m.c_block, cid)
        putconname(m.task,i,name)
    end
end
function set_internal_name(m::MosekModel,c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},D},name::AbstractString) where {D}
    cid = ref2id(c)
    i = getindex(m.c_block, cid)
    putconname(m.task,i,name)
end
function set_internal_name(m::MosekModel, c::MOI.ConstraintIndex{F,D}, name::AbstractString) where {F,D} end


function MOI.set(m    ::MosekModel,
                 ::MOI.ConstraintName,
                 c    ::MOI.ConstraintIndex{F,D},
                 name ::AbstractString) where {F,D}#{F<:MOI.AbstractFunction,D<:AbstractSet}
    if ! haskey(m.constrnames, name)
        m.constrnames[name] = MOI.ConstraintIndex[]
    end
    push!(m.constrnames[name], c)
    set_internal_name(m,c,name)
end

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
    xid = ref2id(func.variable)

    i = getindex(m.c_block,cid)
    j = getindex(m.x_block,xid)

    putaij(m.task,i,j,func.new_coefficient)
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

    putconboundlist(m.task,convert(Vector{Int32},subi),bk,bl,bu)
end

function MOI.modify(m::MosekModel,
                    c::MOI.ConstraintIndex{MOI.VectorAffineFunction{Float64},D},
                    func::MOI.MultirowChange{Float64}) where {D <: MOI.AbstractSet}
    cid = ref2id(c)
    @assert(cid > 0)

    subi = getindexes(m.c_block, cid)[getindex.(func.new_coefficients, 1)]
    xid = ref2id(func.variable)
    j = getindex(m.x_block,xid)

    putaijlist(m.task,convert(Vector{Int32},subi),fill(j,length(subi)),getindex.(func.new_coefficients,2))
end



#MOI.cantransform(m::MosekModel, c::MOI.ConstraintIndex{F,D1}, newdom::D2) where {F <: VariableFunction, D1, D2 } = false
#MOI.cantransform(m::MosekModel, c::MOI.ConstraintIndex{MOI.VectorAffineFunction,D1}, newdom::D2) where {D1 <: VectorLinearDomain, D2 <: VectorLinearDomain} = false
#function MOI.cantransform(m::MosekModel,
#                          cref::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},D1},
#                          newdom::D2) where {D1 <: ScalarLinearDomain,
#                                             D2 <: ScalarLinearDomain}
#    haskey(select(m.constrmap,MOI.ScalarAffineFunction{Float64},D1),cref.value)
#end


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

#MOI.transform(m::MosekModel, c::MOI.ConstraintIndex{F,D1}, newdom::D2) where {F <: MOI.VectorAffineFunction , D1 <: VectorLinearDomain, D2 <: VectorLinearDomain} = false

################################################################################
##  DELETE #####################################################################
################################################################################


#MOI.candelete(
#    m   ::MosekModel,
#    cref::MOI.ConstraintIndex{F,D}) where {F <: Union{MOI.ScalarAffineFunction,
#                                                                       MOI.VectorAffineFunction,
#                                                                       MOI.SingleVariable,
#                                                                       MOI.VectorOfVariables},
#                                                            D <: Union{MOI.LessThan,
#                                                                       MOI.GreaterThan,
#                                                                       MOI.EqualTo,
#                                                                       MOI.Interval,
#                                                                       MOI.Zeros,
#                                                                       MOI.Nonpositives,
#                                                                       MOI.Nonnegatives,
#                                                                       MOI.Reals}} = MOI.isvalid(m,cref)
#
#MOI.candelete(
#    m   ::MosekModel,
#    cref::MOI.ConstraintIndex{F,D}) where {F <: MOI.AbstractFunction,
#                                                            D <: Union{MOI.SecondOrderCone,
#                                                                       MOI.RotatedSecondOrderCone,
#                                                                       MOI.ExponentialCone,
#                                                                       MOI.DualExponentialCone,
#                                                                       MOI.PowerCone,
#                                                                       MOI.DualPowerCone}} = false
#
#MOI.candelete(
#    m   ::MosekModel,
#    cref::MOI.ConstraintIndex{F,D}) where {F <: Union{MOI.SingleVariable,
#                                                                       MOI.VectorOfVariables},
#                                                            D <: Union{MOI.SecondOrderCone,
#                                                                       MOI.RotatedSecondOrderCone,
#                                                                       MOI.ExponentialCone,
#                                                                       MOI.DualExponentialCone,
#                                                                       MOI.PowerCone,
#                                                                       MOI.DualPowerCone}} = true
#

function MOI.delete(
    m::MosekModel,
    cref::MOI.ConstraintIndex{F,D}) where {F <: Union{MOI.ScalarAffineFunction,
                                                                       MOI.VectorAffineFunction},
                                                            D <: Union{MOI.LessThan,
                                                                       MOI.GreaterThan,
                                                                       MOI.EqualTo,
                                                                       MOI.Interval,
                                                                       MOI.Zeros,
                                                                       MOI.Nonpositives,
                                                                       MOI.Nonnegatives,
                                                                       MOI.Reals}}

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






################################################################################
################################################################################
################################################################################

function allocateconstraints(m :: MosekModel,
                             N :: Int)
    numcon = getnumcon(m.task)
    alloced = ensurefree(m.c_block,N)
    id = newblock(m.c_block,N)

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


function allocatevarconstraints(m :: MosekModel,
                                N :: Int)
    nalloc = ensurefree(m.xc_block,N)
    id = newblock(m.xc_block,N)

    M = numblocks(m.xc_block) - length(m.xc_bounds)
    if M > 0
        append!(m.xc_bounds,zeros(Float64,M))
        append!(m.xc_coneid,zeros(Float64,M))
    end
    if nalloc > 0
        append!(m.xc_idxs, zeros(Float64,nalloc))
    end

    return id
end

function allocatevariable(m :: MosekModel, N :: Int)
    @assert(length(m.x_boundflags) == length(m.x_block))
    numvar = getnumvar(m.task)
    alloced = ensurefree(m.x_block, N)
    if alloced > 0
        appendvars(m.task, length(m.x_block) - numvar)
        append!(m.x_boundflags, zeros(Int,length(m.x_block) - numvar))
        append!(m.x_numxc, zeros(Int,length(m.x_block) - numvar))
    end
    return newblock(m.x_block, N)
end

function MOI.is_valid(model::MosekModel,
                      ref::MOI.ConstraintIndex{F, D}) where {F, D}
    return haskey(select(model.constrmap, F, D), ref.value)
end
function MOI.is_valid(model::MosekModel, ref::MOI.VariableIndex)
    return allocated(model.x_block, ref2id(ref))
end


function getvarboundlist(t::Mosek.Task, subj :: Vector{Int32})
    n = length(subj)
    bk = Vector{Boundkey}(undef,n)
    bl = Vector{Float64}(undef,n)
    bu = Vector{Float64}(undef,n)
    for i in 1:n
        bki,bli,bui = getvarbound(t,subj[i])
        bk[i] = bki
        bl[i] = bli
        bu[i] = bui
    end
    bk,bl,bu
end

function getconboundlist(t::Mosek.Task, subj :: Vector{Int32})
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
