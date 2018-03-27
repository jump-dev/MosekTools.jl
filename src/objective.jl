const ObjF = Union{MOI.SingleVariable, MOI.ScalarAffineFunction{Float64}}

MOI.canget(::MosekModel, ::MOI.ObjectiveFunction{<:ObjF}) = true

# TODO get with SingleVariable
function MathOptInterface.get(m::MosekModel, ::MathOptInterface.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}})
    vid = allocatedlist(m.x_block)
    subj = getindexes(m.x_block, vid)
    coeffs = getclist(m.task, subj)
    constant = getcfix(m.task)
    MOI.ScalarAffineFunction(MOI.VariableIndex.(vid), coeffs, constant)
end

MOI.supports(::MosekModel,::MOI.ObjectiveFunction{<:ObjF})  = true
MOI.canset(::MosekModel,::MOI.ObjectiveFunction{<:ObjF})  = true

function MOI.set!(m::MosekModel, ::MOI.ObjectiveFunction, func::MOI.SingleVariable)
    numvar = getnumvar(m.task)
    c = zeros(Float64,numvar)
    vid = ref2id(func.variable)
    subj = getindexes(m.x_block,vid)

    c[subj[1]] = 1.0

    putclist(m.task,Int32[1:numvar...],c)
    putcfix(m.task,0.0)
end

function MOI.set!(m::MosekModel, ::MOI.ObjectiveFunction, func::MOI.ScalarAffineFunction{Float64})
    numvar = getnumvar(m.task)
    c = zeros(Float64,numvar)
    subj = getindexes(m.x_block, ref2id.(func.variables))
    for i in 1:length(subj)
        c[subj[i]] += func.coefficients[i]
    end

    putclist(m.task,Int32[1:numvar...],c)
    putcfix(m.task,func.constant)
end

MOI.canmodifyobjective(m::MosekModel, change :: MOI.ScalarConstantChange) = true
MOI.canmodifyobjective(m::MosekModel, change :: MOI.ScalarCoefficientChange) = true

function MOI.modifyobjective!(m::MosekModel, change :: MOI.ScalarConstantChange)
    putcfix(m.task,change.new_constant)
end

function MOI.modifyobjective!(m::MosekModel, change :: MOI.ScalarCoefficientChange)
    vid = ref2id(change.variable)
    subj = getindexes(m.x_block,vid)
    putcj(m.task,Int32(subj[1]),change.new_coefficient)
end
