

MathOptInterface.canget(::MosekModel, ::MathOptInterface.ObjectiveFunction) = true

function MathOptInterface.set!(m::MosekModel, ::MathOptInterface.ObjectiveFunction, func::MathOptInterface.SingleVariable)
    numvar = getnumvar(m.task)
    c = zeros(Float64,numvar)
    vid = ref2id(func.variable)
    subj = getindexes(m.x_block,vid)

    c[subj[1]] = 1.0

    putclist(m.task,Int32[1:numvar...],c)
    putcfix(m.task,0.0)
end

function MathOptInterface.set!(m::MosekModel, ::MathOptInterface.ObjectiveFunction, func::MathOptInterface.ScalarAffineFunction{Float64})
    numvar = getnumvar(m.task)
    c = zeros(Float64,numvar)
    vids = [ ref2id(vid) for vid in func.variables ]
    subj = Vector{Int}(length(vids))
    for i in 1:length(vids)
        getindexes(m.x_block,vids[i],subj,i)
    end

    for i in 1:length(subj)
        c[subj[i]] += func.coefficients[i]
    end

    putclist(m.task,Int32[1:numvar...],c)
    putcfix(m.task,func.constant)
end

MathOptInterface.canmodifyobjective(m::MosekModel, change :: MathOptInterface.ScalarConstantChange) = true
MathOptInterface.canmodifyobjective(m::MosekModel, change :: MathOptInterface.ScalarCoefficientChange) = true

function MathOptInterface.modifyobjective!(m::MosekModel, change :: MathOptInterface.ScalarConstantChange)
    putcfix(m.task,change.new_constant)
end

function MathOptInterface.modifyobjective!(m::MosekModel, change :: MathOptInterface.ScalarCoefficientChange)
    vid = ref2id(change.variable)
    subj = getindexes(m.x_block,vid)
    putcj(m.task,Int32(subj[1]),change.new_coefficient)
end
