# TODO get with SingleVariable
function MOI.get(m::MosekModel,
                 ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}})
    refs = MOI.get(m, MOI.ListOfVariableIndices())
    cols = columns(m, refs)
    coeffs = getclist(m.task, cols)
    constant = getcfix(m.task)
    @assert length(coeffs) == length(refs)
    terms = MOI.ScalarAffineTerm{Float64}[
        MOI.ScalarAffineTerm(coeffs[i], refs[i]) for i in 1:length(refs)]
    return MOI.ScalarAffineFunction(terms, constant)
end

const ObjF = Union{MOI.SingleVariable, MOI.ScalarAffineFunction{Float64}}
MOI.supports(::MosekModel,::MOI.ObjectiveFunction{<:ObjF})  = true
MOI.supports(::MosekModel,::MOI.ObjectiveSense) = true


function MOI.set(m::MosekModel, ::MOI.ObjectiveFunction,
                 func::MOI.SingleVariable)
    numvar = getnumvar(m.task)
    c = zeros(Float64, numvar)
    col = column(m, func.variable)
    c[col] = 1.0
    putclist(m.task, convert(Vector{Int32}, 1:numvar), c)
    putcfix(m.task, 0.0)
end

function MOI.set(m::MosekModel, ::MOI.ObjectiveFunction,
                 func::MOI.ScalarAffineFunction{Float64})
    numvar = getnumvar(m.task)
    c = zeros(Float64, numvar)
    cols = columns(m, map(t -> t.variable_index, func.terms))
    for i in 1:length(cols)
        c[cols[i]] += func.terms[i].coefficient
    end

    putclist(m.task, convert(Vector{Int32}, 1:numvar), c)
    putcfix(m.task,func.constant)
end

function MOI.modify(m::MosekModel,
                    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
                    change :: MOI.ScalarConstantChange)
    putcfix(m.task,change.new_constant)
end

function MOI.modify(m::MosekModel,
                    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
                    change :: MOI.ScalarCoefficientChange)
    putcj(m.task, column(m, change.variable), change.new_coefficient)
end
