function MOI.get(::Optimizer, ::MOI.ObjectiveFunctionType)
    return MOI.ScalarAffineFunction{Float64}
end
function MOI.get(m::Optimizer, ::MOI.ObjectiveFunction{F}) where F
    obj = MOI.get(m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
    return convert(F, obj)
end
function MOI.get(m::Optimizer,
                 ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}})
    cis = MOI.get(m, MOI.ListOfVariableIndices())
    cols = columns(m, cis).values
    coeffs = getclist(m.task, cols)
    constant = getcfix(m.task)
    @assert length(coeffs) == length(cis)
    terms = MOI.ScalarAffineTerm{Float64}[
        MOI.ScalarAffineTerm(coeffs[i], cis[i]) for i in 1:length(cis)]
    # TODO add matrix terms
    return MOI.ScalarAffineFunction(terms, constant)
end

const ObjF = Union{MOI.VariableIndex, MOI.ScalarAffineFunction{Float64}}
MOI.supports(::Optimizer,::MOI.ObjectiveFunction{<:ObjF})  = true
MOI.supports(::Optimizer,::MOI.ObjectiveSense) = true


function MOI.set(m::Optimizer, ::MOI.ObjectiveFunction,
                 func::MOI.VariableIndex)
    MOI.set(m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
            convert(MOI.ScalarAffineFunction{Float64}, func))
end

function MOI.set(m::Optimizer, ::MOI.ObjectiveFunction,
                 func::MOI.ScalarAffineFunction{Float64})
    cols, values = split_scalar_matrix(m, MOIU.canonical(func).terms,
                                       (j, ids, coefs) -> putbarcj(m.task, j, ids, coefs))
    c = zeros(Float64, getnumvar(m.task))
    for (col, val) in zip(cols, values)
        c[col] += val
    end
    putclist(m.task, convert(Vector{Int32}, 1:length(c)), c)
    putcfix(m.task,func.constant)
    m.has_objective = true
    return
end

function MOI.modify(m::Optimizer,
                    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
                    change :: MOI.ScalarConstantChange)
    putcfix(m.task,change.new_constant)
    m.has_objective = true
    return
end

function MOI.modify(m::Optimizer,
                    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
                    change :: MOI.ScalarCoefficientChange)
    putcj(m.task, column(m, change.variable).value, change.new_coefficient)
    m.has_objective = true
    return
end
