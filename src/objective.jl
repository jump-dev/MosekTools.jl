# Copyright (c) 2017: Ulf Worsøe, Mosek ApS
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

function MOI.get(::Optimizer, ::MOI.ObjectiveFunctionType)
    return MOI.ScalarAffineFunction{Float64}
end

function MOI.get(m::Optimizer, ::MOI.ObjectiveFunction{F}) where {F}
    obj = MOI.get(m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
    return convert(F, obj)
end

function MOI.get(
    m::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
)
    cis = MOI.get(m, MOI.ListOfVariableIndices())
    cols = columns(m, cis).values
    coeffs = getclist(m.task, cols)
    constant = getcfix(m.task)
    @assert length(coeffs) == length(cis)
    terms = MOI.ScalarAffineTerm{Float64}[
        MOI.ScalarAffineTerm(coeffs[i], cis[i]) for i in 1:length(cis)
    ]
    # TODO add matrix terms
    return MOI.ScalarAffineFunction(terms, constant)
end

function MOI.supports(
    ::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
)
    return true
end
MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true

function MOI.set(
    m::Optimizer,
    attr::MOI.ObjectiveFunction,
    func::MOI.ScalarAffineFunction{Float64},
)
    if m.has_psd_in_objective
        throw(
            MOI.SetAttributeNotAllowed(
                attr,
                "Cannot set a different objective if a previous objective was set including the contribution of the entry of a PSD variable.",
            ),
        )
    end
    cols, values = split_scalar_matrix(
        m,
        MOI.Utilities.canonical(func).terms,
        (j, ids, coefs) -> begin
            m.has_psd_in_objective = true
            putbarcj(m.task, j, ids, coefs)
        end,
    )
    c = zeros(Float64, getnumvar(m.task))
    for (col, val) in zip(cols, values)
        c[col] += val
    end
    putclist(m.task, convert(Vector{Int32}, 1:length(c)), c)
    putcfix(m.task, func.constant)
    m.has_objective = true
    return
end

function MOI.modify(
    m::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
    change::MOI.ScalarConstantChange,
)
    putcfix(m.task, change.new_constant)
    m.has_objective = true
    return
end

function MOI.modify(
    m::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
    change::MOI.ScalarCoefficientChange,
)
    if is_matrix(m, change.variable)
        throw(MOI.ModifyObjectiveNotAllowed(change, _MODIFY_PSD_VAR_ERROR))
    end
    putcj(m.task, column(m, change.variable).value, change.new_coefficient)
    m.has_objective = true
    return
end
