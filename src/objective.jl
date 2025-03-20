# Copyright (c) 2017: Ulf Worsøe, Mosek ApS
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# MOI.ObjectiveFunctionType

function MOI.get(::Optimizer, ::MOI.ObjectiveFunctionType)
    return MOI.ScalarAffineFunction{Float64}
end

# MOI.ObjectiveFunction

function MOI.get(m::Optimizer, ::MOI.ObjectiveFunction{F}) where {F}
    obj = MOI.get(m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
    return convert(F, obj)
end

function MOI.get(
    m::Optimizer,
    attr::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
)
    if m.has_psd_in_objective
        msg = "Cannot get objective if it contains the contribution of the entry of a PSD variable."
        throw(MOI.GetAttributeNotAllowed(attr, msg))
    end
    # List of scalar variables only. Exclude matrix variables and those deleted.
    xs = [MOI.VariableIndex(i) for (i, x) in enumerate(m.x_sd) if x.matrix == 0]
    coeffs = Mosek.getclist(m.task, columns(m, xs).values)
    @assert length(coeffs) == length(xs)
    return MOI.ScalarAffineFunction{Float64}(
        MOI.ScalarAffineTerm{Float64}.(coeffs, xs),
        Mosek.getcfix(m.task),
    )
end

function MOI.supports(
    ::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
)
    return true
end

# MOI.ObjectiveSennse

MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true

function MOI.set(
    m::Optimizer,
    attr::MOI.ObjectiveFunction,
    func::MOI.ScalarAffineFunction{Float64},
)
    if m.has_psd_in_objective
        msg = "Cannot set a different objective if a previous objective was set including the contribution of the entry of a PSD variable."
        throw(MOI.SetAttributeNotAllowed(attr, msg))
    end
    cols, values = split_scalar_matrix(
        m,
        MOI.Utilities.canonical(func).terms,
        (j, ids, coefs) -> begin
            m.has_psd_in_objective = true
            Mosek.putbarcj(m.task, j, ids, coefs)
        end,
    )
    c = zeros(Float64, Mosek.getnumvar(m.task))
    for (col, val) in zip(cols, values)
        c[col] += val
    end
    Mosek.putclist(m.task, convert(Vector{Int32}, 1:length(c)), c)
    Mosek.putcfix(m.task, func.constant)
    m.has_objective = true
    return
end

function MOI.modify(
    m::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
    change::MOI.ScalarConstantChange,
)
    Mosek.putcfix(m.task, change.new_constant)
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
    Mosek.putcj(
        m.task,
        column(m, change.variable).value,
        change.new_coefficient,
    )
    m.has_objective = true
    return
end
