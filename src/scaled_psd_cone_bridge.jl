# copy-pased from mosek.jl

"""
    struct ScaledPSDCone <: MOI.AbstractVectorSet
        side_dimension::Int
    end

Compared to MOI.PositiveSemidefiniteConeTriangle, the off-diagonal entries are scaled by `√2`
and it is vectorized with columnwise lower triangular instead of columnwise upper triangular.
"""
struct ScaledPSDCone <: MOI.AbstractVectorSet
    side_dimension::Int
end

function MOI.Utilities.set_with_dimension(::Type{ScaledPSDCone}, dim)
    return ScaledPSDCone(div(-1 + isqrt(1 + 8 * dim), 2))
end

Base.copy(x::ScaledPSDCone) = ScaledPSDCone(x.side_dimension)

MOI.side_dimension(x::ScaledPSDCone) = x.side_dimension

function MOI.dimension(x::ScaledPSDCone)
    return div(x.side_dimension * (x.side_dimension + 1), 2)
end

struct ScaledPSDConeBridge{T,F,G} <: MOI.Bridges.Constraint.SetMapBridge{
    T,
    ScaledPSDCone,
    MOI.PositiveSemidefiniteConeTriangle,
    F,
    G,
}
    constraint::MOI.ConstraintIndex{F,ScaledPSDCone}
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{ScaledPSDConeBridge{T}},
    ::Type{G},
    ::Type{MOI.PositiveSemidefiniteConeTriangle},
) where {T,G<:MOI.AbstractVectorFunction}
    S = MOI.Utilities.scalar_type(G)
    ST = MOI.Utilities.promote_operation(*, T, T, S)
    F = MOI.Utilities.promote_operation(vcat, T, ST)
    return ScaledPSDConeBridge{T,F,G}
end

function MOI.Bridges.map_set(
    ::Type{<:ScaledPSDConeBridge},
    set::MOI.PositiveSemidefiniteConeTriangle,
)
    return ScaledPSDCone(set.side_dimension)
end

function MOI.Bridges.inverse_map_set(
    ::Type{<:ScaledPSDConeBridge},
    set::ScaledPSDCone,
)
    return MOI.PositiveSemidefiniteConeTriangle(set.side_dimension)
end

function _transform_function(
    func::MOI.VectorAffineFunction{T},
    scale::T,
) where {T}
    terms = copy(func.terms)
    for i in eachindex(terms)
        if !MOI.Utilities.is_diagonal_vectorized_index(terms[i].output_index)
            terms[i] = MOI.VectorAffineTerm(
                terms[i].output_index,
                MOI.Utilities.operate_term(*, scale, terms[i].scalar_term),
            )
        end
    end
    return MOI.VectorAffineFunction(
        terms,
        _transform_function(func.constants, scale),
    )
end

function _transform_function(func::MOI.VectorOfVariables, scale::T) where {T}
    return MOI.VectorAffineFunction{T}(
        MOI.VectorAffineTerm{T}[
            MOI.VectorAffineTerm(
                i,
                MOI.ScalarAffineTerm(
                    MOI.Utilities.is_diagonal_vectorized_index(i) ? one(T) : scale,
                    func.variables[i]
                )
            )
            for i in eachindex(func.variables)
        ],
        zeros(T, length(func.variables)),
    )
    new_f = MOI.Utilities.operate(*, T, one(T), func)
    return _transform_function(new_f, scale)
end

function _transform_function(func::Vector{T}, scale::T) where {T}
    func = copy(func)
    for i in eachindex(func)
        if !MOI.Utilities.is_diagonal_vectorized_index(i)
            func[i] *= scale
        end
    end
    return func
end

# Map ConstraintFunction from MOI -> mosek
function MOI.Bridges.map_function(::Type{<:ScaledPSDConeBridge{T}}, f) where {T}
    return _transform_function(f, √convert(T, 2))
end

# Used to map the ConstraintPrimal from mosek -> MOI
function MOI.Bridges.inverse_map_function(::Type{<:ScaledPSDConeBridge{T}}, f) where {T}
    return _transform_function(f, inv(√convert(T, 2)))
end

# Used to map the ConstraintDual from mosek -> MOI
function MOI.Bridges.adjoint_map_function(::Type{<:ScaledPSDConeBridge{T}}, f) where {T}
    return _transform_function(f, inv(√convert(T, 2)))
end

# Used to set ConstraintDualStart
function MOI.Bridges.inverse_adjoint_map_function(
    ::Type{<:ScaledPSDConeBridge{T}},
    f,
) where {T}
    return _transform_function(f, √convert(T, 2))
end
