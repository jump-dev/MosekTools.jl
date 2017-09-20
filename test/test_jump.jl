import MathOptInterfaceMosek
using Mosek.Ext
using Base.Test
using JuMP
const MOI = MathOptInterface

include("jump_sdp.jl")
include("jump_lp.jl")
include("jump_soc.jl")

@testset "JuMP tests" begin
    test_jump_lp(MathOptInterfaceMosek.MosekSolver(LOG=0))
    test_jump_soc(MathOptInterfaceMosek.MosekSolver(LOG=0))
    test_jump_sdp(MathOptInterfaceMosek.MosekSolver(LOG=0))
end
