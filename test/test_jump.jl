using MosekTools
using Test
using JuMP

include("jump_sdp.jl")
include("jump_lp.jl")
include("jump_soc.jl")

@testset "JuMP tests" begin
    test_jump_lp(Mosek.Optimizer(LOG=0))
    test_jump_soc(Mosek.Optimizer(LOG=0))
    test_jump_sdp(Mosek.Optimizer(LOG=0))
end
