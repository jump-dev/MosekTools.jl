using MosekTools
using Test
using JuMP

include("jump_sdp.jl")
include("jump_lp.jl")
include("jump_soc.jl")

@testset "JuMP tests" begin
    o = Mosek.Optimizer()
    MOI.set(o, MOI.RawOptimizerAttribute("LOG"), 0)
    test_jump_lp()
    MOI.empty!(o)
    test_jump_soc(o)
    MOI.empty!(o)
    test_jump_sdp(o)
end
