using MathOptInterfaceMosek


using Base.Test

const MOI = MathOptInterface

include(joinpath(Pkg.dir("MathOptInterface"),"test","contlinear.jl"))
@testset "Continuous linear problems" begin
    contlineartest(MosekSolver(QUIET = true))
end

# include("contquadratic.jl")
# @testset "Continuous quadratic problems" begin
#     # contquadratictest(GurobiSolver())
# end

include(joinpath(Pkg.dir("MathOptInterface"),"test","contconic.jl"))
@testset "Continuous conic problems" begin
    contconictest(MosekSolver(QUIET = true))
end

include(joinpath(Pkg.dir("MathOptInterface"),"test","intlinear.jl"))
@testset "Mixed-integer linear problems" begin
    intlineartest(MosekSolver(QUIET = true))
end

include("test_jump.jl")
