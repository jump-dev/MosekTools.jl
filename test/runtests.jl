using MathOptInterfaceMosek

using Base.Test

using MathOptInterfaceTests
const MOIT = MathOptInterfaceTests

const solver = () -> MosekInstance(QUIET = true)
# 1e-3 needed for rotatedsoc3 test
const config = MOIT.TestConfig(atol=1e-3, rtol=1e-3, query=false)

@testset "Continuous linear problems" begin
    # linear1 is failing because NumberOfConstraints does not take deletion into account
    # linear11 is failing because the following are not implemented:
    # * MOI.cantransformconstraint(instance, c2, MOI.LessThan(2.0))
    # * MOI.get(instance, MathOptInterface.ConstraintFunction())
    MOIT.contlineartest(solver, config, ["linear1", "linear11"])
end

# include("contquadratic.jl")
# @testset "Continuous quadratic problems" begin
#     # contquadratictest(GurobiSolver())
# end

@testset "Continuous conic problems" begin
    # lin1 and soc1 are failing because ListOfConstraints is not implemented
    MOIT.contconictest(solver, config, ["lin1v", "lin1f", "soc1v", "soc1f", "exp", "logdet"])
end

@testset "Mixed-integer linear problems" begin
    MOIT.intlineartest(solver, config, ["int2"])
end

#include("test_jump.jl")
