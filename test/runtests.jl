using MathOptInterfaceMosek

using Base.Test

const MOI = MathOptInterface

using MathOptInterfaceTests
const MOIT = MathOptInterfaceTests

const solver = () -> MosekInstance(QUIET = true)
const config = MOIT.TestConfig(1e-7, 1e-7, true, false, true, true)
const configgeomean = MOIT.TestConfig(1e-4, 1e-4, true, false, true, true)

@testset "Continuous linear problems" begin
    MOIT.contlineartest(solver, config)
end

# include("contquadratic.jl")
# @testset "Continuous quadratic problems" begin
#     # contquadratictest(GurobiSolver())
# end

@testset "Continuous conic problems" begin
    # lin1 and soc1 are failing because ListOfConstraints is not implemented
    # sdp0sv, sdp0sf, sdp1sv, sdp1sf and sdp2 are failing because PositiveSemidefiniteConeScaled is not working yet
    MOIT.contconictest(solver, config, ["lin1", "soc1", "sdp0sv", "sdp0sf", "sdp1sv", "sdp1sf", "sdp2", "geomean", "logdet", "rootdet1qv", "rootdet1qf"])
    MOIT.geomeantest(solver, configgeomean)
end

@testset "Mixed-integer linear problems" begin
    MOIT.intlineartest(solver, config, ["int2"])
end

include("test_jump.jl")
