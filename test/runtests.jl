
# To ensure that tests can run on travis we have to do a little
# hackadoodle here. The tests require a license file. We include
# a license file that is only valid for one day (the day when
# change is submitted).
# If there is no valid license file, we default to that file.


if haskey(ENV,"MOSEKLM_LICENSE_FILE")
    # that's nice
elseif haskey(ENV,"HOME")
    if isfile(joinpath(ENV["HOME"],"mosek","mosek.lic"))
        # our lucky day!
    else
        licfile = joinpath(@__DIR__,"..","test",".dontuse-probablyexpired.lic")
        println("Use license file: $licfile")
        import Mosek
        Mosek.putlicensepath(Mosek.msk_global_env,licfile)
    end
end

using MathOptInterfaceMosek


using Base.Test

using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIB = MOI.Bridges

const optimizer = MosekOptimizer(QUIET = true)
# 1e-3 needed for rotatedsoc3 test
const config = MOIT.TestConfig(atol=1e-3, rtol=1e-3, query=false)

@testset "Continuous linear problems" begin
    # linear1 is failing because NumberOfConstraints does not take deletion into account
    # linear11 is failing because the following are not implemented:
    # * MOI.cantransformconstraint(instance, c2, MOI.LessThan(2.0))
    # * MOI.get(instance, MathOptInterface.ConstraintFunction())
    # linear13 is failing because it is FeasibilitySense
    MOIT.contlineartest(optimizer, config, ["linear1", "linear11", "linear13"])
end

# include("contquadratic.jl")
# @testset "Continuous quadratic problems" begin
#     # contquadratictest(GurobiSolver())
# end

@testset "Continuous conic problems" begin
    # lin1 and soc1 are failing because ListOfConstraints is not implemented
    MOIT.contconictest(MOIB.RootDet{Float64}(MOIB.GeoMean{Float64}(optimizer)), config, ["lin1v", "lin1f", "soc1v", "soc1f", "exp", "psds", "rootdets", "logdet"])
end

@testset "Mixed-integer linear problems" begin
    MOIT.intlineartest(optimizer, config, ["int2"])
end

#include("test_jump.jl")
