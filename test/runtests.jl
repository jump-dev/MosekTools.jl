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


using Test

using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

const optimizer = MosekOptimizer(QUIET = true)
# 1e-3 needed for rotatedsoc3 test
const config = MOIT.TestConfig(atol=1e-3, rtol=1e-3, query=false)

# Mosek does not support names
MOIU.@model(Model,
            (),
            (MOI.EqualTo, MOI.LessThan),
            (MOI.Zeros, MOI.Nonnegatives,),
            (),
            (MOI.SingleVariable,),
            (MOI.ScalarAffineFunction,),
            (MOI.VectorOfVariables,),
            (MOI.VectorAffineFunction,))

@testset "SolverName" begin
    @test MOI.get(optimizer, MOI.SolverName()) == "Mosek"
end

@testset "Copy" begin
    # Currently does not work because get is missing for ConstraintSet
    # and ConstraintFunction, see https://github.com/JuliaOpt/MathOptInterfaceMosek.jl/issues/50
    #MOIT.copytest(optimizer, Model{Float64}())
end

@testset "Continuous linear problems" begin
    # linear1 is failing for two reasons
    # * it does not remove constraints using a variable if this variable is deleted, see https://github.com/JuliaOpt/MathOptInterface.jl/issues/511
    # * it does not support duplicated terms, see https://github.com/JuliaOpt/MathOptInterfaceMosek.jl/issues/41
    MOIT.contlineartest(optimizer, config, ["linear1"])
end

# include("contquadratic.jl")
# @testset "Continuous quadratic problems" begin
#     # contquadratictest(GurobiSolver())
# end

@testset "Continuous conic problems" begin
    MOIT.contconictest(MOIB.SquarePSD{Float64}(MOIB.RootDet{Float64}(MOIB.GeoMean{Float64}(optimizer))),
                       config, ["exp", "rootdets", "logdet"])
end

@testset "Mixed-integer linear problems" begin
    MOIT.intlineartest(optimizer, config, ["int2"])
end

#include("test_jump.jl")
