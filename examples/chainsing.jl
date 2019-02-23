using JuMP
using MosekTools
import MathOptInterface
const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities

"""
Formulate the problem

      min.  sum_{j=0,...,(n-2)/2} s[j] + t[j] + p[j] + q[j]

      s.t.  (1/2, s[j], x[i]+10x[i+1]) \in Qr, 
            (1/2, t[j], 5^{1/2}*(x[i+2]-x[i+3])) \in Qr
            (1/2, r[j], (x[i+1]-2x[i+2])) \in Qr
	    (1/2, u[j], 10^{1/4}*(x[i]-10x[i+3])) \in Qr
	    (1/2, p[j], r[j]) \in Qr
	    (1/2, q[j], u[j]) \in Qr,                       j=0,...,(n-2)/2, i = 2j
        
	    0.1 <= x[i] <= 1.1,                             i=0,2,...,n-2

"""
function chainsing1(M :: Model,n :: Int)
    m = (n-2) >> 1
    
    @variable(M,x[1:n])
    @variable(M,p[1:m])
    @variable(M,q[1:m])
    @variable(M,r[1:m])
    @variable(M,s[1:m])
    @variable(M,t[1:m])
    @variable(M,u[1:m])

    for j in 1:m
        i = ((j - 1) << 1) + 1
        
        # s[j] >= (x[i] + 10*x[i+1])^2
        @constraint(M, [ 0.5, s[j], x[i]+10*x[i+1] ]            in MOI.RotatedSecondOrderCone(3))
        # t[j] >= 5*(x[i+2] - x[i+3])^2
        @constraint(M, [ 0.5, t[j], sqrt(5.0)*(x[i+2] - x[i+3]) ] in MOI.RotatedSecondOrderCone(3))
        # r[j] >= (x[i+1] - 2*x[i+2])^2
        @constraint(M, [ 0.5, r[j], x[i+1] - 2.0*x[i+2] ]         in MOI.RotatedSecondOrderCone(3))
        # 1/10 * u[j] >= (x[i] - 10*x[i+3])^2
        @constraint(M, [ 0.5 * 10^(-1.0), u[j], x[i] - 10.0*x[i+3] ] in MOI.RotatedSecondOrderCone(3))
        # p[j] >= r[j]^2
        @constraint(M, [ 0.5, p[j], r[j] ]                      in MOI.RotatedSecondOrderCone(3))
        # q[j] >= u[j]^2
        @constraint(M, [ 0.5, q[j], u[j] ]                      in MOI.RotatedSecondOrderCone(3))
    end

    #0.1 <= x[i] <= 1.1, i=0,2,...,n-2
    for i in 1:2:n-1
        @constraint(M,x[i] <= 1.1)
        @constraint(M,x[i] >= 0.1)
    end

    @objective(M, Min, sum(s) + sum(t) + sum(p) + sum(q))
end



function chainsing2(M :: Model,n :: Int)
    m = (n-2) >> 1
    
    @variable(M,x[1:n])
    @variable(M,p[1:m])
    @variable(M,q[1:m])
    @variable(M,r[1:m])
    @variable(M,s)
    @variable(M,u[1:m])    


    se = Array{Any}(2 * m + 2)
    

    #@constraint(M, [ 0.5, s,
    #                 @expression(M,    x[i]   + 10*x[i+1]  for i in 1:2:n-3 )...,
    #                 @expression(M, 5*(x[i+2] -    x[i+3]) for i in 1:2:n-3 )... ] in MOI.RotatedSecondOrderCone(2*m+2))
    @constraint(M, @expression(M,  [ 0.5; s; [ x[i]   + 10*x[i+1]  for i in 1:2:n-3 ] ; [ 5*(x[i+2] -    x[i+3]) for i in 1:2:n-3 ] ] ) in MOI.RotatedSecondOrderCone(2*m+2))
    for j in 1:m
        i = ((j - 1) << 1) + 1

        # r[j] >= (x[i+1] - 2*x[i+2])^2        
        @constraint(M, [ 0.5, r[j], x[i+1] - 2*x[i+2] ]         in MOI.RotatedSecondOrderCone(3))
        # u[j] >= sqrt(10)*(x[i] - 10*x[i+3])^2
        @constraint(M, [ 0.5 * 10.0^(-1.0), u[j], x[i] - 10*x[i+3] ] in MOI.RotatedSecondOrderCone(3))
        # p[j] >= r[j]^2
        @constraint(M, [ 0.5, p[j], r[j] ]                      in MOI.RotatedSecondOrderCone(3))
        # q[j] >= u[j]^2
        @constraint(M, [ 0.5, q[j], u[j] ]                      in MOI.RotatedSecondOrderCone(3))
    end

    #0.1 <= x[i] <= 1.1, i=0,2,...,n-2
    for i in 1:2:n-1
        @constraint(M,x[i] <= 1.1)
        @constraint(M,x[i] >= 0.1)
    end

    @objective(M, Min, sum(s) + sum(p) + sum(q))
end

function chainsing3(M :: Model,n :: Int)
    m = (n-2) >> 1
    
    @variable(M,x[1:n])
    @variable(M,p[1:m])
    @variable(M,q[1:m])
    @variable(M,r[1:m])
    @variable(M,s)
    @variable(M,u[1:m])    


    se = Array{Any}(2 * m + 2)
    

    @constraint(M, @expression(M, [ 0.5; s;
                                    [ x[i]   + 10*x[i+1]  for i in 1:2:n-3 ] ;
                                    [ 5*(x[i+2] -    x[i+3]) for i in 1:2:n-3 ] ;
                                    p ; q ]) in MOI.RotatedSecondOrderCone(4*m+2))
    for j in 1:m
        i = ((j - 1) << 1) + 1

        # r[j] >= (x[i+1] - 2*x[i+2])^2
        @constraint(M, [ 0.5, r[j], x[i+1] - 2*x[i+2] ]         in MOI.RotatedSecondOrderCone(3))
        # u[j] >= sqrt(10)*(x[i] - 10*x[i+3])^2
        @constraint(M, [ 0.5 * 10.0^(-1.0), u[j], x[i] - 10*x[i+3] ] in MOI.RotatedSecondOrderCone(3))
    end

    #0.1 <= x[i] <= 1.1, i=0,2,...,n-2
    for i in 1:2:n-1
        @constraint(M,x[i] <= 1.1)
        @constraint(M,x[i] >= 0.1)
    end

    @objective(M, Min, s)
end










function main(argv :: Vector{String})
    n = 8192
    backend = :generic
    method = 1
    timefile = nothing
    taskfile = nothing
    opt = false

    for a in argv
        m = match(r"--([^=]+)(?:=(.*))?",a)
        if m == nothing
            n = parse(Int,argv[1])
        else
            key = m.captures[1]
            val = m.captures[2]

            if key == "optimize"
                opt = true
            elseif key == "method"
                method = parse(Int,val)
            elseif key == "timefile"
                timefile = val
            elseif key == "out"
                taskfile = val
            elseif key == "backend"
                if val     == "generic" val == nothing
                    backend = :generic
                elseif val == "mock"
                    backend = :mock
                elseif val == "mosek"
                    backend = :mosek                    
                end
            elseif key == "help"
                println("Usage: chainsing n [options]")
                println("  n The problem scale, the default being 8192. Must be even.")
                println("Options:")
                println("\t--method=[1|2|3] Which formulation to use")
                println("\t--timefile=filename Which formulation to use")
                println("\t--backend=[generic|mosek|mock] Which backend to use.")
                println("\t\tgeneric Use generic backend, then copy to MOSEK optimizer")
                println("\t\tmosek   Use mosek as backend")
                println("\t\tmock    Use generic backend, then copy to JuMP mock optimizer")
                println("\t--out=filename Write task to this file")
                println("\t--optimize Call optimizer")
                return
            end
        end
    end


    #println("Chainsing")
    #println("\tn        : $n")
    #println("\tmethod   : $method")
    #println("\ttimefile : $timefile")
    #println("\ttaskfile : $taskfile")
    #println("\toptimize : $opt")
    #println("\tbackend  : $backend")    
    
    T = 0.0

    m =
        if  backend == :generic
            solver = Mosek.Optimizer()
            Model( optimizer = solver)
        elseif backend == :mosek
            solver = Mosek.Optimizer()
            Model( mode = JuMP.Direct, backend = solver)
        elseif backend == :mock
            solver = MOIU.MockOptimizer(JuMP.JuMPMOIModel{Float64}())
            Model( optimizer = solver)
        end


    
    if     backend == :generic
        println("Use: Generic backend, then copy to Mosek")
    elseif backend == :mosek
        println("Use: Mosek directly as backend")
    elseif backend == :mock
        println("Use: Generic backend, then copy to mock solver")
    end
    
    T0 = time()
    if method == 1
        chainsing1(m,n)
    elseif method == 2
        chainsing2(m,n)
    elseif method == 3
        chainsing3(m,n)
    end

    if backend != :mosek
        MOIU.attachoptimizer!(m)
    end
    T1 = time()    
    T = T1-T0

    @printf("Model build time: %.2f secs\n",T)
    if timefile != nothing
        open(timefile, "w") do f
            @printf(f,"%.2f", T)
        end
    end

end


main(ARGS)
