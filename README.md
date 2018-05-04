``MathOptInterfaceMosek`` is the
[MathOptInterface.jl](https://github.com/JuliaOpt/MathOptInterface.jl)
implementation for the MOSEK solver. The low-level solver API for MOSEK is
found in the package [Mosek.jl](https://github.com/JuliaOpt/Mosek.jl).

The ``Mosek`` specific model object (used for example with JuMP) is created as 

```
model = MathOptInterfaceMosek.MosekOptimizer()
```

