``MosekTools`` is the
[MathOptInterface.jl](https://github.com/JuliaOpt/MathOptInterface.jl)
implementation for the MOSEK solver. The low-level solver API for MOSEK is
found in the package [Mosek.jl](https://github.com/JuliaOpt/Mosek.jl).
The latest release of this package and the `mosekv8` branch are to be used with
the latest release of Mosek.jl (which uses MOSEK v8). To use MOSEK v9, use
the master branch of both this package and Mosek.jl.

The ``Mosek`` specific model object (used for example with JuMP) is created as
```julia
using MosekTools
model = Mosek.Optimizer()
```
hence to use Mosek in a JuMP model, do, e.g.,
```julia
using JuMP
using MosekTools
model = Model(with_optimizer(Mosek.Optimizer, QUIET=false, INTPNT_CO_TOL_DFEAS=1e-7))
```
The parameter `QUIET` is a special parameter that when set to `true`
disables all Mosek printing output.
All other parameters can be found in the [Mosek doc](https://docs.mosek.com/8.1/capi/param-groups.html#doc-param-groups).
Note that the prefix `MSK_IPAR_` (for integer parameters), `MSK_DPAR_` (for
floating point parameters) or `MSK_SPAR_` (for string parameters) are optional.
If they are not given, they are inferred from the type of the value. For
instance, in the example above, as `1e-7` is a floating point number, the
parameters name used is `MSK_DPAR_INTPNT_CO_TOL_DFEAS`.
