[![Build Status](https://travis-ci.org/jump-dev/MosekTools.jl.svg?branch=master)](https://travis-ci.org/jump-dev/MosekTools.jl)

``MosekTools`` is the
[MathOptInterface.jl](https://github.com/jump-dev/MathOptInterface.jl)
implementation for the MOSEK solver. The low-level solver API for MOSEK is
found in the package [Mosek.jl](https://github.com/MOSEK/Mosek.jl).
The latest release of this package and the `master` branch are to be used with
the latest release of Mosek.jl (which uses MOSEK v10). To use MOSEK v9 (resp. v8), use
the v0.12.x (resp. v0.7.x) releases of this package or the `mosekv9` (resp. `mosekv8`) branch and v1.2.x (resp. v0.9.x) releases of Mosek.jl.

The ``Mosek`` specific model object (used for example with JuMP) is created as
```julia
using MosekTools
model = Mosek.Optimizer()
```
hence to use Mosek in a JuMP model, do, e.g.,
```julia
using JuMP
using MosekTools
model = Model(optimizer_with_attributes(Mosek.Optimizer, "QUIET" => false, "INTPNT_CO_TOL_DFEAS" => 1e-7))
```
The parameter `QUIET` is a special parameter that when set to `true`
disables all Mosek printing output.
All other parameters can be found in the [Mosek doc](https://docs.mosek.com/8.1/capi/param-groups.html#doc-param-groups).
Note that the prefix `MSK_IPAR_` (for integer parameters), `MSK_DPAR_` (for
floating point parameters) or `MSK_SPAR_` (for string parameters) are optional.
If they are not given, they are inferred from the type of the value. For
instance, in the example above, as `1e-7` is a floating point number, the
parameters name used is `MSK_DPAR_INTPNT_CO_TOL_DFEAS`.
