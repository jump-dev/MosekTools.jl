# MosekTools.jl

[MosekTools.jl](https://github.com/jump-dev/MosekTools.jl) is the
[MathOptInterface.jl](https://github.com/jump-dev/MathOptInterface.jl)
implementation for the MOSEK solver.

The low-level solver API for MOSEK is found in the package [Mosek.jl](https://github.com/MOSEK/Mosek.jl).

## Affiliation

MosekTools.jl is maintained by the JuMP community and is not officially
supported by MOSEK. However, Mosek.jl _is_ an officially supported product of
MOSEK.

## License

`MosekTools.jl` is licensed under the [MIT License](https://github.com/jump-dev/MosekTools.jl/blob/master/LICENSE.md).

The underlying solver is a closed-source commercial product for which you must
[obtainn a license](https://www.mosek.com).

## Installation

The latest release of this package and the `master` branch are to be used with
the latest release of Mosek.jl (which uses MOSEK v10).

To use MOSEK v9 (resp. v8), use the v0.12.x (resp. v0.7.x) releases of this
package, and the `mosekv9` (resp. `mosekv8`) branch and v1.2.x (resp. v0.9.x)
releases of Mosek.jl.

See the following table for a summary:

| MOSEK | Mosek.jl | MosekTools.jl release | MosekTools.jl branch |
|-------|----------|-----------------------|----------------------|
| v10   | v10      | v0.13                 | master               |
| v9    | v0.12    | v0.12                 | mosekv9              |
| v8    | v0.9     | v0.7                  | mosekv8              |

## Use with JuMP

```julia
using JuMP
using MosekTools
model = Model(Mosek.Optimizer)
set_attribute(model, "QUIET", true)
set_attribute(model, "INTPNT_CO_TOL_DFEAS", 1e-7)
```

## Options

The parameter `QUIET` is a special parameter that when set to `true`
disables all Mosek printing output.

All other parameters can be found in the [Mosek documentation](https://docs.mosek.com/8.1/capi/param-groups.html#doc-param-groups).

Note that the prefix `MSK_IPAR_` (for integer parameters), `MSK_DPAR_` (for
floating point parameters) or `MSK_SPAR_` (for string parameters) are optional.
If they are not given, they are inferred from the type of the value. For
example, in the example above, as `1e-7` is a floating point number, the
parameters name used is `MSK_DPAR_INTPNT_CO_TOL_DFEAS`.
