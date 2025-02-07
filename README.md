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
[obtain a license](https://www.mosek.com).

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

To use Mosek with JuMP, use `Mosek.Optimizer`:
```julia
using JuMP
import Mosek
import MosekTools
model = Model(Mosek.Optimizer)
set_silent(model)
set_attribute(model, "MSK_DPAR_INTPNT_CO_TOL_DFEAS", 1e-7)
set_attribute(model, "MSK_IPAR_OPTIMIZER", Mosek.MSK_OPTIMIZER_INTPNT)
```
Note that even thought the optimizer is `Mosek.Optimizer`, you must additionally
import `MosekTools`.

## Options

All other parameters can be found in the [Mosek documentation]([https://docs.mosek.com/8.1/capi/param-groups.html#doc-param-groups](https://docs.mosek.com/latest/opt-server/param-groups.html)).

For integer parameters, pass either the value, or the correspondng
constant defined in the `Mosek` package.
```julia
using JuMP
import Mosek
import MosekTools
model = Model(Mosek.Optimizer)
set_attribute(model, "MSK_IPAR_OPTIMIZER", Mosek.MSK_OPTIMIZER_INTPNT)
set_attribute(model, "MSK_IPAR_OPTIMIZER", 4)
set_attribute(model, "MSK_IPAR_CACHE_LICENSE", Mosek.MSK_ON)
set_attribute(model, "MSK_IPAR_CACHE_LICENSE", 1)
set_attribute(model, "MSK_IPAR_CACHE_LICENSE", Mosek.MSK_OFF)
set_attribute(model, "MSK_IPAR_CACHE_LICENSE", 0)
```
