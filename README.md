# MosekTools.jl

[![Build Status](https://github.com/jump-dev/MosekTools.jl/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/jump-dev/MosekTools.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/jump-dev/MosekTools.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jump-dev/MosekTools.jl)

[MosekTools.jl](https://github.com/jump-dev/MosekTools.jl) is the
[MathOptInterface.jl](https://github.com/jump-dev/MathOptInterface.jl)
implementation for the MOSEK solver.

The low-level solver API for MOSEK is found in the package [Mosek.jl](https://github.com/MOSEK/Mosek.jl).

## Affiliation

MosekTools.jl is maintained by the JuMP community and is not officially
supported by MOSEK. However, [Mosek.jl](https://github.com/MOSEK/Mosek.jl) _is_ an officially supported product of
MOSEK.

## Getting help

If you need help, please ask a question on the [JuMP community forum](https://jump.dev/forum).

If you have a reproducible example of a bug, please [open a GitHub issue](https://github.com/jump-dev/MosekTools.jl/issues/new).

## License

`MosekTools.jl` is licensed under the [MIT License](https://github.com/jump-dev/MosekTools.jl/blob/master/LICENSE.md).

The underlying solver is a closed-source commercial product for which you must
[obtain a license](https://www.mosek.com).

## Installation

Install MosekTools as follows:
```julia
import Pkg
Pkg.add("MosekTools")
```

In addition to installing the MosekTools.jl package, this will also download
and install the latest version of Mosek.jl.

Follow the instructions at [Mosek.jl](https://github.com/MOSEK/Mosek.jl) to
obtain and install an appropriate license.

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
Note that even though the optimizer is `Mosek.Optimizer`, you must additionally
import `MosekTools`.

## Options

All other parameters can be found in the [Mosek documentation](https://docs.mosek.com/latest/opt-server/param-groups.html).

For integer parameters, pass either the value, or the corresponding
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
