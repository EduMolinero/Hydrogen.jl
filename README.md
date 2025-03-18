# Hydrogen.jl

The purpose of this `Julia` module is to solve the Time-dependent Schrödinger equation of the one-dimensional Hydrogen atom under the presence of a laser field.

## Usage

The use of the module is rather simple: it only needs an `input_file.ini` file (see the `examples/` folder for some use cases) plus the following `Julia` code:

```julia
include("Hydrogen.jl")
Hydrogen.solve_dynamics("input_file.ini")
