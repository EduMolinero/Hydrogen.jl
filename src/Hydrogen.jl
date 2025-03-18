__precompile__()

module Hydrogen

    using SparseArrays, LinearAlgebra, Arpack
    using DifferentialEquations
    using Interpolations
    using CairoMakie
    using FFTW
    using Parsers
    using Optim
    using MKLSparse
    using Printf
    import IniFile

    include("Constants.jl")
    include("Input.jl")
    include("Model.jl")
    include("Laser.jl")
    include("Timeprop.jl")
    include("Plotting.jl")
    include("Fourier.jl")
    include("Main.jl")
    
end