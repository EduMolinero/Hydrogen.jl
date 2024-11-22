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
    import IniFile

    include("constants.jl")
    include("input.jl")
    include("model.jl")
    include("laser.jl")
    include("timeprop.jl")
    include("plotting.jl")
    include("fourier.jl")
    include("main.jl")
    
end