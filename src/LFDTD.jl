module LFDTD

using StaticArrays
using LinearAlgebra
using DocStringExtensions
using JLD2
using Printf
using CUDA
using Polyester
using Interpolations
using CSV

@template DEFAULT =
    """
    $(SIGNATURES)
    $(DOCSTRING)
    """

const DATA_DIR = joinpath(@__DIR__, "..", "data")

include("constants.jl")
include("device.jl")
include("mesh.jl")
include("fields.jl")
include("fdtd.jl")
include("io.jl")
include("main.jl")

end
