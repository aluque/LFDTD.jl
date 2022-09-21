module LFDTD

using DelimitedFiles
using StaticArrays
using LinearAlgebra
using DocStringExtensions
using JLD2
using Printf
using CUDA
using Polyester
using Interpolations
using CSV
using DSP
using ProgressMeter

ProgressMeter.ijulia_behavior(:append)

export en, tobs, valid, zcenter, zface, rcenter, rface, zgrid, rgrid
export BiGaussian, SoftStep

@template DEFAULT =
    """
    $(SIGNATURES)
    $(DOCSTRING)
    """

const DATA_DIR = joinpath(@__DIR__, "..", "data")
const NOTEBOOK_DIR = joinpath(@__DIR__, "..", "notebooks")

include("tbatch.jl")
include("constants.jl")
include("lookup.jl")
include("device.jl")
include("mesh.jl")
include("fields.jl")
include("fdtd.jl")
include("source.jl")
include("observation.jl")
include("io.jl")
include("electrons.jl")
include("timesteps.jl")
include("main.jl")

end
