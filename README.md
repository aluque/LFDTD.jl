# LFDTD
FDTD simulation of lightning pulses

## Installation
To install the code, first install [julia](https://julialang.org/). Then, from a julia prompt:
```julia
julia> using Pkg; Pkg.add(url="https://github.com/aluque/LFDTD.jl")
```

## Run an example
You can run a jupyter notebook with an example simulation with 
```julia
julia> using LFDTD
julia> using IJulia
julia> notebook(dir=LFDTD.NOTEBOOK_DIR)
```

If you have not installed IJulia, you will be asked whether to install it.  Note that with this procedure the notebook will be read-only so you cannot save your changes.  Select File -> Make a Copy... to start with a writtable notebook.

## Using multple threads
To run simulations using multiple threads, set the environment variable `JULIA_NUM_THREADS` to `auto` before calling julia.

