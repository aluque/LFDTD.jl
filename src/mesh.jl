#=
  Define and manipulate meshes.

  We work on a cell grid of m x n plus additional l ghost cells in all 
  directions.

  Fields may have one or several staggered axes. The dimensions of the arrays
  containing the grids are always 2l + m + 1 (or 2l + n + 1) 
  even if non-staggered axes use only 2l + m values.
=#

struct Mesh{T}
    R::T
    H::T
    m::Int
    n::Int
    l::Int
    dr::T
    dz::T
    dt::T
    
    function Mesh(R, H, m::Int, n::Int, l::Int, dt)
        R, H = promote(R, H)
        T = typeof(R)
        
        dr = R / m
        dz = H / n
        #dt = min(dr, dz) / co.c / sqrt(3)
        
        new{T}(R, H, m, n, l, dr, dz, dt)
    end
end

Base.eltype(mesh::Mesh{T}) where T = T

zcenter(mesh) = LinRange(mesh.dz, mesh.H - mesh.dz, mesh.n)
zface(mesh) = LinRange(0, mesh.H, mesh.n + 1)
rcenter(mesh) = LinRange(mesh.dr, mesh.R - mesh.dr, mesh.m)
rface(mesh) = LinRange(0, mesh.R, mesh.m + 1)

zgrid(mesh, face) = face ? zface(mesh) : zcenter(mesh)
rgrid(mesh, face) = face ? rface(mesh) : rcenter(mesh)

"""
    Return a view of the array `a` defined on `mesh` excluding CPML.

    `rface` and `zface` indicate whether the array is defined on faces
    in the corresponding direction (then they have an extra value).
"""
function valid(a, mesh, rface::Bool=false, zface::Bool=false)
    (;m, n, l) = mesh
    @view a[l + 1:l + m + rface, l + 1:l + n + zface]
end
