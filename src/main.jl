
function main(;outer_steps = 20, inner_steps=nothing)
    device = CPU()

    R::Float64 = 120 * co.kilo
    H::Float64 = 120 * co.kilo

    m = 500
    n = 500
    l = 32
    
    mesh = Mesh(R, H, m, n, l)
    @info "Space/time discretization [m, s]:" dr=mesh.dr dz=mesh.dz dt=mesh.dt

    fields = Fields(Float64, device, mesh)

    load_gas_density!(fields, mesh, joinpath(DATA_DIR, "earth", "gas.dat"))    
    load_electron_density!(fields, mesh, joinpath(DATA_DIR, "earth", "electrons.dat"))
    init_fields(mesh, fields)

    if isnothing(inner_steps)
        inner_steps = div(1e-4, mesh.dt)
        @show inner_steps
    end

    
    t::Float64 = 0
    
    for sout in 1:outer_steps
        @info "" t

        for sin in 1:inner_steps
            t += mesh.dt / 2
            step_1(device, mesh, fields, t)

            t += mesh.dt / 2
            step_2(device, mesh, fields, t)
        end

        save(joinpath("/tmp/", @sprintf("%.5d.jld", sout)), mesh, fields)
    end
    
    return (;mesh, fields)
end

"""
    Set initial conditions for the fields
"""
function init_fields(mesh, fields)
    # set_ne_bc(device, mesh, fields, top=I, bottom=-I, right=I)
    # fields.ne .= 0
end


"""
    Load gas density data
"""
function load_gas_density!(fields, mesh, fname)
    (;l, n) = mesh
    
    f = CSV.File(fname, header=[:z, :ngas])
    interp = LinearInterpolation(f.z .* co.kilo, log.(f.ngas .* co.centi^-3))
    @views fields.ngas[(l + 1):(l + n)] .= exp.(interp.(zcenter(mesh)))
    fields.ngas[l] = fields.ngas[l + 1]
    fields.ngas[l + n + 1] = fields.ngas[l + n]
end

"""
    Load electron density data
"""
function load_electron_density!(fields, mesh, fname)
    (;l, m, n) = mesh
    f = CSV.File(fname, header=[:z, :ne])
    interp = LinearInterpolation(f.z .* co.kilo, log.(f.ne .* co.centi^-3))
    @views fields.ne[(l + 1):(l + m), (l + 1):(l + n)] .= reshape(exp.(interp.(zcenter(mesh))), (1, :))
    @show maximum(fields.ne)

end


"""
   Steps time to integer divisions (i -> i + 1)
"""
function step_1(device, mesh, fields, t)
    update_h(device, mesh, fields)
    cpml_update_h_outer(device, mesh, fields)
    update_ne(device, mesh, fields)
    set_ne_bc(device, mesh, fields, top=I, bottom=I, right=I)

    tau = 1f-5
    i0 = 200 * co.kilo
    i = i0 * exp(-t^2 / tau^2) - exp(-(2 * t)^2 / tau^2)
    width=2f3
    
    source = (r0 = 0, z0 = 10f3, a = i / (π * width^2), width=width, h=10f3)
    gaussian_source(device, mesh, fields, source)    
end

"""
   Steps time to half-integer divisions (i + 1/2 -> i + 3/2)
"""
function step_2(device, mesh, fields, t)
    (;l, dr) = mesh

    update_er_jr(device, mesh, fields)
    cpml_update_er_outer(device, mesh, fields)

    update_ez_jz(device, mesh, fields)
    cpml_update_ez_outer(device, mesh, fields)
    
    update_eabs(device, mesh, fields)
end

