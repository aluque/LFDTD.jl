
function main(;end_t = nothing,
              device = CPU(),

              # Dimensions of the domain
              R = 800 * co.kilo,
              H = 110 * co.kilo,

              # Cells in r- and z-directions
              m = 800,
              n = 110,

              # Cells in CPML layers
              l = 32,

              # Time step
              dt = 1e-7,

              # Location of the observer (r=L, z=K)
              L = 84 * co.kilo,
              K = 400 * co.kilo,

              # Field of view
              #fov = deg2rad(40),
              fov = deg2rad(40),
              
              # Source properties:
              Ipeak = 257.0 * co.kilo,
              z0 = 10.0 * co.kilo,
              z1 = 15.0 * co.kilo,
              τ1 = .4e-4,
              τ2 = .2e-4,

              # Filenames atmospheric profiles
              gas_density_fname = joinpath(DATA_DIR, "earth", "gas.dat"),
              electron_density_fname = joinpath(DATA_DIR, "earth", "electrons.dat"),
              
              output_dt = 1e-4,
              output_folder = "/tmp"
              )

    mesh = Mesh(R, H, m, n, l, dt)
    @info "Space/time discretization [m, s]:" dr=mesh.dr dz=mesh.dz dt=mesh.dt

    fields = Fields(Float64, device, mesh)

    load_gas_density!(fields, mesh, gas_density_fname)
    load_electron_density!(fields, mesh, electron_density_fname)
    init_fields(mesh, fields)

    # Set source properties
    source = BiGaussian(FromPeak; Ipeak, z0, z1, τ1, τ2)
    @info "Total charge transferred by the source is [C]" Q=total_charge(source)

    # Load excitation rate for the N2(a1) state responsible for LBH emissions
    lbhexc = loadtable(joinpath(DATA_DIR, "swarm", "k023.dat"), f=log, g=log, ginv=exp)
    obs = Observer(L, K, fov, LinRange(0., 1e-2, 1000 + 1))
    obssteps = 100
    Δtobs = obssteps * mesh.dt
    
    # Quenching rate for the Ns(a1) state that emits LBH
    # Weightnig quenchings for N2 and O2 (see e.g. Liu 2004)
    kq = (0.79 * 1e-11 * co.centi^3 + 0.21 * 1e-10 * co.centi^3)

    # Einstein coeff (Liu 2004)
    A = 5.5e-5

    nquench = A / kq
    @info "Quenching density" nquench
    
    if isnothing(end_t)
        end_t = 1.1 * max(R, H) / co.c
    end
    nsteps = convert(Int, cld(end_t,  dt))
    
    t::Float64 = 0

    # Write output at intervals output_dt
    output = TimeStepper(0.0, output_dt)

    # Compute emissions at intervals Δtobs
    emissions = TimeStepper(0.0, Δtobs)
    progmeter = Progress(nsteps)
    saved_file = ""
    
    for s in 1:nsteps
        atstep(output, t) do step
            saved_file = joinpath(output_folder, @sprintf("%.5d.jld", step))
            save(saved_file, mesh, fields)
        end
        
        atstep(emissions, t) do step
            update_light_curve!(obs, device, mesh, fields, lbhexc, t,
                                Δtobs, nquench, ntheta=256)
        end
        
        advance!(fields, device, mesh, source, t)
        t += mesh.dt

        ProgressMeter.next!(progmeter, showvalues=[("Sim. time (ms)", t / co.milli),
                                                   (:saved_file, saved_file)])
    end
    
    return (;mesh, fields, obs)
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
    @views fields.ngas[(l + 1):(l + n + 1)] .= exp.(interp.(zface(mesh)))
    fields.ngas[l] = fields.ngas[l + 1]
    fields.ngas[l + n + 1] = fields.ngas[l + n]

    @info "Gas density loaded from $fname"
end

"""
    Load electron density data
"""
function load_electron_density!(fields, mesh, fname)
    (;l, m, n) = mesh
    f = CSV.File(fname, header=[:z, :ne])
    interp = LinearInterpolation(f.z .* co.kilo, log.(f.ne .* co.centi^-3))
    @views fields.ne[(l + 1):(l + m + 1), (l + 1):(l + n + 1)] .= reshape(exp.(interp.(zface(mesh))), (1, :))

    @info "Electron density loaded from $fname"
end


"""
   Steps time
"""
function advance!(fields, device, mesh, source, t)
    (;l, dr, dt) = mesh

    # Steps time to integer divisions (i -> i + 1)

    update_h(device, mesh, fields)
    cpml_update_h_outer(device, mesh, fields)

    update_j(device, mesh, fields)
    update_ne(device, mesh, fields)

    add_source(device, mesh, fields, source, t)
    
    # Steps time to half-integer divisions (i + 1/2 -> i + 3/2)


    update_er(device, mesh, fields)
    cpml_update_er_outer(device, mesh, fields)

    update_ez(device, mesh, fields)
    cpml_update_ez_outer(device, mesh, fields)
    
    update_eabs(device, mesh, fields)
end

