
function main(;end_t = nothing,
              device = CPU(),

              # Dimensions of the domain
              R = 800 * co.kilo,
              H = 90 * co.kilo,

              cells_per_km = 1,
              
              # Cells in r- and z-directions
              m = cells_per_km * ceil(Int, R / co.kilo),
              n = cells_per_km * ceil(Int, H / co.kilo),

              # Cells in CPML layers
              l = 32,

              # Time step
              dt = 2e-7,

              # Altitude of the observer (r=L, z=K)
              K = 400 * co.kilo,

              # Field of view
              #fov = deg2rad(40),
              fov = deg2rad(40),
              
              # Source properties:
              # source = BiGaussian(FromPeak; 
              #                     Ipeak = 466.0 * co.kilo,
              #                     z0 = 10.0 * co.kilo,
              #                     z1 = 11.0 * co.kilo,
              #                     τ1 = 20e-6,
              #                     τ2 = 12.4e-6),

              # ~event 3 in Bjorge-Engeland et al.
              L = 462 * co.kilo,
              source = SoftStep(Ipeak = 466.0 * co.kilo,
                                z0 = 10.0 * co.kilo,
                                z1 = 11.5 * co.kilo,
                                rise = 12.4e-6,
                                decay = 8.2e-6),

              # L = 120 * co.kilo,
              # source = SoftStep(Ipeak = 72.0 * co.kilo,
              #                   z0 = 10.0 * co.kilo,
              #                   z1 = 11.5 * co.kilo,
              #                   rise = 11.6e-6,
              #                   decay = 7e-6),
              
              # Filenames atmospheric profiles
              gas_density_fname = joinpath(DATA_DIR, "earth", "stdatm.csv"),
              electron_density = LogInterpolatedElectronDensity(joinpath(DATA_DIR, "earth", "electrons.dat")),
              
              output_dt = 1e-4,
              output_folder = "/tmp",

              # Min. dielectric relaxation time to prevent instabilities
              # (electron density is cropped at the corresponding value).
              min_relax_time = 10 * dt,

              probe_locations = []
              )

    mesh = Mesh(R, H, m, n, l, dt)
    @info "Space/time discretization [m, s]:" dr=mesh.dr dz=mesh.dz dt=mesh.dt

    fields = Fields(Float64, device, mesh)

    load_gas_density!(fields, mesh, gas_density_fname)
    set_electron_density!(fields, mesh, electron_density)
    crop_electron_density!(fields, mesh, min_relax_time)
    
    init_fields(mesh, fields)

    # Set source properties
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
    A = 1.8e4

    nquench = A / kq
    @info "Quenching density" nquench
    
    if isnothing(end_t)
        end_t = 1.1 * max(R, H) / co.c
    end
    nsteps = convert(Int, cld(end_t,  dt))
    
    t::Float64 = 0

    # Ensure output_folder exists
    isdir(output_folder) || mkpath(output_folder)

    # Write output at intervals output_dt
    output = TimeStepper(0.0, output_dt)

    # Compute emissions at intervals Δtobs
    emissions = TimeStepper(0.0, Δtobs)
    progmeter = Progress(nsteps)
    saved_file = ""

    prob = [eltype(fields)[] for _ in probe_locations]
    probidx = probeinds(probe_locations, mesh)
    
    for s in 1:nsteps
        atstep(output, t) do step
            saved_file = joinpath(output_folder, @sprintf("%.5d.jld", step))
            save(saved_file, mesh, fields)
        end

        for i in eachindex(prob)
            push!(prob[i], fields.eabs[probidx[i]])
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

    saved_file = joinpath(output_folder, "lbh.jld")
    saveobs(saved_file, obs)
    @info "Observer data saved to $(saved_file)"
    
    saved_file = joinpath(output_folder, "probes.jld")
    time = (0:nsteps - 1) .* mesh.dt
    saveprobes(saved_file, probe_locations, time, prob)
    
    @info "Observer data saved to $(saved_file)"
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
    interp = LinearInterpolation(f.z, log.(f.ngas))
    @views fields.ngas[(l + 1):(l + n + 1)] .= exp.(interp.(zface(mesh)))
    fields.ngas[l] = fields.ngas[l + 1]
    fields.ngas[l + n + 1] = fields.ngas[l + n]

    @info "Gas density loaded from $fname"
end

"""
    Load electron density data
"""
function set_electron_density!(fields, mesh, electron_density)
    (;l, m, n) = mesh
    @views fields.ne[(l + 1):(l + m + 1), (l + 1):(l + n + 1)] .= reshape(electron_density.(zface(mesh)), (1, :))

end

"""
    Crops the electron density to ensure that the relaxation time does not
    go below the provided value).
"""
function crop_electron_density!(fields, mesh, min_relax_time)
    @info "Cropping electron density"
    max_ne = @. fields.ngas * co.epsilon_0 / (co.elementary_charge * E_MOBILITY * min_relax_time)
    fields.ne .= min.(fields.ne, reshape(max_ne, (1, :)))
end


function probeinds(locs, mesh)
    (;l, dr, dz) = mesh
    
    inds = map(locs) do (r, z)
        i = l + 1 + round(Int, r / dr)
        j = l + 1 + round(Int, z / dz)
        return CartesianIndex(i, j)
    end
    
    return inds
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

