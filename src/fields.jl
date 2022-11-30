struct Fields{T, AT, VT}
    "electron density"
    ne::AT

    "r-component of the electric field"
    er::AT

    "r-component of current density"
    jr::AT

    "z-component of electric field"
    ez::AT

    "z-component of current density"
    jz::AT

    "phi-component of magnetic field"
    hphi::AT

    "magnitude of electric field"
    eabs::AT

    "Gas density at each altitude"
    ngas::VT

    # Fields of the cpml
    psi_hphi_r::AT
    psi_ez_r::AT
    psi_hphi_p::AT
    psi_er_p::AT
    psi_hphi_l::AT
    psi_er_l::AT

    # 1d fields for the cpml
    ac::VT
    bc::VT
    af::VT
    bf::VT
end

Base.eltype(f::Fields{T}) where {T} = T


"""
    Initialize fields with dimensions given by `mesh` and stored in `device`.
""" 
function Fields(T::Type, device, mesh)
    (;m, n, l) = mesh

    AT = array_type(device, T)
    VT = vector_type(device, T)
    
    dims = (2l + m + 1, 2l + n + 1)

    ne = zeros(device, T, dims...)
    er = zeros(device, T, dims...)
    jr = zeros(device, T, dims...)
    ez = zeros(device, T, dims...)
    jz = zeros(device, T, dims...)
    hphi = zeros(device, T, dims...)
    eabs = zeros(device, T, dims...)

    ngas = zeros(device, T, 2l + n + 1)
    
    # These are defined only on the right CPML
    psi_hphi_r = zeros(device, T, l, 2l + n + 1)
    psi_ez_r = zeros(device, T, l, 2l + n + 1)

    # These are defined on the top/bottom CPML
    psi_hphi_p = zeros(device, T, 2l + m + 1, l)
    psi_er_p = zeros(device, T, 2l + m + 1, l)
    psi_hphi_l = zeros(device, T, 2l + m + 1, l)
    psi_er_l = zeros(device, T, 2l + m + 1, l)

    (;ac, bc, af, bf) = init_cpml(T, device, mesh)
    
    return Fields{T, AT, VT}(ne,
                             er,
                             jr,
                             ez,
                             jz,
                             hphi,
                             eabs,
                             ngas,
                             psi_hphi_r,
                             psi_ez_r,
                             psi_hphi_p,
                             psi_er_p,
                             psi_hphi_l,
                             psi_er_l,
                             ac, bc, af, bf)
end

Fields(mesh, device) = Fields(eltype(mesh), mesh, device)

"""
    Initialize the CPML constant fields (ac, bc, af, bf)
"""
function init_cpml(T::Type, device, mesh)
    (;m, n, l, dr, dz, dt) = mesh

    ac = zeros(device, T, l)
    bc = zeros(device, T, l)
    af = zeros(device, T, l)
    bf = zeros(device, T, l)

    Z = sqrt(co.mu_0 / co.epsilon_0)
    @assert dr == dz "For the cpml we still require dr == dz"

    smax = -(l + 1) * log(1e-6) / (2 * Z * l * dr)
    drcl   = dr * ((1:l) .- 1//2)
    drfl   = dr * (1:l)

    sigmarc = @. smax * (drcl / (l * dr))^4
    sigmarf = @. smax * (drfl / (l * dr))^4

    copyto!(bc, @. exp(-sigmarc * dt / co.epsilon_0))
    copyto!(ac, bc .- 1)

    copyto!(bf, @. exp(-sigmarf * dt / co.epsilon_0))
    copyto!(af, bf .- 1)

    return (;ac, bc, af, bf)
end


"""
    Reduced field E/n in Td
"""
en(fields) = fields.eabs ./ reshape(fields.ngas, 1, :) ./ co.Td


"""
    Charge computed from divergence of electric field.
"""
function charge(fields, mesh)
    (;er, ez) = fields
    (;l, m, n) = mesh
    
    rc = rcenter(mesh)
    rf = rface(mesh)
    der = @views diff(fields.er[l + 1: l + m + 1, l + 1: l + m] .* rf, dims=1)
    dez = @views diff(fields.ez[l + 1: l + m, l + 1: l + m + 1], dims=2)
    @show size(der) size(rc) size(dez)
    q = co.epsilon_0 .* (der ./ rc .+ dez)

    return q
end

function photons(mesh, fields, excrate, nquench, zmin=50 * co.kilo)
    (;R, H, m, n, l, dr, dz, dt) = mesh
    (;eabs, ne, ngas) = fields

    photons = zeros(size(eabs))

        # Ignore emissions below zmin
    jmin = ceil(Int, zmin / dz)
    
    @tbatch for j in jmin:n + 1
        for i in 1:m + 1
            i1 = i + l
            j1 = j + l

            r = (i - 1) * dr
            z = (j - 1) * dz

            # We are here assuming instantaneous de-excitation             
            p = (N2_FRACTION * ne[i1, j1] * ngas[j1] *
                 lookup(excrate, eabs[i1, j1] / ngas[j1] / co.Td))
            p *= nquench / (ngas[j1] + nquench)
            photons[i1, j1] = p
        end
    end

    return photons
end
