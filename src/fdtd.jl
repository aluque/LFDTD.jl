#===============================================================================
  FDTD kernels
===============================================================================#


#===============================================================================
  Constants specific to this file
===============================================================================#

const INV_MU_0 = 1 / co.mu_0
const INV_EPSILON_0 = 1 / co.epsilon_0
const E_CHARGE = co.elementary_charge
const E_CHARGE_MASS = co.elementary_charge / co.electron_mass
const MINUS_INV_E = (-1 / co.elementary_charge)


# The electron transport and ionization parameters are defined for a given
# reference gas density (300 K, 1 bar).
const NGAS_REF = 2.4143235045419996e+25

const E_MOBILITY = 0.0372193 * NGAS_REF
const E_DIFFUSION = 0.18 * NGAS_REF
const IONIZATION_ALPHA = 433200.0 / NGAS_REF
const IONIZATION_FIELD = 2e7 / NGAS_REF
const ATTACHMENT_ALPHA = 2000.0 / NGAS_REF
const ATTACHMENT_FIELD = 3e6 / NGAS_REF
const E_ELECTRON_MOBILITY = E_CHARGE * E_MOBILITY

# Square of (the plasma freq per electron)

# 

#===============================================================================
  Handling Maxwell's equations
===============================================================================#

"""
   Updates the magnetic field.  Assumes that the ghost cells have been set 
   properly.
"""
function update_h(device, mesh, fields)    
    (;R, H, m, n, l, dr, dz, dt) = mesh
    (;er, ez, hphi) = fields
    
    @tbatch for j in 1:n
        for i in 1:m
            i1 = i + l
            j1 = j + l
            
            hphi[i1, j1] +=
                INV_MU_0 * ((dt / dr) * (ez[i1 + 1, j1] - ez[i1, j1]) -
                            (dt / dz) * (er[i1, j1 + 1] - er[i1, j1]))

            if j == 1
                hphi[i1, j1 - 1] = hphi[i1, j1]
            elseif j == n
                hphi[i1, j1 + 1] = hphi[i1, j1]
            end
        end
    end
end


"""
   Updates the r-component of the electric field.
"""
function update_er(device, mesh, fields)
    (;R, H, m, n, l, dr, dz, dt) = mesh
    (;er, hphi, jr) = fields
    
    @tbatch for j in 1:n + 1
        for i in 1:m
            i1 = i + l
            j1 = j + l
            
            er[i1, j1] -= INV_EPSILON_0 * (dt / dz) * (hphi[i1, j1] - hphi[i1, j1 - 1])
            er[i1, j1] -= INV_EPSILON_0 * (dt / 2) * (jr[i1, j1] + jr[i1 + 1, j1])            

            if (i == 1)
                er[i1 - 1, j1] = -er[i1, j1]
            end
        end
    end
end


"""
   Updates the z-component of the electric field.
"""
function update_ez(device, mesh, fields)
    (;R, H, m, n, l, dr, dz, dt) = mesh
    (;ez, hphi, jz, ne, ngas) = fields
    
    @tbatch for j in 1:n
        for i in 2:m + 1
            i1 = i + l
            j1 = j + l
            
            r0    = (2i - 3) * dr / 2
            r1    = r0 + dr
            rmid  = r0 + dr / 2

            ez[i1, j1] += (INV_EPSILON_0 * (dt / dr / rmid) *
                           (r1 * hphi[i1, j1] - r0 * hphi[i1 - 1, j1]))
            ez[i1, j1] -= INV_EPSILON_0 * (dt / 2) * (jz[i1, j1] + jz[i1, j1 + 1])

            if (j == 1)
                ez[i1, j1 - 1] = ez[i1, j1]
            elseif (j == n)
                ez[i1, j1 + 1] = ez[i1, j1]
            end
        end
    end

    @tbatch for j in 1:n
        j1 = j + l
        i1 = l + 1
        
        # See Inan & Marshall, (4.39)
        ez[i1, j1] += 4 * INV_EPSILON_0 * (dt / dr) * hphi[i1, j1]
        ez[i1, j1] -= INV_EPSILON_0 * (dt / 2) * (jz[i1, j1] + jz[i1, j1 + 1])

        if (j == 1)
            ez[i1, j1 - 1] = ez[i1, j1]
        elseif (j == n)
            ez[i1, j1 + 1] = ez[i1, j1]
        end

    end
end


"""
   Updates the r and z-components of the electric current.
"""
function update_j(device, mesh, fields)
    (;R, H, m, n, l, dr, dz, dt) = mesh
    (;ez, er, jz, jr, ne, ngas) = fields
    
    @tbatch for j in 1:n + 1
        for i in 1:m + 1
            i1 = i + l
            j1 = j + l
            
            # Note the missing 1/ϵ0, which we would need to remove later
            ωp2 = (co.elementary_charge^2 / co.electron_mass) * ne[i1, j1]

            # We compute this here to later allow for a field-dependent mobility
            ν = E_CHARGE_MASS * ngas[j1] / E_MOBILITY
            A = exp(-ν * dt)
            K = -expm1(-ν * dt) / ν

            jr[i1, j1] = A * jr[i1, j1] + K * ωp2 * (er[i1, j1] + er[i1 - 1, j1]) / 2
            jz[i1, j1] = A * jz[i1, j1] + K * ωp2 * (ez[i1, j1] + ez[i1, j1 - 1]) / 2
        end
    end
end


"""
   Computes the magnitude of the electric field, eabs.
"""
function update_eabs(device, mesh, fields)
    (;m, n, l) = mesh
    (;ez, er, eabs) = fields
    
    @tbatch for j in 1:n + 1
        for i in 1:m + 1
            i1 = i + l
            j1 = j + l

            er1 = (er[i1, j1] + er[i1 - 1, j1    ]) / 2
            ez1 = (ez[i1, j1] + ez[i1,     j1 - 1]) / 2

            eabs[i1, j1] = sqrt(er1^2 + ez1^2)            
        end
    end
end


"""
   Add a gaussian source to jz
"""
function gaussian_source(device, mesh, fields, source)
    (;R, H, m, n, l, dr, dz, dt) = mesh
    (;jz) = fields
    (;r0, z0, a, width, h) = source

    @tbatch for j in 1:n
        for i in 1:m
            i1 = i + l
            j1 = j + l
            
            lr = (i - 1) * dr - r0;
            lz = (j - 1) * dz - z0;

            jz[i1, j1] += a * exp(-(lr^2 / width^2 + lz^2 / h^2))
            
            ## ((j - 1) * dz > 20e3) && (jz[i1, j1] = 0)
        end    
    end
end




#===============================================================================
  Handling electron density
===============================================================================#

"""
   Update the electron density including reactions and transport
"""
function update_ne(device, mesh, fields, zmin=50 * co.kilo)
    (;R, H, m, n, l, dr, dz, dt) = mesh
    (;jr, jz, eabs, ne, ngas) = fields

    # Ignore ionization below zmin
    jmin = ceil(Int, zmin / dz)
    @tbatch for j in jmin:n
        for i in 1:m
            i1 = i + l
            j1 = j + l

            eabs_ = eabs[i1, j1]
            nu = impact_nu(ngas[j1], eabs_) - attachment_nu(ngas[j1], eabs_)

            # Crank-Nicolson
            ne[i1, j1] = (1 + nu * dt / 2) * ne[i1, j1] / (1 - nu * dt / 2)
        end
    end
end


impact_nu(ngas, eabs) = (E_MOBILITY * eabs * IONIZATION_ALPHA *
                         exp(-ngas * IONIZATION_FIELD / eabs))
attachment_nu(ngas, eabs) = (E_MOBILITY * eabs * ATTACHMENT_ALPHA *
                             exp(-ngas * ATTACHMENT_FIELD / eabs))


#===============================================================================
  CPML layer
===============================================================================#

"""
    Updates the magnetic field in the CPML layer.  
"""
function cpml_update_h_outer(device, mesh, fields)
    (;R, H, m, n, l, dr, dz, dt) = mesh
    (;er, ez, hphi, psi_hphi_r, ac, bc) = fields

    @tbatch for j in 1:n
        for i in 1:l
            # These are indices in the full arrays (er, ez, ...)
            i1 = i + m + l
            j1 = j + l

            # And these are the indices in the psi array
            i2 = i
            j2 = j + l

            dez = (1 / dr) * (ez[i1 + 1, j1] - ez[i1, j1])
            
            psi_hphi_r[i2, j2] = (bc[i] * psi_hphi_r[i2, j2] +
				  ac[i] * INV_MU_0 * dt * dez)
           
            hphi[i1, j1] += (INV_MU_0 * (dt * dez -
		(dt / dz) * (er[i1, j1 + 1] - er[i1, j1]))
                             + psi_hphi_r[i2, j2])

        end
    end
end


"""
    Updates the r-component of the electric field inside the CPML.
"""
function cpml_update_er_outer(device, mesh, fields)
    (;R, H, m, n, l, dr, dz, dt) = mesh
    (;er, hphi) = fields

    @tbatch for j in 1:n + 1
        for i in 1:l
            i1 = i + m + l
            j1 = j + l

            er[i1, j1] -= (INV_EPSILON_0 * (dt / dz) *
                           (hphi[i1, j1] - hphi[i1, j1 - 1]))
        end
    end            
end


"""
    Updates the z-component of the electric field inside the CPML.
"""
function cpml_update_ez_outer(device, mesh, fields)
    (;R, H, m, n, l, dr, dz, dt) = mesh
    (;ez, hphi, psi_ez_r, af, bf) = fields

    @tbatch for j in 1:n
        for i in 1:l
            i1 = i + m + l + 1
            j1 = j + l

            i2 = i
            j2 = j + l

            r0    = dr * (2i + 2m - 1) / 2
            r1    = r0 + dr
            rmid  = r0 + dr / 2

            dhphi = (1 / dr / rmid) *
                (r1 * hphi[i1, j1] - r0 * hphi[i1 - 1, j1])

            psi_ez_r[i2, j2] = (bf[i] * psi_ez_r[i2, j2] +
				af[i] * dt * INV_EPSILON_0 * dhphi)
            ez[i1, j1] += dt * INV_EPSILON_0 * dhphi + psi_ez_r[i2, j2]
        end
    end            
end


