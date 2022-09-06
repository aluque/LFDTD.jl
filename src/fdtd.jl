#===============================================================================
  FDTD kernels
===============================================================================#


#===============================================================================
  Constants specific to this file
===============================================================================#

const INV_MU_0 = 1 / co.mu_0
const INV_EPSILON_0 = 1 / co.epsilon_0
const E_CHARGE = co.elementary_charge
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
    
    @batch for j in 1:n + 1
        for i in 1:m + 1
            i1 = i + l
            j1 = j + l
            
            hphi[i1, j1] +=
                INV_MU_0 * ((dt / dr) * (ez[i1, j1] - ez[i1 - 1, j1     ]) -
                            (dt / dz) * (er[i1, j1] - er[i1,     j1 - 1]))
        end
    end
end


"""
   Updates the r-component of the electric field and current.
"""
function update_er_jr(device, mesh, fields)
    (;R, H, m, n, l, dr, dz, dt) = mesh
    (;er, hphi, jr, ne, ngas) = fields
    
    @batch for j in 1:n
        for i in 1:m + 1
            i1 = i + l
            j1 = j + l
            
            # Diffusion removed because it creates strong fields from the
            # initial electron density.
            
            # jdiff = ((E_CHARGE * E_DIFFUSION / ngas[j1]) *
            #          (ne[i1, j1] - ne[i1 - 1, j1]) / dr)
            
            lhs = er[i1, j1] - dt * INV_EPSILON_0 *
                ((1 / dz) * (hphi[i1, j1 + 1] - hphi[i1, j1]) +
                 jr[i1, j1] / 2) # + jdiff / 2
            
            # The sign of lhs is that of the current at the midpoint.  We use that
            # to know the upwind direction at that time.
            # Note that j contains the electric current.  The electron velocity has
            # the opposite sign (damn Franklin).
            if lhs > 0
                ndown = ne[i1 - 1, j1]
                nup   = ne[i1,     j1]
                nup2  = ne[i1 + 1, j1]
            else
                ndown = ne[i1,     j1]
                nup   = ne[i1 - 1, j1]
                nup2  = ne[i1 - 2, j1]
            end
            
            # This should be generalized to arbitrary field-dependent mobility
            sigma = ((E_ELECTRON_MOBILITY / ngas[j1]) *
                     (nup + safe_koren(nup - nup2, ndown - nup)))
            
            if sigma < 0
                @show i j ngas[j1] nup nup2 ndown
                @error "sigma < 0"
            end
            
                
            er[i1, j1] = lhs / (1 + INV_EPSILON_0 * sigma * dt / 2)
            jr[i1, j1] = er[i1, j1] * sigma #+ jdiff
            
            # Boundary conditions at the conducting surface.
            if j == 1
                er[i1, j1 - 1] = -er[i1, j1]
            elseif j == n
                er[i1, j1 + 1] = er[i1, j1]
            end
        end
    end
end


"""
   Updates the z-component of the electric field and current.
"""
function update_ez_jz(device, mesh, fields)
    (;R, H, m, n, l, dr, dz, dt) = mesh
    (;ez, hphi, jz, ne, ngas) = fields
    
    @batch for j in 1:n + 1
        for i in 1:m
            i1 = i + l
            j1 = j + l
            
            r0    = (i - 1) * dr
            r1    = r0 + dr
            rmid  = r0 + dr / 2

            ngas1 = (ngas[j1] + ngas[j1 - 1]) / 2
            
            # Diffusion is currently removed
            # jdiff = ((E_CHARGE * E_DIFFUSION / ngas1) *
            #          (ne[i1, j1] - ne[i1, j1 - 1]) / dz)
            
            lhs = ez[i1, j1] + dt * INV_EPSILON_0 *
                ((1 / dr / rmid) * (r1 * hphi[i1 + 1, j1] - r0 * hphi[i1, j1]) +
                 jz[i1, j1] / 2) #  + jdiff / 2
            
            if lhs > 0
                ndown = ne[i1, j1 - 1]
                nup   = ne[i1, j1    ]
                nup2  = ne[i1, j1 + 1]
            else
                ndown = ne[i1, j1    ]
                nup   = ne[i1, j1 - 1]
                nup2  = ne[i1, j1 - 2]
            end
            
            sigma = ((E_ELECTRON_MOBILITY / ngas1) *
                     (nup + safe_koren(nup - nup2, ndown - nup)))

            # if INV_EPSILON_0 * sigma * dt > 0.2
            #     @show "Long time-steps" sigma dt (1 / (INV_EPSILON_0 * sigma))
            # end
            
            ez[i1, j1] = lhs / (1 + INV_EPSILON_0 * sigma * dt / 2)
            jz[i1, j1] = ez[i1, j1] * sigma # + jdiff
            
            #= Set the boundary condition for the axis here.
            Performance-wise I am not sure if this is the most efficient way because
            we are leaving many threads inactive.  But on the other hand setting
            it here is so simple perhaps it compensates the overhead of launching
            a new context just to do this. =#
            if i == 1
                ez[i1 - 1, j1] = ez[i1, j1]
            end        
        end
    end
end


"""
   Computes the magnitude of the electric field, eabs.
"""
function update_eabs(device, mesh, fields)
    (;m, n, l) = mesh
    (;ez, er, eabs) = fields
    
    @batch for j in 1:n
        for i in 1:m
            i1 = i + l
            j1 = j + l

            er1 = (er[i1, j1] + er[i1 + 1, j1    ]) / 2
            ez1 = (ez[i1, j1] + ez[i1,     j1 + 1]) / 2

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

    @batch for j in 1:n + 1
        for i in 1:m
            i1 = i + l
            j1 = j + l
            
            lr = (2i - 1) * dr / 2 - r0;
            lz = (j - 1) * dz - z0;
            
            jz[i1, j1] += a * exp(-(lr^2 / width^2 + lz^2 / h^2))
        end    
    end
end


#===============================================================================
  Handling electron density
===============================================================================#

"""
   Update the electron density including reactions and transport
"""
function update_ne(device, mesh, fields)
    (;R, H, m, n, l, dr, dz, dt) = mesh
    (;jr, jz, eabs, ne, ngas) = fields

    @batch for j in 1:n
        for i in 1:m
            i1 = i + l
            j1 = j + l

            rl = (i - 1) * dr
            rr = i * dr
            rc = (rl + rr) / 2
            zc = (2j - 1) / dz
            
            eabs_ = eabs[i1, j1]

            nu = impact_nu(ngas[j1], eabs_) - attachment_nu(ngas[j1], eabs_)

            # Note that j contains the electric current.  The particle current
            # is -j / e.
            adv = (MINUS_INV_E *
		   ((dt / dr / rc) * (rl * jr[i1, j1] -
				      rr * jr[i1 + 1, j1]) + 
		    (dt / dz) * (jz[i1, j1] -
			         jz[i1, j1 + 1])))

            # Crank-Nicolson
            ne[i1, j1] = ((adv +
			   (1 + nu * dt / 2) * ne[i1, j1]) / (1 - nu * dt / 2))
        end
    end
end


impact_nu(ngas, eabs) = (E_MOBILITY * eabs * IONIZATION_ALPHA *
                         exp(-ngas * IONIZATION_FIELD / eabs))
attachment_nu(ngas, eabs) = (E_MOBILITY * eabs * ATTACHMENT_ALPHA *
                             exp(-ngas * ATTACHMENT_FIELD / eabs))

"""
   Update the electron density at the boundaries 
   applying a linear relation given by `top`, `bottom` and `right`.

    For example for the upper bc:

     /          \\           /          \\ 
     | n[j]     |            | n[j - 1] |
     | n[j + 1] |   =  top * | n[j - 2] |
    \\          /           \\          /

    Use `top=I` for Neumann, `-I` for Dirichlet conditions. Same for `bottom`, 
    `right`.
"""
function set_ne_bc(device, mesh, fields; top, bottom, right)
    (;m, n, l) = mesh
    (;ne) = fields

    for i in 1:m
        i1 = i + l

        # Top
        j1 = n + l + 1
        x = SVector(ne[i1, j1 - 1], ne[i1, j1 - 2])
        (ne[i1, j1], ne[i1, j1 + 1]) = top * x

        # Bottom
        j1 = l
        x = SVector(ne[i1, j1 + 1], ne[i1, j1 + 2])
        (ne[i1, j1], ne[i1, j1 - 1]) = bottom * x
    end

    for j in 1:n
        j1 = j + l

        # Right
        i1 = m + l + 1
        x = SVector(ne[i1 - 1, j1], ne[i1 - 2, j1])
        (ne[i1, j1], ne[i1 + 1, j1]) = right * x
    end
end



#===============================================================================
  CPML layer
===============================================================================#

"""
    Updates the magnetic field in the CPML layer.  
"""
function cpml_update_h_outer(device, mesh, fields)
    (;R, H, m, n, l, dr, dz, dt) = mesh
    (;er, ez, hphi, psi_hphi_r, af, bf) = fields

    @batch for j in 1:n + 1
        for i in 1:l
            # These are indices in the full arrays (er, ez, ...)
            i1 = i + m + l + 1
            j1 = j + l

            # And these are the indices in the psi array
            i2 = i
            j2 = j + l

            dez = (1 / dr) * (ez[i1, j1] - ez[i1 - 1, j1])
            
            psi_hphi_r[i2, j2] = (bf[i] * psi_hphi_r[i2, j2] +
				  af[i] * INV_MU_0 * dt * dez)
            
            hphi[i1, j1] += (INV_MU_0 * (dt * dez -
		(dt / dz) * (er[i1, j1] - er[i1, j1 - 1]))
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

    @batch for j in 1:n
        for i in 1:l
            i1 = i + m + l + 1
            j1 = j + l

            er[i1, j1] -= (INV_EPSILON_0 * (dt / dz) *
                           (hphi[i1, j1 + 1] - hphi[i1, j1]))
        end
    end            
end


"""
    Updates the z-component of the electric field inside the CPML.
"""
function cpml_update_ez_outer(device, mesh, fields)
    (;R, H, m, n, l, dr, dz, dt) = mesh
    (;ez, hphi, psi_ez_r, ac, bc) = fields

    @batch for j in 1:n + 1
        for i in 1:l
            i1 = i + m + l
            j1 = j + l

            i2 = i
            j2 = j + l

            r0    = dr * (2i + 2m - 1) / 2
            r1    = r0 + dr
            rmid  = r0 + dr / 2

            dhphi = (1 / dr / rmid) *
                (r1 * hphi[i1 + 1, j1] - r0 * hphi[i1, j1])

            psi_ez_r[i2, j2] = (bc[i] * psi_ez_r[i2, j2] +
				ac[i] * dt * INV_EPSILON_0 * dhphi)
            ez[i1, j1] += dt * INV_EPSILON_0 * dhphi + psi_ez_r[i2, j2]
        end
    end            
end


#===============================================================================
  Accessory functions
===============================================================================#
@inline koren(x) = ifelse(x > 4.0, 1.0, clamp(x, 0., 1.0 / 3.0 + x / 6.0))


"""
    nan-safe Koren limiter function. Returns y * psi(x / y) making sure that
    we never divide by y.
"""
@inline function safe_koren(x, y)
    # theta = x / y
    (x >= 4y) && return y
    (x > 2y / 5) && return y / 3 + x / 6
    (x > 0) && return x
    return 0
end
