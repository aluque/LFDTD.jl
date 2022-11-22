#= 
 Functions to compute the observed signal from a far-away observer (such as ASIM)
=#

N2_FRACTION = 0.79

"""
    Represent an observed located at a distance `L` from the central axes
    and at an altitude `K`.  It collects data at times given in `τ`.
"""
struct Observer{T, R <: LinRange}
    L::T
    K::T
    # Max. half-diagonal angle (in rad)
    fov::T

    τ::R

    # The observed signal
    signal::Matrix{T}

    function Observer(L, K, fov, τ)
        L, K, fov = promote(L, K, fov)
        T = typeof(L)
        
        signal = zeros(eltype(τ), (length(τ), Int(num_threads())))
        new{T, typeof(τ)}(L, K, fov, τ, signal)
    end
end



""" Compute the distance from a point (r, z, θ) to the observer `obs`. """
distance(obs, r, z, θ) = sqrt((z - obs.K)^2 + (obs.L - r * cos(θ))^2 + (r * sin(θ))^2)


""" Observation times substracting the light-travel from the origin. """
tobs(obs, z) = obs.τ .- distance(obs, 0, z, 0) / co.c
tobs(obs) = tobs(obs, 0)


"""
    Given an observer `obs`, cylindrical coordinates `r` and `z` and a time
    difference `Δt`, computes the angle `θ` such that the delay for light 
    to travel to the observer (i.e., L(obs, r, z, θ) / c) matches `Δt`.
"""
function theta_delay(obs, r, z, Δt)
    (;L, K) = obs
    L0 = co.c * Δt

    cosθ = ((K - z)^2 + r^2 + L^2 - L0^2) / (2 * L * r)

    return acos(cosθ)
end


"""
    Check if point with cylindrical coordinates `r`, `z`, `θ` is within
    the field of view of `obs` (square fov).
"""
function inside_fov_rect(obs, r, z, θ)
    (;L, K, fov) = obs
    uz = K - z
    ux = L - r * cos(θ)
    uy = r * sin(θ)
    umax = uz * tan(fov) / sqrt(2)

    return uy < umax && ux < umax
end

function inside_fov(obs, r, z, θ)
    (;L, K, fov) = obs
    uz = K - z
    ux = L - r * cos(θ)
    uy = r * sin(θ)
    umax = uz * tan(fov)

    return ux^2 + uy^2 < umax^2
end

"""
   Updates the r and z-components of the electric current.
"""
function update_light_curve!(obs, device, mesh, fields, excrate, t, Δt,
                             nquench; ntheta=16, zmin=50 * co.kilo)
    (;R, H, m, n, l, dr, dz, dt) = mesh
    (;eabs, ne, ngas) = fields
    (;L, K, τ, signal) = obs
    Δτ = step(τ)
    
    range_θ = LinRange(0, π, ntheta)
    dθ = step(range_θ)

    # Ignore emissions below zmin
    jmin = ceil(Int, zmin / dz)
    
    @tbatch for j in jmin:n + 1
        for i in 1:m + 1
            i1 = i + l
            j1 = j + l
            tid = Threads.threadid()

            r = (i - 1) * dr
            z = (j - 1) * dz

            # Trapezoidal weights
            # Note that i == 1 does not contribute but keeping it
            # anyway to prevent future bugs.
            wi = (i == 1) || (i == m + 1) ? 2 : 1
            wj = (j == 1) || (j == n + 1) ? 2 : 1
            wt = (t == 0) ? 2 : 1
            
            dV = r * Δt * dr * dz * dθ / wi / wj / wt
                        
            # We are here assuming instantaneous de-excitation             
            rW = (N2_FRACTION * ne[i1, j1] * ngas[j1] *
                  lookup(excrate, eabs[i1, j1] / ngas[j1] / co.Td))
            rW *= nquench / (ngas[j1] + nquench)
            
            # Thanks to symmetry, integrate onlt half-circle
            for iθ in eachindex(range_θ)
                # Perhaps this could be improved to reduce no. of comparisons.
                inside_fov(obs, r, z, range_θ[iθ]) || continue
                
                a = (iθ == firstindex(range_θ) || iθ == lastindex(range_θ)) ? 1 : 2

                L = distance(obs, r, z, range_θ[iθ])
                # The measurement point that we contribute to
                k = searchsortedfirst(τ, t + L / co.c)
                f = abs(τ[k] - (t + L / co.c)) / Δτ

                if k <= lastindex(τ) - 1
                    signal[k - 1, tid] += f * rW / (4π * L^2) * a * dV / Δτ
                end                
                if k <= lastindex(τ)
                    signal[k, tid] += (1 - f) * rW / (4π * L^2) * a * dV / Δτ
                end
            end
        end
    end
end

function light_curve(obs)
    return sum(obs.signal, dims=2)
end

function light_curve(obs, decay)
    s = light_curve(obs)
    c = @. exp(-obs.τ[1:end-1] / decay) - exp(-obs.τ[2:end] / decay)
    return conv(s, c)[begin:begin + length(s) - 1]
end
