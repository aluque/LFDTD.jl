#===============================================================================
  Bi-Gaussian
===============================================================================#
"""
    A representation of a bi-gaussian, axial source.
"""
Base.@kwdef struct BiGaussian{T}
    "Lowest altitude bound."
    z0::T

    "Upper altitude bound."
    z1::T

    "Peak current."
    I0::T

    "Approximate duration."
    τ1::T

    "Approximate rise-time."
    τ2::T
end

const FromPeak = Val(:from_peak)

function BiGaussian(::typeof(FromPeak); z0, z1, Ipeak, τ1, τ2)
    k = (τ2 / τ1)^(2 * τ2^2 / (τ1^2 - τ2^2)) - (τ2 / τ1)^(2 * τ1^2 / (τ1^2 - τ2^2))
    I0 = Ipeak / k
    return BiGaussian(;z0, z1, I0, τ1, τ2)
end


"""
    Compute the total charge transferred by a `source`.
"""
function total_charge(source::BiGaussian)
    (;τ1, τ2, I0) = source

    return I0 * sqrt(π) * (τ1 - τ2) / 2
end


"""
    Current (cross-sectional, in A) at time `t`.
"""
function current(source::BiGaussian, t)
    (;τ1, τ2, I0) = source
    
    return t < 0 ? zero(I0) : I0 * (exp(-t^2/τ1^2) - exp(-t^2/τ2^2))
end


"""
    A representation of a bi-exponential, axial source.
"""
Base.@kwdef struct BiExponential{T}
    "Lowest altitude bound."
    z0::T

    "Upper altitude bound."
    z1::T

    "Total transferred charge."
    Q::T

    "Approximate duration."
    τ1::T

    "Approximate rise-time."
    τ2::T
end


"""
    Amplitude of the (cross-sectional) current, `I0` in Ampere.
"""
function current_ampl(source::BiExponential)
    (;τ1, τ2, Q) = source

    return Q / (τ1 - τ2)
end


"""
    Current (cross-sectional, in A) at time `t`.
"""
function current(source::BiExponential, t)
    (;τ1, τ2) = source

    I0 = current_ampl(source)

    return t < 0 ? zero(I0) : I0 * (exp(-t/τ1) - exp(-t/τ2))
end


"""
    A representation of an axial source with soft-step polynomial functions
"""
Base.@kwdef struct SoftStep{T}
    "Lowest altitude bound."
    z0::T

    "Upper altitude bound."
    z1::T

    "Peak current"
    Ipeak::T

    "Rise time"
    rise::T

    "decay"
    decay::T

    "plateau"
    plateau::T = zero(decay)
end

softstep(x) = x^2 * (2 - x)^2

"""
    Current (cross-sectional, in A) at time `t`.
"""
function current(source::SoftStep{T}, t) where T
    (;rise, decay, plateau, Ipeak) = source

    if t < 0
        return zero(Ipeak)
    elseif t < rise
        return Ipeak * softstep(t / rise)
    elseif t < rise + plateau
        return Ipeak
    elseif t < rise + decay + plateau        
        return Ipeak * softstep(1 - (t - rise - plateau) / decay)
    else
        return zero(Ipeak)
    end
end


"""
    Compute the total charge transferred by a `source`.
"""
function total_charge(source::SoftStep)
    (;rise, decay, plateau, Ipeak) = source

    return Ipeak * plateau + 8 * Ipeak * (rise + decay) / 5
end


"""
   Add an on-axis source to jz
"""
function add_source(device, mesh, fields, source, t)
    (;m, n, l, dr, dz) = mesh
    (;jz) = fields
    (;z0, z1) = source
    Ic = current(source, t)
    
    @batch for j in 1:n + 1
        i1 = l + 1
        j1 = j + l
        
        z = (j - 1) * dz
        
        if z0 <= z <= z1
            # The - sign is bc we use the geophysics convention.
            # a positive CG transports positive charge to the ground.
            jz[i1, j1] -= 4Ic / π / dr^2
        end
    end
end

