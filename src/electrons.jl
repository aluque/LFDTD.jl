# Handling of electron density profiles.

"""
    An electron density defined in terms of data linearly interpolated in log
    scale.
"""
struct LogInterpolatedElectronDensity{I}
    interp::I

    function LogInterpolatedElectronDensity(z, ne)
        interp = LinearInterpolation(z, log.(ne), extrapolation_bc=Line())
        new{typeof(interp)}(interp)
    end
end

"""
    Load data from file `fname` and construct a `LogInterpolatedElectronDensity`.
    Optionally, `z_scale` is a multiplicative factor for altitudes and
    `ne_scale` is a multiplicative factor for electron density
    (the resulting values are assumed to be in SI units).
"""
function LogInterpolatedElectronDensity(fname::String; z_scale=co.kilo, ne_scale=co.centi^-3)
    f = CSV.File(fname, header=[:z, :ne])
    @info "Electron density loaded from $fname"
    return LogInterpolatedElectronDensity(f.z * z_scale, f.ne * ne_scale)
end

""" Evaluate electron density at a given altitude. """
(li::LogInterpolatedElectronDensity)(z) = exp(li.interp(z))

"""
    An electron density defined in terms of the Wait-Spies parameterization
    (with the notation of Shao 2013)
"""
struct WaitSpiesElectronDensity{T <: Real}
    q::T
    h::T
end

const WAITSPIES_N0 = 3e7

""" Evaluate electron density at a given altitude. """
(ws::WaitSpiesElectronDensity)(z) = WAITSPIES_N0 * exp(ws.q * (z - ws.h))


"""
    Create a Wait-Spies electron density using the electron density values (`n1`, `n2`)
    at two altitudes (`z1`, `z2`)
"""
function WaitSpiesElectronDensity(z1, n1, z2, n2)
    n0 = WAITSPIES_N0
    q = log(n2 / n1) / (z2 - z1)
    h = (log(n0 / n1) + log(n0 / n2)) / 2q + (z1 + z2) / 2

    return WaitSpiesElectronDensity(q, h)
end


    
