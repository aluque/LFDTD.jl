# Handling of electron density profiles.

"""
    An electron density defined in terms of interpolated data.
"""
struct LogInterpolatedElectronDensity{I}
    interp::I

    function LogInterpolatedElectronDensity(z, ne)
        interp = LinearInterpolation(z, log.(ne))
        new{typeof(interp)}(interp)
    end
end

"""
    Load data from a file and construct a `LogInterpolatedElectronDensity`.
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

WAITSPIES_N0 = 3e7

""" Evaluate electron density at a given altitude. """
(ws::WaitSpiesElectronDensity)(z) = WAITSPIES_N0 * exp(ws.q * (z - ws.h))


    
