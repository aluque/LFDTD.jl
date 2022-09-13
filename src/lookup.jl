# Code for lookup tables with arbitrary transformations in X and Y.

struct LookupTable{TX, TY, F, G, GI}
    # X data (f(x)).  Much faster when this is a uniform range
    fx::TX

    # Y daya (g(y)).
    gy::TY

    # The transform function for x
    f::F

    # The transform function for y
    g::G

    # The inverse of g
    ginv::GI

    function LookupTable(fx, gy; f=identity, g=identity, ginv=identity)
        @assert length(fx) == length(gy)

        new{typeof(fx), typeof(gy), typeof(f), typeof(g), typeof(ginv)}(fx, gy, f, g, ginv)
    end        
end


@inline (tbl::LookupTable)(x) = lookup(tbl, x)


# Unfortunately closures are incompatible with Polyester's @batch in mac
# so we have to define this function instead of the callable.
"""
    Look up and interpolate the value of the table at a given `x`.
"""
@inline function lookup(tbl::LookupTable, x)
    fx = tbl.f(x)
    
    i = searchsortedlast(tbl.fx, fx)

    if i == 0
        # Allow constant-value extrapolation below the range
        return tbl.ginv(tbl.gy[1])
    end
    
    checkbounds(tbl.fx, i)
    checkbounds(tbl.fx, i + 1)
    
    @inbounds w = (fx - tbl.fx[i]) / (tbl.fx[i + 1] - tbl.fx[i])
    @inbounds gy = w * tbl.gy[i + 1] + (1 - w) * tbl.gy[i]

    return tbl.ginv(gy)
end


"""    
    Guess a range from a vector `v`.  Raises a warning if the values are too far 
    from uniform.
"""
function approxrange(v::AbstractVector; atol::Real=0, rtol::Real=0.01)
    @assert issorted(v)
    
    h = diff(v)
    hmean = sum(h) / length(h)
    
    if !all(x -> isapprox(x, hmean; atol, rtol), h)
        @warn "The provided vector is not close enough to uniform range: lookup will be much slower."

        # Note type instability
        return v
        # throw(ArgumentError("The provided vector is not close enough to uniform range"))
    end

    l = length(v)
    
    return LinRange(v[begin], v[end], l)
end


"""
    Load a lookuptable from a delimited file.
"""
function loadtable(fname; f=identity, g=identity, ginv=identity, kw...)
    data = CSV.File(fname, header=[:en, :k], comment="#")
    fx = approxrange(f.(data.en); kw...)
    gy = g.(data.k)

    LookupTable(fx, gy; f, g, ginv)
end
