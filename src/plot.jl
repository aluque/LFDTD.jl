const hν_lbh = co.h * co.c / (200 * co.nano)

function plot_light_curve(obs::Observer; hν=hν_lbh, in_band=0.1, normalize=false, kw...)
    # Liu and Pasko 2004
    decay = 56.3 * co.micro
    y = in_band .* hν .* light_curve(obs, decay) / co.micro
    maybenormalize!(y, normalize)
    
    plt.plot(tobs(obs) / co.micro, y; kw...)
    plt.xlabel("time (μs)")
    plt.ylabel("Irradiance " * (normalize ? "(a. u.)" : "(\$\\mathdefault{μW/m^2}\$)"))
end

function maybenormalize!(y, normalize)
    if normalize
        y ./= maximum(y)
    end

    return y
end
    
function plot_light_curve(fname::String; kw...)
    obs = load(fname, "observer")
    plot_light_curve(obs; kw...)
end
