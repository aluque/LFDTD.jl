using LFDTD: light_curve, tobs, co
using JLD2
using Plots
plotly()

include("event_list.jl")

const decay = 56.3 * co.micro
const in_band = 0.1
const hν_lbh = co.h * co.c / (200 * co.nano)


plot(frame=:box,
     xlabel="time (μs)",
     ylabel="Irradiance (μW/m²)")

for event in EVENT_LIST
    fname = joinpath("/Volumes/data8tb/lfdtd", event.name, "lbh.jld")
    obs = load(fname, "observer")
    plot!(tobs(obs) / co.micro, in_band .* hν_lbh .* light_curve(obs, decay) / co.micro,
          label=event.name)
end

gui()
