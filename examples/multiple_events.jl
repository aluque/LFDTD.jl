using LFDTD

using LFDTD: main, DATA_DIR, Observer, LogInterpolatedElectronDensity, SoftStep, co

include("event_list.jl")


for event in EVENT_LIST    
    LFDTD.main(;
               # Dimensions of the domain
               R=800 * co.kilo,
               H=90 * co.kilo,
               
               cells_per_km=5,
               
               # Time step
               dt=2e-7,
               
               # Altitude of the observer (r=L, z=K)
               K=400 * co.kilo,
               L=event.L,
               
               # Field of view
               #fov = deg2rad(40),
               fov=deg2rad(40),
               
               source=SoftStep(Ipeak=event.Ipeak,
                               z0=10.0 * co.kilo,
                               z1=10.5 * co.kilo,
                               rise=event.rise,
                               decay=event.decay),
               
               # Filenames atmospheric profiles
               gas_density_fname=joinpath(DATA_DIR, "earth", "gas.dat"),
               electron_density=LogInterpolatedElectronDensity(joinpath(DATA_DIR,
                                                                        "earth", "electrons.dat")),
               output_dt=1e-4,

               # Change to corresponding folder
               output_folder=joinpath("/Volumes/data8tb/lfdtd", event.name)
               )
end
    
