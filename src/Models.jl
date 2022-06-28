module Models

export exponential, double_exponential, gaussian, gaussian2d,
       lorentzian, double_lorentzian, sine, squared_errors

export dielectric_real, dielectric_imag, cavity_mode_energy, coupled_energies,
       cavity_transmittance, fsr

include("functions.jl")
include("special_functions.jl")

end # module
