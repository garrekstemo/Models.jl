module Models

export exponential, 
       double_exponential, 
       gaussian, 
       gaussian2d,
       lorentzian, 
       double_lorentzian, 
       sine, 
       damped_sine,
       squared_errors,
       dielectric_real, 
       dielectric_imag, 
       cavity_mode_energy, 
       coupled_energies,
       cavity_transmittance, 
       fsr

include("model_functions.jl")
include("special_functions.jl")

end # module
