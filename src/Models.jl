module Models

export exponential,
       gaussian, 
       gaussian2d,
       lorentzian, 
       sine, 
       damped_sine,
       pseudo_voigt,
       findpeaks,
       squared_errors,
       dielectric_real, 
       dielectric_imag, 
       cavity_mode_energy, 
       coupled_energies,
       cavity_transmittance, 
       fsr,
       transmission,

include("model_functions.jl")
include("special_functions.jl")

end # module
