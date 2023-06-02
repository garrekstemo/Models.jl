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
       pp_spectrum

include("model_functions.jl")
include("special_functions.jl")
include("input-output.jl")

end # module
