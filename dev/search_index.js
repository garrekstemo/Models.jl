var documenterSearchIndex = {"docs":
[{"location":"#Models.jl-Documentation","page":"Models.jl Documentation","title":"Models.jl Documentation","text":"","category":"section"},{"location":"","page":"Models.jl Documentation","title":"Models.jl Documentation","text":"This simple package contains functions and lineshapes commonly used in vibrational polariton spectroscopic experiments. The main purpose is code reusability for my own projects, but on the off chance that it is useful to others, it is available. These functions are best used together with Optim.jl or the global optimization package Optimization.jl for fitting. It would be nice if models and fitting procedures were bundled together in the future, like the lmfit package available for Python (leveraging SciPy).","category":"page"},{"location":"","page":"Models.jl Documentation","title":"Models.jl Documentation","text":"The squared_errors() function is useful for use with Optim.jl for curve fitting. It can be implemented for a set of x- and y-data and  a model, e.g. an exponential curve. (Makie.jl is not a dependency; we only use it to demonstrate the fit.)","category":"page"},{"location":"","page":"Models.jl Documentation","title":"Models.jl Documentation","text":"using Optim\nusing Models\nusing GLMakie\n\nxdata = -30:0.4:1.0\n\nA, τ, y_0 = [0.1, 5.0, 0.3]\nydata = exponential(xdata, [A, τ, y_0]) .+ 1.2 * randn(length(xdata))\n\np0 = [0.1, 1.0, 0.0]\nresult = optimize(b -> squared_errors(b, exponential, xdata, ydata), p0)\nparams = Optim.minimizer(result)\n\nfig = Figure()\ndisplay(fig)\nax = Axis(fig[1, 1])\nlines!(ax, xdata, ydata, label = \"data\")\nlines!(ax, xdata, exponential(xdata, params), label = \"fit\")\naxislegend(ax)","category":"page"},{"location":"","page":"Models.jl Documentation","title":"Models.jl Documentation","text":"(Image: exponential fit image)","category":"page"},{"location":"#General-Functions","page":"Models.jl Documentation","title":"General Functions","text":"","category":"section"},{"location":"","page":"Models.jl Documentation","title":"Models.jl Documentation","text":"Common functions and lineshapes for spectroscopy.","category":"page"},{"location":"","page":"Models.jl Documentation","title":"Models.jl Documentation","text":"Modules = [Models]\nPages = [\"functions.jl\"]","category":"page"},{"location":"#Models.cavity_mode_energy-Tuple{Any, Any}","page":"Models.jl Documentation","title":"Models.cavity_mode_energy","text":"Cavity mode energy as a function of incident angle.\n\nE_textcavity(theta) = E_0 left( 1 - fracsin^2(theta)n^2 right)^-12\n\nwhere E_0 is the energy of the cavity mode at zero degrees incidence angle.\n\n\n\n\n\n","category":"method"},{"location":"#Models.cavity_transmittance-Tuple{Any, Any}","page":"Models.jl Documentation","title":"Models.cavity_transmittance","text":"Cavity transmittance as a function of frequency.\n\nT(nu) = fracT^2 e^-alpha L1 + R^2 e^-2 alpha L - 2 R e^-alpha L cos(4pi n L nu + 2phi)\n\nwhere nu is the frequency, alpha and n are the frequency-dependent absorption and refractive index, L is the cavity length, and phi is the phase accumulated upon reflection against a mirror. T and R are the cavity transmittance and reflectance, respectively, where it is assumed that T = 1 - R.\n\n\n\n\n\n","category":"method"},{"location":"#Models.coupled_energies","page":"Models.jl Documentation","title":"Models.coupled_energies","text":"Coupled energies found by diagonalizing the coupled-oscillator model Hamiltonian:\n\nH =  beginpmatrix E_1(theta)  Omega_R2 Omega_R2  E_2 endpmatrix\n\nDiagonalizing H gives the interaction energies as a function of theta. These exhibit avoided crossing behavior near  where the cavity mode energy equals the material excitation energy.\n\nE_pm(theta) = frac12(E_1(theta) + E_2) pm frac12 sqrt(E_1(theta) - E_2)^2 + Omega_R^2\n\nhttps://en.wikipedia.org/wiki/Avoided_crossing\n\n\n\n\n\n","category":"function"},{"location":"#Models.dielectric_imag-Tuple{Any, Any}","page":"Models.jl Documentation","title":"Models.dielectric_imag","text":"Frequency-dependent imaginary part of the dielectric function.\n\nvarepsilon_2(nu) = fracA Gamma nu(nu^2 - nu_0^2)^2 + (Gammanu)^2\n\n\n\n\n\n","category":"method"},{"location":"#Models.dielectric_real-Tuple{Any, Any}","page":"Models.jl Documentation","title":"Models.dielectric_real","text":"Frequency-dependent real part of the dielectric function in terms of the background index and Lorentzian oscillator,\n\nvarepsilon_1(nu) = n^2 + fracA (nu_0^2 - nu^2)(nu^2 - nu_0^2)^2 + (Gammanu)^2\n\nwhere n is the background real index, nu_0 is the excitation frequency, Gamma is the line width of the oscillator, and A is the  oscillator strength.\n\n\n\n\n\n","category":"method"},{"location":"#Models.double_exponential","page":"Models.jl Documentation","title":"Models.double_exponential","text":"Exponential decay function with two stages represented by two different time constants.\n\nf(x f_1 tau_1 f_2 tau_2 y_0) = f_1 e^-x  tau_1 + f_2 e^-x  tau_2 + y_0\n\n\n\n\n\n","category":"function"},{"location":"#Models.double_lorentzian","page":"Models.jl Documentation","title":"Models.double_lorentzian","text":"Sum of two Lorentzian functions.\n\n\n\n\n\n","category":"function"},{"location":"#Models.exponential","page":"Models.jl Documentation","title":"Models.exponential","text":"Exponential decay function\n\nf(x f_0 tau y_0) = f_0 e^-x  tau + y_0\n\nhttps://en.wikipedia.org/wiki/Exponential_decay\n\n\n\n\n\n","category":"function"},{"location":"#Models.fsr-Tuple{Float64, Float64}","page":"Models.jl Documentation","title":"Models.fsr","text":"Solve for one of three variables in the equation     for free spectral range in terms of the other two.\n\nDelta nu = frac12 L n\n\nwhere Deltanu = nu_2 - nu_2, L is the intracavity length, and n is intracavity index of refraction.\n\n\n\n\n\n","category":"method"},{"location":"#Models.fsr-Tuple{Vector{Float64}}","page":"Models.jl Documentation","title":"Models.fsr","text":"Find the free spectral range for all adjacent peaks. Return the list of FSRs, the average FSR for the range, and the standard deviation.\n\n\n\n\n\n","category":"method"},{"location":"#Models.gaussian","page":"Models.jl Documentation","title":"Models.gaussian","text":"Gaussian function\n\nf(x A mu sigma y_0) = fracAsigma sqrt2pi e^-(x - mu)^2  (2 sigma^2) + y_0\n\nhttps://en.wikipedia.org/wiki/Gaussian_function\n\n\n\n\n\n","category":"function"},{"location":"#Models.gaussian2d","page":"Models.jl Documentation","title":"Models.gaussian2d","text":"Two-dimensional Gaussian function\n\nf(x y A x_0 sigma_x y_0 sigma_y z_0) = A expleft(-left( frac(x - x_0)^22 sigma_x^2 + frac(y - y_0)^22 sigma_y^2 right)right) + z_0\n\n\n\n\n\n","category":"function"},{"location":"#Models.lorentzian","page":"Models.jl Documentation","title":"Models.lorentzian","text":"Lorentzian function with amplitude A, center frequency nu_0, and full width at half maximum (FWHM) 2gamma.\n\nf(nu A nu_0 gamma y_0) = fracApi fracgamma(nu - nu_0)^2 + gamma^2 + y_0\n\nhttps://en.wikipedia.org/wiki/Cauchy_distribution\n\n\n\n\n\n","category":"function"},{"location":"#Models.sine","page":"Models.jl Documentation","title":"Models.sine","text":"Sinusoidal function.\n\n\n\n\n\n","category":"function"},{"location":"#Models.squared_errors-NTuple{4, Any}","page":"Models.jl Documentation","title":"Models.squared_errors","text":"Squared errors\n\nTakes a function, f, and its parameters, p and sums the squared errors given x-data and y-data X and Y, respectively.\n\n\n\n\n\n","category":"method"},{"location":"#Special-Functions","page":"Models.jl Documentation","title":"Special Functions","text":"","category":"section"},{"location":"","page":"Models.jl Documentation","title":"Models.jl Documentation","text":"These function are specific to cavity polariton research, where things like the cavity transmittance, material dielectric function, and cavity free spectral range (FSR) need to be calculated and modeled.","category":"page"},{"location":"","page":"Models.jl Documentation","title":"Models.jl Documentation","text":"Modules = [Models]\nPages = [\"special_functions.jl\"]","category":"page"}]
}
