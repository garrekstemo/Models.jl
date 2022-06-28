# Models.jl Documentation

This simple package contains functions and lineshapes
commonly used in vibrational polariton spectroscopic experiments.
The main purpose is code reusability for my own projects, but
on the off chance that it is useful to others, it is available.
These functions are best used together with [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl) or the global optimization package [Optimization.jl](https://optimization.sciml.ai/stable/) for
fitting. It would be nice if models and fitting procedures were bundled together in the future,
like the [lmfit](https://lmfit.github.io/lmfit-py/index.html) package available for Python (leveraging SciPy).

## General Functions

Common functions and lineshapes for spectroscopy.

```@docs
exponential(x)
```

```@docs
double_exponential(x)
```

```@docs
gaussian(x)
```

```@docs
gaussian2d(x, y)
```

```@docs
lorentzian(ν)
```

```@docs
squared_errors(p)
```

## Special Functions

These function are specific to cavity polariton research,
where things like the cavity transmittance, material dielectric function,
and cavity free spectral range (FSR) need to be calculated and modeled.

```@docs
cavity_transmittance(ν, p)
```

```@docs
dielectric_real(ν, p)
```

```@docs
dielectric_imag(ν, p)
```

```@docs
cavity_mode_energy(θ, p)
```

```@docs
coupled_energies(E_c, E_v, V, branch=0)
```

```@docs
fsr(peak_positions::Array{Float64, 1})
```

```@docs
fsr(x1::Float64, x2::Float64)
```