# Models.jl Documentation

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