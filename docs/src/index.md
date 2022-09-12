# Models.jl Documentation

This simple package contains functions and [spectral line shapes](https://en.wikipedia.org/wiki/Spectral_line_shape)
commonly used in mid-infrared spectroscopy.
The main purpose is code reusability for my own projects, but
on the off chance that it is useful to others, it is available to use freely.
These functions are best used together with [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl) or the global optimization package [Optimization.jl](https://optimization.sciml.ai/stable/) for
fitting. It would be nice if models and fitting procedures were bundled together in the future,
like the [lmfit](https://lmfit.github.io/lmfit-py/index.html) package available for Python (leveraging SciPy).

The `squared_errors()` function is useful for use with `Optim.jl` for
curve fitting. It can be implemented for a set of x- and y-data and 
a model, e.g. an exponential curve. (Makie.jl is not a dependency; we only use it to demonstrate the fit.)

```@setup 1
using Pkg
Pkg.add(["Optim", "CairoMakie"])
Pkg.add(url="https://github.com/garrekstemo/Models.jl")
```

```@example 1
using Optim
using Models
using CairoMakie

xdata = -30:0.4:1.0

A, τ = [0.1, 5.0]
ydata = exponential(xdata, [A, τ]) .+ 1.2 * randn(length(xdata))

p0 = [0.1, 1.0]
result = optimize(b -> squared_errors(b, exponential, xdata, ydata), p0)
params = Optim.minimizer(result)

f = Figure()
ax = Axis(f[1, 1])

lines!(ax, xdata, ydata, label = "data")
lines!(ax, xdata, exponential(xdata, params), label = "fit")
axislegend(ax)

f
```
