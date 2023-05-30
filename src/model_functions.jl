"""
    exponential(x, p = [f_0, τ])

Exponential decay function with amplitude ``f_0``, and decay constant τ.

```math
\\begin{aligned}
    f(x; f_0, \\tau) = f_0 e^{-x / \\tau}
\\end{aligned}
```

[Exponential decay](https://en.wikipedia.org/wiki/Exponential_decay)
"""
function exponential(x, p = [1.0, 1.0])
    f_0, τ = p
    return f_0 * exp(-x / τ)
end

"""
    gaussian(x, p = [A, μ, σ])

Gaussian function with amplitude ``A``, center ``μ``, and width ``σ``.

```math
\\begin{aligned}
    f(x; A, \\mu, \\sigma) = \\frac{A}{\\sigma \\sqrt{2\\pi}} e^{-(x - \\mu)^2 / (2 \\sigma^2)}
\\end{aligned}
```

[Gaussian function](https://en.wikipedia.org/wiki/Gaussian_function)
"""
function gaussian(x, p = [1.0, 0.0, 0.1])
    A, μ, σ = p
    return A * exp( -(x - μ)^2 / (2 * σ^2) )
end

"""
    gaussian2d(x, y, p = [A, x_0, σ_x, y_0, σ_y])

Two-dimensional Gaussian function centered at ``(x_0, y_0)`` and x-width
``σ_x`` and y-width ``σ_y``, and amplitude ``A``.

```math
\\begin{aligned}
    f(x, y; A, x_0, \\sigma_x, y_0, \\sigma_y) = A \\exp\\left(-\\left( \\frac{(x - x_0)^2}{2 \\sigma_x^2} + \\frac{(y - y_0)^2}{2 \\sigma_y^2} \\right)\\right)
\\end{aligned}
```

[Two-dimensional Gaussian function](https://en.wikipedia.org/wiki/Gaussian_function#Two-dimensional_Gaussian_function)
"""
function gaussian2d(x, y, p = [1.0, 0.0, 1.0, 0.0, 1.0])
    G = Array{Float64}(undef, length(x), length(y))
    A, x_0, σ_x, y_0, σ_y = p

    for (i, x_i) in enumerate(x)
        for (j, y_j) in enumerate(y)
            G[i, j] = A / (σ_x * σ_y * 2π) * exp( -(x_i - x_0)^2 / (2 * σ_x^2) ) * exp( -(y_j - y_0)^2 / (2 * σ_y^2) )
        end
    end
    return G
end

"""
    lorentzian(ν, p = [A, ν_0, γ])

Lorentzian function with amplitude `A`, center frequency `ν_0`, and
full width at half maximum (FWHM) `2γ`.

```math
\\begin{aligned}
    f(\\nu; A, \\nu_0, \\gamma) = \\frac{A}{\\pi} \\frac{\\gamma}{(\\nu - \\nu_0)^2 + \\gamma^2}
\\end{aligned}
```

[https://en.wikipedia.org/wiki/Cauchy_distribution](https://en.wikipedia.org/wiki/Cauchy_distribution)
"""
function lorentzian(ν, p = [1.0, 0.0, 0.1])
    A, ν_0, γ  = p
    return A / π * γ / ( (ν - ν_0)^2 + γ^2 )
end
    
"""
    sine(x, p = [A, ω, ϕ])

Sinusoidal function.

```math
\\begin{aligned}
    \\sin(\\omega t + \\phi)
\\end{aligned}
```

[https://en.wikipedia.org/wiki/Sine_wave](https://en.wikipedia.org/wiki/Sine_wave)
"""
function sine(t, p = [1.0, 10.0, 0.0])
    A, ω, ϕ = p
    return A * sin(t * ω + ϕ)
end

"""
    damped_sine(x, p = [A, ω, ϕ, τ])

Damped sine function.

```math
\\begin{aligned}
    f(t) = A e^{-\\frac{t}{\\tau}} \\sin(\\omega t + \\phi)
\\end{aligned}
```

[https://en.wikipedia.org/wiki/Damping](https://en.wikipedia.org/wiki/Damping)
"""
function damped_sine(t, p = [1.0, 10.0, 0.0, 1.0])
    A, ω, ϕ, τ = p
    return A * exp(-t / τ) * sin(t * ω + ϕ)
end

"""
    pseudo_voigt(ω, p)

Weighted sum of a Lorentzian and a Gaussian function
with the same center `ω_0` and amplitude `f_0`.

p = [f_0, ω_0, σ, α]
"""
function pseudo_voigt(ω, p)
    f_0, ω_0, σ, α = p
    σ_g = σ / sqrt(2 * log(2))
    return (1 - α) * f_0 * exp(-(ω - ω_0)^2 / (2 * σ_g^2)) / (σ_g * sqrt(2 * π)) + α * f_0 * σ / ((ω - ω_0)^2 + σ^2) / π
            
end

"""
    localmax(y::AbstractVector; height = 0.0)

Find the local maxima of a vector and return their indices.
Optionally, a subset of peaks can be selected by specifying peak properties.
"""
function findpeaks(y::AbstractVector; height = 0.0)
    maxima = Real[]
    for i in 2:length(y) - 1
        if y[i] > y[i - 1] && y[i] > y[i + 1]
            if y[i] > height
                push!(maxima, i)
            end
        end
    end
    return maxima
end

"""
    squared_errors(p, f, X, Y)

Takes a function, `f`, and its parameters, `p`
and sums the squared errors given x-data and y-data
`X` and `Y`, respectively.

```math
\\begin{aligned}
    \\text{err} = \\sum_{i=1}^{n} \\left( Y_i(X) - \\hat{Y}_i \\right)^2
\\end{aligned}
```
"""
function squared_errors(p, f::F, X, Y) where F <: Function
    error = 0.0
    for i in eachindex(X)
        y_i = f(X[i], p)
        error += (Y[i] - y_i)^2
    end
    return error
end

