"""
    exponential(x, p = [f_0, τ])

Exponential decay function with amplitude ``f_0``, and decay constant τ.

```math
\\begin{aligned}
    f(x; f_0, \\tau) = f_0 e^{-x / \\tau}
\\end{aligned}
```

[https://en.wikipedia.org/wiki/Exponential_decay](https://en.wikipedia.org/wiki/Exponential_decay)
"""
function exponential(x, p = [1.0, 1.0])
    f_0, τ = p
    return @. f_0 * exp(- x / τ)
end

"""
    double_exponential(x, p = [f_1, τ_1, f_2, τ_2])

Exponential decay function with two stages represented by
two different time constants, ``\\tau_1`` and ``\\tau_2``.

```math
\\begin{aligned}
    f(x; f_1, \\tau_1, f_2, \\tau_2) = f_1 e^{-x / \\tau_1} + f_2 e^{-x / \\tau_2}
\\end{aligned}
```

"""
function double_exponential(x, p = [1.0, 1.0, 1.0, 1.0])
    f_1, τ_1, f_2, τ_2 = p
    return @. f_1 * exp(- x / τ_1) + f_2 * exp(- x / τ_2)
end

"""
    gaussian(x, p = [A, μ, σ])

Gaussian function with amplitude ``A``, center ``μ``, and width ``σ``.

```math
\\begin{aligned}
    f(x; A, \\mu, \\sigma) = \\frac{A}{\\sigma \\sqrt{2\\pi}} e^{-(x - \\mu)^2 / (2 \\sigma^2)}
\\end{aligned}
```

[https://en.wikipedia.org/wiki/Gaussian_function](https://en.wikipedia.org/wiki/Gaussian_function)
"""
function gaussian(x, p = [1.0, 0.0, 0.1])

    A, μ, σ = p
    return @. A * exp( -(x - μ)^2 / (2 * σ^2) )
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
    return @. A / π * γ / ( (ν - ν_0)^2 + γ^2 )
end

"""
    double_lorentzian(ν, p = [A_1, σ_1, ν_1, A_2, σ_2, ν_2])

Sum of two Lorentzian functions.
"""
function double_lorentzian(ν, p = [1.0, -0.5, 0.1, 1.0, 0.5, 0.1])

    A_1, ν_1, σ_1, A_2, ν_2, σ_2 = p
    return @. A_1 / π * σ_1 / ( (ν - ν_1)^2 + σ_1^2 ) + A_2 / π * σ_2 / ( (ν - ν_2)^2 + σ_2^2 )
end
    
"""
    sine(x, p = [A, ω, ϕ])

Sinusoidal function.
"""
function sine(x, p = [1.0, 10.0, 0.0])
    A, ω, ϕ = p
    return @. A * sin(x * ω + ϕ)
end

"""
    squared_errors(p, f, X, Y)

Takes a function, `f`, and its parameters, `p`
and sums the squared errors given x-data and y-data
`X` and `Y`, respectively.
"""
function squared_errors(p, f, X, Y)
    error = 0.0

    for i in eachindex(X)
        y_i = f(X[i], p)
        error += (Y[i] - y_i)^2
    end
    return error
end