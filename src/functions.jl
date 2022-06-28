"""Squared errors

Takes a function, `f`, and its parameters, `p`
and sums the squared errors given x-data and y-data
`X` and `Y`, respectively.
"""
function squared_errors(p, f::Function, X, Y)
    error = 0.0

    for i in 1:length(X)
        y_i = f(X[i], p)
        error += (Y[i] - y_i)^2
    end
    return error
end

"""Exponential decay function.

``f(x; f_0, \\tau, y_0) = f_0 e^{-x / \\tau} + y_0``

[https://en.wikipedia.org/wiki/Exponential_decay](https://en.wikipedia.org/wiki/Exponential_decay)
"""
function exponential(x, p = [1.0, 1.0, 0.0])
    f_0, τ, y_0 = p
    return @. f_0 * exp(- x / τ) + y_0
end

"""
Exponential decay function with two stages represented by
two different time constants.

``f(x; f_1, \\tau_1, f_2, \\tau_2, y_0) = f_1 e^{-x / \\tau_1} + f_2 e^{-x / \\tau_2} + y_0``

"""
function double_exponential(x, p = [1.0, 1.0, 1.0, 1.0, 0.0])
    f_1, τ_1, f_2, τ_2, y_0 = p
    return @. f_1 * exp(- x / τ_1) + f_2 * exp(- x / τ_2) + y_0
end

"""
Gaussian function.
``f(x; A, \\mu, \\sigma, y_0) = \\frac{A}{\\sigma \\sqrt{2\\pi}} e^{-(x - \\mu)^2 / (2 \\sigma^2)} + y_0``

[https://en.wikipedia.org/wiki/Gaussian_function](https://en.wikipedia.org/wiki/Gaussian_function)
"""
function gaussian(x, p = [1.0, 0.0, 1.0, 0.0])

    A, μ, σ, y_0 = p
    return @. A * exp( -(x - μ)^2 / (2 * σ^2) ) + y_0
end

"""
Two-dimensional Gaussian lineshape.

``f(x, y; A, x_0, \\sigma_x, y_0, \\sigma_y, z_0) = A \\exp\\left(-\\left( \\frac{(x - x_0)^2}{2 \\sigma_x^2} + \\frac{(y - y_0)^2}{2 \\sigma_y^2} \\right)\\right) + z_0``
"""
function gaussian2d(x, y, p = [1.0, 0.0, 1.0, 0.0, 1.0, 0.0])

    G = Array{Float64}(undef, length(x), length(y))
    A, x_0, σ_x, y_0, σ_y, z_0 = p

    for (i, x_i) in enumerate(x)
        for (j, y_j) in enumerate(y)
            G[i, j] = A / (σ_x * σ_y * 2π) * exp( -(x_i - x_0)^2 / (2 * σ_x^2) ) * exp( -(y_j - y_0)^2 / (2 * σ_y^2) ) + z_0
        end
    end
    return G
end

"""
Lorentzian function with amplitude ``A``, center frequency ``\\nu_0``, and full width at half maximum (FWHM) ``2\\gamma``.

``f(\\nu; A, \\nu_0, \\gamma, y_0) = \\frac{A}{\\pi} \\frac{\\gamma}{(\\nu - \\nu_0)^2 + \\gamma^2} + y_0``

[https://en.wikipedia.org/wiki/Cauchy_distribution](https://en.wikipedia.org/wiki/Cauchy_distribution)
"""
function lorentzian(ν, p = [1.0, 1.0, 0.0, 0.0])

    A, γ, ν_0, y_0 = p
    return @. A / π * γ / ( (ν - ν_0)^2 + γ^2 ) + y_0
end

"""
Sum of two Lorentzian functions.
"""
function double_lorentzian(ν, p = [1.0, 1.0, -1.0, 1.0, 1.0, 1.0, 0.0])

    A_1, σ_1, ν_1, A_2, σ_2, ν_2, y_0 = p
    return @. A_1 / π * σ_1 / ( (ν - ν_1)^2 + σ_1^2 ) + A_2 / π * σ_2 / ( (ν - ν_2)^2 + σ_2^2 ) + y_0
end
    
"""
Sinusoidal function.
"""
function sine(x, p = [1.0, 1.0, 0.0, 0.0])
    A, ω, ϕ, y_0 = p
    return @. A * sin(x * ω + ϕ) + y_0
end