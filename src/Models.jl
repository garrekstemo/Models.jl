module Models

function exponential(Float64::x, Vector{Float64}::p)
    f0, τ, y0 = p
    return @. f0 * exp(- x / τ) + y0
end

function double_exponential(Float64::x, Vector{Float64}::p)
    f1, τ1, f2, τ2, y0 = p
    return @. f1 * exp(- x / τ1) + f2 * exp(- x / τ2) + y0
end

end # module
