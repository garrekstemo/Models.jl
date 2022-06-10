module Models

function exponential(x, f0, τ)
    return f0 * exp(- x / τ)
end

function double_exponential(x, f1, τ1, f2, τ2)
    return exponential(x, f1, τ1) + exponential(x, f2, τ2)
end

end # module
