"""
    dielectric_real(ν, p)

`p = [ν_0, A, Γ, n_eff]`

Frequency-dependent real part of the dielectric function in terms
of the background index and Lorentzian oscillator,

```math
\\begin{aligned}
\\varepsilon_1(\\nu) = n^2 + \\frac{A (\\nu_0^2 - \\nu^2)}{(\\nu^2 - \\nu_0^2)^2 + (\\Gamma\\nu)^2},
\\end{aligned}
```

where ``n`` is the background real index, ``\\nu_0`` is the excitation
frequency, ``\\Gamma`` is the line width of the oscillator, and ``A`` is the 
oscillator strength.
"""
function dielectric_real(ν, p)
    ν_0, A, Γ, n_eff = p
    return n_eff^2 + A * (ν_0^2 - ν^2) / ((ν^2 - ν_0^2)^2 + Γ^2 * ν^2)
end

"""
    dielectric_imag(ν, p)

`p = [ν_0, A, Γ, y_0]`

Frequency-dependent imaginary part of the dielectric function.

```math
\\begin{aligned}
\\varepsilon_2(\\nu) = \\frac{A \\Gamma \\nu}{(\\nu^2 - \\nu_0^2)^2 + (\\Gamma\\nu)^2}
\\end{aligned}
```
"""
function dielectric_imag(ν, p)
    ν_0, A, Γ, y_0 = p
    return A * Γ * ν / ((ν^2 - ν_0^2)^2 + Γ^2 * ν^2) + y_0
end

"""
    cavity_mode_energy(θ::Real, p)

`p = [E_0, n_eff]`

Cavity mode energy as a function of incident angle.

```math
\\begin{aligned}
    E_\\text{cavity}(\\theta) = E_0 \\left( 1 - \\frac{\\sin^2(\\theta)}{n^2} \\right)^{-1/2},
\\end{aligned}
```

where ``E_0`` is the energy of the cavity mode at zero degrees incidence angle.

"""
function cavity_mode_energy(θ::Real, p)
    E_0, n_eff = p
    return E_0 / sqrt(1 - (sin(θ) / n_eff)^2)
end

"""
    coupled_energies(E_c::Vector{Real}, E_v::Real, V::Real, branch=0)

Coupled energies found by diagonalizing the coupled-oscillator model Hamiltonian
with principle energies `E_c` (cavity photon) and `E_v` (material excitation),
and interaction energy `V`. The lower and upper branches are called with `branch=0`
and `branch=1`, respectively.

```math
\\begin{aligned}
H = 
\\begin{pmatrix}
    E_1(\\theta) & \\Omega_R/2 \\\\\
    \\Omega_R/2 & E_2
\\end{pmatrix}
\\end{aligned}
```

Diagonalizing ``H`` gives the interaction energies as a function
of ``\\theta``. These exhibit avoided crossing behavior near 
where the cavity mode energy equals the material excitation energy.

```math
\\begin{aligned}
    E_{\\pm}(\\theta) = \\frac{1}{2}(E_1(\\theta) + E_2) \\pm \\frac{1}{2} \\sqrt{(E_1(\\theta) - E_2)^2 + \\Omega_R^2}
\\end{aligned}
```

[https://en.wikipedia.org/wiki/Avoided_crossing](https://en.wikipedia.org/wiki/Avoided_crossing)
"""
function coupled_energies(E_c::Vector{T1}, E_v::T2, V::T3, branch=0) where {T1, T2, T3 <: Real}

    if branch == 0
        return 0.5 * ( (E_v + E_c) - sqrt(V^2 + (E_v - E_c)^2) )
    elseif branch == 1
        return 0.5 * ( (E_v + E_c) + sqrt(V^2 + (E_v - E_c)^2) )
    end
end


"""
    cavity_transmittance(ν, p)
    
    p = [n, α, L, T, R, ϕ]

Cavity transmittance as a function of frequency.

```math
\\begin{aligned}
    T(\\nu) = \\frac{T^2 e^{-\\alpha L}}{1 + R^2 e^{-2 \\alpha L} - 2 R e^{-\\alpha L} \\cos(4\\pi n L \\nu + 2\\phi)},
\\end{aligned}
```

where ``\\nu`` is the frequency, ``\\alpha`` and ``n`` are the frequency-dependent absorption and refractive index,
``L`` is the cavity length, and ``\\phi`` is the phase accumulated upon reflection against a mirror.
``T`` and ``R`` are the cavity transmittance and reflectance, respectively, where it is assumed that
``T = 1 - R``.
"""
function cavity_transmittance(ν, p)
    n, α, L, T, R, ϕ = p
    return T^2 * exp(-α * L) / (1 + R^2 * exp(-2 * α * L) - 2 * R * exp(-α * L) * cos(4*π * n * L * ν + 2*ϕ))
end

"""
    fsr(peak_positions::Vector{Real)
    
Find the free spectral range for all adjacent peaks.
Return the list of FSRs, the average FSR for the range, and the
standard deviation.
"""
function fsr(peak_positions::Vector{Real})
	fsrs = Float64[]
	i = 1
	while i < length(peak_positions)
		fsr_i = peak_positions[i + 1] - peak_positions[i]
		push!(fsrs, fsr_i)
		i += 1
	end
	mean(fsrs), std(fsrs), fsrs
end

"""
    fsr(x1::Real, x2::Real)

Solve for one of three variables in the equation
    for free spectral range in terms of the other two.

```math
\\begin{aligned}
    \\Delta \\nu = \\frac{1}{2 L n},
\\end{aligned}
```

where ``\\Delta\\nu = |\\nu_2 - \\nu_1|``, ``L`` is the intracavity length,
and ``n`` is intracavity index of refraction.
"""
function fsr(x1::T1, x2::T2) where {T1, T2 <: Real}
    return 1 / (2 * x1 * x2)
end
