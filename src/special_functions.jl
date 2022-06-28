"""
Frequency-dependent real part of the dielectric function in terms
of the background index and Lorentzian oscillator,

``\\varepsilon_1(\\nu) = n^2 + \\frac{A (\\nu_0^2 - \\nu^2)}{(\\nu^2 - \\nu_0^2)^2 + (\\Gamma\\nu)^2},``

where ``n`` is the background real index, ``\\nu_0`` is the excitation
frequency, ``\\Gamma`` is the line width of the oscillator, and ``A`` is the 
oscillator strength.
"""
function dielectric_real(ν, p)
    ν_0, A, Γ, n_eff = p
    return @. n_eff^2 + A * (ν_0^2 - ν^2) / ((ν^2 - ν_0^2)^2 + Γ^2 * ν^2)
end

"""
Frequency-dependent imaginary part of the dielectric function.

``\\varepsilon_2(\\nu) = \\frac{A \\Gamma \\nu}{(\\nu^2 - \\nu_0^2)^2 + (\\Gamma\\nu)^2}``
"""
function dielectric_imag(ν, p)
    ν_0, A, Γ, y_0 = p
    return @. A * Γ * ν / ((ν^2 - ν_0^2)^2 + Γ^2 * ν^2) + y_0
end

"""
Cavity mode energy as a function of incident angle.

``E_\\text{cavity}(\\theta) = E_0 \\left( 1 - \\frac{\\sin^2(\\theta)}{n^2} \\right)^{-1/2},``

where ``E_0`` is the energy of the cavity mode at zero degrees incidence angle.

"""
function cavity_mode_energy(θ, p)
    E_0, n_eff = p
    @. E_0 / sqrt(1 - (sin(θ) / n_eff)^2)
end

"""
Coupled energies found by diagonalizing the coupled-oscillator model Hamiltonian:

``H = 
\\begin{pmatrix}
E_1(\\theta) & \\Omega_R/2 \\\\\
\\Omega_R/2 & E_2
\\end{pmatrix}.``

Diagonalizing ``H`` gives the interaction energies as a function
of ``\\theta``. These exhibit avoided crossing behavior near 
where the cavity mode energy equals the material excitation energy.

``
E_{\\pm}(\\theta) = \\frac{1}{2}(E_1(\\theta) + E_2) \\pm \\frac{1}{2} \\sqrt{(E_1(\\theta) - E_2)^2 + \\Omega_R^2}
``

[https://en.wikipedia.org/wiki/Avoided_crossing](https://en.wikipedia.org/wiki/Avoided_crossing)
"""
function coupled_energies(E_c, E_v, V, branch=0)

    if branch == 0
        return @. 0.5 * ( (E_v + E_c) - sqrt(V^2 + (E_v - E_c)^2) )
    elseif branch == 1
        return @. 0.5 * ( (E_v + E_c) + sqrt(V^2 + (E_v - E_c)^2) )
    end
end

"""
Cavity transmittance as a function of frequency.

``T(\\nu) = \\frac{T^2 e^{-\\alpha L}}{1 + R^2 e^{-2 \\alpha L} - 2 R e^{-\\alpha L} \\cos(4\\pi n L \\nu + 2\\phi)},``

where ``\\nu`` is the frequency, ``\\alpha`` and ``n`` are the frequency-dependent absorption and refractive index,
``L`` is the cavity length, and ``\\phi`` is the phase accumulated upon reflection against a mirror.
``T`` and ``R`` are the cavity transmittance and reflectance, respectively, where it is assumed that
``T = 1 - R``.
"""
function cavity_transmittance(ν, p)
    n, α, L, T, R, ϕ = p
    return @. T^2 * exp(-α * L) / (1 + R^2 * exp(-2 * α * L) - 2 * R * exp(-α * L) * cos(4*π * n * L * ν + 2*ϕ))
end

"""
Find the free spectral range for all adjacent peaks.
Return the list of FSRs, the average FSR for the range, and the
standard deviation.
"""
function fsr(peak_positions::Array{Float64, 1})
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
Solve for one of three variables in the equation
    for free spectral range in terms of the other two.

``\\Delta \\nu = \\frac{1}{2 L n},``

where ``\\Delta\\nu = |\\nu_2 - \\nu_2|``, ``L`` is the intracavity length,
and ``n`` is intracavity index of refraction.
"""
function fsr(x1::Float64, x2::Float64)
    return 1 / (2 * x1 * x2)
end
