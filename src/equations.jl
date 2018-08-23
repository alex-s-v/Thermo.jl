
using Polynomials: Poly, roots

"""
The Universal gas constant.
"""
const R = 8.3144598484848485

"""
The phase which determines a computational algorithm.
"""
@enum Phase LIQUID VAPOR

"""
    EquilibriumData(T::Float64, P::Float64, x::Float64, y::Float64, Vv::Float64,
                    Vl::Float64, ϕl::Array{Float64, 1}, ϕv::Array{Float64, 1})

Contains the data about an equilibrium point.
"""
type EquilibriumData
    T::Float64
    P::Float64
    x::Array{Float64, 1}
    y::Array{Float64, 1}
    ϕv::Array{Float64, 1}
    Vv::Float64
    ϕl::Array{Float64, 1}
    Vl::Float64
end

"""
    Mixture(name::AbstractString, x::Array{Float64, 1},
            y::Array{Float64, 1}, Tb::Array{Float64, 1},
            Tcr::Array{Float64, 1}, Pcr::Array{Float64, 1},
            ω::Array{Float64, 1}, k::Array{Float64, 2} = [])

Construct an instance of a mixture from the data about its individual
components (Where `ω` is a Pitzer acentric factors and `k` is a Binary interaction
parameters).

Mole fractions `x` and `y` must not include the last components (For example,
if a number of components is 4 then a length of `x` and `y` must be equals 3).
Sum of `x` and sum of `y` must be less or equal to 1.
"""
type Mixture
    name::AbstractString
    x::Array{Float64, 1}
    y::Array{Float64, 1}
    Tb::Array{Float64, 1}
    Tcr::Array{Float64, 1}
    Pcr::Array{Float64, 1}
    ω::Array{Float64, 1}
    k::Array{Float64, 2}
    RTcr::Array{Float64, 1}
    RTcr²::Array{Float64, 1}
    RTcr²⁵::Array{Float64, 1}
    function Mixture(name::AbstractString, x::Array{Float64, 1},
                     y::Array{Float64, 1}, Tb::Array{Float64, 1},
                     Tcr::Array{Float64, 1}, Pcr::Array{Float64, 1},
                     ω::Array{Float64, 1}, k = nothing)
        append!(x, 1 - sum(x))
        append!(y, 1 - sum(y))
        if !(length(x) == length(y) == length(Tb) ==
             length(Tcr) == length(Pcr) == length(ω))
            error("Length of all vectors must be the same.")
        end
        RTcr = R * Tcr
        RTcr² = RTcr.^2
        RTcr²⁵ = R^2 * Tcr.^2.5
        if k == nothing
            l = length(x)
            k = zeros((l, l))
        end
        return new(name, x, y, copy(Tb), copy(Tcr), copy(Pcr), copy(ω), copy(k),
                   RTcr, RTcr², RTcr²⁵)
    end
end

################################################################################
## Equations of state
################################################################################


abstract type Equation end

"""
    SRK(mix::Mixture)

Contains the data that necessary for solving
Soave-Redlich-Kwong equation of state.
"""
type SRK <: Equation
    mix::Mixture
    Ωa::Float64
    Ωb::Float64
    m̄::Array{Float64, 1}
    function SRK(mix::Mixture)
        m̄ = 0.48508 + 1.55171 * mix.ω - 0.15613 * mix.ω.^2
        return new(mix, 0.42747, 0.08664, m̄)
    end
end

"""
    PR(mix::Mixture)

Contains the data that necessary for solving
Peng-Robinson equation of state.
"""
type PR <: Equation
    mix::Mixture
    Ωa::Float64
    Ωb::Float64
    m̄::Array{Float64, 1}
    function PR(mix::Mixture)
        m̄ = 0.37464 + 1.54226 * mix.ω - 0.26992 * mix.ω.^2
        return new(mix, 0.457235, 0.077796, m̄)
    end
end

"""
    RK(mix::Mixture)

Contains the data that necessary for solving
Redlich-Kwong equation of state.
"""
type RK <: Equation
    mix::Mixture
    Ωa::Float64
    Ωb::Float64
    function RK(mix::Mixture)
        return new(mix, 0.42747, 0.08664)
    end
end

"""
    VdW(mix::Mixture)

Contains the data that necessary for solving
Van-der-Waals equation of state.
"""
type VdW <: Equation
    mix::Mixture
    Ωa::Float64
    Ωb::Float64
    function VdW(mix::Mixture)
        return new(mix, 0.421875, 0.125000)
    end
end

################################################################################
## Generic Functions
################################################################################


"""
    calc_compres_factor(coefs::Array{Float64, 1})

Solves the equation of state in a form of a cubic polynomial and return an
array of 3 roots if all of them is real or an array of 1 real root otherwise.
"""
function calc_compres_factor(coefs::Array{Float64, 1})

    z = roots(Poly(coefs))
    if any(!isreal, z)
        return nothing
    end
    return z
end

################################################################################
## Specific Functions
################################################################################


"""
    calc_params(eq::Equation, T::Float64, P::Float64, phase::Phase)

Computes the parameters for futher calculations of the fugacity coefficient.
"""
function calc_params end

function calc_params(eq::Union{SRK,PR}, T::Float64, P::Float64,
                     phase::Phase)

    frac = phase == LIQUID ? eq.mix.x : eq.mix.y

    ᾱ = (1 + eq.m̄ .* (1 - sqrt.(T ./ eq.mix.Tcr))).^2
    ā = (eq.Ωa * eq.mix.RTcr² ./ eq.mix.Pcr) .* ᾱ
    b̄ = eq.Ωb * eq.mix.RTcr ./ eq.mix.Pcr

    aa = (1 - eq.mix.k) .* sqrt.(ā*ā')
    top = 2 * aa * frac
    a = sum(aa .* (frac * frac'))
    b = sum(frac .* b̄)

    RT = R * T
    Ā = top ./ a
    B̄ = b̄ * P / RT

    A = a * P / RT^2
    B = b * P / RT

    return (A, B, Ā, B̄)
end

function calc_params(eq::RK, T::Float64, P::Float64, phase::Phase)

    frac = phase == LIQUID ? eq.mix.x : eq.mix.y

    ā = eq.Ωa * eq.mix.RTcr²⁵ ./ eq.mix.Pcr
    b̄ = eq.Ωb * eq.mix.RTcr ./ eq.mix.Pcr

    Ā = ā * P / (R^2 * T^2.5)
    B̄ = b̄ * P / (R * T)

    A = sum(frac .* sqrt.(Ā))^2
    B = sum(frac .* B̄)

    return (A, B, Ā, B̄)
end

function calc_params(eq::VdW, T::Float64, P::Float64, phase::Phase)

    frac = phase == LIQUID ? eq.mix.x : eq.mix.y

    Ā = eq.Ωa * eq.mix.RTcr² ./ eq.mix.Pcr
    B̄ = eq.Ωb * eq.mix.RTcr ./ eq.mix.Pcr

    A = sum(frac .* sqrt.(Ā))^2
    B = sum(frac .* B̄)

    return (A, B, Ā, B̄)
end

"""
    calc_fugacity_coef(eq::Equation, T::Float64, P::Float64, phase::Phase)

Computes the fugacity coefficients for a components in a given mixture by
a given equation in a given phase.

Returns tuple of the fugacity coefficients and molar volume of the phase.
"""
function calc_fugacity_coef(eq::Equation, T::Float64, P::Float64, phase::Phase)

    A, B, Ā, B̄ = calc_params(eq, T, P, phase)
    RT::Float64 = R*T
    eq_type = typeof(eq)

    coefs::Array{Float64, 1} = []
    if eq_type == SRK || eq_type == RK
        coefs = [-A * B, A - B - B^2, -1.0, 1.0]
    elseif eq_type == PR
        coefs = [-A * B + B^2 + B^3, A - 2 * B - 3 * B^2, -1.0 + B, 1.0]
    elseif eq_type == VdW
        coefs = [-A * B * P^2 / RT^3, A * P / RT^2, -(B * P / RT + 1), 1.0]
    end

    zs = calc_compres_factor(coefs)
    if zs == nothing
        return nothing
    end
    Z::Float64 = phase == LIQUID ? min(zs...) : max(zs...)
    V::Float64 = Z * RT / P

    if eq_type == SRK
        h1 = (Z - 1) * B̄ / B
        h2 = Z - B
        h3 = (A / B) * (B̄ / B - Ā)
        h4 = 1 + B / Z
    elseif eq_type == PR
        h1 = B̄ * (Z - 1) / B
        h2 = Z - B
        h3 = (A / (2.8284 * B)) * (B̄ / B - Ā)
        h4 = (Z + 2.4142 * B) / (Z - 0.4142 * B)
    elseif eq_type == RK
        h1 = (Z - 1) * B̄ / B
        h2 = Z - B
        h3 = (A / B) * (B̄ / B - 2 * sqrt.(Ā / A))
        h4 = 1 + B / Z
    elseif eq_type == VdW
        h1 = B̄ / (V - B)
        h2 = Z * (1 - B / V)
        h3 = -2 * sqrt.(Ā * A) / (RT * V)
        h4 = e
    end
    ϕ::Array{Float64, 1} = exp.(h1 - log(h2) + h3 * log(h4))

    return (ϕ, V)
end

"""
    solve(eq::Equation, T::Float64, P::Float64)

Solves the equation of state and returns thermodynamic data for given point.
"""
function solve(eq::Equation, T::Float64, P::Float64)

    v = calc_fugacity_coef(eq, T, P, VAPOR)
    l = calc_fugacity_coef(eq, T, P, LIQUID)

    if v == nothing || l == nothing
        return nothing
    end
    return EquilibriumData(T, P, copy(eq.mix.x), copy(eq.mix.y), v..., l...)
end
