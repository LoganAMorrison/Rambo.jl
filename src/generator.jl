include("four_vector.jl")
using SpecialFunctions

mutable struct PhaseSpacePoint
    four_momenta::Array{FourMomentum{Float64}, 1}
    weight::Float64
end

"""
    PhaseSpacePoint(nfsp)

Basis constuctor of a phase-space-point.

# Arguments
- `nfsp::Int64`: number of final-state particles
"""
function PhaseSpacePoint(nfsp::Int64)
    four_momenta::Array{FourMomentum{Float64}, 1} = [
        FourMomentum(0.0, 0.0, 0.0, 0.0) for _ in 1:nfsp]
    PhaseSpacePoint(four_momenta, 1.0)
end


"""
    get_scale_factor(ps_point, cme, masses)

Function for finding the scaling parameter to turn massless four-vectors
into four-vectors with the correct masses.

# Arguments
- `ps_point::PhaseSpacePoint`: event four_momenta and weight
- `cme::Float64`: center-of-mass energy
- `masses::Array{Float64, 1}`: masses of final-state particles
"""
function get_scale_factor(ps_point::PhaseSpacePoint, cme::Float64, masses::Array{Float64,1})

    # Initial guess
    ξ0::Float64 = sqrt(1.0 - (sum(masses) / cme)^2)

    is_done::Bool = false
    iter_count::Int64 = 0
    ξ::Float64 = ξ0

    tol::Float64=1e-4
    max_iter::Int64=50
    while is_done == false
        # Compute residual and derivative of residual
        f::Float64, ∂f::Float64 = -cme, 0.0
        for (m, p) in zip(masses, ps_point.four_momenta)
            Δf::Float64 = sqrt(m^2 + ξ^2 * p.e^2)
            f += Δf
            ∂f += ξ * p.e^2 / Δf
        end

        # Newton correction
        Δξ = - (f / ∂f)
        ξ += Δξ

        iter_count += iter_count
        if abs(Δξ) < tol || iter_count > max_iter
            is_done = true
        end
    end
    ξ
end

"""
    initialize_four_momenta!(ps_point)

Initialize the four-momenta with isotropic, random four-momenta with energies,
q₀, distributed according to q₀ * exp(-q₀).

# Arguments
- `ps_point::PhaseSpacePoint`: event four_momenta and weight
"""
function initialize_four_momenta!(ps_point::PhaseSpacePoint)
    for p in ps_point.four_momenta
        ρ1::Float64 = rand(Float64)
        ρ2::Float64 = rand(Float64)
        ρ3::Float64 = rand(Float64)
        ρ4::Float64 = rand(Float64)

        c::Float64 = 2ρ1 - 1.0
        ϕ::Float64 = 2π * ρ2

        p.e = -log(ρ3 * ρ4)
        p.x = p.e * sqrt(1.0 - c^2) * cos(ϕ)
        p.y = p.e * sqrt(1.0 - c^2) * sin(ϕ)
        p.z = p.e * c
    end
end

"""
    boost_four_momenta!(ps_point, cme)

Boost the four-momenta into the center-of-mass frame and compute the inital
weight of the event.

# Arguments
- `ps_point::PhaseSpacePoint`: event four_momenta and weight
- `cme::Float64`: center-of-mass energy
"""
function boost_four_momenta!(ps_point::PhaseSpacePoint, cme::Float64)

    # Total momentum
    sum_qs::FourMomentum{Float64} = sum(ps_point.four_momenta)
    mass_Q::Float64 = mass(sum_qs)

    # Boost three-vector
    bx::Float64 = -sum_qs.x / mass_Q
    by::Float64 = -sum_qs.y / mass_Q
    bz::Float64 = -sum_qs.z / mass_Q
    # Boost factors
    x::Float64 = cme / mass_Q
    γ::Float64 = sum_qs.e / mass_Q
    a::Float64 = 1.0 / (1.0 + γ)

    for p in ps_point.four_momenta
        b_dot_q::Float64 = bx * p.x + by * p.y + bz * p.z

        pe::Float64 = x * (γ * p.e + b_dot_q)
        px::Float64 = x * (p.x + bx * p.e + a * b_dot_q * bx)
        py::Float64 = x * (p.y + by * p.e + a * b_dot_q * by)
        pz::Float64 = x * (p.z + bz * p.e + a * b_dot_q * bz)

        p.e = pe
        p.x = px
        p.y = py
        p.z = pz
    end

    n::Int64 = length(ps_point.four_momenta) # Number of final-state particles
    ps_point.weight = (π/2)^(n-1)*cme^(2n-4)/gamma(n)/gamma(n-1)*(2π)^(4-3n)
end

"""
    correct_masses!(ps_point, masses, cme)

Correct the masses of the four-momenta and correct the weight of the event.

# Arguments
- `ps_point::PhaseSpacePoint`: event four_momenta and weight
- `cme::Float64`: center-of-mass energy
- `masses::Array{Float64, 1}`: masses of final-state particles
"""
function correct_masses!(ps_point::PhaseSpacePoint, cme::Float64, masses::Array{Float64,1})

    term1::Float64 = 0.0
    term2::Float64 = 0.0
    term3::Float64 = 1.0

    ξ::Float64 = get_scale_factor(ps_point, cme, masses)

    for (m, p) in zip(masses, ps_point.four_momenta)
        p.e = sqrt(m^2 + ξ^2 * p.e^2)
        p.x = ξ * p.x
        p.y = ξ * p.y
        p.z = ξ * p.z

        modulus::Float64 = sqrt(p.x^2 + p.y^2 + p.z^2)

        term1 += modulus / cme
        term2 += modulus^2 / p.e
        term3 *= modulus / p.e
    end

    n::Int64 = length(ps_point.four_momenta)
    term1 = term1^(2n-3)
    term2 = 1.0 / term2

    ps_point.weight *= term1 * term2 * term3 * cme
end

"""
    generate_point(cme, masses)

Generate a single phase-space point.

# Arguments
- `cme::Float64`: center-of-mass energy.
- `masses::Array{Float64, 1}`: masses of the final-state particles.
- `msqrd::Function`: squared matrix element (default is flat)
"""
function generate_phase_space_point(cme::Float64, masses::Array{Float64, 1}; msqrd::Function=fms->1.0)
    point = PhaseSpacePoint(length(masses))

    initialize_four_momenta!(point)
    boost_four_momenta!(point, cme)
    correct_masses!(point, cme, masses)
    point.weight *= msqrd(point.four_momenta)
    point
end

"""
    generate_space(cme, masses, num_fsp)

Generate a single phase-space point.

# Arguments
- `cme::Float64`: center-of-mass energy
- `masses::Array{Float64, 1}`: masses of the final-state particles
- `nevents::Int64=10000`: number of events to generate
- `msqrd::Function=fms->1.0`: squared matrix element
"""
function generate_phase_space(cme::Float64, masses::Array{Float64, 1}; nevents::Int64=10000, msqrd::Function=fms->1.0)
    points::Array{PhaseSpacePoint, 1} = [PhaseSpacePoint(length(masses)) for _ in 1:nevents]

    for point in points
        initialize_four_momenta!(point)
        boost_four_momenta!(point, cme)
        correct_masses!(point, cme, masses)
        point.weight *= msqrd(point.four_momenta)
    end
    points
end
