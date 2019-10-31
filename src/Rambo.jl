module Rambo

include("four_vector.jl")
include("generator.jl")

"""
    integrate_phase_space(cme, masses, <keyword arguments>)

Compute integral over phase-space and the estimated error.

# Arguments
- `cme::Float64`: center-of-mass energy
- `masses::Array{Float64, 1}`: masses of the final-state particles
- `nevents::Int64=10000`: number of events to generate
- `msqrd::Function=fms->1.0`: squared matrix element
"""
function integrate_phase_space(cme::Float64, masses::Array{Float64, 1};
                               nevents::Int64=10000, msqrd::Function=fms->1.0)
    points::Array{PhaseSpacePoint, 1} = generate_phase_space(cme, masses; nevents=nevents, msqrd=msqrd)
    weights::Array{Float64, 1} = [point.weight for point in points]

    avg::Float64 = sum(weights) / nevents
    var::Float64 = sum((weights .- avg) .^2) / (nevents - 1)

    (avg, sqrt(var / nevents))
end

"""
    cross_section(cme, isp_masses, fsp_masses[nevents=10000,msqrd=fms->1.0])

Compute the cross section for a center-of-mass energy `cme`, initial-state
particle masses `isp_masses`, final-state particle masses `fsp_masses`. Default
keyword arguments assume a flat matrix element.

# Arguments
- `cme::Float64`: center-of-mass energy
- `isp_masses::Array{Float64, 1}`: initial-state particle masses
- `fsp_masses::Array{Float64, 1}`: final-state particle masses
- `nevents::Int64=10000`: number of phase-space points to generate
- `msqrd::Function=fms->1.0`: function for squared matrix element
"""
function cross_section(cme::Float64, isp_masses::Array{Float64, 1},
                       fsp_masses::Array{Float64, 1};
                       nevents::Int64=10000, msqrd::Function=fms->1.0)
    ps_integral, ps_error = integrate_phase_space(cme, fsp_masses;
                                                  nevents=nevents, msqrd=msqrd)

    m1::Float64 = isp_masses[1]
    m2::Float64 = isp_masses[2]
    # Compute the energies of the isp in the CM frame
    eng1::Float64 = (cme^2 + m1^2 - m2^2) / (2cme)
    eng2::Float64 = (cme^2 + m2^2 - m1^2) / (2cme)

    # Compute the magnitude of the 3 momenta of the isp in CM frame
    p::Float64 = sqrt((m1 - m2 - cme) * (m1 + m2 - cme) *
                      (m1 - m2 + cme) * (m1 + m2 + cme)) / (2cme)

    # Absolute and relative velocities in CM frame
    v1::Float64 = p / eng1
    v2::Float64 = p / eng2
    vrel::Float64 = v1 + v2

    cs_prefactor::Float64 = 1.0 / (2eng1 * 2eng2 * vrel)

    (ps_integral * cs_prefactor, ps_error * cs_prefactor)
end

"""
    decay_width(cme, fsp_masses, <keyword arguments>)

Compute the decay width for a center-of-mass energy `cme`, final-state particle
masses `fsp_masses`. Default keyword arguments assume a flat matrix element.

# Arguments
- `cme::Float64`: center-of-mass energy
- `isp_masses::Array{Float64, 1}`: initial-state particle masses
- `nevents::Int64=10000`: number of phase-space points to generate
- `msqrd::Function=fms->1.0`: function for squared matrix element
"""
function decay_width(cme::Float64, fsp_masses::Array{Float64, 1};
                     nevents::Int64=10000, msqrd::Function=fms->1.0)
    ps_integral, ps_error = integrate_phase_space(cme, fsp_masses;
                                                  nevents=nevents, msqrd=msqrd)

    (ps_integral / (2cme), ps_error / (2cme))
end

export FourMomentum
export PhaseSpacePoint
export generate_phase_space_point
export generate_phase_space
export integrate_phase_space
export cross_section
export decay_width
export scalar_product



end # module
