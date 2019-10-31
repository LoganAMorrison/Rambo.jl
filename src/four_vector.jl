mutable struct FourMomentum{T<:Real}
    e::T
    x::T
    y::T
    z::T

    function FourMomentum(e::T, x::T, y::T, z::T) where T <: Real
        new{T}(e, x, y, z)
    end
end

"""
    scalar_product(fm1, fm2)

Compute the scalar product between to four-momenta.

# Arguments
- `fm1::FourMomentum`: first four-momentum
- `fm1::FourMomentum`: second four-momentum
"""
function scalar_product(
    fm1::FourMomentum{T},
    fm2::FourMomentum{T}
) where T <: Real
    fm1.e * fm2.e - (fm1.x * fm2.x + fm1.y * fm2.y + fm1.z * fm2.z)
end

"""
    mass(fm)

Compute the mass of a four-momentum vector, i.e. the sqrt of the magnitude.

# Arguments
- `fm::FourMomentum`: four-momentum to compute mass of
"""
function mass(fm::FourMomentum{T}) where T <: Real
    sqrt(abs(scalar_product(fm, fm)))
end

function Base.:+(fm1::FourMomentum{T}, fm2::FourMomentum{T}) where T <: Real
    FourMomentum(fm1.e + fm2.e, fm1.x + fm2.x, fm1.y + fm2.y, fm1.z + fm2.z)
end

function Base.:-(fm1::FourMomentum{T}, fm2::FourMomentum{T}) where T <: Real
    FourMomentum(fm1.e - fm2.e, fm1.x - fm2.x, fm1.y - fm2.y, fm1.z - fm2.z)
end

function Base.:*(fm::FourMomentum{T}, x::T) where T <: Real
    FourMomentum(fm.e * x, fm.x * x, fm.y * x, fm.z * x)
end

function Base.:*(x::T, fm::FourMomentum{T}) where T <: Real
    fm * x
end
