const Z = @SMatrix [0 0]

# 1×4 matrix that produces the z-component of the E-field when applied to the coefficients.
Eₗ4coeff(ν::Integer, q::Number, r::Real) = hcat(Fₗ4coeff(ν, q, r), Z)

# Matrix that produces the z-component of the H-field when applied to the coefficients.
Hₗ4coeff(ν::Integer, q::Number, r::Real) = hcat(Z, Fₗ4coeff(ν, q, r))
function Fₗ4coeff(ν::Integer, q::Number, r::Real)
    z = q * r
    M = @SMatrix [Jfun(ν,z)  Yfun(ν,z)]

    return M
end

# 1×4 matrix that produces the φ-component of the E-field when applied to the coefficients.
function Eᵩ4coeff(ν::Integer, 
                  μ::Real,  # relative permeability
                  k::Real,  # vacuum wavenumber
                  β::Real,  # propagation constant
                  q²::Real, q::Number, r::Real)
    M₁ = (β/r) * (im*ν) * Fₗ4coeff(ν,q,r)  # Fₗ4coeff = Eₗ4coeff[1:2]
    M₂ = -(k*μ*η₀) * deriv(r->Fₗ4coeff(ν,q,r), r)  # Fₗ4coeff = Hₗ4coeff[3:4]
    M = -(im/q²) * hcat(M₁, M₂)

    return M
end

# 1×4 matrix that produces the φ-component of the H-field when applied to the coefficients.
function Hᵩ4coeff(ν::Integer, 
                  ε::Real,  # relative permittivity (a.k.a. dielectric constant)
                  k::Real,  # vacuum wavenumber
                  β::Real,  # propagation constant
                  q²::Real, q::Number, r::Real)
    M₁ = (k*ε/η₀) * deriv(r->Fₗ4coeff(ν,q,r), r)  # Fₗ4coeff = Eₗ4coeff[1:2]
    M₂ = (β/r) * (im*ν) * Fₗ4coeff(ν,q,r)  # Fₗ4coeff = Hₗ4coeff[3:4]
    M = -(im/q²) * hcat(M₁, M₂)

    return M
end

# 1×4 matrix that produces the r-component of the E-field when applied to the coefficients.
function Eᵣ4coeff(ν::Integer, 
                  μ::Real,  # relative permeability
                  k::Real,  # vacuum wavenumber
                  β::Real,  # propagation constant
                  q²::Real, q::Number, r::Real)
    M₁ = β * deriv(r->Fₗ4coeff(ν,q,r), r)  # Fₗ4coeff = Eₗ4coeff[1:2]
    M₂ = (k*μ*η₀ / r) * (im*ν) * Fₗ4coeff(ν,q,r)  # Fₗ4coeff = Hₗ4coeff[3:4]
    M = -(im/q²) * hcat(M₁, M₂)

    return M
end

# 1×4 matrix that produces the r-component of the H-field when applied to the coefficients.
function Hᵣ4coeff(ν::Integer, 
                  ε::Real,  # relative permittivity (a.k.a. dielectric constant)
                  k::Real,  # vacuum wavenumber
                  β::Real,  # propagation constant
                  q²::Real, q::Number, r::Real)
    M₁ = -(k*ε/η₀) / r * (im*ν) * Fₗ4coeff(ν,q,r)  # Fₗ4coeff = Eₗ4coeff[1:2]
    M₂ = β * deriv(r->Fₗ4coeff(ν,q,r), r)  # Fₗ4coeff = Hₗ4coeff[3:4]
    M = -(im/q²) * hcat(M₁, M₂)

    return M
end
