export dispfun

# 4×4 matrix that produces the tangential components of the E- and H-fields (Eₗ, Hₗ, Eᵩ, Hᵩ)
# when applied to the coefficients
function Fₜ4coeff(ν::Integer, 
                  ε::Real,  # relative permittivity (a.k.a. dielectric constant)
                  μ::Real,  # relative permeability
                  k::Real,  # vacuum wavenumber
                  β::Real,  # propagation constant
                  r::Real)
    q² = k^2 * μ * ε - β^2
    q = √complex(q²)

    M = vcat(Eₗ4coeff(ν, q, r), 
             Hₗ4coeff(ν, q, r), 
             Eᵩ4coeff(ν, μ, k, β, q², q, r), 
             Hᵩ4coeff(ν, ε, k, β, q², q, r))

    return M
end

# Transfer matrix that transforms the coefficients of an inner region to the coefficients of
# the adjacent outer region.  The interface between the two regions is at radius r.
function Tcoeff(ν::Integer, 
                r::Real,  # radius of interface between two regions
                ε::SReal{2},  # relative permittivities of [outer, inner] regions
                μ::SReal{2},  # relative permeabilities of [outer, inner] regions
                k::Real,  # vacuum wavenumber
                β::Real)  # propagation constant
    εₒ, εᵢ = ε
    μₒ, μᵢ = μ

    Tᵢ = Fₜ4coeff(ν, εᵢ, μᵢ, k, β, r)
    Tₒ = Fₜ4coeff(ν, εₒ, μₒ, k, β, r)
    # @show Tₒ; println()
    # @show ν, εₒ, μₒ, k, β, r

    T = Tₒ \ Tᵢ

    return T
end

const coeff_in = @SMatrix [1 0; 0 0; 0 1; 0 0]  # generate coefficients of innermost region when applied to first half of four unknowns
const coeff_out = @SMatrix [1 0; im 0; 0 1; 0 im]  # generate coefficients of outermost region when applied to last half of four unknowns

transfer_matrix(ν::Integer,
                i::Integer,  # index of interface
                r::AbsVecReal,  # radius of interface between regions
                ε::AbsVecReal,  # relative permittivities of regions
                μ::AbsVecReal,  # relative permeabilities of regions
                k::Real,  # vacuum wavenumber
                β::Real) =  # propagation constant
    Tcoeff(ν, r[i], SVec{2}(ε[i+1],ε[i]), SVec{2}(μ[i+1],μ[i]), k, β)

function dispfun(ν::Integer, fd::FiberDesign, k::Real, β::Real)
    @unpack r, μmdl, εmdl = fd

    λ = 2π / k
    ε = [εmdlᵢ(λ) for εmdlᵢ = εmdl]
    μ = [μmdlᵢ(λ) for μmdlᵢ = μmdl]

    disp = dispfun(ν, r, ε, μ, k, β)

    return disp
end

function dispfun(ν::Integer, 
                 r::AbsVecReal,  # radius of interface regions
                 ε::AbsVecReal,  # relative permittivities of regions
                 μ::AbsVecReal,  # relative permeabilities of regions
                 k::Real,  # vacuum wavenumber
                 β::Real)  # propagation constant
    M₁ = coeff_in
    for i = eachindex(r)
        M₁ = transfer_matrix(ν, i, r, ε, μ, k, β) * M₁
    end

    M₂ = -coeff_out
    M = hcat(M₁, M₂)

    disp = det(M)
    # @show M; println()

    return disp
end
