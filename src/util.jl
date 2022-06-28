# The argument z of Jfun and Yfun are complex in general, because q² used to construct z = qr
# can be either positive or negative.

Jfun(ν::Integer, z::Real) = besselj(ν,z)
function Jfun(ν::Integer, z::Number)
    ℜz, ℑz = real(z), imag(z)
    if iszero(ℑz)  # z is real
        J = complex(Jfun(ν,ℜz))
    else
        @argcheck iszero(ℜz)  # z should be purely imaginary
        J = im^ν * besseli(ν,ℑz)
    end

    return J
end

Yfun(ν::Integer, z::Real) = bessely(ν,z)
function Yfun(ν::Integer, z::Number)
    ℜz, ℑz = real(z), imag(z)
    if iszero(ℑz)  # z is real
        Y = complex(Yfun(ν,ℜz))
    else
        @argcheck iszero(ℜz)  # z should be purely imaginary
        Y = im^(ν+1) * besseli(ν,ℑz) - 2/π * (-im)^ν * besselk(ν,ℑz)
    end

    return Y
end

deriv(Ffun::Any, z::Real) = derivative(Ffun, z)
function deriv(Ffun::Any, z::Number)
    ℜz, ℑz = real(z), imag(z)
    if iszero(ℑz)  # z is real
        F′ = complex(deriv(Ffun,z))
    else
        @argcheck iszero(ℜz)  # z should be purely imaginary

        dualz = im * Dual(ℑz,1)  # treated as purely imaginary inside Ffun
        dualF = Ffun(dualz)

        F′ = real(dualF).partials[1] + im * imag(dualF).partials[1]
    end

    return F′
end
