using Polynomials
using LinearAlgebra

export ChebyshevPolynomial
export chebyshev_fit, eval_fit
export mateval, inner

struct ChebyshevPolynomial{T}
    coef::Vector{T}
    order::Int
end

chebyshev_fit(f::Function, npoints::Int; order::Int=5) = chebyshev_fit(f, collect(-1:2/(npoints-1):1), order=order)
chebyshev_fit(f::Function, xvals::Vector; order::Int=5) = chebyshev_fit(xvals, [f(x) for x in xvals], order=order)

function chebyshev_fit(x::Vector, y::Vector; order::Int=4)

    if order < 1
        throw(ArgumentError("Chebyshev degree (order) has to be greater than 0."))
    end

    # x-values need to be scaled and shifted such that the fitting domain
    # becomes -1:1, a requirement from the Chebyshev procedure
    if maximum(x) > 1 || minimum(x) < -1
        throw(BoundsError("Chebyshev fitting can only be performed within the [-1,+1] domain."))
    end

    # Create a matrix T composed of Chebyshev degree Tn(x) along rows (i.e. each row is a degree [T₀(x), T₁(x), T₂(x)...])
    # with different values of x along columns. For example, if we have three values of x (x₀, x₁, x₂) and a degree 3 
    # Chebyshev polynomial the T arrays would look like
    # +--------+-------+-------+
    # | T₀(x₀) | T₁(x₀)| T₂(x₀)|
    # +--------+-------+-------+
    # | T₀(x₁) | T₁(x₁)| T₂(x₁)|
    # +--------+-------+-------+
    # | T₀(x₂) | T₁(x₂)| T₂(x₂)|
    # +--------+-------+-------+

    T = zeros(length(x), order+1)

    # Since T₀(x) = 1, the first column is simply
    T[:,1] .= 1

    # T₁(x) = x, hence
    T[:,2] .= x

    # The remaining terms are populated using the recursive relation
    # Tₙ(x) = 2xTₙ₋₁ - Tₙ

    # Loop through columns (n)
    for n in 3:(order+1)
        for i in eachindex(x)
            T[i,n] = 2*x[i]*T[i,n-1] - T[i,n-2]
        end
    end

    # Solve the matrix equation Tc = y 
    c = T \ y

    return ChebyshevPolynomial(c, order)
end

function eval_fit(x, cb::ChebyshevPolynomial{T}) where T

    # Get shifted value of x
    xs = cb.a * x + cb.b

    Tn = zeros(T, cb.order+1)

    Tn[1] = 1
    if cb.order == 0
        return cb.coef[1]
    end

    Tn[2] = xs

    for n in 3:(cb.order+1)
        Tn[n] = 2xs*Tn[n-1] - Tn[n-2]
    end

    return dot(Tn, cb.coef)
end

# Auxiliary function to evalute mattrix functions
function mateval(f::Function, A::Matrix)
    λ, U = eigen(A)

    fvals = diagm(f.(λ))

    return U * fvals * U'
end

# Auxiliary function to evalute mattrix functions using Chebyshev expansion
function mateval(cb::ChebyshevPolynomial, A::Matrix)

    # Zeroth order is just an identity matrix
    T0 = zeros(size(A))
    T0[diagind(T0)] .= 1.0 

    if cb.order == 0
        return cb.coef[1] .* T0
    end

    # First order is just A
    T1 = similar(A)
    T1 .= A

    out = cb.coef[1] .* T0 + cb.coef[2] .* T1

    T2 = similar(A)

    # Forst third-order and beyond use regular recursive relation
    # Tn = 2xTₙ₋₁ - Tₙ₋₂
    for n in 3:(cb.order+1)
        T2 .= 2 .*A*T1 - T0 
        out += cb.coef[n] .* T2

        T0 .= T1
        T1 .= T2
    end

    return out
end

"""
    inner(A, z0, cb)

Compute the inner product ⟨z₀|f(A)|z₀⟩ where f(A) is approximated by a Chebyshev expansion (cb).
The algorithm implemented here is described by Hallman (https://arxiv.org/abs/2101.00325v1 - Algorithm 3.1)
"""
function inner(A, z0, cb::ChebyshevPolynomial)

    # Alias to Chebyshev coefficients (so we can use 0-indexing)
    α(n) = cb.coef[n+1]

    z1 = A*z0
    ζ0 = z0⋅z0
    ζ1 = z0⋅z1
    s = α(0)*ζ0 + α(1)*ζ1 + α(2)*(2 * (z1⋅z1) - ζ0)

    zⱼ₋₂ = z0
    zⱼ₋₁ = z1
    for j = 2:ceil(Int, cb.order/2)
        zⱼ = 2 * (A*zⱼ₋₁) - zⱼ₋₂
        s = s + α(2j-1) * (2 * (zⱼ₋₁⋅zⱼ) - ζ1) 
        if 2j-1 == cb.order
            break
        end
        s = s + α(2j) * (2 * (zⱼ⋅zⱼ) - ζ0) 

        zⱼ₋₂ = zⱼ₋₁
        zⱼ₋₁ = zⱼ
    end

    return s
end