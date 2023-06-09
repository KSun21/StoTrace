using Polynomials
using LinearAlgebra

export ChebyshevPolynomial
export chebyshev_fit, eval_fit

struct ChebyshevPolynomial{T,P}
    coef::Vector{T}
    order::Int
    domain::Tuple{P,P}
    a::P
    b::P
end

function get_fitting_domain(domain::Tuple, npoints::Int)
    if npoints < 2
        throw(ArgumentError("Number of fitting points (npoints) cannot be less than 2. (Though 2 would not be great either, but you do you)."))
    end

    x0 = domain[1]
    xf = domain[2]
    if x0 ≥ xf
        throw(ArgumentError("Invalid domain. Please input an ordered tuple, not $domain."))
    end

    Δx = (domain[2] - domain[1]) / (npoints-1) 
    return collect(x0:Δx:xf)
end

chebyshev_fit(f::Function, domain::Tuple, npoints::Int; order::Int=5) = chebyshev_fit(f, get_fitting_domain(domain, npoints), order=order)
chebyshev_fit(f::Function, xvals::Vector; order::Int=5) = chebyshev_fit(xvals, [f(x) for x in xvals], order=order)
chebyshev_fit(yvals::Vector, domain::Tuple, npoints::Int; order::Int=5) = chebyshev_fit(get_fitting_domain(domain, npoints), yvals, order=order)

function chebyshev_fit(xvals::Vector, yvals::Vector; order::Int=5)

    if order < 2
        throw(ArgumentError("Chebyshev degree (order) has to be greater than 1."))
    end

    # Get domain
    domain = (xvals[1], xvals[end])

    # x-values need to be scaled and shifted such that the fitting domain
    # becomes -1:1, a requirement from the Chebyshev procedure
    # Using a*x₀ + b = -1 and a*xf + b = 1 we get:
    a = 2/(domain[2] - domain[1])
    b = (domain[1] + domain[2]) / (domain[1] - domain[2])

    xs = a .* xvals .+ b

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

    T = zeros(length(xs), order)

    # Since T₀(x) = 1, the first column is simply
    T[:,1] .= 1

    # T₁(x) = x, hence
    T[:,2] .= xs

    # The remaining terms are populated using the recursive relation
    # Tₙ(x) = 2xTₙ₋₁ - Tₙ

    # Loop through columns (n)
    for n in 3:order
        for i in eachindex(xs)
            T[i,n] = 2*xs[i]*T[i,n-1] - T[i,n-2]
        end
    end

    # Solve the matrix equation Tc = y 
    c = T \ yvals

    return ChebyshevPolynomial(c, order, domain, a, b)
end

function eval_fit(x, cb::ChebyshevPolynomial{T}) where T

    # Get shifted value of x
    xs = cb.a * x + cb.b

    Tn = zeros(T, cb.order)

    Tn[1] = 1
    if cb.order == 1
        return cb.coef[1]
    end

    Tn[2] = xs

    for n in 3:cb.order
        Tn[n] = 2xs*Tn[n-1] - Tn[n-2]
    end

    return dot(Tn, cb.coef)
end