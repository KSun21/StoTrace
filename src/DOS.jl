export DOS

function gaussian(x, x0, σ)
    pref = 1/(σ√(2π))
    return pref * exp(-0.5*( (x-x0)/σ )^2)
end


function DOS(H, Emin, Emax, δE; cheb_order=300, cheb_npoints=600, k=30, ϵ=0.05, real=true, jackson=true)

    # Get range of values for density of states evaluation.
    # A delta function (approximated as a Gaussian) will be placed at each
    # point of this array. The important parameter here is δE, which gives the resolution of 
    # the DOS. It is related to the number of bins in a histogram.

    rE = (Emin+δE/2):δE:(Emax-δE/2)
    println(rE)
    Ne = length(rE)

    b = 0.5*(Emax + Emin)
    a = (Emax - Emin)/(2-ϵ)

    # Rescale the energy range to [-1,1]
    Ek = (rE .- b) ./ a

    # Compute the σ value for the gaussian used to approximate the delta function
    # We want the Gaussian to mostly vanish for (x-x0) values greater than δE.
    # where x0 is the center of each gaussian (previous array we just created)
    # This can be achieved if σ = δE/3. Since 99.8% of the normal distribution lies 
    # within ± 3σ = ± δE. Note tha we will take the scaled δE value.
    σ = step(Ek)/3

    # For each value of Ek in (∑ₖδ(E-Ek)) we will fit a chebyshev function...
    # Note that this will result in many fit evaluations which includes solving systems
    # of cheb_npoints x chebo_order size MANY times... While this is awful for small matrices (H)
    # It should be a fixed cost that get diluted as H gets really big. It will be hard to assess the impact of 
    # this step on the overall efficience untill we run tests including very large matrices.
    # However, for order cheb_order = 300 and cheb_npoints = 600) it takes around a second

    cbs = ChebyshevPolynomial{Float64}[]
    y = zeros(Float64, cheb_npoints)
    κ = 0:(cheb_npoints-1)
    xvals = collect(cos.(π * (κ .+ 0.5) / cheb_npoints))
    for x0 in Ek
        f(x) = gaussian(a*x+b, x0, σ)
        y .= f.(xvals)
        cb = chebyshev_fit(xvals, y, order=cheb_order, jackson=jackson)
        push!(cbs, cb)
    end

    # Perform the multiple cb Hutchinson method. Note that the output is an array, not a number
    sH = (H - b*I)/a
    s = zeros(Float64, Ne)

    l = size(sH, 1)
    for i = 1:k
        v = real ? rademacher(l) : complex_unit(l)
        inner!(s, sH, v, cbs)
    end

    return s ./ k
end

function inner!(s::Vector{Float64}, A, z0, cbs::Vector{ChebyshevPolynomial{Float64}})

    O = cbs[1].order
    # Alias to Chebyshev coefficients (so we can use 0-indexing)
    α(i,n) = cbs[i].coef[n+1]

    z1 = A*z0
    ζ0 = z0⋅z0
    ζ1 = z0⋅z1
    for i in eachindex(s)
        s[i] += α(i,0)*ζ0 + α(i,1)*ζ1 + α(i,2)*(2 * (z1⋅z1) - ζ0)
    end

    zⱼ₋₂ = z0
    zⱼ₋₁ = z1
    for j = 2:ceil(Int, O/2)
        zⱼ = 2 * (A*zⱼ₋₁) - zⱼ₋₂

        for i in eachindex(s)
            s[i] += α(i, 2j-1) * (2 * (zⱼ₋₁⋅zⱼ) - ζ1) 
        end

        if 2j-1 == O
            break
        end

        for i in eachindex(s)
            s[i] += α(i, 2j) * (2 * (zⱼ⋅zⱼ) - ζ0) 
        end

        zⱼ₋₂ = zⱼ₋₁
        zⱼ₋₁ = zⱼ
    end
end
