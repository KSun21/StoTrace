using Makie
using WGLMakie
using StoTrace
using LinearAlgebra


function plotexpconv(;N=100, αvals=[-50, -30, -20, -10, -1, 1, 10, 30, 50], order=[10, 20, 30, 40, 50])

    fig = Figure()
    ax = Axis(fig[1,1], ylabel="Relative Error", xlabel="Chebyshev Order")

    # Get hermitian matrix
    H = rand(N,N)
    H .*= rand([-1,1], N,N)
    H .+= H'
    cheb_scale!(H)

    z = rand(size(H,1))

    errs = zeros(length(order))

    for α in αvals

        # Compute exact value
        exact = z ⋅ (exp(α .* H) * z)
        println("\nExact: $exact")

        errs .= 0.0
        for i in eachindex(order)
            cb = chebyshev_fit(x-> exp(α .* x), (-1, 1), 1000, order=order[i])
            approx = inner(H, z, cb)
            errs[i] = abs(approx - exact)/exact
            println("Order $(order[i]) -> $approx -> $(errs[i])")
        end
        lines!(ax, order, errs, label=L"\alpha = %$α")
        scatter!(ax, order, errs, marker= α > 0 ? :utriangle : :dtriangle, label=L"\alpha = %$α")
    end

    axislegend(ax, position=:rt, merge=true)

    ylims!(ax, -0.1, 1)

    fig
end