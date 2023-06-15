using WGLMakie
using Makie
using StoTrace
using LinearAlgebra
using Statistics

function get_rand_herm(N; normalize=false)

    H = rand(N,N)
    H = H + H'

    if normalize
        H = H ./ maximum(abs, H)
    end

    return H
end

function compare_mat(f; order=[5,10,15,20,30,35,40,45,50], Nvals=[5,50,100,500,1000], NR=10)

    fontsize_theme = Theme(fontsize = 23)
    set_theme!(fontsize_theme)

    fig = Figure()
    ax = Axis(fig[1,1], title="Error in the matrix function evaluation",
    xlabel="Chebyshev Order", ylabel="Relative Error")

    # Loop thorugh matrix sizes
    for N in Nvals

        # Loop number of realizations
        errors = zeros(NR, length(order))
        for j in 1:NR
            H = get_rand_herm(N, normalize=true)
            exact = StoTrace.mateval(f, H) |> tr

            # Loop through approximation orders
            for i in eachindex(order)
                cb = chebyshev_fit(f, (-1,1), 100, order=order[i])
                approx = StoTrace.mateval(H, cb) |> tr
                errors[j,i] = abs(exact-approx)/abs(exact)
            end

        end
        RMS = [mean(errors[:,i]) for i = eachindex(order)]
        STD = [std(errors[:,i]) for i = eachindex(order)]

        scatter!(ax, order, RMS, label="$N")
        errorbars!(ax, order, RMS, STD)
    end

    Legend(fig[1,2], ax, "Matrix Size")
    fig
end