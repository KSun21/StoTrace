using WGLMakie
using Makie
using StoTrace
using LinearAlgebra
using Statistics

function get_rand_herm(N; normalize=false, eignorm=false, colnorm=false)

    H = rand(N,N)
    H = H + H'

    if normalize
        H = H ./ maximum(abs, H)
    end

    if eignorm
        e,_ = eigen(H)
        m = maximum(abs.(e))
        H = H ./ m
    elseif colnorm
        colsums = [sum(H[:,i]) for i in 1:N]
        m = maximum(colsums)
        H = H ./ m
    end

    return H
end

function compare_mat(f; order=[5,10,15,20,30,35,40,45,50], Nvals=[5,50,100,500,1000], NR=10, normalize=false, eignorm=false,colnorm=false)

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
            H = get_rand_herm(N, normalize=normalize, eignorm=eignorm,colnorm=colnorm)
            exact = StoTrace.mateval(f, H) |> tr
            println("Exact: $exact")

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

function compare_sparse_mat(f; order=[5,10,15,20,30,35,40,45,50], N=50, dvals=[0.2, 0.4, 0.6, 0.8], NR=10)

    fontsize_theme = Theme(fontsize = 23)
    set_theme!(fontsize_theme)

    fig = Figure()
    ax = Axis(fig[1,1], title="Error in the matrix function evaluation",
    xlabel="Chebyshev Order", ylabel="Relative Error")

    # Loop thorugh matrix sizes
    for d in dvals

        # Loop number of realizations
        errors = zeros(NR, length(order))
        for j in 1:NR
            H = Array(sprand(N,N,d/2))
            H = H + H'
            colsums = [sum(H[:,i]) for i in 1:N]
            m = maximum(colsums)
            H = H ./ m
            exact = mateval(f, H) |> tr

            # Loop through approximation orders
            for i in eachindex(order)
                cb = chebyshev_fit(f, (-1,1), 100, order=order[i])
                approx = mateval(H, cb) |> tr
                errors[j,i] = abs(exact-approx)/abs(exact)
            end

        end
        RMS = [mean(errors[:,i]) for i = eachindex(order)]
        STD = [std(errors[:,i]) for i = eachindex(order)]

        scatter!(ax, order, RMS, label="$d")
        errorbars!(ax, order, RMS, STD)
    end

    Legend(fig[1,2], ax, "Matrix density")
    fig
end