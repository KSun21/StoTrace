using Makie
using WGLMakie
using StoTrace
using LinearAlgebra
using Statistics

function get_error(f, N, o, k; verb=false)

    pln(x) = verb ? println(x) : nothing

    # Get matrix and its eigenvalues
    H = rand(N,N) + diagm(repeat([20], N))
    H = H + H'

    evals, evec = eigen(H)

    exact = sum(f.(evals))
    pln("Exact value: $exact")
    Emin = evals[1]
    Emax = evals[end]
    pln("Spectrum bounds $Emin - $Emax")

    approx = trace_estimate(f, H, Emin, Emax, cheb_order=o, k=k, verb=verb)

    pln("Approx value: $approx")

    return abs(approx-exact)/exact
end


function plot_error(f; Nvals=[50, 500, 100, 1000], realizations=5, order=[10,25,50,75,100,150,200], kvals=[5,10,15])

    fig = Figure()

    @assert iseven(length(Nvals))

    nrows = length(Nvals) ÷ 2

    axs = Axis[]
    for i = 1:nrows
        push!(axs, Axis(fig[i,1], title="Matrix size = $(Nvals[2*i-1])"))
        push!(axs, Axis(fig[i,2], title="Matrix size = $(Nvals[2*i])"))
    end
    linkaxes!(axs...)

    for n in eachindex(axs)
        ax = axs[n]
        N = Nvals[n]
        for k in kvals 
            errs = zeros(realizations, length(order))
            for (i,o) in enumerate(order)
                for r in 1:realizations
                    errs[r,i] = get_error(f, N, o, k)
                end
            end
            avgerr = [mean(errs[:,o]) for o = eachindex(order)]
            stderr = [std(errs[:,o]) for o = eachindex(order)]
            scatter!(ax, order, avgerr, label="$k")
            lines!(ax, order, avgerr, label="$k")
            errorbars!(ax, order, avgerr, stderr)
        end
    end

    Label(fig[nrows+1, 1:end], "Chebyshev Order")
    Label(fig[1:end, 0], "Relavive Absolute Error", rotation=π/2)
    Legend(fig[1:nrows,3], axs[1], "# of vectors", merge=true)
    fig
end
