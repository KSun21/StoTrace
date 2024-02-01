module StoTrace

export trace_estimate

using LinearAlgebra

include("Chebyshev.jl")
include("Hutchinson.jl")
include("DOS.jl")

function trace_estimate(f::Function, A, λmin, λmax; order=4, npoints=1000, k=10, ϵ=0.01, verb=false, real=true, chebyshev=false, hermitian=true)

    if chebyshev
        return chebyshev_trace_estimate(f, A, λmin, λmax; order=order, npoints=npoints, k=k, ϵ=ϵ, verb=verb, real=real)
    end

    pln(x) = verb ? println(x) : nothing

    pln("Getting Polynomial fit using order = $order")
    pl = polynomial_fit(f, λmin, λmax, npoints, order=order)

    pln("Done.")
    pln("Last Polynomial coefficient: $(pl.coef[end])")

    pln("Starting Hutchinson iterations")
    out = hutchinson(k, A, pl, real=real, hermitian=hermitian)

    pln("Final estimated trace: $out")

    return out
end

function chebyshev_trace_estimate(f::Function, A, λmin, λmax; order=4, npoints=1000, k=10, ϵ=0.01, verb=false, real=true)

    pln(x) = verb ? println(x) : nothing

    b = 0.5*(λmax + λmin)
    a = (λmax - λmin)/(2-ϵ)

    pln("Getting Chebyshev fit using order = $order")
    cb = chebyshev_fit(x->f(a*x+b), npoints, order=order)

    pln("Done.")
    pln("Last Chebyshev coefficient: $(cb.coef[end])")

    pln("Starting Hutchinson iterations")
    sA = (A - b*I)/a
    out = hutchinson(k, sA, cb, real=real)

    pln("Final estimated trace: $out")

    return out
end

end # module
