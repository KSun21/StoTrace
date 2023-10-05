module StoTrace

export trace_estimate

using LinearAlgebra

include("Chebyshev.jl")
include("Hutchinson.jl")
include("DOS.jl")

function trace_estimate(f::Function, A, λmin, λmax; cheb_order=4, cheb_npoints=1000, k=10, ϵ=0.01, verb=false, real=true, jackson=false)

    pln(x) = verb ? println(x) : nothing

    b = 0.5*(λmax + λmin)
    a = (λmax - λmin)/(2-ϵ)

    pln("Getting Chebyshev fit using order = $cheb_order")
    cb = chebyshev_fit(x->f(a*x+b), cheb_npoints, order=cheb_order, jackson=jackson)

    pln("Done.")
    pln("Last Chebyshev coefficient: $(cb.coef[end])")

    pln("Starting Hutchinson iterations")
    sA = (A - b*I)/a
    out = hutchinson(k, sA, cb, real=real)

    pln("Final estimated trace: $out")

    return out
end

end # module
