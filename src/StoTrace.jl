module StoTrace

using LinearAlgebra


function rademacher(n)
    v = rand([-1.0, 1.0], n)
    return v
end

function hutchinson(m, M)

    out = 0.0

    l = size(M, 1)
    temp = zeros(l)
    invm = 1/m

    for i = 1:m
        v = rademacher(l)
        mul!(temp, M, v)

        out += invm * dot(v, temp)
    end

    return out
end

include("chebyshev.jl")

end # module
