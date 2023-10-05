using StoTrace
using Distributions

export hutchinson

function rademacher(n)
    return rand([-1.0, 1.0], n)
end

function complex_unit(n)
    dist = Uniform(0, 2π)
    return [exp(-im*θ) for θ = rand(dist, n)]
end

function hutchinson(m, M, cb; real=true)

    out = 0.0

    l = size(M, 1)
    for i = 1:m
        v = real ? rademacher(l) : complex_unit(l)
        out += inner(M, v, cb)
    end

    return out / m
end