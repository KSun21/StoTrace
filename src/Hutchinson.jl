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

function rademacher!(v::Array{R}) where R <: Real
    for i in eachindex(v)
        v[i] = rand([-one(R), one(R)])
    end
end

function rademacher!(v::Array{C}) where C <: Complex
    dist = Uniform(0, 2π)
    for i in eachindex(v)
        θ = rand(dist)
        v[i] = exp(-im*θ)
    end
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