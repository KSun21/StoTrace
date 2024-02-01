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

function hutchinson(m, M, pl; real=true, hermitian=true)

    out = 0.0

    l = size(M, 1)
    v = real ? zeros(l) : zeros(ComplexF64, l)
    for i = 1:m
        rademacher!(v)
        out += hermitian ? inner(M, v, pl) : non_hermitian_inner(M, v, pl)
    end

    return out / m
end