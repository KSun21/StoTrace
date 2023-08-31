using StoTrace

export hutchinson

function rademacher(n)
    return rand([-1.0, 1.0], n)
end

function hutchinson(m, M, cb)

    out = 0.0

    l = size(M, 1)
    for i = 1:m
        v = rademacher(l)
        out += inner(M, v, cb)
    end

    return out / m
end