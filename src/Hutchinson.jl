using StoTrace

export trace1

function rademacher(n)
    return rand([-1.0, 1.0], n)
end

function trace1(A, cb, m)

    out = 0.0

    w1 = zeros(size(A,1))
    w2 = zeros(size(A,1))
    w3 = zeros(size(A,1))
    u = zeros(size(A,1))
    for i = 1:m
        v = rademacher(size(A,1))
        w1 .= v
        w2 .= cb.a .* (A*v) - cb.b .* v
        #mul!(w2, A, v, cb.a, cb.b)
        u .= cb.coef[1]*w1 + cb.coef[2]*w2
        for j in 3:cb.order
            w3 .= 2cb.a .* A * w2 + 2cb.b .* w2 - w1
            u .+= cb.coef[j] .* w3
            w1 .= w2
            w2 .= w3
        end
        out += dot(v, u)
    end

    return out / m
end