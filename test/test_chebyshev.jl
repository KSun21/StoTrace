@testset "Chebyshev Fit" begin
    
    # Test error throws
    @test_throws ArgumentError chebyshev_fit(sin, (-1,1), 1)
    @test_throws ArgumentError chebyshev_fit(sin, (-1,1), 2, order=1)
    @test_throws ArgumentError chebyshev_fit(sin, (1,-1), 2, order=2)

    # Trivial test for chebyshev_fit, use a linear function
    f(x) = 0.5 * x - 2
    @test chebyshev_fit(f, (-1,1), 10, order=2).coef ≈ [-2, 0.5]

    # Higher order coefficients need to return 0
    @test chebyshev_fit(f, (-1,1), 10, order=5).coef ≈ [-2, 0.5, 0.0, 0.0, 0.0]

    # If we try fitting a cubic polynomial the analysis gets more trick
    # because the polinomial coefficients pₙ are combinations of the Chebyshev 
    # coefficients cₙ
    (p1, p2, p3, p4) = (0.7, 0.3, -0.5, 1.1)
    f(x) = p1 + p2*x + p3*x^2 + p4*x^3
    c = chebyshev_fit(f, (-1,1), 4, order=4).coef
    @test (p4 ≈ 4c[4]) & (p3 ≈ 2c[3]) & (p2 ≈ c[2] - 3c[4]) & (p1 ≈ c[1] - c[3])

    # Numerical tests. Now we look at some numerical fit and set some RMS tolerance
    tols = [1e-1, 1e-3, 1e-5, 1e-7]
    ord = [5, 10, 15, 20]
    xvals = 1:0.01:5
    f(x) = 3/x^2 - 10*log(x)
    exact = f.(xvals)
    for (t,o) in zip(tols, ord) 
        cb = chebyshev_fit(f, (1,5), 100, order=o)
        fitted = [eval_fit(x, cb) for x in xvals]
        rms = sum(((fitted .- exact).^2) / length(xvals)) |> sqrt
        @test rms < t
    end

    # Compare coefficients with numpy
    x = collect(0:0.01:5)
    y = exp.(x)
    cb = chebyshev_fit(x, y, order=5)
    @test cb.coef ≈ [40.04581473, 60.95298089, 31.03200979, 11.12392713,  3.27698755]

    x = collect(1:0.01:5)
    y = log.(x)
    cb = chebyshev_fit(x, y, order=7)
    @test cb.coef ≈ [0.96244263,  0.76380517, -0.14585858,  0.03701245, -0.01059841, 0.00307569, -0.00097732]

    x = collect(10:1.0:100)
    y = [(x^2 - 1/x)*sin(x) for x in x]
    cb = chebyshev_fit(x, y, order=10)
    @test cb.coef ≈ [-487.36553758,-966.80521864,-963.92179365,-947.14866806,-918.00999684,-883.81423028,-798.09985261,-734.84252216,-543.63284474,-457.03974943]
end