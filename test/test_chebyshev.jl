using LinearAlgebra

@testset "Chebyshev Fit" begin
    
    # Test error throws
    @test_throws ArgumentError chebyshev_fit(sin, 1000, order=0)
    @test_throws BoundsError chebyshev_fit(sin, -2:0.01:1, order=1)

    # Trivial test for chebyshev_fit, use a linear function
    f(x) = 0.5 * x - 2
    @test chebyshev_fit(f, 100, order=1).coef ≈ [-2, 0.5]

    # Higher order coefficients need to return 0
    @test chebyshev_fit(f, 100, order=4).coef ≈ [-2, 0.5, 0.0, 0.0, 0.0]

    # If we try fitting a cubic polynomial the analysis gets more trick
    # because the polinomial coefficients pₙ are combinations of the Chebyshev 
    # coefficients cₙ
    (p1, p2, p3, p4) = (0.7, 0.3, -0.5, 1.1)
    f(x) = p1 + p2*x + p3*x^2 + p4*x^3
    c = chebyshev_fit(f, 1000, order=3).coef
    @test (p4 ≈ 4c[4]) & (p3 ≈ 2c[3]) & (p2 ≈ c[2] - 3c[4]) & (p1 ≈ c[1] - c[3])

    # Numerical tests. Now we look at some numerical fit and set some RMS tolerance
    tols = [1e-1, 1e-3, 1e-5, 1e-7]
    ord = [4, 9, 14, 19]
    xvals = -1:0.01:1
    f(x) = 3*x^2 - 10*log(x+2)
    exact = f.(xvals)
    for (t,o) in zip(tols, ord) 
        cb = chebyshev_fit(f, 100, order=o)
        fitted = [eval_fit(x, cb) for x in xvals]
        rms = sum(((fitted .- exact).^2) / length(xvals)) |> sqrt
        @test rms < t
    end
end