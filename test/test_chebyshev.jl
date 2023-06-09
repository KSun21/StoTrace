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

end