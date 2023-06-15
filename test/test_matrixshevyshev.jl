@testset "Chebyshev Matrix Fit" begin

    # Start out with a small matrix
    H = rand(5,5)
    H = H + H'
    H = H ./ maximum(abs.(H))

    cb = chebyshev_fit(exp, (-1,1), 200, order=20)

    exact = StoTrace.mateval(exp, H) |> tr 
    approx = StoTrace.mateval(H, cb) |> tr

    @test abs(exact - approx)/exact < 1e-5


end