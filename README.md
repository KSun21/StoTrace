# StoTrace

## Setup

1. Clone the code
`git clone git@github.com:gustavojra/StoTrace.git`
2. Activate the package, you can do that when you initialize julia as `julia --project=PATH_TO_CODE`
or in the julia terminal you can run
```
using Pkg
Pkg.activate('PATH_TO_CODE')
```

## Basically usage.

To compute $\text{tr}f(A)$ stochastic use the `trace_estimate` function. This function takes four mandatory positional arguments

`trace_estimate(f::Function, A, λmin, λmax)`

where `f` is the function to by applied on the matrix `A`. `λmin` and `λmax` give the bounds to the spectrum of `A`, in julia they can be obtained, for example, using the package [`KrylovKit`](https://github.com/Jutho/KrylovKit.jl) where one can get the largest and the smallest eigenvalues, for example:

```
using KrylovKit

# Get random hermitian matrix
H = rand(50,50)
H = H + H'

# Compute maximum and minimum eigenvalues
λmax = eigsolve(H, 1, :LR)[1][1]
λmin = eigsolve(H, 1, :SR)[1][1]
```

## Minimal Example
```
using StoTrace

# Example matrix
H = [1 0;
     0 2]

# For which we want tr(log(H)). Since this is a diagonal matrix the exact result
# is log(1) + log(2) = 0.6931

# Get the stochastic approximation
trace_estimate(log, H, 1, 2)
```
This returns `0.6931819586364876`

# Approximation parameters

For a finer control, the keyword arguments can be played with. The complete function call is

`trace_estimate(f::Function, A, λmin, λmax; cheb_order=4, cheb_npoints=1000, k=10, ϵ=0.01, verb=false)`

| Keyword Argument      | Description |
| ----------- | ----------- |
| cheb_order      | Approximation order of the Chebyshev polynomial.  **(Important!)**     |
| cheb_points      | Number of grid points used for Chebyshev fitting.    (Not very important)  |
| k   | Number of vectors used in the stochastic evaluation (Hutchinson step).    (Not very important)    |
|  ϵ  | Off-set used in the scaling of the matrix A. There is not much need to play with this.    |
|  verb  | Flag to make the code spit a bunch of information as the computation carries on. Only useful for debugging.    |

