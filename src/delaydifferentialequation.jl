export dde_stability

"""
    dde_stability(A, B, τ; N=10)

Calculate (a finite number of) the characteristic roots of the linear delay differential
equation

```math
x'(t) = Ax(t) + Bx(t-τ)
```

where `A` and `B` can be (square) matrices. `N` controls the order of the Chebyshev
polynomial approximation to the infinitesimal generator; higher values of `N` result in a
more accurate approximation. A property of the algorithm used is that characteristic roots
closest to the origin are approximated best.

## Reference

Pseudospectral Differencing Methods for Characteristic Roots of Delay Differential
Equations; D. Breda, S. Maset, and R. Vermiglio, SIAM Journal on Scientific Computing 2005
27:2, 482-495
"""
function dde_stability(A, B, τ; N=10)
    @assert size(A) == size(B)
    @assert size(A, 1) == size(A, 2)
    T = eltype(A) == Int ? Float64 : eltype(A)
    n = size(A, 1)
    (x, D) = cheb_diff(T, N)
    M = kron(-D, Matrix(I / convert(T, τ), n, n))
    M[1:n, 1:n] .= A
    M[1:n, (n + 1):(end - n)] .= 0
    M[1:n, (end - n + 1):end] .= B
    return eigen(M)
end

# Tests
# 1. x'(t) = 2x(t) - exp(1)*x(t - 1) has a double root at λ = 1
# 2. x'(t) = -x(t - τ) goes unstable at τ=π/2
# 3. x'(t) = α*A*x(t) + (1-α)*A*x(t-τ)
#    A = [
#          50  284   41   23   50   32
#        -280  -46  -19  -37  -10  -28
#          35   -1   26  143   35   17
#          5   -31 -139  -22    5  -13
#          20  -16   11   -7 -115  137
#         -10  -46  -19  -37 -145 -163
#    ]
#    with α=0.9 is stable for 3π/(270*(2*0.9 - 1)) ≤ τ ≤ 4π/270
#    that is 0.04363323129985824 ≤ τ ≤ 0.046542113386515455
