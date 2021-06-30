export cheb_diff, cheb_mesh

"""
    cheb_mesh(T::Type, N::Integer)
    cheb_mesh(N::Integer)

Return the Chebyshev nodes of degree `N` on the interval [0, 1].
"""
cheb_mesh(T::Type, N::Integer) = (-cos.((T(π) * (0:N)) ./ N) .+ 1) ./ 2
cheb_mesh(N::Integer) = cheb_mesh(Float64, N)

"""
    cheb_diff(T::Type, N::Integer)
    cheb_diff(N::Integer)

Return a tuple `(x, D)` where `x` contains the Chebyshev nodes of degree `N` on the interval
[0, 1] using the numerical type `T` and `D` is the Chebyshev differentiation matrix on the
same interval.

## Example - differentiate a sinusoid

```julia
julia> (x, D) = cheb_diff(10) ;

julia> D * sinpi.(x) ≈ π * cospi.(x)
true
```
"""
function cheb_diff(T::Type, N::Integer)
    if N == 0
        return (ones(T, 1), zeros(T, 1, 1))
    else
        x = cheb_mesh(T, N)
        c = [2; ones(T, N - 1, 1); 2] .* (-1) .^ (0:N)
        dx = x .- x'
        D = (c * (1 ./ c)') ./ (dx + I)
        D[diagind(D)] .-= vec(sum(D; dims=2))
        return (x, D)
    end
end
cheb_diff(N::Integer) = cheb_diff(Float64, N)
