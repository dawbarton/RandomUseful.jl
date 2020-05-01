export fourier_diff

"""
    fourier_diff([T=Float64,] N; order=1)

Create a Fourier differentiation matrix with numerical type T on the domain
`x = range(0, 2π, length=N+1)[1:end-1]`.
"""
function fourier_diff(T::Type{<:Number}, N::Integer; order=1)
    D = zeros(T, N, N)
    n1 = (N - 1) ÷ 2
    n2 = N ÷ 2
    x = LinRange{T}(0, π, N+1)
    if order == 1
        for i in 2:N
            sgn = (one(T)/2 - iseven(i))
            D[i, 1] = iseven(N) ? sgn*cot(x[i]) : sgn*csc(x[i])
        end
    elseif order == 2
        D[1, 1] = iseven(N) ? -N^2*one(T)/12 - one(T)/6 : -N^2*one(T)/12 + one(T)/12
        for i in 2:N
            sgn = -(one(T)/2 - iseven(i))
            D[i, 1] = iseven(N) ? sgn*csc(x[i]).^2 : sgn*cot(x[i])*csc(x[i])
        end
    else
        error("Not implemented")
    end
    for j in 2:N
        D[1, j] = D[N, j-1]
        D[2:N, j] .= D[1:N-1, j-1]
    end
    return D
end
fourier_diff(N::Integer; kwargs...) = fourier_diff(Float64, N; kwargs...)
