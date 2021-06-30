@testset "delaydifferentialequation" begin
# Tests
# 1. x'(t) = 2x(t) - exp(1)*x(t - 1) has a double root at λ = 1
    ev = dde_stability(2, -exp(1), 1)
    @test count((abs.(real.(ev.values) .- 1) .< 1e-6) .& (abs.(imag.(ev.values)) .< 1e-6)) == 2
# 2. x'(t) = -x(t - τ) goes unstable at τ=π/2
    ev = dde_stability(0, -1, π/2 - 0.01)
    @test all(real.(ev.values) .< 0)
    ev = dde_stability(0, -1, π/2 + 0.01)
    @test any(real.(ev.values) .> 0)
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
    A = [
          50  284   41   23   50   32
        -280  -46  -19  -37  -10  -28
          35   -1   26  143   35   17
          5   -31 -139  -22    5  -13
          20  -16   11   -7 -115  137
         -10  -46  -19  -37 -145 -163
    ]
    α = 0.9
    ev = dde_stability(α * A, (1 - α)*A, 0.043)
    @test any(real.(ev.values) .> 0)
    ev = dde_stability(α * A, (1 - α)*A, 0.044)
    @test all(real.(ev.values) .< 0)
    ev = dde_stability(α * A, (1 - α)*A, 0.047)
    @test any(real.(ev.values) .> 0)
end
