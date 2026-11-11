export β, Δt, V, ∇V, langevin_step

const β = 2.5
const Δt = 5.0e-3

const σ = sqrt(2Δt / β)

function V(x::Point)
    a, b = x
    return 3exp(-a^2) * (exp(-(b - 1 / 3)^2) - exp(-(b - 5 / 3)^2)) -
           5exp(-b^2) * (exp(-(a - 1)^2) + exp(-(a + 1)^2)) + 0.2a^4 + 0.2(b - 1 / 3)^4 # three-well "entropic switch" potential
end

function ∇V(x::Point)
    a, b = x
    e_ga = exp(-a^2)
    e_l = exp(-(a + 1)^2 - b^2)
    e_r = exp(-(a - 1)^2 - b^2)
    e_1 = exp(-(b - 1 / 3)^2)
    e_5 = exp(-(b - 5 / 3)^2)
    da = 0.8a^3 + 10 * ((a - 1) * e_r + (a + 1) * e_l) - 6a * e_ga * (e_1 - e_5)
    db = 10b * (e_r + e_l) + 3e_ga * (2(b - 5 / 3) * e_5 - 2(b - 1 / 3) * e_1) +
         0.8(b - 1 / 3)^3
    return Point(da, db)
end

langevin_step(x::Point, rng) = x - Δt * ∇V(x) + σ * randn(rng, Point)