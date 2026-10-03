# Julia helpers for the propulsion models. Loaded into `DyadShip.Propulsion` by the
# generated module, so Dyad components in this submodule call them unqualified.

"""
    table_hold(xs, ys, x)

Piecewise-linear interpolation of the knots `(xs[i], ys[i])` with `xs` strictly
increasing, holding the end values outside `[xs[1], xs[end]]` (the Modelica
`CombiTable1Ds` `HoldLastPoint` extrapolation). Written as a sum of clipped ramps,

    y = ys[1] + Σᵢ sᵢ (min(max(x, xs[i]), xs[i+1]) - xs[i]),   sᵢ = Δyᵢ / Δxᵢ,

so it stays a closed-form expression when `x` is symbolic.
"""
function table_hold(xs::AbstractVector, ys::AbstractVector, x)
    acc = ys[1]
    for i in 1:(size(xs, 1) - 1)
        slope = (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i])
        acc = acc + slope * (min(max(x, xs[i]), xs[i + 1]) - xs[i])
    end
    return acc
end

"""
    table_in_range(xs, x)

`1` when `xs[1] <= x <= xs[end]` (the argument lies inside the tabulated range of
`table_hold`), `0` when the held end value is being used instead.
"""
function table_in_range(xs::AbstractVector, x)
    return ifelse(x < xs[1], 0, ifelse(x > xs[size(xs, 1)], 0, 1))
end
