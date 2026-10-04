# Julia helpers for the propulsion models. Loaded into `DyadShip.Propulsion` by the
# generated module, so Dyad components in this submodule call them unqualified.

"""
    positive_part(x)

`max(x, 0)`. A plain function registered as symbolic, so it is evaluated pointwise and
is never turned into an event by discontinuity detection.
"""
positive_part(x::Real) = max(x, zero(x))

"""
    table_hold(xs, ys, x)

Piecewise-linear interpolation of the knots `(xs[i], ys[i])` with `xs` strictly
increasing, holding the end values outside `[xs[1], xs[end]]` (the Modelica
`CombiTable1Ds` `HoldLastPoint` extrapolation). Written as a sum of clipped ramps,

    y = ys[1] + Σᵢ sᵢ (min(max(x, xs[i]), xs[i+1]) - xs[i]),   sᵢ = Δyᵢ / Δxᵢ.

Registered as symbolic, so it is evaluated pointwise.
"""
function table_hold(xs::AbstractVector, ys::AbstractVector, x::Real)
    acc = ys[1] + zero(x)
    for i in 1:(length(xs) - 1)
        slope = (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i])
        acc = acc + slope * (min(max(x, xs[i]), xs[i + 1]) - xs[i])
    end
    return acc
end

"""
    table_covers(xs, x)

Coverage flag of a table used with `table_hold` for a nonnegative load `x`: `1.0` when
`x <= 0` (nothing is looked up) or `xs[1] <= x <= xs[end]`, `0.0` when a held end value
is in use. A plain function registered as symbolic, so it is evaluated pointwise and is
never turned into an event.
"""
table_covers(xs::AbstractVector, x::Real) = (x <= 0 || first(xs) <= x <= last(xs)) ? 1.0 : 0.0

"""
    table_knots_valid(xs)

`1.0` when the knots `xs` are all finite, nonnegative and strictly increasing, else
`0.0`. Usable on plain numbers before a run, and registered as a symbolic function so a
component can `assert` it on a parameter array.
"""
function table_knots_valid(xs::AbstractVector)
    all(isfinite, xs) || return 0.0
    first(xs) >= 0 || return 0.0
    all(xs[i + 1] > xs[i] for i in 1:(length(xs) - 1)) || return 0.0
    return 1.0
end

"""
    table_values_valid(ys)

`1.0` when the table values `ys` are all finite and positive, else `0.0`.
"""
table_values_valid(ys::AbstractVector) = all(y -> isfinite(y) && y > 0, ys) ? 1.0 : 0.0

"""
    rate_valid(x)

`1.0` when the rate `x` is finite and nonnegative, else `0.0`.
"""
rate_valid(x::Real) = isfinite(x) && x >= 0 ? 1.0 : 0.0

@register_symbolic positive_part(x::Real)
@register_symbolic table_hold(xs::AbstractVector, ys::AbstractVector, x::Real)
@register_symbolic table_covers(xs::AbstractVector, x::Real)
@register_symbolic table_knots_valid(xs::AbstractVector)
@register_symbolic table_values_valid(ys::AbstractVector)
@register_symbolic rate_valid(x::Real)
