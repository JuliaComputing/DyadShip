# Rainflow cycle counting and Miner damage accumulation.
#
# `ShipSIM.Components.DataProcessing.RainflowCounter` and `FatigueCounter` run the count
# inside the model, on a `when u <> pre(u)` event with an `algorithm` section that keeps
# two stacks and walks them with `while` loops. Dyad's clocked sublanguage could express
# a bounded version of that — a fixed-size stack updated one-hot, one extraction per tick
# — but it is the wrong home for it. Counting cycles is post-processing over a series of
# turning points, not a dynamic system: it has no state the rest of the model can see, it
# needs an unbounded stack to be exact, and clocked array state does not cross-compile
# (SynchToolkit #212), so a model carrying one could not be deployed.
#
# The model's job is to produce the turning points, which `Machinery.EventPeakSampler`
# already does exactly — it clocks on the zero crossings of the signal's derivative, so
# every sample it emits is a reversal. These functions consume that series.

"""
    turning_points(x) -> Vector{Float64}

Reversals of `x`: the first and last points and every local extremum, with runs of equal
or monotonically ordered values collapsed. A series already produced by
`EventPeakSampler` passes through unchanged apart from its endpoints.
"""
function turning_points(x::AbstractVector{<:Real})
    n = length(x)
    n <= 2 && return collect(float.(x))
    tp = Float64[float(x[1])]
    for i in 2:(n - 1)
        d1 = x[i] - x[i - 1]
        d2 = x[i + 1] - x[i]
        if d1 * d2 < 0
            push!(tp, float(x[i]))
        end
    end
    push!(tp, float(x[n]))
    return tp
end

"""
    rainflow_count(x) -> Vector{@NamedTuple{range::Float64, mean::Float64, count::Float64}}

Rainflow cycle count of `x` by the three-point method of ASTM E1049-85.

The series is reduced to its turning points and pushed onto a stack. Whenever the last
three points `a, b, c` satisfy `|b - a| <= |c - b|`, the inner range `(a, b)` is a closed
cycle: it is recorded with `count = 1.0` and both points are removed, leaving `c` to be
tested against what was under them. What cannot be closed is the residue, reported as
half cycles (`count = 0.5`) — the standard convention, and what lets the counted ranges
account for the whole signal rather than silently dropping its largest excursion.
"""
function rainflow_count(x::AbstractVector{<:Real})
    cycles = @NamedTuple{range::Float64, mean::Float64, count::Float64}[]
    stack = Float64[]
    for p in turning_points(x)
        push!(stack, p)
        while length(stack) >= 3
            a, b, c = stack[end - 2], stack[end - 1], stack[end]
            if abs(b - a) <= abs(c - b)
                push!(cycles, (range = abs(b - a), mean = (a + b) / 2, count = 1.0))
                deleteat!(stack, (length(stack) - 2, length(stack) - 1))
            else
                break
            end
        end
    end
    for i in 1:(length(stack) - 1)
        push!(cycles, (range = abs(stack[i + 1] - stack[i]),
                       mean = (stack[i] + stack[i + 1]) / 2,
                       count = 0.5))
    end
    return cycles
end

"""
    rainflow_matrix(cycles, range_edges, mean_edges) -> Matrix{Float64}

Bin `cycles` into the range × mean histogram `RainflowCounter` reports, counting each
cycle by its `count` so half cycles contribute a half. Values outside the edges are
clamped into the end bins; `range_edges` and `mean_edges` must be sorted and have at
least two entries each.
"""
function rainflow_matrix(cycles, range_edges::AbstractVector{<:Real},
                         mean_edges::AbstractVector{<:Real})
    nr = length(range_edges) - 1
    nm = length(mean_edges) - 1
    (nr >= 1 && nm >= 1) || throw(ArgumentError("need at least two edges in each direction"))
    M = zeros(Float64, nr, nm)
    bin(v, edges) = clamp(searchsortedlast(edges, v), 1, length(edges) - 1)
    for c in cycles
        M[bin(c.range, range_edges), bin(c.mean, mean_edges)] += c.count
    end
    return M
end

"""
    miner_damage(cycles; C, m) -> Float64

Palmgren-Miner damage sum for `cycles` against the S-N line `N = C * S^(-m)`, the
accumulation `FatigueCounter` performs. Failure is conventionally taken at 1.0.
Zero-range cycles contribute nothing.
"""
function miner_damage(cycles; C::Real, m::Real)
    d = 0.0
    for c in cycles
        c.range > 0 || continue
        d += c.count / (C * c.range^(-m))
    end
    return d
end
