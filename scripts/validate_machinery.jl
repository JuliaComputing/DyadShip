#!/usr/bin/env julia
# Check the electrical power load analysis: independently scheduled consumers on seeded
# clocked random starts, and the bus load they add up to.
#
# Kept out of validate_heattransfer.jl because six periodic clocks over six hours is a
# few minutes of solve on its own.
#
# Usage:  ~/dyad-fleet/heavy -m 8G ../julia-dyad.sh scripts/validate_machinery.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using DyadShip, Printf, Test
using SciMLBase: successful_retcode
using DyadInterface: symbolic_container

const MC = DyadShip.Machinery

@testset "Machinery" begin

@testset "electrical load balance" begin
    res = MC.ShipLoadBalanceTransient(); sol = res.sol; m = symbolic_container(res)
    @test successful_retcode(sol)
    ts = collect(0.0:60.0:21600.0)
    load = [sol(t; idxs = m.bank.y) for t in ts]
    avg  = [sol(t; idxs = m.bank.average) for t in ts]
    @test all(isfinite, vcat(load, avg))

    installed = 37000 + 15000 + 11000 + 6600 + 5500 + 4000
    Kr = 0.85

    # Every consumer starts off, so the bus starts unloaded.
    @test load[1] == 0
    @test all(load .>= 0)
    # and no instant can exceed every consumer running at once.
    @test maximum(load) <= installed * Kr * 1.001
    @test maximum(load) > 0.3 * installed

    # Each consumer follows its own duty: the mean power is its nominal times Kr times
    # the fraction of the run it should be running for.
    duties = [0.3, 0.4, 0.6, 0.5, 0.2, 0.7]
    powers = [37000, 15000, 11000, 6600, 5500, 4000]
    for i in 1:6
        yi = sol[getproperty(m.bank, Symbol("loads⸺$i")).y]
        @test maximum(yi) ≈ powers[i] * Kr rtol=1e-6
        @test 0.5 * duties[i] < (sum(yi) / length(yi)) / (powers[i] * Kr) < 1.5 * duties[i]
    end

    # The consumers are on independent seeds, so the bus takes many distinct levels
    # rather than switching as a block.
    levels = length(unique(round.(load; digits = 0)))
    @test levels > 50
    @printf("  bus load 0 .. %.0f W over %d distinct levels\n", maximum(load), levels)

    # The running average is what an EPLA is actually asked for; it should settle well
    # before the end of the run rather than still drifting.
    half = avg[findfirst(t -> t >= 10800, ts)]
    @test abs(avg[end] - half) < 0.1 * avg[end]
    diversity = avg[end] / installed
    @test 0.15 < diversity < 0.6
    @printf("  average %.0f W against %d W installed: diversity factor %.3f\n",
            avg[end], installed, diversity)

    # The schedulers are clocked, so their own signal is readable one entry per tick even
    # where the held continuous output is not.
    sched1 = getproperty(m.bank, Symbol("scheds⸺1"))
    w = sol[sched1.logic.Work]
    @test all(x -> x ≈ 0 || x ≈ 1, w)
    @test 0 < sum(w) < length(w)
    starts = count(i -> w[i] > 0.5 && w[i-1] < 0.5, 2:length(w))
    @printf("  consumer 1: %d starts in 6 h on a nominal 3/h (%d ticks)\n", starts, length(w))
    @test starts >= 2
end


@testset "rainflow counting" begin
    # A nested sequence that returns to where it started: two small reversals close as
    # full cycles of range 2, and the outer swing closes as well, because the signal comes
    # back to 0. Nothing is left over.
    x = [0.0, 3.0, 1.0, 3.0, 1.0, 5.0, 0.0]
    cyc = DyadShip.rainflow_count(x)
    full = filter(c -> c.count == 1.0, cyc)
    @test length(full) == 3
    @test length(filter(c -> c.count == 0.5, cyc)) == 0
    @test sort([c.range for c in full]) ≈ [2.0, 2.0, 5.0]
    @test all(c -> c.mean ≈ 2.0, filter(c -> c.range ≈ 2.0, full))
    @test only(filter(c -> c.range ≈ 5.0, full)).mean ≈ 2.5
    @printf("  nested sequence: ranges %s, no residue\n",
            string(sort([c.range for c in full])))

    # Turning-point reduction drops everything that is not a reversal and leaves the count
    # unchanged.
    dense = [0.0, 1.0, 2.0, 3.0, 2.0, 1.0, 3.0, 2.0, 1.0, 3.0, 4.0, 5.0, 2.0, 0.0]
    @test DyadShip.turning_points(dense) == x
    @test DyadShip.rainflow_count(dense) == cyc

    # A sequence that does not return closes nothing: all three ranges are residue half
    # cycles, which is how the largest excursion stays counted.
    res = DyadShip.rainflow_count([0.0, 5.0, 1.0, 4.0])
    @test all(c -> c.count == 0.5, res)
    @test sort([c.range for c in res]) ≈ [3.0, 4.0, 5.0]

    # No counted range can exceed the signal's own peak-to-peak.
    @test maximum(c -> c.range, cyc) ≈ maximum(x) - minimum(x)

    # Four periods of a sine: three closed cycles at the full peak-to-peak, plus the
    # opening quarter that pairs with the closing one.
    t = range(0, 4; length = 4001)
    cs = DyadShip.rainflow_count(10 .* sin.(2π .* t))
    big = filter(c -> c.count == 1.0 && c.range > 15, cs)
    @test length(big) == 3
    @test all(c -> isapprox(c.range, 20.0; rtol = 1e-3), big)
    @test all(c -> isapprox(c.mean, 0.0; atol = 1e-6), big)
    @test sum(c -> c.count, cs) ≈ 4.5
    @printf("  four periods of a 10-amplitude sine: %d closed cycles of range %.3f\n",
            length(big), big[1].range)

    # The matrix bins by range and mean and conserves the total count.
    M = DyadShip.rainflow_matrix(cyc, [0.0, 1.0, 3.0, 6.0], [-1.0, 1.0, 3.0])
    @test sum(M) ≈ sum(c -> c.count, cyc)
    @test M[2, 2] ≈ 2.0          # the two range-2, mean-2 cycles
    @test M[3, 2] ≈ 1.0          # the range-5, mean-2.5 cycle

    # Miner damage: n cycles at one range against N = C·S^-m is exactly n/N, and a
    # zero-range cycle does no damage.
    @test DyadShip.miner_damage([(range = 4.0, mean = 0.0, count = 3.0)]; C = 1e12, m = 3) ≈
          3.0 / (1e12 * 4.0^-3)
    @test DyadShip.miner_damage([(range = 0.0, mean = 0.0, count = 5.0)]; C = 1e12, m = 3) == 0.0
    @printf("  Miner damage of the nested sequence at C=1e12, m=3: %.4e\n",
            DyadShip.miner_damage(cyc; C = 1e12, m = 3))
end
end
