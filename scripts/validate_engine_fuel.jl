#!/usr/bin/env julia
# Focused checks of the SimpleDieselEngine fuel metering: sign and operating-mode
# boundary, SFOC validity flag, idle fuel parameter, exact rate/counter consistency, and
# a regression of the forward ramp against the values the engine produced before fuel
# was restricted to the brake power (identical by construction for forward operation).
#
# Runs the three engine analyses only; never the generated test suite.
#
# Usage:  julia --project=<environment with DyadShip> scripts/validate_engine_fuel.jl

using DyadShip, Printf, Test
using DyadInterface: symbolic_container

const PR = DyadShip.Propulsion

# Avoids a direct SciMLBase dependency in the running environment.
successful_retcode(sol) = string(sol.retcode) == "Success"

# Trapezoid integral of a sampled rate, the integration policy used to difference the
# cumulative counter against the rate samples.
trapz(ts, ys) = sum((ts[i + 1] - ts[i]) * (ys[i] + ys[i + 1]) / 2 for i in 1:length(ts) - 1)

@testset "SimpleDieselEngine fuel metering" begin

@testset "forward ramp (regression of normal operation)" begin
    res = PR.SimpleDieselEngineTransient(); sol = res.sol; m = symbolic_container(res)
    @test successful_retcode(sol)
    e = m.engine
    # Values recorded from the previous engine (signed-power fuel) at the same
    # operating points. Forward operation must be unchanged.
    reference = [
        # t      rpm         P_W           inst_fuel        Fuel          KWh
        (10.0, 1800.0000, 355305.724, 1.82587664e-02, 0.18444315, 1.00222023),
        (20.0, 1800.0000, 355305.758, 1.82587681e-02, 0.36703083, 1.98918067),
        (30.0, 1800.0000, 355305.758, 1.82587681e-02, 0.54961851, 2.97614111),
    ]
    for (t, rpm, P, mdot, F, E) in reference
        @test sol(t; idxs = e.rpm) ≈ rpm rtol = 1e-4
        @test sol(t; idxs = e.ShaftPower) ≈ P rtol = 1e-4
        @test sol(t; idxs = e.Inst_Fuel) ≈ mdot rtol = 1e-4
        @test sol(t; idxs = e.Fuel) ≈ F rtol = 1e-3
        @test sol(t; idxs = e.KWh) ≈ E rtol = 1e-3
    end
    # Units: kg/s = g/kWh · kW / 3.6e6 at the steady state, with the SFOC held at the
    # table's first value because 355 kW is below its 605 kW lower knot.
    P_kW = sol(30.0; idxs = e.ShaftPower) / 1000
    @test sol(30.0; idxs = e.sfoc) ≈ 185.0
    @test sol(30.0; idxs = e.Inst_Fuel) ≈ 185.0 * P_kW / 3.6e6 rtol = 1e-9
    @test sol(30.0; idxs = e.SFOC_valid) == 0
    # The table range is entered only during the 2 s torque burst (864 kW at t = 2).
    @test sol(2.0; idxs = e.SFOC_valid) == 1
    @test 178 < sol(2.0; idxs = e.sfoc) < 185
    ts = collect(0.0:0.05:30.0)
    F = [sol(t; idxs = e.Fuel) for t in ts]
    @test all(diff(F) .>= 0)
    @test all(t -> sol(t; idxs = e.Inst_Fuel) >= 0, ts)
    @printf("  ramp: Fuel(30) = %.5f kg, KWh(30) = %.5f, SFOC held at %.0f g/kWh\n",
        F[end], sol(30.0; idxs = e.KWh), sol(30.0; idxs = e.sfoc))
end

@testset "reversing shaft (sign and mode boundary)" begin
    res = PR.SimpleDieselEngineReversingTransient(); sol = res.sol; m = symbolic_container(res)
    @test successful_retcode(sol)
    e = m.engine
    ts = collect(0.0:0.02:30.0)
    rpm  = [sol(t; idxs = e.rpm) for t in ts]
    P    = [sol(t; idxs = e.ShaftPower) for t in ts]
    mdot = [sol(t; idxs = e.Inst_Fuel) for t in ts]
    F    = [sol(t; idxs = e.Fuel) for t in ts]
    E    = [sol(t; idxs = e.KWh) for t in ts]
    tau  = [sol(t; idxs = e.tau_cmd) for t in ts]
    q    = [sol(t; idxs = e.SFOC_valid) for t in ts]
    @test minimum(rpm) < -400 && maximum(rpm) > 1400
    # Engine torque is never negative; the engine never generates.
    @test all(tau .>= 0)
    # Driven backwards: negative shaft power, zero fuel rate.
    back = rpm .< -1
    @test count(back) > 100
    @test all(P[back] .< 0)
    @test all(mdot[back] .== 0)
    # Signed net work does decrease while driven, as documented.
    @test E[findlast(back)] < E[findfirst(back)]
    # Firing: positive power, positive fuel.
    fwd = (rpm .> 1) .& (tau .> 1)
    @test all(P[fwd] .> 0)
    @test all(mdot[fwd] .> 0)
    # Everywhere: rate nonnegative, quality flag is 0/1 and 1 whenever no
    # load-dependent fuel is metered.
    @test all(mdot .>= 0)
    @test all(x -> x == 0 || x == 1, q)
    @test all(q[back] .== 1)
    # The counter is nondecreasing at every solver step and exactly flat over the
    # reverse-drive steps. The dense interpolant inside the one step that straddles the
    # zero-power crossing can overshoot by the solver tolerance; bound that wobble.
    Fs = sol[e.Fuel]; rpms = sol[e.rpm]
    @test all(diff(Fs) .>= 0)
    backs = rpms .< -1
    @test count(backs) >= 3
    @test maximum(Fs[backs]) == minimum(Fs[backs])
    @test all(diff(F) .>= -1e-4 * F[end])
    @test maximum(F[back]) - minimum(F[back]) <= 1e-3 * F[end]
    # The counter is exactly the integral of the rate.
    @test F[end] ≈ trapz(ts, mdot) rtol = 2e-3
    @test F[end] > 0
    @printf("  reversing: rpm %.0f .. %.0f, Fuel(30) = %.5f kg, counter flat for %.2f s of reverse drive\n",
        minimum(rpm), maximum(rpm), F[end], count(back) * 0.02)
end

@testset "stopped shaft with idle fuel" begin
    res = PR.SimpleDieselEngineIdleTransient(); sol = res.sol; m = symbolic_container(res)
    @test successful_retcode(sol)
    e = m.engine
    for t in (0.0, 25.0, 50.0, 100.0)
        @test abs(sol(t; idxs = e.rpm)) < 1e-9
        @test sol(t; idxs = e.ShaftPower) == 0
        @test sol(t; idxs = e.KWh) == 0
        @test sol(t; idxs = e.Inst_Fuel) ≈ 0.002 rtol = 1e-12
        @test sol(t; idxs = e.Fuel) ≈ 0.002 * t rtol = 1e-9 atol = 1e-12
        @test sol(t; idxs = e.SFOC_valid) == 1
    end
    @printf("  idle: Fuel(100) = %.4f kg at m_dot_idle = 0.002 kg/s\n", sol(100.0; idxs = e.Fuel))
end

@testset "default idle fuel is zero" begin
    res = PR.SimpleDieselEngineTransient(); sol = res.sol; m = symbolic_container(res)
    @test sol(0.0; idxs = m.engine.Inst_Fuel) == 0
    @test sol(0.0; idxs = m.engine.Fuel) == 0
end

end
