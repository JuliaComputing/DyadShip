#!/usr/bin/env julia
# Focused checks of the SimpleDieselEngine fuel metering: sign and operating-mode
# boundary, SFOC table-coverage flag and interval counter, idle fuel parameter, exact
# rate/counter consistency, counter sampling (accepted steps, dense output, `saveat`),
# rejected parameter values, and a regression of the forward ramp against the values the
# engine produced before fuel was restricted to the brake power (identical by
# construction for forward operation).
#
# Runs the four engine analyses only; never the generated test suite.
#
# Usage:  julia --project=<environment with DyadShip> scripts/validate_engine_fuel.jl

using DyadShip, Printf, Test
using DyadInterface: symbolic_container
using ModelingToolkit: assertions

const PR = DyadShip.Propulsion

# Avoids a direct SciMLBase dependency in the running environment.
successful_retcode(sol) = string(sol.retcode) == "Success"

# Trapezoid integral of a sampled rate, the integration policy used to difference the
# cumulative counter against the rate samples.
trapz(ts, ys) = sum((ts[i + 1] - ts[i]) * (ys[i] + ys[i + 1]) / 2 for i in 1:length(ts) - 1)

# Dense (interpolated) samples of one variable.
dense(sol, var, ts) = [sol(t; idxs = var) for t in ts]

# Outcome of an analysis that is expected to refuse its parameters: `:threw` when the
# model or problem cannot be built, `:failed` when the solve ends without success, and
# `:accepted` when it runs to a successful return code (the case that must not happen).
function outcome(run)
    res = try
        run()
    catch
        return :threw
    end
    return successful_retcode(res.sol) ? :accepted : :failed
end

ramp_with(; kw...) = PR.SimpleDieselEngineTransient(model = PR.DieselEngineRamp(; name = :DieselEngineRamp, kw...))

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
    F = dense(sol, e.Fuel, ts)
    @test all(diff(F) .>= 0)
    @test all(t -> sol(t; idxs = e.Inst_Fuel) >= 0, ts)
    # Interval quality. The accepted steps of the torque burst (2 s to 2.55 s) are all
    # inside the table and the extrapolation counter does not move over them.
    burst = findall(t -> 2.001 < t <= 2.55, sol.t)
    @test length(burst) >= 5
    @test all(sol[e.SFOC_valid][burst] .== 1)
    Fxb = sol[e.Fuel_extrapolated][burst]
    @test maximum(Fxb) == minimum(Fxb)
    @test sol[e.Fuel][last(burst)] > sol[e.Fuel][first(burst)]
    Fx = dense(sol, e.Fuel_extrapolated, ts)
    # At the 355 kW steady state every kilogram is metered with the held end value.
    @test sol(30.0; idxs = e.Fuel_extrapolated) - sol(20.0; idxs = e.Fuel_extrapolated) ≈
          sol(30.0; idxs = e.Fuel) - sol(20.0; idxs = e.Fuel) rtol = 1e-6
    @test 0 < Fx[end] < F[end]
    @printf("  ramp: %.1f %% of Fuel(30) metered outside the SFOC table (Fuel_extrapolated = %.5f kg)\n",
        100 * Fx[end] / F[end], Fx[end])
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
    # The analysis declares its zero-power crossings as tstops, so steps end there: the counter is nondecreasing and exactly flat over the reverse drive,
    # at the accepted steps and on the dense grid alike, with no tolerance.
    Fs = sol[e.Fuel]; rpms = sol[e.rpm]
    @test all(diff(Fs) .>= 0)
    backs = rpms .< -1
    @test count(backs) >= 3
    @test maximum(Fs[backs]) == minimum(Fs[backs])
    @test all(diff(F) .>= 0)
    @test maximum(F[back]) == minimum(F[back])
    # Accepted steps end at the two crossings (10 s and 20 s).
    @test minimum(abs.(sol.t .- 10)) < 1e-6
    @test minimum(abs.(sol.t .- 20)) < 1e-6
    # Total against a tight-tolerance run of the same analysis.
    ref = PR.SimpleDieselEngineReversingTransient(abstol = 1e-10, reltol = 1e-10)
    @test successful_retcode(ref.sol)
    Fref = ref.sol(30.0; idxs = symbolic_container(ref).engine.Fuel)
    @test F[end] ≈ Fref rtol = 1e-4
    # Interval quality: the coverage flag is 1 at t = 12 (driven, no fuel) and at
    # t = 28 (firing inside the table), yet the engine crossed the low-load range in
    # between. The extrapolation counter records it; the endpoint flags cannot.
    # The table-end crossings are not step-aligned, so the extrapolation counter is
    # asserted at accepted steps and its dense fall-back is reported, not bounded.
    Fx = dense(sol, e.Fuel_extrapolated, ts)
    Fxs = sol[e.Fuel_extrapolated]
    @test all(diff(Fxs) .>= 0)
    @test maximum(Fxs[backs]) == minimum(Fxs[backs])
    @printf("  reversing: Fuel_extrapolated dense grid %d decreases (largest fall-back %.3e kg)\n",
        count(diff(Fx) .< 0), maximum(accumulate(max, Fx) .- Fx))
    @test sol(12.0; idxs = e.SFOC_valid) == 1 && sol(28.0; idxs = e.SFOC_valid) == 1
    @test sol(28.0; idxs = e.Fuel_extrapolated) > sol(12.0; idxs = e.Fuel_extrapolated)
    @test Fx[end] ≈ trapz(ts, mdot .* (1 .- q)) rtol = 2e-2
    @test all(Fxs .<= Fs)
    @printf("  reversing: %d accepted steps, Fuel(30) = %.8f kg (tight-tolerance run %.8f), Fuel_extrapolated(30) = %.5f kg\n",
        length(sol.t), F[end], Fref, Fx[end])
    # The counter is exactly the integral of the rate.
    @test F[end] ≈ trapz(ts, mdot) rtol = 2e-3
    @test F[end] > 0
    @printf("  reversing: rpm %.0f .. %.0f, Fuel(30) = %.5f kg, counter flat for %.2f s of reverse drive\n",
        minimum(rpm), maximum(rpm), F[end], count(back) * 0.02)
end

@testset "counter sampling: accepted steps, dense output and saveat" begin
    # Same analysis without its tstops, to show what they are for. The accepted steps
    # stay nondecreasing; the interpolant does not.
    ts = collect(0.0:0.02:30.0)
    res = PR.SimpleDieselEngineReversingTransient(tstops = Float64[])
    sol = res.sol; e = symbolic_container(res).engine
    @test successful_retcode(sol)
    @test all(diff(sol[e.Fuel]) .>= 0)
    F = dense(sol, e.Fuel, ts)
    # `saveat` stores the interpolant of the same accepted steps; it does not make the
    # solver step to the requested times.
    ressa = PR.SimpleDieselEngineReversingTransient(tstops = Float64[], saveat = 0.02)
    solsa = ressa.sol; esa = symbolic_container(ressa).engine
    @test successful_retcode(solsa)
    @test length(solsa.t) == length(ts)
    @test length(sol.t) < length(ts) ÷ 10
    Fsa = solsa[esa.Fuel]
    @test maximum(abs.(Fsa .- F)) <= 1e-12
    @test count(diff(Fsa) .< 0) == count(diff(F) .< 0)
    @printf("  no tstops: %d accepted steps with %d decreases; dense grid %d decreases (largest fall-back %.3e kg); saveat grid %d decreases\n",
        length(sol.t), count(diff(sol[e.Fuel]) .< 0), count(diff(F) .< 0),
        maximum(accumulate(max, F) .- F), count(diff(Fsa) .< 0))
    # With the tstops (the analysis default) a saveat grid is nondecreasing too.
    reson = PR.SimpleDieselEngineReversingTransient(saveat = 0.02)
    Fon = reson.sol[symbolic_container(reson).engine.Fuel]
    @test successful_retcode(reson.sol)
    @test all(diff(Fon) .>= 0)
    @printf("  with tstops: saveat grid %d decreases\n", count(diff(Fon) .< 0))
end

@testset "whole run inside the SFOC table" begin
    res = PR.SimpleDieselEngineInRangeTransient(); sol = res.sol; e = symbolic_container(res).engine
    @test successful_retcode(sol)
    # Everything below is read at the accepted steps.
    ts = sol.t; P = sol[e.ShaftPower] ./ 1000; mdot = sol[e.Inst_Fuel]; F = sol[e.Fuel]
    # Table coverage for the whole interval: flag, power range, and the interval counter.
    @test all(sol[e.SFOC_valid] .== 1)
    @test all(605 .<= P .<= 1210)
    @test all(sol[e.Fuel_extrapolated] .== 0)
    @test all(diff(F) .>= 0)
    @test F[1] == 0
    # Steady state: 1800 rpm, rate = SFOC · kW / 3.6e6 with the SFOC inside the table.
    @test sol[e.rpm][end] ≈ 1800 rtol = 1e-6
    @test 178 < sol[e.sfoc][end] < 185
    @test mdot[end] ≈ sol[e.sfoc][end] * P[end] / 3.6e6 rtol = 1e-12
    # Counter against a tight-tolerance run, at the tolerance the analysis requests.
    ref = PR.SimpleDieselEngineInRangeTransient(abstol = 1e-10, reltol = 1e-10)
    @test successful_retcode(ref.sol)
    Fref = ref.sol[symbolic_container(ref).engine.Fuel][end]
    @test F[end] ≈ Fref rtol = 1e-6
    # Rate against counter. The trapezoid of the rate samples differs from the counter
    # by the quadrature error of the trapezoid rule on these steps, which is estimated
    # by comparing it with Simpson's rule on the same steps (midpoint rate from the dense
    # output). The gap must stay within twice that estimate plus the solver tolerance on
    # the counter; no bound is fitted to the result.
    T = trapz(ts, mdot)
    S = sum((ts[i + 1] - ts[i]) * (mdot[i] + 4 * sol((ts[i] + ts[i + 1]) / 2; idxs = e.Inst_Fuel) + mdot[i + 1]) / 6
            for i in 1:length(ts) - 1)
    @test abs(T - F[end]) <= 2 * abs(T - S) + 1e-6 * F[end]
    @test abs(S - F[end]) <= abs(T - F[end])
    @printf("  in range: %d accepted steps, power %.1f .. %.1f kW, Fuel(30) = %.8f kg (tight-tolerance run %.8f)\n",
        length(ts), minimum(P), maximum(P), F[end], Fref)
    @printf("  in range: trapezoid - counter = %.3e kg, Simpson - counter = %.3e kg, quadrature estimate |T - S| = %.3e kg\n",
        T - F[end], S - F[end], abs(T - S))
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
        # Coverage flag only: it is 1 although the idle rate is a synthetic number.
        @test sol(t; idxs = e.SFOC_valid) == 1
        @test sol(t; idxs = e.Fuel_extrapolated) == 0
    end
    @printf("  idle: Fuel(100) = %.4f kg at m_dot_idle = 0.002 kg/s\n", sol(100.0; idxs = e.Fuel))
end

@testset "stopped shaft without idle fuel (zero case)" begin
    res = PR.SimpleDieselEngineIdleTransient(model = PR.DieselEngineIdle(; name = :DieselEngineIdle, engine__m_dot_idle = 0.0))
    sol = res.sol; e = symbolic_container(res).engine
    @test successful_retcode(sol)
    @test all(sol[e.Inst_Fuel] .== 0)
    @test all(sol[e.Fuel] .== 0)
    @test all(sol[e.Fuel_extrapolated] .== 0)
    @test all(sol[e.ShaftPower] .== 0)
    @test all(sol[e.SFOC_valid] .== 1)
end

@testset "default idle fuel is zero" begin
    res = PR.SimpleDieselEngineTransient(); sol = res.sol; m = symbolic_container(res)
    @test sol(0.0; idxs = m.engine.Inst_Fuel) == 0
    @test sol(0.0; idxs = m.engine.Fuel) == 0
end

@testset "parameter checks" begin
    # The plain-Julia checks.
    @test PR.table_knots_valid([605.0, 907.5, 1028.5, 1210.0]) == 1
    @test PR.table_knots_valid([0.0, 1210.0]) == 1
    @test PR.table_knots_valid([605.0, 605.0, 1028.5, 1210.0]) == 0      # duplicate knot
    @test PR.table_knots_valid([1210.0, 1028.5, 907.5, 605.0]) == 0      # decreasing
    @test PR.table_knots_valid([-1.0, 907.5, 1028.5, 1210.0]) == 0       # negative power
    @test PR.table_knots_valid([605.0, NaN, 1028.5, 1210.0]) == 0
    @test PR.table_knots_valid([605.0, 907.5, 1028.5, Inf]) == 0
    @test PR.table_values_valid([185.0, 179.0, 178.0, 182.0]) == 1
    @test PR.table_values_valid([185.0, 0.0, 178.0, 182.0]) == 0
    @test PR.table_values_valid([185.0, -1.0, 178.0, 182.0]) == 0
    @test PR.table_values_valid([185.0, NaN, 178.0, 182.0]) == 0
    @test PR.table_values_valid([185.0, Inf, 178.0, 182.0]) == 0
    @test PR.rate_valid(0.0) == 1 && PR.rate_valid(0.002) == 1
    @test PR.rate_valid(-0.001) == 0 && PR.rate_valid(NaN) == 0 && PR.rate_valid(Inf) == 0

    # The component carries the three assertions.
    eng = PR.SimpleDieselEngine(; name = :engine)
    @test length(assertions(eng)) == 3

    # Accepted: the default table passed explicitly, and a longer table with its length.
    base = PR.SimpleDieselEngineTransient()
    Fbase = base.sol(30.0; idxs = symbolic_container(base).engine.Fuel)
    same = ramp_with(engine__SFOC_P = [605.0, 907.5, 1028.5, 1210.0], engine__SFOC_g = [185.0, 179.0, 178.0, 182.0])
    @test successful_retcode(same.sol)
    @test same.sol(30.0; idxs = symbolic_container(same).engine.Fuel) == Fbase
    # A table of another length is declared with its length (construction only: the
    # structural `n_sfoc` cannot be overridden through an enclosing test component).
    five = PR.SimpleDieselEngine(; name = :engine, n_sfoc = 5,
        SFOC_P = [300.0, 605.0, 907.5, 1028.5, 1210.0], SFOC_g = [185.0, 185.0, 179.0, 178.0, 182.0])
    @test length(assertions(five)) == 3

    # The metering does not depend on event detection: switching it on changes nothing
    # beyond solver accuracy, although the engine starts at exactly zero power.
    lifted = PR.SimpleDieselEngineTransient(automatic_discontinuity_detection = true)
    @test successful_retcode(lifted.sol)
    el = symbolic_container(lifted).engine
    @test lifted.sol(30.0; idxs = el.Fuel) ≈ Fbase rtol = 1e-5
    @test lifted.sol(30.0; idxs = el.Fuel_extrapolated) ≈ base.sol(30.0; idxs = symbolic_container(base).engine.Fuel_extrapolated) rtol = 1e-5
    @test lifted.sol[el.SFOC_valid][end] == 0

    # Rejected: every case must end without a successful return code or refuse to build.
    bad = [
        "duplicate SFOC_P knot"      => (; engine__SFOC_P = [605.0, 605.0, 1028.5, 1210.0]),
        "decreasing SFOC_P"          => (; engine__SFOC_P = [1210.0, 1028.5, 907.5, 605.0]),
        "negative SFOC_P knot"       => (; engine__SFOC_P = [-1.0, 907.5, 1028.5, 1210.0]),
        "NaN in SFOC_P"              => (; engine__SFOC_P = [605.0, NaN, 1028.5, 1210.0]),
        "infinite SFOC_P knot"       => (; engine__SFOC_P = [605.0, 907.5, 1028.5, Inf]),
        "negative SFOC_g"            => (; engine__SFOC_g = [-1.0, -1.0, -1.0, -1.0]),
        "zero SFOC_g"                => (; engine__SFOC_g = [185.0, 0.0, 178.0, 182.0]),
        "NaN in SFOC_g"              => (; engine__SFOC_g = [185.0, NaN, 178.0, 182.0]),
        "infinite SFOC_g"            => (; engine__SFOC_g = [185.0, Inf, 178.0, 182.0]),
        "negative m_dot_idle"        => (; engine__m_dot_idle = -0.001),
        "NaN m_dot_idle"             => (; engine__m_dot_idle = NaN),
        "5 knots with n_sfoc = 4"    => (; engine__SFOC_P = [300.0, 605.0, 907.5, 1028.5, 1210.0],
                                           engine__SFOC_g = [185.0, 185.0, 179.0, 178.0, 182.0]),
        "SFOC_P longer than SFOC_g"  => (; engine__SFOC_P = [300.0, 605.0, 907.5, 1028.5, 1210.0]),
    ]
    for (label, kw) in bad
        o = outcome(() -> ramp_with(; kw...))
        @printf("  rejected input, %-26s: %s\n", label, o)
        @test o != :accepted
    end
end

end
