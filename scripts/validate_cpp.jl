#!/usr/bin/env julia
# Controllable pitch propeller: the fixed-pitch propeller it shares its equations with is
# unchanged; at three pitch settings and a constant shaft speed the thrust, torque and
# power are those of the Wageningen B-series regression at the pitch in use; the ship
# slows as the pitch is reduced. Also: the KVLCC2 vessel held as a subcomponent.
#
# Runs three analyses only; never the generated test suite.
#
# Usage:  julia --project=<environment with DyadShip> scripts/validate_cpp.jl

using DyadShip, Printf, Test
using DyadInterface: symbolic_container

const S = DyadShip.Ship6DOF
successful_retcode(sol) = string(sol.retcode) == "Success"

@testset "controllable pitch propeller" begin

trial = S.SpeedTrialTransient(); tm = symbolic_container(trial)
@testset "fixed-pitch propeller unchanged by the split into a base and two variants" begin
    @test successful_retcode(trial.sol)
    @test trial.sol[tm.ship.Surge][end] ≈ 6.053806940 rtol = 1e-9
    @test trial.sol[tm.prop.Thrust][end] ≈ 129102.465380 rtol = 1e-9
    @test trial.sol[tm.prop.ShaftPower][end] ≈ 1110404.213262 rtol = 1e-9
end

@testset "three pitch settings at a constant shaft speed" begin
    res = S.CPPShipPitchStepsTransient(); sol = res.sol; m = symbolic_container(res); p = m.prop
    @test successful_retcode(sol)
    D = sol.ps[p.Diameter]; rho = sol.ps[p.SeaDensity]; td = sol.ps[p.ThrustDeduction]
    rows = []
    for (t, pd) in ((599.0, 1.0), (1199.0, 0.8), (1799.0, 0.6))
        at(x) = sol(t; idxs = x)
        n = at(p.rpm) / 60
        @test at(p.Pitch_ratio) ≈ pd rtol = 1e-6            # the blades have reached the order
        @test at(p.rpm) ≈ 100 rtol = 1e-6                   # the shaft speed is held
        # Open-water coefficients: the series regression at the pitch in use.
        kt = S.wageningen_kt(at(p.J), at(p.Pitch_ratio), 0.55, 4)
        kq = S.wageningen_kq(at(p.J), at(p.Pitch_ratio), 0.55, 4)
        # (values between solver steps are interpolated one by one, hence 1e-6)
        @test at(p.Kt) ≈ kt rtol = 1e-6
        @test at(p.Kq) ≈ kq rtol = 1e-6
        # Thrust, torque and power from them.
        torque = kq * rho * n^2 * D^5
        @test at(p.Thrust) ≈ kt * rho * n^2 * D^4 * (1 - td) rtol = 1e-6
        @test at(p.ShaftPower) ≈ 2pi * n * torque rtol = 1e-6
        # The ship has settled, and thrust balances resistance.
        @test abs(at(m.ship.Surge) - sol(t - 100; idxs = m.ship.Surge)) < 2e-2
        @test at(p.Thrust) ≈ -at(m.hydro.Resistance) rtol = 2e-2
        push!(rows, (pd = pd, u = at(m.ship.Surge), J = at(p.J), thrust = at(p.Thrust), torque = torque, power = at(p.ShaftPower)))
    end
    # At the design pitch the ship is the fixed-pitch ship.
    @test rows[1].u ≈ trial.sol[tm.ship.Surge][end] rtol = 1e-3
    @test rows[1].power ≈ trial.sol[tm.prop.ShaftPower][end] rtol = 1e-3
    # Less pitch: less speed, thrust, torque and power, at the same shaft speed.
    for f in (:u, :thrust, :torque, :power)
        @test getproperty(rows[1], f) > getproperty(rows[2], f) > getproperty(rows[3], f)
    end
    # The pitch ratio stays inside the range of the regression.
    @test all(0.5 .<= sol[p.Pitch_ratio] .<= 1.4)
    @printf("  %6s %10s %8s %11s %12s %10s\n", "P/D", "speed m/s", "J", "thrust kN", "torque kNm", "power kW")
    for r in rows
        @printf("  %6.2f %10.3f %8.3f %11.1f %12.1f %10.1f\n", r.pd, r.u, r.J, r.thrust / 1e3, r.torque / 1e3, r.power / 1e3)
    end
end

@testset "KVLCC2 vessel held as a subcomponent" begin
    res = S.KVLCC2VesselWithPlantTransient(); sol = res.sol; m = symbolic_container(res)
    @test successful_retcode(sol)
    @test sol[m.vessel.prop.rpm][end] ≈ 73.7 rtol = 1e-6
    @test sol[m.plant.shaft_power][end] ≈ sol[m.vessel.prop.ShaftPower][end] rtol = 1e-6
    @test sol[m.vessel.speed][end] == sol[m.vessel.ship.Surge][end]
    @test sol[m.vessel.speed][end] ≈ 15.5 * 1852 / 3600 rtol = 1e-2
    @test maximum(abs.(sol[m.vessel.ship.Yaw])) < 1e-2
    @printf("  vessel as subcomponent: %.2f kn at %.1f rev/min, %.2f MW from the plant\n", sol[m.vessel.speed][end] * 3600 / 1852, sol[m.vessel.prop.rpm][end], sol[m.plant.shaft_power][end] / 1e6)
end

end
