#!/usr/bin/env julia
# KVLCC2 at full scale: the hull preset reaches the solved model unchanged, the ship holds
# its approach speed, and the 35° turning circles are compared with the values printed in
# Yasukawa and Yoshimura, J. Mar. Sci. Technol. 20 (2015) 37-52, Table 6 (their full-scale
# simulation) and Table 4 (their free-running test of a 7 m model). Nothing is tuned: the
# deviations are printed, and only a coarse 25 % band is asserted.
#
# Also: the powered variant with a plant in the plant slot, and the powering table from the
# published open-water curves with the full-scale resistance estimate.
#
# Runs a handful of analyses only; never the generated test suite.
#
# Usage:  julia --project=<environment with DyadShip> scripts/validate_kvlcc2.jl

using DyadShip, Printf, Test, TOML
using DyadInterface: symbolic_container

const S = DyadShip.Ship6DOF
const LPP = 320.0
successful_retcode(sol) = string(sol.retcode) == "Success"

# Published turning indices in ship lengths. The source's rudder angle is positive to
# starboard, this library's to port, so the source's +35° is a starboard turn.
const PUBLISHED = Dict(
    :starboard => (full_scale_simulation = (AD = 3.62, DT = 3.71), model_simulation = (AD = 3.31, DT = 3.36), model_test = (AD = 3.25, DT = 3.34)),
    :port      => (full_scale_simulation = (AD = 3.56, DT = 3.59), model_simulation = (AD = 3.26, DT = 3.26), model_test = (AD = 3.11, DT = 3.08)),
)

function turning_circle(rudder_deg; kw...)
    res = S.KVLCC2TurningCircleTransient(model = S.KVLCC2TurningCircle(; name = :KVLCC2TurningCircle, rudder_deg = rudder_deg, kw...))
    sol = res.sol; m = symbolic_container(res)
    t = sol.t; x = sol[m.ship.pos_x]; y = sol[m.ship.pos_y]; psi = sol[m.ship.Yaw]
    i0 = findfirst(>=(100.0), t)
    turned = abs.(psi .- psi[i0])
    function at(angle)                     # track position, relative to rudder execute, at a change of heading
        k = findfirst(i -> i > i0 && turned[i] >= angle, eachindex(t))
        f = (angle - turned[k - 1]) / (turned[k] - turned[k - 1])
        return (x[k - 1] + f * (x[k] - x[k - 1]) - x[i0], y[k - 1] + f * (y[k] - y[k - 1]) - y[i0])
    end
    return (; sol, m, i0, AD = at(pi / 2)[1] / LPP, DT = abs(at(pi)[2]) / LPP, side = at(pi)[2])
end

# Open-water point of the published curves for a given resistance coefficient and speed:
# the shaft speed at which (1 - t) T = R2 U², then torque and power. Plain arithmetic,
# independent of the model.
const KT = (0.2931, -0.2753, -0.1385); const KQ = (0.030710, -0.018559, -0.020451)
const DP = 9.86; const TP = 0.220; const WP = 0.35; const RHO = 1025.0
function open_water_point(R2, U)
    T = R2 * U^2 / (1 - TP)
    f(n) = (KT[1] + KT[2] * (U * (1 - WP) / (n * DP)) + KT[3] * (U * (1 - WP) / (n * DP))^2) * RHO * n^2 * DP^4 - T
    a, b = 0.2, 6.0
    for _ in 1:100
        mid = (a + b) / 2
        f(a) * f(mid) <= 0 ? (b = mid) : (a = mid)
    end
    J = U * (1 - WP) / (a * DP)
    Q = (KQ[1] + KQ[2] * J + KQ[3] * J^2) * RHO * a^2 * DP^5
    return (rpm = 60a, J = J, thrust = T, torque = Q, power = 2pi * a * Q)
end

@testset "KVLCC2" begin

@testset "the optional thrust polynomial leaves the sample ship unchanged" begin
    res = S.SpeedTrialTransient(); sol = res.sol; m = symbolic_container(res)
    @test sol[m.ship.Surge][end] ≈ 6.053806940 rtol = 1e-9
    @test sol[m.prop.Thrust][end] ≈ 129102.465380 rtol = 1e-9
    @test sol[m.prop.ShaftPower][end] ≈ 1110404.213262 rtol = 1e-9
end

port = turning_circle(35.0)
starboard = turning_circle(-35.0)

@testset "hull preset and approach" begin
    sol, m, i0 = port.sol, port.m, port.i0
    @test successful_retcode(sol)
    @test successful_retcode(starboard.sol)
    # Every value of the file is the value in the solved problem.
    file = TOML.parsefile(joinpath(pkgdir(DyadShip), "assets", "presets", "Hull", "kvlcc2_full_scale.toml"))
    @test length(file) == 27
    for (key, value) in file
        @test sol.ps[getproperty(m.hydro, Symbol(key))] == value
    end
    # Arithmetic of the converted values, from the printed ones.
    mnd = 312600 / (0.5 * 320^2 * 20.8)
    @test file["mx"] ≈ 0.022 / mnd rtol = 1e-6
    @test file["my"] ≈ 0.223 / mnd rtol = 1e-6
    @test file["Jz"] ≈ 0.011 * 0.5 * 320^4 * 20.8 / 312600 rtol = 1e-6
    @test file["R2"] ≈ 0.5 * 1025 * 320 * 20.8 * 0.022 rtol = 1e-9
    # Straight approach at 15.5 kn: speed held, thrust balances resistance, upright.
    U0 = 15.5 * 1852 / 3600
    @test sol[m.ship.Surge][i0] ≈ U0 rtol = 2e-3
    @test sol[m.prop.Thrust][i0] ≈ -sol[m.hydro.Resistance][i0] rtol = 5e-3
    # The reaction to the propeller torque gives a small heel and a slight yaw.
    @test abs(sol[m.ship.Yaw][i0]) < 1e-3
    @test sol[m.ship.Draft][i0] ≈ 20.8 rtol = 1e-3
    @test abs(sol[m.ship.Heel][i0]) < 2e-3
    @test sol[m.ship.Displacement][i0] ≈ 1025 * 312600 rtol = 1e-3
end

@testset "turning circles against the published values" begin
    # Rudder to port turns the ship to port (positive y), and conversely.
    @test port.side > 0
    @test starboard.side < 0
    # Flow straightening at the rudder in the steady turn: the source's value for a port
    # turn (its beta_R < 0) is 0.395 and for a starboard turn 0.640.
    @test port.sol[port.m.rudder.Gamma_R][end] == 0.395
    @test starboard.sol[starboard.m.rudder.Gamma_R][end] == 0.640
    # The turn costs speed and settles to a steady rate.
    @test 0.4 < port.sol[port.m.ship.Surge][end] / port.sol[port.m.ship.Surge][port.i0] < 0.7
    for (name, run) in ((:port, port), (:starboard, starboard))
        ref = PUBLISHED[name]
        for (label, r) in pairs(ref)
            @printf("  %-9s turn: advance %.3f Lpp, tactical diameter %.3f Lpp; %-21s %.2f, %.2f; deviation %+5.1f %%, %+5.1f %%\n",
                name, run.AD, run.DT, replace(string(label), "_" => " "), r.AD, r.DT, 100 * (run.AD / r.AD - 1), 100 * (run.DT / r.DT - 1))
        end
        # Coarse band against the source's own full-scale simulation; not a tuned bound.
        @test abs(run.AD / ref.full_scale_simulation.AD - 1) < 0.25
        @test abs(run.DT / ref.full_scale_simulation.DT - 1) < 0.25
    end
end


@testset "the earlier turning circles are unchanged by the hull split and the torque curve" begin
    @test port.AD ≈ 4.158 atol = 5e-4
    @test port.DT ≈ 3.983 atol = 5e-4
    @test starboard.AD ≈ 4.042 atol = 5e-4
    @test starboard.DT ≈ 3.752 atol = 5e-4
    # The propeller now loads the shaft with the published open-water torque.
    ow = open_water_point(75046.4, port.sol[port.m.ship.Surge][port.i0])
    @test port.sol[port.m.prop.ShaftPower][port.i0] ≈ ow.power rtol = 5e-3
end

@testset "powered variant: powering table with the full-scale resistance estimate" begin
    R2 = 26908.9
    @printf("  %-34s %8s %9s %8s %10s %11s %9s\n", "resistance law", "speed kn", "shaft rpm", "J", "thrust kN", "torque kNm", "power MW")
    for kn in (10.0, 12.5, 14.0, 15.5)
        U = kn * 1852 / 3600
        ow = open_water_point(R2, U)
        res = S.KVLCC2PoweringTransient(model = S.KVLCC2Powering(; name = :KVLCC2Powering, rpm = ow.rpm, U0 = U))
        sol = res.sol; m = symbolic_container(res); p = m.plant
        @test successful_retcode(sol)
        @test sol.ps[m.hydro.R2] == R2
        # The ship settles at the speed the arithmetic gives for this shaft speed.
        @test sol[m.ship.Surge][end] ≈ U rtol = 5e-3
        @test abs(sol[m.ship.Surge][end] - sol(1000.0; idxs = m.ship.Surge)) < 5e-3
        # The rudder holds the course.
        @test maximum(abs.(sol[m.ship.Yaw])) < 1e-2
        @test sol[m.prop.rpm][end] ≈ ow.rpm rtol = 1e-6
        # Net thrust balances the hull resistance plus the drag of the rudder in the
        # slipstream (about 1 % of it); torque and power are the open-water values.
        @test sol[m.prop.Thrust][end] > -sol[m.hydro.Resistance][end]
        @test sol[m.prop.Thrust][end] ≈ -sol[m.hydro.Resistance][end] rtol = 2e-2
        @test sol[m.prop.Thrust][end] / (1 - TP) ≈ ow.thrust rtol = 1e-2
        @test sol[m.prop.ShaftPower][end] ≈ ow.power rtol = 1e-2
        # The plant in the slot delivers exactly what the propeller absorbs.
        @test sol[p.shaft_power][end] ≈ sol[m.prop.ShaftPower][end] rtol = 1e-6
        @test sol[p.shaft_power][end] / (sol[m.prop.rpm][end] * pi / 30) ≈ ow.torque rtol = 1e-2
        @test all(diff(sol[p.shaft_energy]) .>= 0)
        @test all(diff(sol[m.service_distance]) .>= 0)
        @test sol[p.main_engine_fuel_consumed_mass][end] == 0
        @printf("  %-34s %8.2f %9.2f %8.3f %10.0f %11.0f %9.2f\n", "full-scale estimate (R2 = 26 909)", sol[m.ship.Surge][end] * 3600 / 1852, sol[m.prop.rpm][end],
            ow.J, sol[m.prop.Thrust][end] / (1 - TP) / 1e3, sol[p.shaft_power][end] / (sol[m.prop.rpm][end] * pi / 30) / 1e3, sol[p.shaft_power][end] / 1e6)
    end
    # For comparison, by arithmetic only: the model resistance coefficient at full scale.
    for kn in (10.0, 12.5, 14.0, 15.5)
        ow = open_water_point(75046.4, kn * 1852 / 3600)
        @printf("  %-34s %8.2f %9.2f %8.3f %10.0f %11.0f %9.2f\n", "model R'0 at full scale (75 046)", kn, ow.rpm, ow.J, ow.thrust / 1e3, ow.torque / 1e3, ow.power / 1e6)
    end
end

@testset "turning circles with the full-scale resistance estimate (printed, not judged)" begin
    U0 = 15.5 * 1852 / 3600
    ow = open_water_point(26908.9, U0)
    for (name, deg) in ((:port, 35.0), (:starboard, -35.0))
        run = turning_circle(deg; R2_resistance = 26908.9, rpm = ow.rpm)
        @test successful_retcode(run.sol)
        ref = PUBLISHED[name].full_scale_simulation
        @printf("  %-9s turn at %.1f rpm: advance %.3f Lpp, tactical diameter %.3f Lpp; source full-scale simulation %.2f, %.2f; deviation %+5.1f %%, %+5.1f %%\n",
            name, ow.rpm, run.AD, run.DT, ref.AD, ref.DT, 100 * (run.AD / ref.AD - 1), 100 * (run.DT / ref.DT - 1))
    end
end

end
