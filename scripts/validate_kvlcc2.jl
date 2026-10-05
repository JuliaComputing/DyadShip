#!/usr/bin/env julia
# KVLCC2 at full scale: the hull preset reaches the solved model unchanged, the ship holds
# its approach speed, and the 35° turning circles are compared with the values printed in
# Yasukawa and Yoshimura, J. Mar. Sci. Technol. 20 (2015) 37-52, Table 6 (their full-scale
# simulation) and Table 4 (their free-running test of a 7 m model). Nothing is tuned: the
# deviations are printed, and only a coarse 25 % band is asserted.
#
# Runs three analyses only; never the generated test suite.
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

function turning_circle(rudder_deg)
    res = S.KVLCC2TurningCircleTransient(model = S.KVLCC2TurningCircle(; name = :KVLCC2TurningCircle, rudder_deg = rudder_deg))
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
    @test abs(sol[m.ship.Yaw][i0]) < 1e-6
    @test sol[m.ship.Draft][i0] ≈ 20.8 rtol = 1e-3
    @test abs(sol[m.ship.Heel][i0]) < 1e-4
    @test sol[m.ship.Displacement][i0] ≈ 1025 * 312600 rtol = 1e-3
end

@testset "turning circles against the published values" begin
    # Rudder to port turns the ship to port (positive y), and conversely.
    @test port.side > 0
    @test starboard.side < 0
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

end
