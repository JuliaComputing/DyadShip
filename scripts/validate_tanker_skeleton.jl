#!/usr/bin/env julia
# Focused checks of the tanker skeleton: the hull parameters added for it leave the
# existing speed trial unchanged; mass bookkeeping of the two load conditions and of a
# ballast transfer; hydrostatic response; power bookkeeping of the plant on this ship.
#
# Runs four analyses only; never the generated test suite.
#
# Usage:  julia --project=<environment with DyadShip> scripts/validate_tanker_skeleton.jl

using DyadShip, Printf, Test
using DyadInterface: symbolic_container

const S = DyadShip.Ship6DOF
successful_retcode(sol) = string(sol.retcode) == "Success"

# Draft at which the sample hull's displacement polynomial (ShipBody default
# `Disp_Table`: 1.2514e6 T + 4.1292e4 T² kg) carries a given mass.
draft_for(mass) = (-1.2514e6 + sqrt(1.2514e6^2 + 4 * 4.1292e4 * mass)) / (2 * 4.1292e4)

# Masses of the two conditions, added up here from the documented numbers [kg].
const WING = 2 * 0.45 * 100 * 1025
const LADEN = 2500e3 + 3 * 900e3 + 2 * 0 + 150e3 + WING
const BALLAST = 2500e3 + 3 * 20e3 + 2 * 700e3 + 150e3 + WING

@testset "tanker skeleton" begin

@testset "hull parameters leave the governor-driven ship unchanged" begin
    res = S.SpeedTrialTransient(); sol = res.sol; m = symbolic_container(res)
    @test successful_retcode(sol)
    @test sol[m.ship.Surge][end] ≈ 6.053806940 rtol = 1e-9
    @test sol[m.prop.ShaftPower][end] ≈ 1110404.213262 rtol = 1e-9
    @test sol[m.ship.pos_x][end] ≈ 3548.759866 rtol = 1e-9
end

results = Dict{String, Any}()

for (label, analysis, total, cargo, peak) in (("laden", S.TankerLadenTrialTransient, LADEN, 900e3, 0.0),
                                               ("ballast", S.TankerBallastTrialTransient, BALLAST, 20e3, 700e3))
@testset "speed trial, $label" begin
    res = analysis(); sol = res.sol; m = symbolic_container(res); p = m.plant
    @test successful_retcode(sol)
    val(x) = sol[x]; last_of(x) = sol[x][end]
    burnt = val(p.main_engine_fuel_consumed_mass) .+ val(p.genset_fuel_consumed_mass)

    @testset "mass bookkeeping" begin
        @test val(m.total_mass)[1] == total
        @test val(m.liquid_mass)[1] == total - 2500e3
        # Cargo and peak tanks have no flow; the wing tanks are not pumped here.
        for tank in (m.cargo_1, m.cargo_2, m.cargo_3)
            @test all(val(tank.mass) .== cargo)
        end
        for tank in (m.fore_peak, m.aft_peak)
            @test all(val(tank.mass) .== peak)
        end
        @test maximum(abs.(val(m.mass_a.mass) .+ val(m.mass_b.mass) .- WING)) <= 1e-9 * WING
        # The bunker tank loses exactly what the two engines meter, at every step.
        @test maximum(abs.(val(m.bunker.mass) .- (150e3 .- burnt))) <= 1e-6
        @test maximum(abs.(val(m.total_mass) .- (total .- burnt))) <= 1e-6
        @test last_of(m.bunker.mass) < 150e3
    end

    @testset "hydrostatics" begin
        # The hull carries the total mass, at the draft its displacement polynomial gives
        # (under way the ship also squats and trims a little).
        @test last_of(m.ship.Displacement) ≈ last_of(m.total_mass) rtol = 5e-3
        @test last_of(m.ship.Draft) ≈ draft_for(total) rtol = 1e-2
        @test abs(last_of(m.ship.Heel)) < 1e-2
    end

    @testset "power bookkeeping" begin
        @test last_of(m.prop.rpm) ≈ 80 rtol = 1e-6
        @test last_of(p.main_engine.ShaftPower) ≈ last_of(p.shaft_power) rtol = 1e-5
        @test last_of(p.shaft_power) ≈ last_of(m.prop.ShaftPower) rtol = 1e-5
        @test 0 < last_of(m.prop.ThrustPower) < last_of(m.prop.ShaftPower)
        J = sol.ps[p.J_engine]; w = val(p.engine_inertia.w)
        @test last_of(p.main_engine.KWh) * 3.6e6 ≈ last_of(p.shaft_energy) + J * (w[end]^2 - w[1]^2) / 2 rtol = 1e-5
        @test maximum(abs.(val(p.genset_electric_power) .- val(p.hotel_electric_power))) <= 1e-6 * last_of(p.hotel_electric_power)
        @test all(diff(val(p.main_engine_fuel_consumed_mass)) .>= 0)
        @test all(diff(val(p.genset_fuel_consumed_mass)) .>= 0)
        @test all(diff(val(m.service_distance)) .>= 0)
    end
    results[label] = (draft = last_of(m.ship.Draft), surge = last_of(m.ship.Surge), power = last_of(p.shaft_power))
    @printf("  %-7s: total %.1f t, draft %.3f m (polynomial %.3f m), %.3f m/s at 80 rpm, shaft %.1f kW, fuel burnt %.2f kg, %.1f m\n",
        label, last_of(m.total_mass) / 1e3, last_of(m.ship.Draft), draft_for(total), last_of(m.ship.Surge), last_of(p.shaft_power) / 1e3, burnt[end], last_of(m.service_distance))
end
end

@testset "load conditions differ as they should" begin
    @test results["laden"].draft > results["ballast"].draft + 0.5
    @test results["ballast"].surge >= results["laden"].surge
end

@testset "ballast transfer under way" begin
    res = S.TankerBallastTransferTransient(); sol = res.sol; m = symbolic_container(res); p = m.plant
    @test successful_retcode(sol)
    t = sol.t; val(x) = sol[x]
    a = val(m.mass_a.mass); b = val(m.mass_b.mass)
    before = findall(<=(100.0), t); during = findall(x -> 100.0 <= x <= 400.0, t); after = findall(>=(600.0), t)
    burnt = val(p.main_engine_fuel_consumed_mass) .+ val(p.genset_fuel_consumed_mass)

    # Nothing moves before the pump starts.
    @test maximum(abs.(a[before] .- a[1])) <= 1e-6 * a[1]
    # Water goes from the port tank to the starboard tank and none is lost.
    @test all(diff(a[during]) .<= 0)
    @test all(diff(b[during]) .>= 0)
    @test maximum(abs.(a .+ b .- WING)) <= 1e-9 * WING
    moved = a[1] - a[end]
    @test moved > 10e3
    @test b[end] - b[1] ≈ moved rtol = 1e-9
    # The transfer has stopped and the tanks stay within their fill limits.
    @test maximum(abs.(a[after] .- a[end])) <= 1e-3 * moved
    @test 0.05 * 100 * 1025 <= a[end]
    @test b[end] <= 0.85 * 100 * 1025
    # The only change of the ship's liquid mass is the fuel burnt.
    @test maximum(abs.(val(m.liquid_mass) .- (val(m.liquid_mass)[1] .- burnt))) <= 1e-6
    # The ship heels as the water moves, then holds the heel; draft and speed stay.
    heel = val(m.ship.Heel)
    @test abs(heel[last(before)]) < 1e-3
    @test abs(heel[end]) > 5e-3
    @test abs(heel[end] - heel[first(after)]) < 1e-3
    @test val(m.ship.Draft)[end] ≈ val(m.ship.Draft)[last(before)] rtol = 1e-3
    @test val(p.shaft_power)[end] ≈ val(m.prop.ShaftPower)[end] rtol = 1e-5
    @printf("  transfer: %.2f t moved port to starboard (%.1f -> %.1f t and %.1f -> %.1f t), heel %.4f rad, liquid mass change %.2f kg = fuel burnt %.2f kg\n",
        moved / 1e3, a[1] / 1e3, a[end] / 1e3, b[1] / 1e3, b[end] / 1e3, heel[end], val(m.liquid_mass)[end] - val(m.liquid_mass)[1], burnt[end])
end

end
