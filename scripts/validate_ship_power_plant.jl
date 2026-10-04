#!/usr/bin/env julia
# Focused checks of the ship power plant on the sample hull: that the hull split left the
# existing speed trial unchanged, signs, the power balance from the engines through the
# gear and the bus to the propeller and the hotel load, and the observation counters.
#
# Runs four analyses only; never the generated test suite.
#
# Usage:  julia --project=<environment with DyadShip> scripts/validate_ship_power_plant.jl

using DyadShip, Printf, Test, TOML
using DyadInterface: symbolic_container

const S = DyadShip.Ship6DOF
const PRESETS = joinpath(pkgdir(DyadShip), "assets", "presets", "Synthetic")

successful_retcode(sol) = string(sol.retcode) == "Success"
nondecreasing(v) = all(diff(v) .>= 0)

# An out-of-table fuel counter has a rate that jumps at a table end, and the one accepted
# step that contains the crossing can end slightly below its start. The allowance is the
# local error the analyses request (abstol = reltol = 1e-6), not a figure taken from a run.
almost_nondecreasing(v) = all(diff(v) .>= -(1e-6 .+ 1e-6 .* abs.(v[1:(end - 1)])))

# Independent piecewise-linear lookup of an engine preset file (callers stay inside it).
function preset_rate(file, P_W)
    a = TOML.parsefile(joinpath(PRESETS, file))
    xs, ys = Float64.(a["SFOC_P"]), Float64.(a["SFOC_g"])
    x = P_W / 1000
    i = clamp(searchsortedlast(xs, x), 1, length(xs) - 1)
    sfoc = ys[i] + (ys[i + 1] - ys[i]) * (x - xs[i]) / (xs[i + 1] - xs[i])
    return (rate = sfoc * x / 3.6e6 + a["m_dot_idle"], inside = xs[1] <= x <= xs[end])
end

@testset "ship power plant" begin

@testset "hull split leaves the governor-driven ship unchanged" begin
    res = S.SpeedTrialTransient(); sol = res.sol; m = symbolic_container(res)
    @test successful_retcode(sol)
    # Values of the same analysis before `StandardShip` was split into
    # `SingleScrewHull` plus its governor.
    @test sol[m.ship.Surge][end] ≈ 6.053806940 rtol = 1e-9
    @test sol[m.prop.Thrust][end] ≈ 129102.465380 rtol = 1e-9
    @test sol[m.prop.ShaftPower][end] ≈ 1110404.213262 rtol = 1e-9
    @test sol[m.ship.pos_x][end] ≈ 3548.759866 rtol = 1e-9
end

@testset "transit without a shaft machine" begin
    res = S.DieselShipTransitTransient(); sol = res.sol; m = symbolic_container(res); p = m.plant
    @test successful_retcode(sol)
    t = sol.t; late = findall(>=(100.0), t)
    val(x) = sol[x]
    last_of(x) = sol[x][end]

    @testset "signs" begin
        @test all(val(m.ship.Surge) .> 0)
        @test all(val(m.prop.Thrust)[late] .> 0)
        @test all(val(p.shaft_power) .> 0)                    # plant drives the shaft ahead
        @test all(val(p.main_engine.ShaftPower) .>= 0)
        @test all(val(p.main_engine_fuel_mass_rate) .>= 0) && all(val(p.genset_fuel_mass_rate) .>= 0)
        @test all(val(p.hotel_electric_power) .> 0) && all(val(p.genset_electric_power) .> 0)
        @test all(val(p.pto_electric_power) .== 0) && last_of(p.pto_electric_energy) == 0
    end

    @testset "propulsion line at the steady state" begin
        # Gear kinematics, exact at every step.
        @test maximum(abs.(val(p.main_engine.core.rpm) .- sol.ps[p.gear_ratio] .* val(m.prop.rpm))) <= 1e-6
        @test last_of(m.prop.rpm) ≈ 80 rtol = 1e-6             # the order is met
        # Engine brake power = power at the propeller flange (ideal gear) = power the
        # propeller absorbs, once the inertias no longer accelerate.
        @test last_of(p.main_engine.ShaftPower) ≈ last_of(p.shaft_power) rtol = 1e-5
        @test last_of(p.shaft_power) ≈ last_of(m.prop.ShaftPower) rtol = 1e-5
        # Hull: thrust power is below shaft power, and the speed has settled.
        @test 0 < last_of(m.prop.ThrustPower) < last_of(m.prop.ShaftPower)
        @test abs(last_of(m.ship.Surge) - sol(550.0; idxs = m.ship.Surge)) < 1e-3
        # The engine-driven ship ends where the ideal-governor ship ends at the same order.
        ref = S.SpeedTrialTransient(model = S.SpeedTrial(; name = :SpeedTrial, rpm = 80.0)); rm = symbolic_container(ref)
        @test last_of(m.ship.Surge) ≈ ref.sol[rm.ship.Surge][end] rtol = 1e-4
        @test last_of(m.prop.ShaftPower) ≈ ref.sol[rm.prop.ShaftPower][end] rtol = 1e-3
        @test last_of(m.prop.Thrust) ≈ ref.sol[rm.prop.Thrust][end] rtol = 1e-3
    end

    @testset "energy over the run" begin
        # Work done by the engine = work delivered to the propeller shaft + change of the
        # engine rotor's kinetic energy (the gear is ideal). Both sides are solver states.
        J = sol.ps[p.J_engine]; w = val(p.engine_inertia.w)
        engine_work = last_of(p.main_engine.KWh) * 3.6e6
        @test engine_work ≈ last_of(p.shaft_energy) + J * (w[end]^2 - w[1]^2) / 2 rtol = 1e-5
        @test last_of(p.shaft_energy) > 0
    end

    @testset "bus" begin
        # One source, one load: what the set delivers is what the hotel load receives.
        @test maximum(abs.(val(p.genset_electric_power) .- val(p.hotel_electric_power))) <= 1e-6 * last_of(p.hotel_electric_power)
        # Generating-set engine power = electrical power + loss in the series resistance.
        i = last_of(p.genset_machine.emf.i); R = sol.ps[p.R_genset]
        @test last_of(p.genset_engine.ShaftPower) ≈ last_of(p.genset_electric_power) + R * i^2 rtol = 1e-6
        # The conductance load receives the demand at nominal voltage, slightly less here.
        demand = 400e3
        @test 0.98 * demand < last_of(p.hotel_electric_power) < demand
        @test last_of(p.genset_engine.core.rpm) ≈ sol.ps[p.genset_rpm] rtol = 1e-6
    end

    @testset "fuel per consumer, coverage and counters" begin
        # Rates equal an independent lookup of the two preset files at the engine powers.
        a = preset_rate("diesel_engine_a.toml", last_of(p.main_engine.ShaftPower))
        b = preset_rate("diesel_engine_b.toml", last_of(p.genset_engine.ShaftPower))
        @test a.inside && b.inside
        @test last_of(p.main_engine_fuel_mass_rate) ≈ a.rate rtol = 1e-9
        @test last_of(p.genset_fuel_mass_rate) ≈ b.rate rtol = 1e-9
        # After the start both engines stay inside their tables.
        @test all(val(p.main_engine_sfoc_covered)[late] .== 1) && all(val(p.genset_sfoc_covered)[late] .== 1)
        for x in (p.main_engine_fuel_extrapolated_mass, p.genset_fuel_extrapolated_mass)
            @test maximum(val(x)[late]) == minimum(val(x)[late])
        end
        # Counters start at zero and do not decrease at the accepted steps.
        for x in (p.main_engine_fuel_consumed_mass, p.genset_fuel_consumed_mass, p.hotel_electric_energy,
                  p.genset_electric_energy, p.shaft_energy, p.main_engine_running_time, p.genset_running_time, m.service_distance)
            @test val(x)[1] == 0
            @test nondecreasing(val(x))
        end
        for x in (p.main_engine_fuel_extrapolated_mass, p.genset_fuel_extrapolated_mass)
            @test val(x)[1] == 0
            @test almost_nondecreasing(val(x))
            @printf("  transit: out-of-table counter, largest fall at an accepted step %.3e kg\n", -min(0.0, minimum(diff(val(x)))))
        end
        @test all(val(p.main_engine_fuel_extrapolated_mass) .<= val(p.main_engine_fuel_consumed_mass))
        # No stopped state: running time is the elapsed time. Straight course: the
        # distance is the advance along world x up to the hull's pitch and heave motion.
        @test last_of(p.main_engine_running_time) ≈ t[end] rtol = 1e-9
        @test last_of(p.genset_running_time) ≈ t[end] rtol = 1e-9
        @test last_of(m.service_distance) ≈ last_of(m.ship.pos_x) rtol = 1e-3
        @printf("  transit: %.2f rpm, %.3f m/s, shaft %.1f kW, thrust %.1f kN; main engine %.5f kg/s, set %.5f kg/s (%.1f kW electrical); fuel %.3f + %.3f kg over %.0f s, %.1f m\n",
            last_of(m.prop.rpm), last_of(m.ship.Surge), last_of(p.shaft_power) / 1e3, last_of(m.prop.Thrust) / 1e3,
            last_of(p.main_engine_fuel_mass_rate), last_of(p.genset_fuel_mass_rate), last_of(p.genset_electric_power) / 1e3,
            last_of(p.main_engine_fuel_consumed_mass), last_of(p.genset_fuel_consumed_mass), t[end], last_of(m.service_distance))
    end
end

@testset "transit with a shaft machine: take-in, then take-off" begin
    res = S.DieselShipShaftMachineTransient(); sol = res.sol; m = symbolic_container(res); p = m.plant
    @test successful_retcode(sol)
    t = sol.t; val(x) = sol[x]
    lo = findall(x -> 100.0 <= x <= 300.0, t)        # low order, engine below the matching speed
    i_lo = last(lo); i_hi = length(t)
    @test minimum(abs.(t .- 300.0)) < 1e-9           # a step ends at the order change

    @testset "bus balance at every accepted step" begin
        resid = val(p.genset_electric_power) .+ val(p.pto_electric_power) .- val(p.hotel_electric_power)
        @test maximum(abs.(resid)) <= 1e-6 * maximum(val(p.hotel_electric_power))
    end

    @testset "signs of the two modes" begin
        # Take-in: the machine draws from the bus and helps drive the shaft.
        @test all(val(p.pto_electric_power)[lo] .< 0)
        @test val(p.genset_electric_power)[i_lo] > val(p.hotel_electric_power)[i_lo]
        @test val(p.main_engine.ShaftPower)[i_lo] < val(p.shaft_power)[i_lo]
        # Take-off: the machine feeds the bus and loads the engine.
        @test val(p.pto_electric_power)[i_hi] > 0
        @test val(p.genset_electric_power)[i_hi] < val(p.hotel_electric_power)[i_hi]
        @test val(p.main_engine.ShaftPower)[i_hi] > val(p.shaft_power)[i_hi]
        # The signed energy falls during take-in and rises during take-off.
        E = val(p.pto_electric_energy)
        @test E[i_lo] < 0 && E[i_hi] > E[i_lo]
        @test all(val(p.shaft_power) .> 0)
    end

    @testset "engine shaft power balance at the end of each phase" begin
        R = sol.ps[p.R_pto]
        for i in (i_lo, i_hi)
            ip = val(p.pto.emf.i)[i]
            # Engine brake power = propeller shaft power + machine shaft power, the
            # latter being its electrical output plus the loss in its series resistance.
            @test val(p.main_engine.ShaftPower)[i] ≈ val(p.shaft_power)[i] + val(p.pto_electric_power)[i] + R * ip^2 rtol = 1e-4
            @test val(p.shaft_power)[i] ≈ val(m.prop.ShaftPower)[i] rtol = 1e-4
        end
        @test val(m.prop.rpm)[i_lo] ≈ 75 rtol = 1e-5
        @test val(m.prop.rpm)[i_hi] ≈ 86 rtol = 1e-5
    end

    @testset "fuel moves between the consumers; counters" begin
        late = findall(>=(100.0), t)
        @test all(val(p.main_engine_sfoc_covered)[late] .== 1) && all(val(p.genset_sfoc_covered)[late] .== 1)
        # Take-off relieves the generating set and loads the main engine.
        @test val(p.genset_fuel_mass_rate)[i_hi] < val(p.genset_fuel_mass_rate)[i_lo]
        @test val(p.main_engine_fuel_mass_rate)[i_hi] > val(p.main_engine_fuel_mass_rate)[i_lo]
        for x in (p.main_engine_fuel_consumed_mass, p.genset_fuel_consumed_mass, p.hotel_electric_energy,
                  p.genset_electric_energy, p.main_engine_running_time, p.genset_running_time, m.service_distance)
            @test nondecreasing(val(x))
        end
        for x in (p.main_engine_fuel_extrapolated_mass, p.genset_fuel_extrapolated_mass)
            @test almost_nondecreasing(val(x))
            @test maximum(val(x)[late]) == minimum(val(x)[late])
            @printf("  shaft machine: out-of-table counter, largest fall at an accepted step %.3e kg\n", -min(0.0, minimum(diff(val(x)))))
        end
        @printf("  shaft machine: take-in %.1f kW at %.1f rpm (set %.1f kW, main engine %.1f kW for %.1f kW on the shaft); take-off %.1f kW at %.1f rpm (set %.1f kW, main engine %.1f kW for %.1f kW on the shaft)\n",
            -val(p.pto_electric_power)[i_lo] / 1e3, val(m.prop.rpm)[i_lo], val(p.genset_electric_power)[i_lo] / 1e3, val(p.main_engine.ShaftPower)[i_lo] / 1e3, val(p.shaft_power)[i_lo] / 1e3,
            val(p.pto_electric_power)[i_hi] / 1e3, val(m.prop.rpm)[i_hi], val(p.genset_electric_power)[i_hi] / 1e3, val(p.main_engine.ShaftPower)[i_hi] / 1e3, val(p.shaft_power)[i_hi] / 1e3)
    end
end

end
