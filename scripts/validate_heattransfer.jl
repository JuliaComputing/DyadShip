#!/usr/bin/env julia
# Check the ported heat-transfer components against the correlations recomputed
# independently in Julia: the four ConvectionFactors and their branch structures
# (the cylinder's natural/forced switchover, the internal-convection sign flip
# through ΔT = 0), and the ConvRadSunWall composition over a full day.
#
# Usage:  ~/dyad-fleet/heavy -m 8G ../julia-dyad.sh scripts/validate_heattransfer.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using DyadShip, Printf, Test
using SciMLBase: successful_retcode
using DyadInterface: symbolic_container

const TH = DyadShip.Thermal

# Reference correlations, written from the literature rather than from the port.
churchill_chu(Ra, Pr) = (0.6 + 0.387 * Ra^(1/6) / (1 + (0.559/Pr)^(9/16))^(8/27))^2
function churchill_bernstein(Re, Pr)
    base = 0.3 + (0.62 * sqrt(Re) * Pr^(1/3)) / ((1 + (0.4/Pr)^(2/3))^(1/4))
    Re < 4000    ? base :
    Re < 20000   ? base * (1 + (Re/282000)^(5/8))^(4/5) :
    Re < 400000  ? base * (1 + sqrt(Re/282000)) :
                   base * (1 + (Re/282000)^(5/8))^(4/5)
end
petukhov(Re, Pr) = (0.037 * Re^0.8 * Pr) / (1 + 2.443 / Re^0.1 * (Pr^(2/3) - 1))
# Matches ExternalConvection including its rest regularisation (see its docstring).
function nusselt_jurges(u; u_reg = 1e-4)
    ur = sqrt(u^2 + u_reg^2)
    7.13 * ur^0.78 + 5.35 * exp(-0.6ur)
end

# Air properties, the defaults both components carry.
const K_AIR, NU_AIR, RHO_AIR, CP_AIR, BETA_AIR = 0.024, 1.562e-5, 1.184, 1007.0, 0.003354016
const ALPHA = K_AIR / (RHO_AIR * CP_AIR)
const PR_CYL = NU_AIR / ALPHA
const PR_PLATE = (NU_AIR * RHO_AIR) * CP_AIR / K_AIR

@testset "Thermal" begin

@testset "speed sweep" begin
    res = TH.ConvectionSpeedSweepAnalysis(); sol = res.sol; m = symbolic_container(res)
    @test successful_retcode(sol)
    ts = collect(0.0:1.0:200.0)
    grab(sym) = [sol(t; idxs = sym) for t in ts]

    h_cyl, Re_cyl, Ra_cyl = grab(m.cyl.h), grab(m.cyl.Re_D), grab(m.cyl.Ra_D)
    h_ww, h_lw = grab(m.windward.h), grab(m.leeward.h)
    h_plate, Re_plate = grab(m.plate.h), grab(m.plate.Re)
    T_body = grab(m.body.T)
    @test all(isfinite, vcat(h_cyl, Re_cyl, Ra_cyl, h_ww, h_lw, h_plate, Re_plate, T_body))

    # At rest the cylinder is pure natural convection.
    D = 0.151
    ΔT0 = T_body[1] - 298.15
    Ra0 = 9.80665 * BETA_AIR * abs(ΔT0) * D^3 / (NU_AIR * ALPHA)
    @test Ra_cyl[1] ≈ Ra0 rtol=1e-9
    @test h_cyl[1] ≈ churchill_chu(Ra0, PR_CYL) * K_AIR / D rtol=1e-6
    @printf("  cylinder at rest: h = %.2f W/(m2.K), Ra_D = %.3e\n", h_cyl[1], Ra_cyl[1])

    # At the end of the ramp forced convection dominates, in the 2e4..4e5 band.
    ΔTend = T_body[end] - 298.15
    Reend = 20.0 * D / NU_AIR
    @test Re_cyl[end] ≈ Reend rtol=1e-6
    @test 20000 < Reend < 400000
    Raend = 9.80665 * BETA_AIR * abs(ΔTend) * D^3 / (NU_AIR * ALPHA)
    href = max(churchill_bernstein(Reend, PR_CYL), churchill_chu(Raend, PR_CYL)) * K_AIR / D
    @test h_cyl[end] ≈ href rtol=1e-6
    @test h_cyl[end] > 5 * h_cyl[1]
    @printf("  cylinder at 20 m/s: h = %.2f W/(m2.K), Re_D = %.3e (forced/natural = %.1f)\n",
            h_cyl[end], Re_cyl[end], h_cyl[end] / h_cyl[1])
    @test T_body[end] < T_body[1]

    # At rest neither wall is in anyone's lee, so both take W = 1; once the wind
    # moves, the leeward wall is exactly the windward wall halved.
    u_wind = grab(m.wind.y)
    @test h_ww[1] ≈ nusselt_jurges(0.0) rtol=1e-9
    @test h_ww[1] ≈ 5.35 rtol=2e-3          # the regularisation costs 0.1% at rest
    @test h_lw[1] ≈ h_ww[1] rtol=1e-12
    moving = findall(>(0.0), u_wind)
    @test length(moving) == length(ts) - 1
    @test all(h_lw[moving] .≈ 0.5 .* h_ww[moving])
    @test h_ww[end] ≈ nusselt_jurges(20.0) rtol=1e-6
    @printf("  wall: windward h = %.2f -> %.2f, leeward half of it whenever the wind moves\n",
            h_ww[1], h_ww[end])

    # Flat plate: the run crosses out of the laminar range into the fitted one.
    @test Re_plate[1] ≈ 0 atol=1e-3
    @test Re_plate[end] ≈ 20.0 * 2.0 / NU_AIR rtol=1e-6
    @test Re_plate[end] > 5e5
    @test h_plate[end] ≈ petukhov(Re_plate[end], PR_PLATE) * K_AIR / 2.0 rtol=1e-6
    icross = findfirst(>(5e5), Re_plate)
    @printf("  plate: Re = %.2e at 20 m/s, crosses 5e5 at t = %.0f s, h = %.2f W/(m2.K)\n",
            Re_plate[end], ts[icross], h_plate[end])
end

@testset "internal convection sign flip" begin
    res = TH.InternalConvectionCrossoverAnalysis(); sol = res.sol; m = symbolic_container(res)
    @test successful_retcode(sol)
    ts = collect(0.0:1.0:200.0)
    grab(sym) = [sol(t; idxs = sym) for t in ts]

    ΔT = grab(m.floor_surf.ΔT)
    h_floor, h_ceil, h_wall = grab(m.floor_surf.h), grab(m.ceil_surf.h), grab(m.wall_surf.h)
    @test all(isfinite, vcat(ΔT, h_floor, h_ceil, h_wall))

    # Air ramps 293.15 -> 303.15 past a 298.15 K surface, so ΔT falls through zero.
    @test ΔT[1] > 0 && ΔT[end] < 0
    icross = argmin(abs.(ΔT))
    @printf("  ΔT crosses zero at t = %.0f s (ΔT = %+.3f K)\n", ts[icross], ΔT[icross])

    # Away from the crossing, the vertical wall follows 1.3*|ΔT|^(1/3) exactly and
    # the horizontal pair swaps which of 1.51 / 0.76 it uses.
    for i in (1, length(ts))
        @test h_wall[i] ≈ 1.3 * abs(ΔT[i])^(1/3) rtol=1e-4
    end
    @test h_floor[1] ≈ 1.51 * abs(ΔT[1])^(1/3) rtol=1e-4     # surface warmer: plume
    @test h_ceil[1] ≈ 0.76 * abs(ΔT[1])^(1/3) rtol=1e-4      # surface warmer: stratified
    @test h_floor[end] ≈ 0.76 * abs(ΔT[end])^(1/3) rtol=1e-4
    @test h_ceil[end] ≈ 1.51 * abs(ΔT[end])^(1/3) rtol=1e-4
    @test h_floor[1] > h_ceil[1] && h_floor[end] < h_ceil[end]
    @printf("  floor h: %.2f -> %.2f, ceiling h: %.2f -> %.2f W/(m2.K)\n",
            h_floor[1], h_floor[end], h_ceil[1], h_ceil[end])

    # The rounded-off cube root stays small and finite at the crossing.
    @test h_wall[icross] < 1.0
    @printf("  wall h at the crossing: %.4f W/(m2.K) (upstream's |ΔT|^(1/3) is non-differentiable here)\n",
            h_wall[icross])
end


@testset "external wall, winter solstice" begin
    res = TH.SunWallWinterAnalysis(); sol = res.sol; m = symbolic_container(res)
    @test successful_retcode(sol)
    ts = collect(0.0:300.0:86400.0)
    grab(sym) = [sol(t; idxs = sym) for t in ts]
    hours = ts ./ 3600

    T_s, T_n, T_x = grab(m.south.port_Wall.T), grab(m.north.port_Wall.T), grab(m.shaded.port_Wall.T)
    S_s, S_n, S_x = grab(m.south.S_eff), grab(m.north.S_eff), grab(m.shaded.S_eff)
    sunh = grab(m.sun.SunHeight)
    @test all(isfinite, vcat(T_s, T_n, T_x, S_s, S_n, S_x, sunh))

    T_air, T_sky, SIGMA = 298.15, 278.15, 5.670374419e-8
    Surf, eps_w, vf = 30.0, 0.9, 0.5

    # The two long-wave paths carry exactly the view-factor split of one emissivity-area.
    # This is what the flattened version got wrong: it radiated the whole surface to sky.
    # Read off the saved steps rather than through `sol(t)`: the dense output interpolates
    # an observed T^4 series independently of the temperature series, which puts a few
    # tenths of a percent between them between steps and is ambiguous at the sun gate.
    Tn_s = sol[m.south.port_Wall.T]
    q_sky = sol[m.south.rad_sky.Q_flow]
    q_air = sol[m.south.rad_air.Q_flow]
    # The tolerance is the nonlinear surface balance's, not floating point's: the outer
    # face carries no capacitance, so its temperature comes out of a Newton solve of
    # radiation (in T^4) against convection, solar gain and conduction, and Q_flow agrees
    # with its own defining equation to about 3e-6.
    @test all(isapprox.(q_sky, vf * eps_w * Surf * SIGMA .* (Tn_s .^ 4 .- T_sky^4); rtol = 1e-4))
    @test all(isapprox.(q_air, (1 - vf) * eps_w * Surf * SIGMA .* (Tn_s .^ 4 .- T_air^4);
                        rtol = 1e-4, atol = 1e-3))
    Gr_sky = q_sky[1] / (SIGMA * (Tn_s[1]^4 - T_sky^4))
    Gr_air = q_air[1] / (SIGMA * (Tn_s[1]^4 - T_air^4))
    @test Gr_sky + Gr_air ≈ eps_w * Surf rtol=1e-9
    @printf("  radiation split: Gr_sky = %.2f, Gr_air = %.2f, sum = %.2f (= e*Surf = %.2f)\n",
            Gr_sky, Gr_air, Gr_sky + Gr_air, eps_w * Surf)

    # At midnight there is no sun and the sky path cools the wall.
    @test S_s[1] == 0
    @test q_sky[1] > 0
    @test T_s[1] < T_air

    # The sun-height gate: below minSunHeight nothing is counted.
    below = findall(<(10.0), sunh)
    @test !isempty(below)
    @test all(grab(m.south.S_gated)[below] .== 0)

    # In winter the sun never gets behind the building: the north wall is dark all day and
    # the south wall takes a near-square-on 23 deg sun at solar noon.
    @test maximum(S_n) == 0
    @test maximum(S_s) > 800
    @test maximum(T_s) > maximum(T_n) + 5
    @test abs(hours[argmax(S_s)] - 13.55) < 1.0
    @printf("  winter: south peak S_eff %.0f W/m2 at %.1f h (T %.1f -> %.1f C), north %.0f W/m2\n",
            maximum(S_s), hours[argmax(S_s)], T_s[1] - 273.15, maximum(T_s) - 273.15, maximum(S_n))

    # The overhang is sized for a high sun and admits the whole of a 23 deg one, which is
    # the point of the device: shade in summer, gain in winter. Its summer effect is
    # checked in the next testset.
    @test all(S_x .≈ S_s)
    @test maximum(T_x) ≈ maximum(T_s)
    @printf("  overhang admits the low winter sun in full: peak S_eff %.0f W/m2\n", maximum(S_x))

    # Transmission through the glazed fraction, and none where tau_trans = 0.
    q_tr_s, q_tr_n = grab(m.south.transmitted.Q_flow), grab(m.north.transmitted.Q_flow)
    @test all(q_tr_s .≈ 0.1 * Surf * 0.75 .* S_s)
    @test all(q_tr_n .== 0)
    @printf("  transmitted through the glazing: peak %.0f W\n", maximum(q_tr_s))

    # Convection is identical on all three walls: the breeze is square-on to every normal.
    h_s, h_n = grab(m.south.conv.h), grab(m.north.conv.h)
    @test all(h_s .≈ h_n)
    @printf("  film coefficient, 3 m/s square-on: %.2f W/(m2.K) on every wall\n", h_s[1])
end

@testset "external wall, summer solstice" begin
    res = TH.SunWallSummerAnalysis(); sol = res.sol; m = symbolic_container(res)
    @test successful_retcode(sol)
    ts = collect(0.0:300.0:86400.0)
    grab(sym) = [sol(t; idxs = sym) for t in ts]
    hours = ts ./ 3600

    S_s, S_n, S_x = grab(m.south.S_eff), grab(m.north.S_eff), grab(m.shaded.S_eff)
    T_s, T_x = grab(m.south.port_Wall.T), grab(m.shaded.port_Wall.T)
    @test all(isfinite, vcat(S_s, S_n, S_x, T_s, T_x))

    # The orientations invert: the sun rises north of east and sets north of west, so the
    # north wall is lit at dawn and dusk while a vertical south wall catches only the
    # cosine of a 70 deg midday sun. This is the pair that pins the solar azimuth down.
    @test maximum(S_n) > 300
    @test maximum(S_s) < 400
    @test maximum(S_n) > maximum(S_s)
    @test abs(hours[argmax(S_s)] - 13.55) < 1.0

    # The north wall's gain is two lobes either side of midday, not one peak at it.
    lit_n = findall(>(1.0), S_n)
    @test !isempty(lit_n)
    @test minimum(hours[lit_n]) < 9
    @test maximum(hours[lit_n]) > 18
    @test all(S_n[findall(h -> 11 < h < 16, hours)] .== 0)
    @printf("  summer: north peak S_eff %.0f W/m2, lit %.1f-%.1f h and dark at midday; south peak %.0f W/m2 at %.1f h\n",
            maximum(S_n), minimum(hours[lit_n]), maximum(hours[lit_n]),
            maximum(S_s), hours[argmax(S_s)])

    # Now the overhang earns its keep: a high sun is cut back, and it can only ever
    # subtract.
    @test all(S_x .<= S_s .+ 1e-9)
    @test any(S_x .< S_s .- 1e-6)
    @test maximum(T_x) < maximum(T_s)
    @printf("  overhang: peak S_eff %.0f W/m2 (%.0f%% of the unshaded wall), peak T %.1f C vs %.1f C\n",
            maximum(S_x), 100 * maximum(S_x) / maximum(S_s),
            maximum(T_x) - 273.15, maximum(T_s) - 273.15)
end

@testset "ship compartment at sea" begin
    res = TH.DeckhouseAtSeaAnalysis(); sol = res.sol; m = symbolic_container(res)
    @test successful_retcode(sol)
    ts = collect(0.0:300.0:18000.0)
    grab(sym) = [sol(t; idxs = sym) for t in ts]
    hours = ts ./ 3600

    h_out, h_in = grab(m.live.h_outer), grab(m.live.h_inner)
    T_live, T_fix, T_steel = grab(m.live.T_air), grab(m.fixed.volume.medium.T), grab(m.steel.T_air)
    T_face, T_dew = grab(m.live.wall.Ta[1]), grab(m.dew.Tdew)
    @test all(isfinite, vcat(h_out, h_in, T_live, T_fix, T_steel, T_face, T_dew))

    # The weather-side coefficient is the wind's, alongside and at sea.
    @test h_out[1] ≈ nusselt_jurges(0.0) rtol=1e-9
    @test h_out[end] ≈ nusselt_jurges(20.0) rtol=1e-9
    @test h_out[end] / h_out[1] > 13

    # The room-side coefficient stays in the natural-convection band for a bulkhead.
    @test all(0.5 .< h_in .< 4.0)

    # Alongside, the constant U_outer = 5 is very nearly the true still-air 5.36, so the
    # two models agree: this is the control, not a coincidence.
    i2h = findfirst(==(2.0), hours)
    @test abs(T_live[i2h] - T_fix[i2h]) < 0.5
    @printf("  alongside at 2 h: live %.2f C, fixed %.2f C (h_out %.2f vs U_outer 5.0)\n",
            T_live[i2h] - 273.15, T_fix[i2h] - 273.15, h_out[i2h])

    # At sea the live compartment sheds heat the fixed-U one never notices; the fixed-U
    # air temperature is still climbing while the live one has turned over.
    @test T_live[end] < T_fix[end] - 2.0
    @test T_live[end] < T_live[i2h]
    @test T_fix[end] > T_fix[i2h]
    @printf("  at sea (5 h):     live %.2f C, fixed %.2f C — fixed U overstates by %.2f K\n",
            T_live[end] - 273.15, T_fix[end] - 273.15, T_fix[end] - T_live[end])

    # Behind insulation the films are a small part of the resistance; on bare steel they
    # are nearly all of it, so the same wind costs three times as much.
    drop_live = T_live[i2h] - T_live[end]
    drop_steel = T_steel[i2h] - T_steel[end]
    @test drop_live < 3.0
    @test drop_steel > 8.0
    @test drop_steel > 3 * drop_live
    @printf("  wind costs the insulated deckhouse %.2f K and the bare steel one %.2f K\n",
            drop_live, drop_steel)

    # The bulkhead's inner face stays clear of the dew point, so it does not sweat.
    @test all(T_face .> T_dew .+ 5)
    @printf("  inner face %.1f C against a %.1f C dew point: no condensation\n",
            T_face[end] - 273.15, T_dew[end] - 273.15)
end


@testset "environment: a day of weather" begin
    res = TH.WeatherDayTransient(); sol = res.sol; m = symbolic_container(res)
    @test successful_retcode(sol)
    ts = collect(0.0:600.0:86400.0)
    grab(sym) = [sol(t; idxs = sym) for t in ts]
    hours = ts ./ 3600

    cloud, T_sky, T_air = grab(m.env.CloudRatio), grab(m.env.T_sky), grab(m.env.T_air)
    Wx, Wy, Wz = grab(m.env.WindVector_x), grab(m.env.WindVector_y), grab(m.env.WindVector_z)
    Tdew = grab(m.env.Tdew)
    @test all(isfinite, vcat(cloud, T_sky, T_air, Wx, Wy, Wz, Tdew))

    # The clocked signal is the authority: one entry per tick, 1440 ticks in a day at 60 s.
    cl = sol[m.env.latch.cloud]
    @test length(cl) == 1440
    @test all(0 .<= cl .<= 1)
    @test cl[1] == 0

    # Clear morning, then a step to a 0.35 transmission at 13:00 that the component has to
    # recover as a 0.65 cover from the irradiance shortfall alone.
    morning = findall(h -> 8 < h < 12, hours)
    evening = findall(h -> 16 < h < 22, hours)
    @test all(cloud[morning] .< 0.02)
    @test all(abs.(cloud[evening] .- 0.65) .< 0.02)
    # and it is held through the night rather than decaying.
    @test abs(cloud[end] - 0.65) < 0.02
    @printf("  cloud ratio: %.3f through the clear morning, %.3f after the 13:00 step, %.3f at midnight\n",
            maximum(cloud[morning]), cloud[evening[1]], cloud[end])

    # No spike at sunrise: the cover is only measured where there is sun to measure it
    # against, so the ratio never leaves [0, 0.7] on this trace.
    @test maximum(cl) < 0.7

    # The sky is always colder than the air, and much colder under a clear sky than under
    # cloud — that gap is what drives night-time radiative cooling.
    @test all(T_sky .< T_air)
    clear_gap = T_air[findfirst(h -> h ≈ 12.0, hours)] - T_sky[findfirst(h -> h ≈ 12.0, hours)]
    cloudy_gap = T_air[findfirst(h -> h ≈ 16.0, hours)] - T_sky[findfirst(h -> h ≈ 16.0, hours)]
    @test clear_gap > cloudy_gap + 3
    @printf("  sky depression: %.1f K under clear sky, %.1f K under 0.65 cover\n",
            clear_gap, cloudy_gap)

    # A wind from due south travels due north.
    @test all(abs.(Wx) .< 1e-9)
    @test all(Wy .≈ 5.0)
    @test all(Wz .== 0)

    # Dew point below the air temperature by the right amount for 70 % RH.
    @test all(Tdew .< T_air)
    @test all(3.5 .< (T_air .- Tdew) .< 6.0)
    @printf("  dew point %.2f K below air at 70%% RH\n", sum(T_air .- Tdew) / length(ts))

    # The air heat port carries the air temperature.
    @test all(grab(m.env.Air_Temp.T) .≈ T_air)
end

end
