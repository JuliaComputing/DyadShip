#!/usr/bin/env julia
# Run the Ship6DOF validation analyses and summarise them the way a
# manoeuvring report would: roll period and damping, speed-trial balance,
# IMO turning-circle metrics, rudder-return yaw-rate decay, zig-zag
# overshoots, the crash stop, the wing-sail polar and sailing ship, the pod
# turning circle, the crane operation and the autopilot transit. Writes PNG
# figures into assets/.
#
# Usage:  ~/dyad-fleet/heavy -m 6G ../julia-dyad.sh scripts/validate_6dof.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using DyadShip, Printf, Statistics
using DyadInterface: symbolic_container
using Plots
gr()

const ASSETS = abspath(joinpath(@__DIR__, "..", "assets"))
const S6 = DyadShip.Ship6DOF
mkpath(ASSETS)

sample(sol, sym, ts) = [sol(t; idxs = sym) for t in ts]

println("== Roll decay")
res = S6.RollDecayTransient(); sol = res.sol; m = symbolic_container(res)
ts = collect(0:0.05:120); heel = rad2deg.(sample(sol, m.ship.Heel, ts))
zc = [ts[i] for i in 2:length(ts) if heel[i-1] > 0 && heel[i] <= 0]
period = mean(diff(zc)); peaks = [maximum(abs.(heel[(ts .>= zc[i]) .& (ts .< zc[i+1])])) for i in 1:length(zc)-1]
decrement = mean(log.(peaks[1:end-1] ./ peaks[2:end])); zeta = decrement / sqrt(4pi^2 + decrement^2)
@printf("  retcode %s, roll period %.2f s (linear estimate 2π k_xx/sqrt(g GM) = %.2f s), damping ratio %.3f\n", sol.retcode, period, 2pi * 0.3557 * 20 / sqrt(9.80665 * 2.89), zeta)
p = plot(ts, heel, xlabel = "t [s]", ylabel = "heel [deg]", label = "heel", title = "Roll decay from 10 deg", lw = 2)
plot!(p, ts, rad2deg.(sample(sol, m.ship.Trim, ts)) .* 100, label = "trim x100")
savefig(p, joinpath(ASSETS, "ship6dof_roll_decay.png"))

println("== Speed trial")
res = S6.SpeedTrialTransient(); sol = res.sol; m = symbolic_container(res)
u_end = sol(600; idxs = m.ship.Surge)
@printf("  retcode %s, 100 rpm: %.2f m/s (%.1f kn), thrust %.0f kN, resistance %.0f kN, shaft power %.0f kW, J = %.3f\n", sol.retcode, u_end, u_end / 0.5144, sol(600; idxs = m.prop.Thrust) / 1e3, -sol(600; idxs = m.hydro.Resistance) / 1e3, sol(600; idxs = m.prop.ShaftPower) / 1e3, sol(600; idxs = m.prop.J))

println("== Turning circle (35 deg to port)")
res = S6.TurningCircleTransient(); sol = res.sol; m = symbolic_container(res)
ts = collect(0:0.5:900); x = sample(sol, m.ship.pos_x, ts); y = sample(sol, m.ship.pos_y, ts); psi = sample(sol, m.ship.Yaw, ts)
i0 = findfirst(>=(20.0), ts); x0 = x[i0]; y0 = y[i0]
i90 = findfirst(>=(pi / 2), psi); i180 = findfirst(>=(pi), psi)
advance = x[i90] - x0; transfer = y[i90] - y0; tactical = y[i180] - y0
U = hypot(sol(900; idxs = m.ship.Surge), sol(900; idxs = m.ship.Sway)); r = sol(900; idxs = m.ship.YawRate)
drift = rad2deg(atan(-sol(900; idxs = m.ship.Sway), sol(900; idxs = m.ship.Surge)))
@printf("  retcode %s, advance %.0f m (%.2f L), transfer %.0f m (%.2f L), tactical diameter %.0f m (%.2f L)\n", sol.retcode, advance, advance / 100, transfer, transfer / 100, tactical, tactical / 100)
@printf("  steady turn: speed %.2f m/s (%.0f%% of approach), yaw rate %.2f deg/s, radius %.0f m, drift %.1f deg, heel %.2f deg\n", U, 100U / 6.69, rad2deg(r), U / r, drift, rad2deg(sol(900; idxs = m.ship.Heel)))
p = plot(x .- x0, y .- y0, aspect_ratio = 1, xlabel = "advance [m]", ylabel = "transfer [m]", label = "track", title = "Turning circle, 35 deg rudder", lw = 2)
hline!(p, [tactical], ls = :dash, label = @sprintf("tactical diameter %.2f L", tactical / 100))
savefig(p, joinpath(ASSETS, "ship6dof_turning_circle.png"))

println("== Rudder return (20 deg pulse, 50-150 s)")
res = S6.RudderReturnTransient(); sol = res.sol; m = symbolic_container(res)
ts = collect(0:0.5:600); r6 = rad2deg.(sample(sol, m.ship.YawRate, ts))
r_off = sol(150; idxs = m.ship.YawRate); t_half = ts[findfirst(i -> ts[i] > 150 && abs(sample(sol, m.ship.YawRate, [ts[i]])[1]) < 0.5abs(r_off), 1:length(ts))]
@printf("  retcode %s, yaw rate at rudder centring %.2f deg/s, halved after %.0f s, %.3f deg/s at 600 s\n", sol.retcode, rad2deg(r_off), t_half - 150, r6[end])
p = plot(ts, r6, xlabel = "t [s]", ylabel = "yaw rate [deg/s]", label = "Ship6DOF", title = "Rudder 20 deg from 50 s, centred at 150 s", lw = 2)
res2 = DyadShip.Ship.RudderReturn2DTransient(); sol2 = res2.sol; m2 = symbolic_container(res2)
plot!(p, ts, rad2deg.(sample(sol2, m2.hull.r, ts)), label = "planar HullMMG", lw = 2)
savefig(p, joinpath(ASSETS, "ship6dof_rudder_return.png"))

println("== Zig-zag 20/20")
res = S6.ZigZagTransient(); sol = res.sol; m = symbolic_container(res)
ts = collect(0:0.25:600); psi = rad2deg.(sample(sol, m.ship.Yaw, ts)); rud = sample(sol, m.rudder.Rudder_position, ts)
# Heading extremes beyond the switching angle; the clocked relay state itself
# is only available per tick through sol[m.controller.relay.s]
ext = [i for i in 2:length(ts)-1 if (psi[i] - psi[i-1]) * (psi[i+1] - psi[i]) < 0 && abs(psi[i]) > 20]
overshoots = [abs(psi[i]) - 20 for i in ext[1:min(end, 4)]]
nrev = count(i -> sol[m.controller.relay.s][i] != sol[m.controller.relay.s][i-1], 2:length(sol[m.controller.relay.s]))
@printf("  retcode %s, %d relay reversals, heading extremes at %s s, overshoot angles %s deg\n", sol.retcode, nrev, join(round.(ts[ext[1:min(end, 4)]], digits = 0), ", "), join(round.(overshoots, digits = 1), ", "))
p = plot(ts, psi, xlabel = "t [s]", ylabel = "[deg]", label = "heading", title = "20/20 zig-zag", lw = 2)
plot!(p, ts, rud, label = "rudder", lw = 2)
savefig(p, joinpath(ASSETS, "ship6dof_zigzag.png"))

println("== Crash stop (four-quadrant propeller, 110 rpm ahead to 80 rpm astern)")
res = S6.CrashStopTransient(); sol = res.sol; m = symbolic_container(res)
ts = collect(0:0.5:300); u = sample(sol, m.ship.Surge, ts); istop = findfirst(<=(0), u)
@printf("  retcode %s, %s\n", sol.retcode, istop === nothing ? "not stopped within 300 s" : @sprintf("stopped after %.0f s, head reach %.0f m (%.2f L)", ts[istop], sol(ts[istop]; idxs = m.ship.pos_x), sol(ts[istop]; idxs = m.ship.pos_x) / 100))
p = plot(ts, u, xlabel = "t [s]", ylabel = "[m/s], [rpm/10]", label = "surge speed", title = "Crash stop", lw = 2)
plot!(p, ts, sample(sol, m.prop.rpm, ts) ./ 10, label = "shaft rpm / 10", lw = 2)
savefig(p, joinpath(ASSETS, "ship6dof_crash_stop.png"))

println("== Wing sail sweep, 10 m/s beam wind")
res = S6.WingSailSweepTransient(); sol = res.sol; m = symbolic_container(res)
ts = collect(0:0.5:185); ang = sample(sol, m.sail.Sail_position, ts); fx = sample(sol, m.mount.frame_b.f[1], ts) ./ 1e3; fy = sample(sol, m.mount.frame_b.f[2], ts) ./ 1e3
ib = argmax(fx)
@printf("  retcode %s, peak forward thrust %.1f kN at sail angle %.0f deg (attack angle %.1f deg), lateral force there %.1f kN\n", sol.retcode, fx[ib], ang[ib], rad2deg(sample(sol, m.sail.AttackAngle, [ts[ib]])[1]), fy[ib])
p = plot(ang, fx, xlabel = "sail angle [deg]", ylabel = "force on mount [kN]", label = "forward", title = "Wing sail polar, 10 m/s beam wind from port", lw = 2)
plot!(p, ang, fy, label = "lateral (+ port)", lw = 2)
savefig(p, joinpath(ASSETS, "ship6dof_wingsail_sweep.png"))

println("== Four wing sails, 15 m/s wind from the port beam")
res = S6.FourWingSailsTransient(); sol = res.sol; m = symbolic_container(res)
fx = sum(sol(600; idxs = getproperty(m, Symbol("sail$i")).force.force_x) for i in 1:4)
@printf("  retcode %s, speed %.2f m/s at 100 rpm (6.05 m/s without sails), sail thrust %.0f kN, attack angle %.1f deg, heel %.2f deg\n", sol.retcode, sol(600; idxs = m.ship.Surge), fx / 1e3, rad2deg(sol(600; idxs = m.sail1.AttackAngle)), rad2deg(sol(600; idxs = m.ship.Heel)))
res = S6.FourWingSailsAHTransient(); sol = res.sol; m = symbolic_container(res)
ts = collect(0:1:900)
@printf("  with anti-heeling (pump enabled at 300 s): heel %.2f deg at 290 s, %.2f deg at 900 s; tank levels %.0f%% / %.0f%%\n", rad2deg(sol(290; idxs = m.ship.Heel)), rad2deg(sol(900; idxs = m.ship.Heel)), sol(900; idxs = m.tank_port.fill_level), sol(900; idxs = m.tank_starboard.fill_level))
p = plot(ts, rad2deg.(sample(sol, m.ship.Heel, ts)), xlabel = "t [s]", ylabel = "[deg], [m³/h / 100]", label = "heel", title = "Four wing sails with anti-heeling from 300 s", lw = 2)
plot!(p, ts, sample(sol, m.antiheeling.pump_flow, ts) ./ 100, label = "pump flow / 100", lw = 2)
savefig(p, joinpath(ASSETS, "ship6dof_sails_antiheeling.png"))

res = S6.FourWingSailsAHTanksTransient(); sol = res.sol; m = symbolic_container(res)
@printf("  with BallastTank masses instead of an applied torque: heel %.2f deg at 290 s, %.2f deg at 900 s; tank masses %.0f t / %.0f t, draft +%.3f m\n", rad2deg(sol(290; idxs = m.ship.Heel)), rad2deg(sol(900; idxs = m.ship.Heel)), sol(900; idxs = m.tank_port.mass) / 1e3, sol(900; idxs = m.tank_starboard.mass) / 1e3, sol(900; idxs = m.ship.Draft) - 4)
ts = collect(0:1:900)
p = plot(ts, rad2deg.(sample(sol, m.ship.Heel, ts)), xlabel = "t [s]", ylabel = "[deg], [t / 100]", label = "heel", title = "Anti-heeling with moving ballast masses", lw = 2)
plot!(p, ts, sample(sol, m.tank_port.mass, ts) ./ 1e5, label = "port tank mass / 100 t", lw = 2)
plot!(p, ts, sample(sol, m.tank_starboard.mass, ts) ./ 1e5, label = "starboard tank mass / 100 t", lw = 2)
savefig(p, joinpath(ASSETS, "ship6dof_ballast_tanks.png"))

println("== Moist air: dew-point sensor on a ventilated space (HVACComponents)")
res = DyadShip.MoistAir.MoistAirDewPointTransient(); sol = res.sol; m = symbolic_container(res)
Ps = 610.78 * exp(17.558 * 30 / (241.88 + 30)); Pv = 0.8 * Ps; Tref = 241.88 * log(Pv / 610.78) / (17.558 - log(Pv / 610.78))
@printf("  retcode %s, outlet dew point %.2f C at 900 s (Magnus reference for 30 C / 80 %%: %.2f C), outlet RH %.3f\n", sol.retcode, sol(900; idxs = m.dew.Tdew) - 273.15, Tref, sol(900; idxs = m.dew.RH))

println("== Pod turning circle (35 deg azimuth)")
res = S6.PodTurningCircleTransient(); sol = res.sol; m = symbolic_container(res)
U = hypot(sol(600; idxs = m.ship.Surge), sol(600; idxs = m.ship.Sway)); r = sol(600; idxs = m.ship.YawRate)
@printf("  retcode %s, steady radius %.0f m (%.2f L), speed %.2f m/s, pod thrust %.0f kN, cavitation ratio %.2f\n", sol.retcode, abs(U / r), abs(U / r) / 100, U, sol(600; idxs = m.pod.Thrust) / 1e3, sol(600; idxs = m.pod.prop.Cavitation_Warning))

println("== Crane operation, 50 t load")
res = S6.CraneOperationTransient(); sol = res.sol; m = symbolic_container(res)
ts = collect(0:0.5:400); heel = rad2deg.(sample(sol, m.ship.Heel, ts)); tension = sample(sol, m.crane.CableTension, ts) ./ 1e3
@printf("  retcode %s, heel from %.2f to %.2f deg after slewing to port, cable tension %.0f-%.0f kN\n", sol.retcode, heel[1], heel[end], minimum(tension[tension .> 0]), maximum(tension))
p = plot(ts, heel, xlabel = "t [s]", ylabel = "[deg], [m]", label = "heel", title = "Crane: luff 20-40 s, slew 45-100 s, pay out 100-300 s", lw = 2)
plot!(p, ts, sample(sol, m.load.r_0[3], ts), label = "load height", lw = 2)
savefig(p, joinpath(ASSETS, "ship6dof_crane.png"))

println("== Multi-waypoint transit (clocked waypoint sequencer), 10 m/s wind from NE")
res = S6.WaypointTransitTransient(); sol = res.sol; m = symbolic_container(res)
idx = sol[m.route.counter.i]; sw = findall(i -> idx[i] != idx[i-1], 2:length(idx)) .+ 1
ts = collect(0:2:3600); x = sample(sol, m.ship.pos_x, ts); y = sample(sol, m.ship.pos_y, ts)
@printf("  retcode %s, waypoint switches at samples %s of the 1 s clock, track ends at (%.0f, %.0f)\n", sol.retcode, string(sw), x[end], y[end])
p = plot(x, y, aspect_ratio = 1, xlabel = "x [m]", ylabel = "y [m]", label = "track", title = "Three-waypoint transit", lw = 2)
scatter!(p, [5000, 8000, 12000], [0, 3000, 3000], label = "waypoints", ms = 6)
savefig(p, joinpath(ASSETS, "ship6dof_waypoints.png"))

println("== Autopilot transit, 10 m/s wind from NE")
res = S6.FullShip6DOFTransient(); sol = res.sol; m = symbolic_container(res)
ts = collect(0:2:2400); d = sample(sol, m.pilot.distance_to_target, ts)
ia = findfirst(<=(100), d)
@printf("  retcode %s, %s, steady rudder offset %.1f deg, steady heel %.2f deg\n", sol.retcode, ia === nothing ? "not within 100 m" : @sprintf("within 100 m of the waypoint at %.0f s", ts[ia]), sol(900; idxs = m.rudder.Rudder_position), rad2deg(sol(900; idxs = m.ship.Heel)))
p = plot(sample(sol, m.ship.pos_x, ts), sample(sol, m.ship.pos_y, ts), aspect_ratio = 1, xlabel = "x [m]", ylabel = "y [m]", label = "track", title = "Waypoint transit, wind 10 m/s from NE", lw = 2)
scatter!(p, [10000], [1000], label = "waypoint", ms = 6)
savefig(p, joinpath(ASSETS, "ship6dof_transit.png"))
println("figures written to ", ASSETS)
