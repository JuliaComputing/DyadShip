# DyadShip

A [Dyad](https://help.juliahub.com/dyad/) port of the Modelica naval-architecture
library [ShipSIM](https://github.com/BasilioPV/ShipSIM) by Basilio Puente and
M Dolores Fernandez, built on the 3D `MultibodyComponents` library, plus a
WaterLily.jl CFD-driven Flettner rotor propulsor.

> **This is a rewrite, not a translation.** The Modelica components have been
> reimplemented in Dyad on top of `MultibodyComponents`, `RotationalComponents`,
> `BlockComponents` and the rest of the Dyad standard libraries. Per-component
> docstrings state what was simplified or corrected relative to upstream;
> `PORTING_NOTES.md` is the map of what exists, what changed and what was left out.

## Six-degree-of-freedom ship (`dyad/Ship6DOF`)

The primary stack mirrors ShipSIM's architecture on `Frame3D` connectors:

| Component | What it does |
|---|---|
| `ShipBody` | Rigid body (yaw-pitch-roll Euler angles) with draft-polynomial hydrostatics: displacement, centre of buoyancy and metacentric radii as functions of the instantaneous draft, heel and trim. |
| `HydrodynamicXYY` | MMG surge/sway/yaw forces with the empirical derivative estimates of Clarke, Smitt, Khattab, Lee & Shin, Kijima and Yoshimura, resistance curve, and added mass with the coupling terms. Forces act at the centre of forces, so a turn heels the hull. |
| `HydrodynamicZRP` | Heave/roll/pitch damping and added mass, with defaults from damping ratios. |
| `Propeller1Q` / `Propeller4Q` | Wageningen B-series open-water characteristics: the full Oosterveld & van Oossanen polynomial, or the 14 four-quadrant Fourier data sets for astern and crash-stop work. Both feed the rudder with Brix's slipstream model. |
| `Rudder` | Steering-gear angle and rate limits, NACA 0012/0015 `Cl/Cd/Cm(α, Re)` tables (`assets/naca*.csv`), flow straightening, Söding slipstream and hull-interaction factors. |
| `ShipWind` | Fujiwara superstructure wind loads at the centre of the lateral area. |
| `WingSail` | Rigid wing sail on a servo-driven mast revolute, NACA tables, forces at the quarter chord of the rotated sail. |
| `POD4Q` | Azimuthing pod: servo revolute, strut and an internal `Propeller4Q` that reads its own advance speed. |
| `AntiHeeling` / `BallastTank` / `VariableMass` | Clocked hysteresis pump controller with ramped flow and a latched transfer direction; the tank liquids as variable-mass points on the hull, so the righting moment, the added weight and the change of the ship's centre of gravity and roll inertia come from where the water is. |
| `AntiHeelingCircuit` / `CentrifugalPump` / `TankLiquidMass` | The same tanks as `IncompressibleFlowComponents` open tanks (their liquid a `TankLiquidMass` on the hull) joined by two antiparallel affinity-law pump and linear-valve branches: the transfer rate follows the pump characteristic against the level difference and the valve losses. A partial the ship assembly extends, because fluid ports cannot pass through a wrapper. |
| `FuelTank` | Consumable liquid on `VariableMass`: burn or bunker fuel and the hull's draft, trim and inertia follow. |
| `RollTunedMassDamper` | The building tuned mass damper for roll: a mass on a transverse `Prismatic` slide high in the ship with a `SpringDamper` tuned to the roll period (Den Hartog); damps roll, does not correct static heel. |
| `Crane` / `Cable` | Slewing and luffing revolutes on position servos, boom, tension-only cable with a clocked break latch; loads react into the hull. |
| `ApparentSpeedXY` | Frame-based apparent wind / current sensor. |
| `WaypointAutopilot` / `WaypointSequencer` | `LimPID` heading autopilot with throttle ramp, and a clocked waypoint table that advances its index on arrival (`DiscreteComponents`). |
| `StandardShip` | The ShipSIM sample hull (100 m, 5681 t) with propeller, rudder and hydrodynamics wired up; the manoeuvring analyses extend it. |

### Validation

`scripts/validate_6dof.jl` runs the standard manoeuvring tests and writes the
figures in `assets/ship6dof_*.png`. Results for the sample hull (rudder rate
2.5 °/s, approach 6.69 m/s, 100 rpm unless noted):

| Analysis | Result |
|---|---|
| `RollDecayTransient` | roll period 9.1 s (linear estimate 8.4 s; the difference is the sway added-mass coupling of a hull rolling about a CoG 5 m above its hydrodynamic centre), damping ratio 0.047 for the 0.05 setting |
| `RollDecayTMDTransient` | same release with a 100 t tuned mass damper 15 m above the roll centre: first peaks 4.6°, 1.6°, 0.7° against 8.5°, 6.2°, 4.6°, heel 0.1° at 60 s; the mass strokes 6.4 m |
| `SpeedTrialTransient` | 6.05 m/s at 100 rpm, thrust 129 kN against 128 kN resistance, 1.1 MW shaft power |
| `TurningCircleTransient` (35°) | advance 3.9 L, transfer 0.7 L, tactical diameter 3.2 L; steady turn at 44 % of approach speed, 0.95 °/s, 35° drift, 0.35° outward heel |
| `RudderReturnTransient` | yaw rate halves 24 s after the rudder is centred and decays to zero |
| `ZigZagTransient` (20/20) | overshoot angles 32°, 31°, 27° |
| `CrashStopTransient` | 110 rpm ahead to 80 rpm astern: stopped after 287 s, head reach 9.9 L |
| `FullShip6DOFTransient` | 10 km waypoint transit in a 10 m/s wind from the north-east: arrives at 1768 s holding a 3.8° rudder offset and 0.2° heel |
| `WingSailSweepTransient` | one sail in a 10 m/s beam wind: peak forward thrust 12 kN at a 15° attack angle |
| `FourWingSailsTransient` | four sails, 15 m/s from the port beam, autopilot holding course: 7.17 m/s at 100 rpm against 6.05 m/s without sails, 29 kN sail thrust, 1.2° heel |
| `FourWingSailsAHTransient` | same with the anti-heeling system enabled at 300 s: heel back to zero by 750 s, tank levels 65 % port / 25 % starboard |
| `FourWingSailsAHTanksTransient` | same with the tanks as moving masses instead of an applied torque: heel back to zero by 700 s with 67 t / 25 t in the tanks, ship 7 cm deeper from the ballast |
| `FourWingSailsAHCircuitTransient` | same with the tanks as `IncompressibleFlowComponents` open tanks and a pump/valve circuit: the pump runs at 200 m³/h (57 kg/s) against a 3.0 to 4.7 m head, heel 1.15° to 0.005° by 700 s, tank fills 0.45 / 0.45 to 0.65 / 0.25 |
| `BunkeringTransient` | 300 t of fuel loaded forward in one hour at rest: draft +0.20 m, trim 0.52° by the bow, matching the hydrostatic estimate |
| `MoistAir.MoistAirDewPointTransient` | 30 °C, 80 % air through a ventilated space on `HVACComponents` moist-air media: outlet dew point 26.20 °C (Magnus reference 26.17 °C) |
| `PodTurningCircleTransient` | 35° pod azimuth: steady radius 1.0 L at 2.1 m/s, no cavitation |
| `CraneOperationTransient` | 50 t load luffed, slewed 90° to port and lowered: ship heels to 2.1°, cable tension 490 kN |
| `WaypointTransitTransient` | three-leg route with a 45° dog-leg: the clocked sequencer switches waypoints at 851 s and 1635 s, arrival at 2400 s |

![Turning circle](assets/ship6dof_turning_circle.png)

The zig-zag overshoots are large because the upstream Khattab estimate of the
yaw damping leaves the bare hull linearly course-unstable
(`HydrodynamicXYY.CourseStability < 0`); the rudder's fin effect holds the
course. Override `N_r` for a stiffer hull.

The zig-zag switch, the anti-heeling on/off controller, the cable break and
the waypoint table are clocked components built on `DiscreteComponents`
(sampled every 0.05–1 s on a `PeriodicClock`, with the state held between
samples), not continuous relay approximations.

`ManualShip6DOFTransient` exposes shaft rpm and rudder angle as tunable
parameters for interactive (WASM) use.

## Planar stack and Flettner rotor

`dyad/Ship` and `dyad/Propulsion` hold the earlier `PlanarMechanics`
(surge/sway/yaw) port: `HullMMG`, `Rudder`, `Propeller1Q/4Q`, `ShipWind`,
`HeadingAutoPilot`, `ManualShip`, the `FullShip*` transits and the
Flettner-rotor propulsor. The planar wing sail, pod, crane, cable and
anti-heeling components were replaced by the 3D versions. Two sign errors in it were fixed
in this pass (rudder inflow angle, wind lateral force and moment); its
analyses still use a hull mass well below the sample ship's displacement, so
prefer `Ship6DOF` for manoeuvring studies.

`Propulsion.SimpleDieselEngine` meters fuel on the brake power only: the
engine torque is never negative, so a shaft driven backwards against it is a
motored engine that reports negative `ShaftPower` and no fuel, never a fuel
credit (upstream integrates the signed sensed power and goes negative). The
SFOC table is a parameter pair (`SFOC_P`, `SFOC_g`) whose end values are held
outside the tabulated 605–1210 kW. `SFOC_valid` is a table-coverage flag (0
while a held end value is in use), not a calibration statement, and
`Fuel_extrapolated` [kg] accumulates the fuel metered outside the table, so
an interval of the `Fuel` counter is inside the table exactly when
`Fuel_extrapolated` did not change over it. `m_dot_idle` adds a no-load rate
and defaults to zero as upstream. `Inst_Fuel` [kg/s] is nonnegative, `Fuel`
[kg] is exactly its integral and nondecreasing, `KWh` stays the signed net
work. The default tables are the upstream ShipSIM values, for which upstream
cites no engine or test: they are reference values, not an OEM calibration
(sources and manufacturer documents for context are linked in the component
docstring). `assert`s reject unordered, non-finite or non-positive tables and
a negative idle rate.

The fuel counters have been nondecreasing at every accepted solver step in
every run checked. For `Fuel_extrapolated`, whose rate jumps to zero where
the brake power enters the table, the engine ends a solver step at each table
entry (`Propulsion.StepAtRisingCrossing`), so the steps that follow lie
inside the table and leave the counter exactly constant.
Interpolated values, including a `saveat` grid, which stores the interpolant
rather than extra steps, can fall back inside a step that spans a kink of the
rate (4.6e-5 kg on 0.41 kg in `DieselEngineReversing` without step
alignment; the size follows the step length, not the solver tolerance).
`DieselEngineReversing` declares its zero-power crossings as `tstops`, which
makes its dense fuel counter nondecreasing too; where crossing times are not
known, export counters at the accepted steps. The metering is evaluated
pointwise and does not depend on `automatic_discontinuity_detection`, whose
events miss a condition that starts exactly on its threshold (an engine
starting at zero power).
`scripts/validate_engine_fuel.jl` checks the forward ramp against the
previous numbers, a reversing shaft imposed by a velocity source, a run held
inside the SFOC table from start to end (`DieselEngineInRange`, synthetic
load, rate samples reconciled with the counter), a stopped shaft with and
without idle fuel, the sampling behaviour and the rejected
parameter values.

### How a powered ship comes together

`Ship6DOF.PoweredSingleScrewShip` is the sample hull with its propeller,
rudder and a free propeller shaft (`SingleScrewHull`, which is `StandardShip`
without its ideal speed governor) plus a distance counter. A ship is that
partial with one `plant` subcomponent connected to `shaft.spline_a`.

A plant is any component that extends `Ship6DOF.ShipPowerPlantPorts`, which
fixes what the hull needs and what the machinery reports: the
`propeller_flange`, the orders `shaft_speed_order` [rpm] and
`hotel_power_demand` [W], and the observation outputs in SI units (fuel mass
rate, fuel counter, out-of-table fuel counter and coverage flag for the main
engine and for the generating set; shaft, hotel, generating-set and
shaft-machine power with their energy integrals; running states and running
times). Dyad has no replaceable components, so the slots inside a plant are
conventional subcomponent names (`main_engine`, `gearbox`, `pto`,
`genset_engine`, `genset_machine`, `hotel_load`) with a partial component for
the connectors of each kind: `Propulsion.DieselEnginePorts` for an engine,
`Ship6DOF.ShaftMachinePorts` for an electrical machine on a shaft.

Two plants are built from this library and the standard component libraries
only: `DieselMechanicalPlant` (synthetic engine preset A through an ideal
gear to the propeller, synthetic preset B driving a DC machine that feeds a
conductance hotel load) and `DieselShaftMachinePlant`, the same with a
`DCShaftMachine` on the engine shaft and the bus. That machine is an ideal
EMF behind a resistance with no control, so it takes power from the bus
below its matching engine speed and feeds the bus above it. Everything in
them is synthetic or ideal; they show the structure, not an installation.
`scripts/validate_ship_power_plant.jl` checks signs, the power balance from
engines to propeller and hotel load, and the counters.

A plant written in another package extends these partials directly. The
generated code of the extending package names the libraries that the
inherited hull composes, so its module must import them as well:
`MultibodyComponents`, `RotationalComponents`, `BlockComponents` and
`ElectricalComponents`, not only the libraries its own components name.

### A tanker skeleton

`Ship6DOF.TankerSkeleton` puts the pieces of a ship in one assembly: the
powered hull above with its body reduced to a lightship mass, three cargo
tanks and two peak ballast tanks (`FuelTank`, used as a liquid mass with no
flow), two wing ballast tanks joined by the pump-and-valve circuit
`AntiHeelingCircuit`, a bunker tank drained by the fuel rates of the plant's
engines, and `DieselMechanicalPlant`. `TankerLaden` and `TankerBallast` are
two load conditions; the analyses are a speed trial in each and a ballast
transfer under way. All numbers are synthetic round values on the library's
100 m sample hull, which is not a tanker hull.

Each tank is a point mass at half the liquid height on the multibody hull,
so draft, trim and heel follow the loading. The tank models have no free
surface, no sloshing, no cargo piping and no flooding; the docstring lists
what each piece can and cannot represent.
`scripts/validate_tanker_skeleton.jl` checks the mass bookkeeping (the bunker
tank loses exactly what the engines meter; a ballast transfer conserves
water), the hydrostatic response and the power bookkeeping.

### A published tanker hull: KVLCC2

`Ship6DOF.KVLCC2Ship` is the KVLCC2 research tanker at full scale (320 m,
312 600 m³). Its hull forces come from the file
`assets/presets/Hull/kvlcc2_full_scale.toml`: the principal particulars,
resistance coefficient, added masses and sixteen hull derivatives of
Yasukawa and Yoshimura, "Introduction of MMG standard method for ship
maneuvering predictions", J. Mar. Sci. Technol. 20 (2015) 37–52 (open
access). The sidecar next to the file gives the table and page of every
value and the arithmetic of the converted ones. `Ship6DOF.Propeller1Q` takes
the published open-water thrust polynomial through its new optional
polynomial parameters.

The source's model has three degrees of freedom and its own rudder and wake
models. This one has six degrees of freedom, the library's `Rudder` and a
constant wake fraction; the vertical hydrostatics, the height of the centre
of gravity and the positions of propeller and rudder are stated assumptions.
With nothing tuned, the 35° turning circles from 15.5 kn come out with an
advance of 4.04 (port) and 4.16 (starboard) ship lengths and a tactical
diameter of 3.75 and 3.98, against 3.56 / 3.62 and 3.59 / 3.71 in the
source's own full-scale simulation: the advance is 13–15 % larger and the
tactical diameter 4–7 % larger. Against the source's free-running test of a
7 m model, which has no full-scale counterpart, they are 28–30 % and
19–22 % larger. `scripts/validate_kvlcc2.jl` prints the comparison.

`Propulsion.SyntheticDieselEngineA` and `SyntheticDieselEngineB` are parameter
presets of the one `SimpleDieselEngine`: a wrapper without equations holds the
engine as `core` and loads `assets/presets/Synthetic/diesel_engine_{a,b}.toml`
with an `apply` clause, which is the path Dyad 3.4 supports for a table whose
length differs from the default (six and three knots here; the file sets the
structural `n_sfoc` together with the arrays). Both files are **synthetic**:
round numbers invented to exercise the mechanism, not measurements,
manufacturer data or a calibration, as their `*.provenance.toml` sidecars
record. `scripts/validate_engine_presets.jl` solves both and checks the
metered fuel against the files; `test/invalid_engine_presets` holds files that
must be refused and shows where each kind is caught (compiler diagnostic with
exit status 0, construction, or the engine's assertions). A table whose
length disagrees with `n_sfoc` is refused when the engine is constructed.

The `FlettnerRotor` component reads `Cl(ξ)`, `Cd(ξ)` from
`assets/flettner_coeffs.csv`, produced offline by
`scripts/run_waterlily_flettner.jl` from WaterLily.jl simulations of a
rotating cylinder in cross-flow (G. D. Weymouth's SpinCyl pattern). Three
rendered transits compare diesel-only, rotor with bow-quarter wind and rotor
with beam wind:

| Analysis | Rotor | Mean wind from | Time to target |
|---|---|---|---|
| `ShipRenderTransient` | none | NE (45°) | not reached by 2200 s |
| `ShipFlettnerRenderTransient` | yes | NE (45°), bow quarter | not reached by 1500 s |
| `ShipFlettnerFavorableRenderTransient` | yes | N (0°), beam | reaches at t ≈ 1200 s |

<video src="assets/ship_flettner_favorable_animation.mp4" controls width="720">
  <a href="assets/ship_flettner_favorable_animation.mp4">assets/ship_flettner_favorable_animation.mp4</a>
</video>

`FlettnerRotorOnline` and `scripts/run_flettner_cosim_verification.jl` run
the same transit with WaterLily stepped live in a `PeriodicCallback`
(package extension `FlettnerCFDLiveExt`, GPU-capable through CUDA) to
verify the table approach; `assets/flettner_cosim_*` hold the comparison
plots and animations. Table and live CFD agree on shape, sign and
trajectory, with live magnitudes 20–30 % lower during the ramp-up.

## Getting started

Dyad models live in `dyad/`; the Dyad compiler emits Julia into `generated/`
(never edit those files). The wrappers `../julia-dyad.sh` and `../dyad.sh`
select the `dyad-3.4.0` JuliaUp channel and `dyad-cli@3.4.0`; run heavy
commands through `~/dyad-fleet/heavy` on this machine.

```sh
# Compile Dyad -> Julia
../dyad.sh compile

# Run an analysis from Julia
../julia-dyad.sh -e 'using DyadShip; res = DyadShip.Ship6DOF.TurningCircleTransient(); println(res.sol.retcode)'

# Manoeuvring validation report + figures
../julia-dyad.sh scripts/validate_6dof.jl

# Convection factors, external wall, ship compartment and the weather boundary
../julia-dyad.sh scripts/validate_heattransfer.jl

# Electrical load analysis and rainflow counting
../julia-dyad.sh scripts/validate_machinery.jl

# Planar transit animations
../julia-dyad.sh scripts/render_all.jl

# Re-characterise the Flettner rotor with WaterLily (scripts/Project.toml)
../julia-dyad.sh --project=scripts scripts/run_waterlily_flettner.jl
```

Two dependencies are not open source. `HVACComponents` (moist air, used by
`dyad/MoistAir`) and `IncompressibleFlowComponents` (tanks, pumps and valves,
used by `Ship6DOF.AntiHeelingCircuit`) are © JuliaHub, all rights reserved,
and are provided under the JuliaHub end user license agreement
(https://juliahub.com/company/eula) from the private `DyadHVACRegistry` and
`DyadThermoFluidRegistry` registries; see `PORTING_NOTES.md` for how they are
installed. They are not redistributed here, and using the models that
depend on them requires access to those registries under that agreement.

Accessing results follows the Dyad convention:

```julia
using DyadShip
using DyadInterface: symbolic_container
res = DyadShip.Ship6DOF.ZigZagTransient()
m = symbolic_container(res)
heading_deg = rad2deg.(res.sol[m.ship.Yaw])
```

## Layout

- `dyad/Ship6DOF/` — 6-DOF ship stack, analyses, `definitions.jl` (Wageningen
  polynomials, four-quadrant Fourier sets, draft polynomials).
- `dyad/Ship/`, `dyad/Propulsion/` — planar stack, Flettner rotor, transits.
- `dyad/Machinery/` — on-off consumer, the seeded clocked `RandomStart` scheduler
  and the `ElectricalLoad` bank that make up a ship's electrical power load
  analysis, plus the continuous peak sampler and the event-driven
  `EventPeakSampler` built on `DiscreteComponents` clocks.
- `dyad/Thermal/` — solar irradiation, sun screen, plate and cylinder
  transients, temperature dataset, air exchanger, the four convection factors
  (horizontal cylinder, forced flat plate, internal and external wall surfaces)
  on `ThermalComponents.Interfaces.ConvectiveElement1D`, and `ConvRadSunWall`,
  an external wall composing those with two view-factor-split
  `ThermalComponents.BodyRadiation` paths, `IrradiationOnPlane` and `SunScreen`.
  `dyad/Thermal/ShipCompartment.dyad` is a moist-air compartment behind
  weather-exposed bulkheads whose film coefficients are those convection
  factors rather than fixed `U` values, so its heat loss follows the ship's
  apparent wind.
- `dyad/MoistAir/` — `SourceMoistAir` and `DewTemperature` on the
  `HVACComponents` moist-air medium (`MoistAirFluidPort`).
- `dyad/Environment.dyad`, `VariableEnvironment.dyad`, `ApparentSpeedXY.dyad` —
  signal-level environment helpers.
- `assets/` — wing-profile and Flettner coefficient tables, temperature CSV,
  validation figures, animations.
- `scripts/` — validation, rendering and WaterLily characterisation scripts
  (`scripts/Project.toml` carries the WaterLily/CUDA dependencies).
- `src/DyadShip.jl`, `src/Rainflow.jl`, `ext/FlettnerCFDLiveExt.jl` — Julia
  module wrapper, rainflow cycle counting and Miner damage over the turning
  points `EventPeakSampler` emits, and the live-CFD extension.
- `PORTING_NOTES.md` — port status, conventions, corrections relative to upstream.
- `AGENTS.md` — toolchain notes and Dyad/MTK learnings for agents.

## License

The Dyad rewrite in this repository is © 2025 JuliaHub and contributors, Panagiotis Georgakopoulos. The
upstream Modelica `ShipSIM` library is © Basilio Puente and M Dolores
Fernandez, distributed under the 3-clause BSD license.

The `HVACComponents` and `IncompressibleFlowComponents` libraries this
repository depends on are © JuliaHub, all rights reserved, and are subject to
the JuliaHub end user license agreement rather than the licenses above; the
`MoistAir` module and the `AntiHeelingCircuit` models are unusable without
them. The Dyad toolchain itself is provided by JuliaHub for educational and
personal use, with commercial use requiring a license.
