# Porting notes: ShipSIM → DyadShip

Status of the port of [ShipSIM](https://github.com/BasilioPV/ShipSIM) (Modelica)
to Dyad, consolidated from the earlier `REPORT.md` / `SUCCESSES.md` /
`PORT_PLAN.md` / `HARD.md` / `STILL_HARD.md` / `INVESTIGATION_yaw_resistance.md`
files. Per-component assumptions live in each component's docstring; this file
is the map.

## Two stacks

| | `dyad/Ship6DOF` (primary) | `dyad/Ship`, `dyad/Propulsion` (planar, legacy) |
|---|---|---|
| Multibody library | `MultibodyComponents` 3D (`Frame3D`, `Body`, `WorldForce`/`WorldTorque`, `FixedTranslation`) | `MultibodyComponents.PlanarMechanics` (`Frame2D`) |
| Degrees of freedom | surge, sway, heave, roll, pitch, yaw | surge, sway, yaw |
| Hydrostatics | draft-polynomial displacement, CoB, BM (upstream `ShipModelTh`) | none |
| Hull forces | `HydrodynamicXYY` (MMG + added mass), `HydrodynamicZRP` | `HullMMG` (MMG, cross terms only) |
| Propellers | `Propeller1Q` (full Wageningen B polynomial), `Propeller4Q` (14 four-quadrant Fourier sets, Burrill cavitation check), `POD4Q` (servo revolute + `Propeller4Q`) | quadratic-in-J fits |
| Rudder | servo revolute, NACA 0012/0015 `Cl/Cd/Cm(α, Re)` tables, Söding slipstream and hull factors | quadratic fits |
| Sails, deck gear | `WingSail` (servo revolute, tables), `AntiHeeling` (hysteresis relay, ramp, roll torque), `Crane` + `Cable` (two servo revolutes, tension-only cable with latching break) | — |
| Validated by | roll decay, speed trial, turning circle, rudder return, zig-zag, crash stop, sail polar, four-sail ship with anti-heeling, pod turning circle, crane operation, autopilot transit (`scripts/validate_6dof.jl`) | autopilot transit, Flettner-rotor renders |

The planar stack is kept because the Flettner-rotor work (CFD table, live
WaterLily co-simulation, the three rendered transits) is built on it. Its
`FullShip*` analyses use `mass = 1e6`, `Iz = 1e8`, far below the 5.7e6 kg /
4.3e9 kg m² of the hull the hydrodynamic estimates describe, so their
transients are faster than a real 100 m ship. New maneuvering work should use
`Ship6DOF`.

## Coordinate conventions (both stacks)

World: x east, y north, z up, free surface at z = 0
(`MultibodyComponents.World(n = [0, 0, -1])` in 3D). Ship: x forward, y port,
z up, origin at the aft perpendicular / centreline / keel in 3D (upstream
convention; CoG at `[50.364, 0, 7]` for the sample hull). Positive yaw is bow
to port, positive rudder order turns the ship to port, positive heel is
starboard down, the `Trim` output is positive bow up (upstream formula) while
the pitch Euler angle is positive bow down. Wind and current directions are
"coming from", 0° = north, 90° = east.

## Upstream component map

| ShipSIM | DyadShip | Notes |
|---|---|---|
| `Components.Environment`, `VariableEnvironment` | root `Environment`, `VariableEnvironment` | signal-only; no `inner`/`outer` in Dyad |
| `SubComponents.ApparentSpeedXY` | `Ship6DOF.ApparentSpeedXY` (frame-based), root `ApparentSpeedXY` (signal) | |
| `Ship.ShipModelTh` | `Ship6DOF.ShipBody` | yaw-pitch-roll Euler body, `draft_poly` hydrostatics |
| `Ship.HidrodynamicXYY` | `Ship6DOF.HydrodynamicXYY`, planar `HullMMG` | see corrections below |
| `Ship.HidrodynamicZRP` | `Ship6DOF.HydrodynamicZRP` | upstream is a zero-coefficient placeholder; defaults here from damping ratios |
| `Ship.ShipWind` | `Ship6DOF.ShipWind`, planar `Ship.ShipWind` | Fujiwara regression |
| `Propulsion.Propeller1Q` | `Ship6DOF.Propeller1Q`, planar `Propulsion.Propeller1Q` | 3D: exact Oosterveld & van Oossanen polynomial |
| `Propulsion.Propeller4Q` | `Ship6DOF.Propeller4Q`, planar `Propulsion.Propeller4Q` | 3D: all 14 Fourier data sets; planar: mirrored 1Q |
| `Propulsion.POD4Q` | `Ship6DOF.POD4Q` | servo revolute, strut, internal `Propeller4Q` |
| `Propulsion.Rudder` | `Ship6DOF.Rudder`, planar `Propulsion.Rudder` | 3D: 2D tables from `assets/naca*.csv` |
| `AlternativePropulsion.WingSail` | `Ship6DOF.WingSail` | servo revolute, NACA tables, forces at the rotated quarter chord |
| (new) | planar `Propulsion.FlettnerRotor`, `FlettnerRotorOnline` | WaterLily-derived Magnus coefficients |
| `AutoPilot.SimpleAutoPilot` | `Ship.HeadingAutoPilot` (PI + throttle ramp), `Ship6DOF.WaypointAutopilot` (`LimPID`) + `WaypointSequencer` (clocked waypoint table) | |
| `AntiHeelingSystem.AntiHeeling`, `Tank` | `Ship6DOF.AntiHeeling`, `Ship.Tank`, `Ship6DOF.BallastTank` / `AntiHeelingCircuit` (+ `VariableMass`, `TankLiquidMass`) | clocked on/off hysteresis with latched direction and slew-rate ramp; tanks as moving masses on the hull, optionally with a pump/valve circuit |
| (new) | `Ship6DOF.FuelTank` | consumable liquid mass on the hull |
| (new) | `Ship6DOF.RollTunedMassDamper` | inertial (tuned-mass) roll damper on a prismatic slide |
| (new) | `Ship6DOF.AntiHeelingCircuit`, `CentrifugalPump`, `TankLiquidMass` | pump/valve ballast transfer on `IncompressibleFlowComponents` |
| `Machines.SimpleDieselEngine` | `Propulsion.SimpleDieselEngine` | tables inlined as `ifelse` |
| `Machines.Crane`, `SubComponents.Cable` | `Ship6DOF.Crane`, `Ship6DOF.Cable` | 3D; cable break as a clocked latch (`CableBreakLatch`) |
| `Electrical.OnOffConsumer` | `Machinery.OnOffConsumer` | driven by a work signal instead of a random schedule |
| `DataProcessing.PeakSampler` | `Machinery.EventPeakSampler`, `Machinery.PeakSampler` | event-driven (`DiscreteComponents.ZeroCrossingClock`) and continuous peak-hold |
| `Others.Solar.*`, `Others.HeatTransfer.*` | `Thermal.*` | see docstrings |
| `Others.MoistAir.SourceMoistAir`, `DewTemperature` | `MoistAir.SourceMoistAir`, `MoistAir.DewTemperature` | on `HVACComponents` moist-air media and `MoistAirFluidPort` |

Not ported: `RainflowCounter` / `FatigueCounter` (needs `algorithm` + `pre`;
belongs in a Julia post-processor), `TriggerConsumer` / `StartGenerator`
(event-driven schedules), `SunIrradianceMultibody`, `EnvironmentHeatTransfer`, `ConvectionFactors.*`,
`SubComponents.Ikeda` (partial roll-damping stub upstream), `VariableTranslation`
(the 3D components apply forces at variable points algebraically instead).

## Corrections relative to upstream (documented in the docstrings)

- `HidrodynamicXYY.mx` (Brix surge added mass) divides by the water density
  inside the square root, which makes the fraction negative; the port uses
  `1 / (π sqrt(Lpp³/∇) - 14)`.
- `HidrodynamicZRP` applies the roll damping coefficient and roll rate to the
  pitch axis; the port uses the pitch coefficient and pitch rate.
- `ShipWind` normalises the lateral force with `A_T`; the port uses `A_L`
  (Fujiwara's definition). The lateral force and yaw moment keep the upstream
  minus signs.
- `Rudder`: the effective angle of attack is `α = δ - γ_R β_R`, so a drift
  angle reduces the effective rudder angle in a turn and a centred rudder
  produces the restoring fin force. Both the previous planar port (and, read
  literally, the upstream sign chain) had the inflow term with the opposite
  sign; on the planar stack that made a centred rudder a destabilising fin
  and left the hull in a steady turn after the rudder was returned to zero.
  The 6-DOF `RudderReturn` and planar `RudderReturn2D` analyses are the
  regression checks.
- Planar `HullMMG`: the added-mass cross terms are `+Δ my v r` (surge) and
  `-Δ mx u r` (sway) as fractions of displacement; the earlier port had both
  signs flipped and divided by 100.
- `Propeller4Q` data set `B4_70_14` is assigned twice upstream (the second
  block zeroes rows 8–31); the port keeps the non-zero rows.
- `Tank` drains the tank on the side the anti-heeling system moves ballast to;
  the port fills it, consistent with the righting moment `AntiHeeling` applies.
- `FourWingSails` upstream leaves the rudder fixed, so the ship weathercocks
  into the wind and the sails stall; the port steers with `WaypointAutopilot`
  and sheets the sails to 50°.

## Discrete events are available after all

`DiscreteComponents` (registry, not bundled with the kernel; 0.2.0 resolves on
dyad-3.3.0) brings a clocked sublanguage to Dyad: components and connectors
tagged `@[clk]`, `x@(clk-1)` for the previous sample, `PeriodicClock`,
`ZeroCrossingClock`, `Sampler`, `ZeroOrderHold`, delays, discrete PID and
binary (hysteresis) controllers, seeded noise. `EventPeakSampler`, the
`ZigZagRelay`, `AntiHeelingRelay`, `CableBreakLatch` and `WaypointIndex`
clocked components use it (each wrapped as Sampler → clocked logic →
ZeroOrderHold on a `PeriodicClock`), and the zig-zag, sailing, crane and
waypoint analyses validate them. The random schedules of `TriggerConsumer` /
`StartGenerator` can be written the same way with the seeded noise sources.
`RainflowCounter` needs loops over a stack and remains a Julia post-processor.
Two caveats: a periodic clock forces a solver step per tick (a 0.1 s clock
over 600 s is 6000 events, about ten times the step count of the continuous
relay), and the `ZeroOrderHold` continuous output reads back as its final
value through symbolic indexing of the solution, so post-process the clocked
signal itself (`sol[m.path.relay.s]`, one entry per tick) or the physical
response.
Installing it precompiles a Lustre toolchain (`Heptagon_jll`,
`Clang_unified_jll`): use `JULIA_NUM_PRECOMPILE_TASKS=1` and a 12 GB cap.

## Other registry libraries checked (September 2026)

- `DyadFMUGeneration` 0.2.1 installs on dyad-3.3.0 and provides
  `FMUAnalysis` (FMI 2/3, model-exchange and co-simulation, named inputs and
  outputs): the upstream `FMU` package can be reproduced by extending it from a
  ship component.
- `DyadControlSystems` 2.3.0 installs and provides linear analysis, PID
  autotuning, closed-loop and LQG analyses; a candidate for tuning
  `WaypointAutopilot` from a linearised `StandardShip`.
- `BlockComponents` 4.6 (dyad-3.4.0-rc1) adds no discrete or hysteresis
  blocks; nothing there we lack.
- `MultibodyComponents` mesh visualisers load STL/OBJ through `MeshIO`, not
  DXF; the upstream `Ship.dxf` hull would need converting.
- `BuildingsHeatTransfer` (convection factors), `PrimitiveComponents` and
  `TranslatedComponents` are pinned to older `DyadEcosystemDependencies` and
  do not resolve on the 3.3.0 sysimage.
- Two further registries exist on GitHub (`JuliaComputing/DyadHVACRegistry`,
  `JuliaComputing/DyadThermoFluidRegistry`) and are not served by the JuliaHub
  package server: clone them into `~/.julia/registries/` and set
  `JULIA_PKG_USE_CLI_GIT=true` so Pkg fetches the private package repos with
  the system git. Both libraries below target kernel 3.3.0 and install on the
  3.3.0 sysimage (verified 2026-09-07):
  - `HVACComponents` 0.3.0: moist-air, dry-air and refrigerant media
    (`BaseProperties_MoistAir` with p, T, Xi states and relative humidity),
    `MoistAirFluidPort`, `Boundary_pTPhi` / `MassFlowSource_TPhi` sources,
    mixing volumes, lumped rooms with 1D walls, heat exchangers, compressors.
    `MoistAir.SourceMoistAir` / `DewTemperature` are built on it (project
    dependency); `EnvironmentHeatTransfer`'s air side can follow.
  - `IncompressibleFlowComponents` 0.1.0: liquid media (constant or
    polynomial in T), `FluidPort` with `path medium`, pumps and fans with
    characteristics, valves, open tanks, volumes, ε-NTU and UA heat
    exchangers, pT and mass-flow boundaries. The anti-heeling
    circuit (`Ship6DOF.AntiHeelingCircuit`) is built on it (project
    dependency). Lessons: the library `Pump` is a flow source that needs a
    speed floor and then chatters at its zero-flow kink whenever the head
    across an idle branch changes sign, so the circuit uses an affinity-law
    `CentrifugalPump` (pressure-flow relation with a check valve and a small
    laminar term, non-singular at rest); `FluidPort`s must not pass through
    wrapper components (the stream expansion turns the port enthalpy into a
    runaway state), hence the circuit is a partial the ship test extends; an
    enum-variant field such as `LinearFlowCoefficient.K(value = …)` must be a
    literal, not a parameter. A clocked model (any `DiscreteComponents`
    clock) runs a callback at every step and every tick whose default
    `CheckInit` re-checks the algebraic residual without re-solving it, so
    a discrete output (the latched `direction`, the on/off state) must not
    enter the fluid loop's algebraic equations directly: at the tick where
    it jumps the pump flow equation is left with a residual of the size of
    the jump ("DAE initialization failed … normresid = 0.039", always at an
    integer second) and the run aborts. Pump speeds and valve openings
    therefore pass through `SlewRateLimiter` states (`T_pump`, `T_valve`).
    `ODEAlg.FBDF()` is no cure: the per-step callback resets the BDF history
    each step, so the step size never grows beyond its initial 6e-7 s.
  - `MediaComponents` / `FluidComponents` (general media and distributed
    pipes) exist but are earlier-stage (`kernel = 3.3.0-rc4`, water only).

## Known limitations

- With the default empirical derivatives the linear course-stability index of
  the sample hull is slightly negative (`HydrodynamicXYY.CourseStability`); the
  hull holds a straight course through the rudder's fin effect, which shows as
  a slow yaw-rate decay after rudder centring and 20/20 zig-zag overshoots
  around 30°. Override `N_r` (e.g. Clarke's estimate, about twice Khattab's)
  for a stiffer hull.
- No wave excitation, no shallow-water effects, no propeller transverse
  thrust in the crash stop.
- Clocked logic is sampled (0.05–1 s periods), not event-exact; use a
  `ZeroCrossingClock` where the crossing instant matters (as in
  `EventPeakSampler`), at the cost of the `dtmax` caveat in DiscreteComponents.
- Pods, sails and the crane have no aerodynamic/hydrodynamic interaction
  between units (no sail-sail interference, no pod-hull interaction).
- Small-angle hydrostatics (about ±20° heel/trim), as upstream.
