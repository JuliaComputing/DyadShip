# DyadShip

This is intended to be a rewrite of the modelica library https://github.com/BasilioPV/ShipSIM into dyad.

## Learning about dyad / julia

See the ./agent_resources folder for information. Follow what is there very carefully, unless there are guidelines in this document about how to run dyad/julia (for these, see below).

This will be an autonomous session so you don't need to make plots, but do validate your models.

## Heavy commands

Every compile, simulation, test or package operation goes through
`~/dyad-fleet/heavy -m 8G <command>` (memory-capped cgroup; see the global
`CLAUDE.md`). Never run the whole generated test suite; run single analyses
from a script instead.

## Running Julia

**Always** use `../julia-dyad.sh` (i.e. `/home/pgeorgakopoulos/dyad-ship/julia-dyad.sh`) to invoke Julia in this repo. Do **not** call `julia` directly.

The wrapper sets the environment this project requires:

- `JULIAUP_DEPOT_PATH=~/.julia/juliaup-depots/juliahub.com`
- `JULIAUP_SERVER=https://juliahub.com/juliabin`
- `JULIAUP_CHANNEL=dyad-3.3.0`
- `JULIA_PKG_SERVER=juliahub.com`

It launches the `dyad-3.3.0` channel and uses `--project=@.` so the active directory's project is used. Any extra args are forwarded.

Examples:

```sh
../julia-dyad.sh                              # REPL in current project
../julia-dyad.sh -e 'using Pkg; Pkg.status()' # one-off command
../julia-dyad.sh script.jl arg1 arg2          # run a script
```

## Running Dyad

Use `../dyad.sh` (i.e. `/home/pgeorgakopoulos/dyad-ship/dyad.sh`) to invoke the Dyad CLI. It runs `npx --yes @juliacomputing/dyad-cli@3.3.0` (GitHub Packages registry, token in `~/.npmrc`) and forwards all arguments verbatim.

Examples:

```sh
../dyad.sh --help                  # show CLI help
../dyad.sh compile                 # compile ./dyad in the current package
../dyad.sh render <component>      # render a model
```

## Toolchain pins (dyad-3.3.0)

The `dyad-3.3.0` channel's sysimage bakes in BlockComponents 4.5.1,
RotationalComponents 2.5.4, TranslationalComponents 2.5.0, ThermalComponents
2.0.5, ElectricalComponents 2.2.1 and DyadInterface 7.2.1; `Project.toml`
must pin exactly those (Pkg reports "package in sysimage!" otherwise).
`MultibodyComponents` is not in the sysimage and resolves freely (0.2.4 at
the time of writing, with the 3D library). Bundled library sources live under
`<juliaup>/julia-1.12.7+dyad-3x3x0…/share/julia/stdlib/v1.12/<Lib>/dyad/`.
`RotationalComponents.Sources.SpeedSource` is deprecated in favour of
`VelocitySource` (same ports).

## TASK (DONE)

For each component defined in ShipSIM, make one component in the ./dyad folder. 

1. Start with an investigation step; Make a list of what exists and plan to move them over.
2. Investigate which libraries you're going to use. BlockComponents, RotationalComponents, MultibodyComponents, TranslationalComponents are all available.
3. loop through the plan, translating one step at the time, starting from the simpler models moving up. Use composition to make complex models out of simple ones.
    1. Each model should be in a file
    2. If a model is too hard, skip it, but place it in HARD.md. Prefer to not skip components, try at least three times.
    3. For each component, make an analysis (in its own file) for one or more simulations that show how this works.
    4. Keep every success under a SUCCESSES.md, along with assumptions that were potentially used
    5. A component must compile cleanly with dyad
    6. data should be put in /assets - unless there is a better practice in agent_resources
4. Loop until can no longer make progress for 3 turns
5. When you stop make a thorough report on how you did

## TASK 2

1. Work similarly

## Learnings from a previous session

A prior autonomous session ported 18 ShipSIM components into Dyad and validated 5 transient analyses end-to-end. Read `REPORT.md`, `SUCCESSES.md`, and `HARD.md` for the full state. The notes below are observations that apply to *any* ShipSIM → Dyad work, not just to task-specific decisions.

### What went well

- **Tier the work simple → complex.** Doing parameter-only / signal-only components first (Environment, DewTemperature, IrradiationOnPlane) built up confidence with the toolchain before tackling propulsion and hull dynamics.
- **Validate end-to-end with composed analyses, not just per-component compile.** Compile-passes-but-runs-broken happened more than once. The closed-loop `PropellerOnHullTransient` exposed sign and inertia-allocation bugs that single-component compiles didn't.
- **Keep simplifications explicit in the docstring of each ported component.** The next agent can see at a glance what was dropped vs upstream (e.g. "multibody → signal", "Wageningen poly → quadratic-in-J fit", "OnOffController hysteresis → smooth saturation").
- **Run `dyad compile` after every file.** It's fast and catches syntax issues immediately. Filter the noisy JuliaHub banner with `grep -v -E "(Important Note|Commercial use|JuliaHub|public GitHub|To report|Error output \(0\))"`. If `generated/` updates after a compile call, it succeeded — even if stderr is loud.

### What didn't / pitfalls to avoid

- **`MultibodyComponents` ≥ 0.2 ships both the 2D `PlanarMechanics` sub-library and a full 3D library** (`Body`, `BodyShape`, `FixedTranslation`, `FixedRotation`, `Revolute`/`Prismatic`/`Spherical`…, `WorldForce`/`WorldTorque`/`Force`/`Torque`, `AbsoluteVelocity`/`AbsoluteAngularVelocity`/`TransformAbsoluteVector` sensors, `Cable`, visualizers). Sources are under `~/.julia/packages/MultibodyComponents/<hash>/dyad/`; `examples/floating_wind_turbine.dyad` (`BuoyantBody`) is the template for a free-floating 6-DOF body driven by algebraic forces. Namespacing: `MultibodyComponents.Body(...)`, enums as `MultibodyComponents.OrientationState.Euler()` / `MultibodyComponents.ResolveInFrame.FrameB()`, helper functions as `MultibodyComponents.resolve1/resolve2/cross/angular_velocity2(...)`. The `Ship6DOF` submodule is built on it; see the "6-DOF pass" learnings below.
- **`Modelica.Fluid` has no Dyad equivalent.** `HydraulicComponents` covers liquid-only hydraulics, not moist-air media. Components built on `Modelica.Fluid.Sources.MassFlowSource_T` with `Modelica.Media.Air.MoistAir` (e.g. `SourceMoistAir`) cannot be ported without first porting the medium model. The downstream signal-only sensors (e.g. `DewTemperature`) are fine.
- **No `when` / `discrete` events in Dyad as documented.** Anything using Modelica's `when ZeroCrossing(...) then`, `discrete Integer i + pre(i)`, `OnOffController` with hysteresis, or `TriggeredTrapezoid` cannot be ported 1:1 today. Workarounds:
  - For threshold-crossing-with-hysteresis: replace with a smooth `clamp((|x| - threshold) / band, 0, 1)` saturating activation. Loses crisp on/off but works in equation form.
  - For peak/zero-crossing samplers: skip and put in `HARD.md`. A `PeriodicSample`-style block (if/when one lands) might cover a subset.
- **Algebraic `if/elseif/else` chains in Modelica equations must become nested `ifelse(...)`.** Dyad has no algebraic `if/elseif`, only the `ifelse(cond, a, b)` function. For a 4-branch piecewise (e.g. SunScreen, max-torque-vs-RPM, SFOC) the nesting can get deep — write the data points inline as parameters and order branches by ascending threshold to keep it readable.
- **`Modelica.Blocks.Tables.CombiTable*` doesn't have a default-data shortcut in Dyad.** `BlockComponents.Tables.Interpolation` requires a `DyadData.DyadTimeseries` or 2D table from a CSV. For small fixed-data tables (< ~10 points) it's much simpler to inline the table as nested `ifelse` piecewise-linear, rather than ship a CSV. Reserve CSV-driven tables for ≥ 50-row datasets.
- **`dyad.sh` now exports the four `JULIAUP_*` / `JULIA_PKG_SERVER` env vars internally** so `../dyad.sh compile` runs cleanly. (Earlier the wrapper only `exec`'d node, and the dyad CLI's spawned `julia` failed with `ERROR: Invalid Juliaup channel \`dyad-3.0.0-rc5\`` unless the caller pre-exported them.) If you ever need to debug or override, the env block lives in `/home/pgeorgakopoulos/dyad-ship/dyad.sh` next to the `exec npx --yes @juliacomputing/dyad-cli@3.3.0 "$@"` line — keep it in sync with `julia-dyad.sh`. Compile output is still noisy on stderr (the JuliaHub banner); `generated/*.jl` mtime is the authoritative success signal.
- **MTK index reduction can choke on `der(...)` inside an algebraic floor like `sqrt(der(x)^2 + ε^2)`.** The first `Propeller1Q` had `(I+Ia)·der(w_eff) = flange.tau - Torque_Kq` with `w_eff = sqrt(der(flange.phi)^2 + w_floor^2)`. MTK saw two derivative orders of `phi` and reported `ExtraVariablesSystemException: 3 highest order derivative variables and 2 equations`. Fix: introduce an explicit state `w` with `w = der(flange.phi)` and apply the floor only inside non-derivative algebraic uses; never put a `der(...)` *inside* a `sqrt`.
- **Connector flow rules surprise: a single connection with a forced flow makes the system over-determined.** Setting `prop.flange.tau = 50000` directly with no other connection on `flange` produced `ExtraEquationsSystemException` because the connector also implies "sum of flows = 0 ⇒ flange.tau = 0" at an unconnected port. Use a proper source: `RotationalComponents.Sources.TorqueSource` + a `Constant` on its input + a `Fixed` grounding the support spline.
- **`RotationalComponents.Sources.TorqueSource` needs its support spline grounded.** It extends `PartialElementaryOneSplineAndSupport`, which has `support` and `phi_support`. Always connect the support to a `RotationalComponents.Components.Fixed` or another grounded element.
- **External libraries must be added to `Project.toml` AND `using`-imported in `src/<ModuleName>.jl`.** The Dyad compiler's auto-imports cover ModelingToolkit, Markdown, Moshi, OrdinaryDiffEqDefault, RuntimeGeneratedFunctions — that's it. For BlockComponents/Rotational/Translational/Thermal you need `Pkg.add(...)` *and* a manual `using ...` line in `src/DyadShip.jl`.
- **Don't lump shaft inertia inside a propulsion component.** Originally I had `Propeller1Q` include its own propeller + added-water inertia internally with `der(w)` on a one-sided `Spline`. MTK didn't index-reduce cleanly. Keeping the propeller purely *algebraic* on `flange.tau = Torque_Kq` and forcing the user to attach an external `RotationalComponents.Components.Inertia(J = ...)` worked first try. Same shape applies to the diesel engine: end the component in a `Spline` and let the caller put the inertia in series.
- **`PartialCompliant` (translational) and `PartialTwoSplines` (rotational) are the right base classes for spring-like components.** Don't reinvent the flange wiring; `extends TranslationalComponents.Interfaces.PartialCompliant` saves writing the `s_rel`, `f`, and force-distribution equations.
- **Array-typed connectors are awkward in Dyad.** `Modelica.Blocks.Interfaces.RealInput[3]` for a vector input becomes either three scalar `RealInput()`s or an internal `variable v::Real[3]` plus three scalar inputs that copy into it. There's no clean array-RealInput pattern at this Dyad version.
- **`agent_resources/stdlib_reference/` is documentation, not a library.** Quoted from `libraries.md` but worth re-emphasizing — don't try to `using ./agent_resources/...`. The reference `.dyad` files are excerpts to read, never to copy.
- **Regenerated analysis specs now default `optimize` to `OptimizationLevel.Aggressive()`** (was `OptimizationLevel.None()`). This behavior change arrived with the kernel 3.2.0-rc2 codegen; validated on ShipTransient / Plate / Buoyancy / SimpleAutoPilot. Scripts comparing against stored 3.1.1 outputs may see small numerical deltas.
- **`unit = "..."` annotations on Modelica `Real` variables don't all map to a Dyad type alias.** Many ShipSIM variables are `Real` with a unit annotation rather than a typed quantity (e.g. `Real Q(unit = "m^3/h")`, `Real rpm(unit = "rpm")`). Don't blindly translate to a typed Dyad alias — `m³/h` and `rpm` aren't in the standard `types.md` list. Keep them as plain `Real` (or define a project-local `type` alias) so the unit checker doesn't reject them.

### Submodule organization (third-pass learnings)

The library is organized into 4 submodules under `dyad/`, plus a few foundational
helpers at the root:

- **Top level**: `Environment.dyad`, `VariableEnvironment.dyad`, `ApparentSpeedXY.dyad`
  — small, foundational, no domain-specific dependencies.
- **`Ship6DOF/`**: the primary, 3D-multibody stack: `ShipBody`, `HydrodynamicXYY`,
  `HydrodynamicZRP`, `ApparentSpeedXY`, `Propeller1Q`, `Propeller4Q`, `Rudder`,
  `ShipWind`, `WingSail`, `POD4Q`, `AntiHeeling`, `Crane`, `Cable`,
  `WaypointAutopilot`, the `StandardShip` partial assembly and the validation
  analyses (speed trial, turning circle, rudder return, zig-zag, roll decay,
  crash stop, sail sweep, four wing sails, pod turning circle, crane operation,
  autopilot transit, manual ship). Julia helpers (Wageningen
  polynomials, 4Q Fourier sets, draft polynomials) in `Ship6DOF/definitions.jl`.
- **`MoistAir/`**: `SourceMoistAir`, `DewTemperature` on `HVACComponents` (`MoistAirFluidPort`, `medium_data::HVACComponents.MediaModels.AbstractMoistAirMedium = HVACComponents.moistair`, property functions `HVACComponents.h_pTX / T_phX / d_pTX / massFraction_pTϕ / saturationPressure / relativeHumidity_pTX`, `instream(port.h_outflow)` for stream variables at a sensor).
- **`Ship/`** (planar): `HullMMG`, `ShipWind`, `Tank`, `HeadingAutoPilot`,
  `ManualShip`, the `FullShip*` assemblies (autopilot, disturbed, render, Flettner).
- **`Propulsion/`** (planar): `Propeller1Q`, `Propeller4Q`, `Rudder`,
  `FlettnerRotor`, `FlettnerRotorOnline`, `SimpleDieselEngine`.
- **`Machinery/`**: `OnOffConsumer`, `PeakSampler` (the crane and cable moved to `Ship6DOF`).
- **`Thermal/`**: heat transfer + solar + moist air — `PlateTransient`,
  `CylinderTransient`, `SimpleAirExchanger`, `ConvRadSunWall`, `TemperatureDataset`,
  `DewTemperature`, `SolarIrradiation`, `IrradiationOnPlane`, `SunScreen`. Also
  carries a `definitions.jl` with the Julia helpers (`sun_height_deg`,
  `sun_vector_world`, `discretize_cylinder_*`).

**Why some files live at the root:**
- **Submodule name vs. component name collision crashes precompile.** A `dyad/Foo/` directory
  containing `Foo.dyad` (component `Foo`) makes Julia generate `module Foo ... function
  Foo(...) ... end end`, which raises `cannot define function Foo; it already has a
  value`. The `Environment.dyad` component therefore lives at the root rather than in
  an `Environment/` submodule. Either rename the submodule or rename the component to
  resolve this.
- **`dyad/<Sub>/definitions.jl` vs. `dyad/definitions.jl`.** Julia helpers in
  `dyad/<Sub>/definitions.jl` are visible only inside `DyadShip.<Sub>`. Helpers used
  by components in the same submodule should live there. Helpers used at the top
  level go in `dyad/definitions.jl`. The `discretize_cylinder_*` and `sun_*` helpers
  moved into `dyad/Thermal/definitions.jl` because their consumers (`CylinderTransient`,
  `SolarIrradiation`) are both in `Thermal/`.
- **Cross-submodule references must be fully qualified from the root** (per
  `library_namespacing.md`). E.g. `Ship/Ship_analysis.dyad` references the propeller
  as `DyadShip.Propulsion.Propeller1Q(...)`, not `Propeller1Q(...)`. Within the same
  submodule, the bare name works.
- **Julia callers must use `DyadShip.Ship6DOF.TurningCircleTransient()` / `DyadShip.Ship.ShipTransient()`**, not `ShipTransient()`,
  because the analyses are exported only inside their submodule (not re-exported at
  the root).

**Tip:** when any rename / move breaks compile, delete the matching `generated/<Sub>/`
folder before recompiling — stale auto-generated files (especially `Sub_definition.jl`)
can collide with new code.

### Second-pass learnings (multibody + array discretization + Julia helpers)

- **Component-array comprehensions stumble on outer parameters being substituted into element constructors.** In `[HeatCapacitor(C = C_node, T0 = T_init) for i in 1:N]`, MTK errors with `Could not evaluate value of parameter caps⸺i₊T0. Missing values for variables in expression caps⸺i₊T_init` — the outer `T_init` symbol gets re-scoped to `caps[i].T_init` (which doesn't exist) during substitution. **Fix:** declare every outer parameter that flows into the comprehension as a `structural parameter`, and likewise for derived params (`final parameter dx = ...` becomes `structural parameter dx = ...`). Then they're literals at compile time and the substitution loop terminates. Pattern landed `PlateTransient` and `CylinderTransient` cleanly.
- **Component-array comprehensions can't size on derived expressions like `nNodes - 1`.** Compile error: `unimplemented - Unable to size arrays based on complex range expressions`. Workaround: declare a separate `structural parameter nInterior::Integer` for the size and let the user keep it consistent with `nNodes` when overriding.
- **Component arrays accessed from Julia use the `caps⸺i` (with the special `⸺` separator) symbol, not normal indexing.** `model.plate.caps[1]` fails; use `getproperty(model.plate, Symbol("caps⸺", 1))`. The `⸺` is U+2E3A (TWO-EM DASH) — copy from a printed `propertynames(...)` result.
- **`HeatCapacitor.T0` is a *required* parameter — no default.** Forgot once and got "missing value" on every node. Must always pass an initial temperature, even if you also write `initial T = ...` in `relations` (and in fact, doing both raises an extra-equations exception).
- **Julia helpers that get called from Dyad parameters or relations must avoid Julia-only `if/else` ternaries on `Symbolics.Num` inputs.** `sin_alt > 0 ? a : b` raises `non-boolean (Symbolics.Num) used in boolean context` once the helper sees a symbolic time argument. Rewrite as plain arithmetic and gate on the Dyad side with `ifelse(...)`. Same for `max(-1, min(1, x))` clamps — replace with `sqrt(x^2 + ε)` or similar smooth analogues.
- **Type the helper signatures as `Real`, not concrete types.** `function f(t::Real, ...)` (or even untyped) plays nicer with MTK Symbolic propagation than `function f(t::Float64, ...)`.
- **`HullMMG` pattern for planar hulls:** wrap a `PlanarMechanics.Body` and a `WorldForceTorque(resolve_in_frame = FrameB)` together. The `WorldForceTorque` carries body-frame drag + `Fx_extra`/`Fy_extra`/`Mz_extra` real inputs. External components (propeller, rudder, wind) connect their `frame_a` to the wrapper's `frame_a`, and forces sum at the body node automatically. Saves writing connector-level flow equations by hand.
- **`PlanarMechanics.World` is required even when "you don't need gravity".** Set `g = 0` for ship/horizontal models; without any `World`, mtkcompile errors with `Could not evaluate value of parameter g`.
- **`Revolute.phi` is a state, not a constraint variable.** Forcing `boom_pivot.phi = command` causes an extra-equations error. Drive it via a `RotationalComponents.Sources.Position` connected to `Revolute.flange_a`, with the position source's `support` grounded by `RotationalComponents.Components.Fixed`.
- **Dyad `if … end` blocks inside `relations` work for *enum-keyed* equation inclusion** (see `world_force_torque.dyad`), but the cases must use literal enum values like `MultibodyComponents.ResolveInFrame.World()` — not arbitrary boolean expressions. For value-level branching, use algebraic `ifelse(cond, a, b)`.
- **`dyad/definitions.jl` is auto-included by the generated module** (the generated `definitions.jl` does `if isfile(joinpath((@__DIR__) |> Base.dirname, "dyad", "definitions.jl"))` then include it). Drop Julia helpers there. They get the project's full `using ...` scope from `src/<ModuleName>.jl`.
- **Saving CSV assets to `assets/` and referencing them as `dyad://<PackageName>/<filename>.csv` works without any registration.** `BlockComponents.Tables.Interpolation` reads them through `DyadData.DyadTimeseries`. Keep the asset small (< 100 KB) for compile-time embedding.
- **Tension-only `ifelse(...)` cables can create algebraic loops with passive bodies** in PlanarMechanics. The 2D `Crane` couldn't close with the tension-only `Cable` + a free `Body` (cf. `STILL_HARD.md`). Use a `PlanarMechanics.Spring` or `SpringDamper` to break the loop at the cost of allowing compression.

### 6-DOF pass learnings (MultibodyComponents 3D, dyad-3.3.0)

- **Free-floating body:** `MultibodyComponents.Body(orientation_state = MultibodyComponents.OrientationState.Euler(), sequence = [3, 2, 1], statePriority = 100, linearStatePriority = 100, r_0(initial = …), v_0(initial = …), phi(initial = [yaw, pitch, roll]))` with nothing connecting it to the world. `sequence = [3, 2, 1]` is yaw-pitch-roll (`phi[1]` = yaw, unbounded). `frame_a.f` / `frame_a.tau` are resolved in the body frame; `body.v_0` is the world velocity, `body.w_a` the body-frame angular velocity.
- **Apply algebraic forces with `WorldForce(resolve_in_frame = ResolveInFrame.FrameB())` + `WorldTorque(FrameB)` connected to the same frame,** setting `force_x/y/z`, `torque_x/y/z` from equations (the `BuoyantBody` example pattern). A force acting away from the frame origin is applied at the origin plus `MultibodyComponents.cross(r, F)` as torque; no `VariableTranslation` needed.
- **Reading kinematics inside a force component works without sensors:** `MultibodyComponents.resolve2(frame_a.R, der(frame_a.r_0))` and `MultibodyComponents.angular_velocity2(frame_a.R)`. MTK also accepts `der()` of those algebraic velocities for added-mass terms (`-Δ mx der(u)`); the resulting acceleration loop is handled by structural simplification (turning circle: ~800 steps for 900 s).
- **Array parameters with structural row counts** (`structural parameter n::Integer = 7; parameter T::Real[n, 2] = [[…], …]`) can be passed to a Julia helper (`draft_poly(T, x)`); `structural parameter PropModel::String` selects a Julia `Dict` entry (`wageningen_4q_ct(PropModel, beta)`).
- **`BlockComponents.Tables.InterpolatedTable` with `DyadData.DyadInterpolationTable2D`** wants the Modelica `CombiTable2D` CSV layout (header row = axis-2 values, first column = axis-1 values) **with every cell parseable as Float64** — an integer-valued axis column makes the ND interpolator fail with `no method matching validate_size_u` (mixed `Vector{Int64}` / `Vector{Float64}` axes). Clamp the inputs to the tabulated range yourself; there is no extrapolation option.
- **Driving a 3D `Revolute` with `RotationalComponents.Sources.Position`:** connect `drive.spline → joint.axis` and `drive.support → joint.support`, and use `ref_type = ReferenceType.Exact()` with a differentiable reference (e.g. a `SlewRateLimiter` output). `Filtered(f_crit = …)` on a revolute is structurally unbalanced (`ExtraVariablesSystemException … drive₊wˍt`). This is the servo pattern of `WingSail`, `POD4Q`, `Rudder` and `Crane`.
- **Tension-only cable:** `extends MultibodyComponents.PartialLineForce` and set the scalar `f` (positive = tension); instantiate with `fixed_rotation_at_frame_b = true` so a `Body` with no orientation state can hang from it (see `Ship6DOF.Cable`, `CraneOperation`).
- **Two same-named conditional subcomponents are not supported** (`prop = Propeller1Q() if !flag` / `prop = Propeller4Q() if flag`, the pattern MultibodyComponents' `Spring` uses internally): the compiler emits only the last declaration, and the other branch fails at build time with `prop not defined`. Write a separate assembly (`CrashStop` vs `StandardShip`) instead.
- **Private registries (HVAC, thermo-fluid):** the JuliaHub package server only serves `DyadRegistry`. For `HVACComponents` / `IncompressibleFlowComponents`, `git clone https://github.com/JuliaComputing/DyadHVACRegistry.git ~/.julia/registries/DyadHVACRegistry` (same for `DyadThermoFluidRegistry`), then run Pkg with `JULIA_PKG_USE_CLI_GIT=true` (Pkg's LibGit2 fails with "Unable to parse URL" on the private package repos; `Pkg.Registry.add(url = …)` fails the same way). Precompile serially under `heavy -m 12G`.
- **Discrete events exist: `DiscreteComponents`** (registry package, `Pkg.add`, needs `using DiscreteComponents` in `src/DyadShip.jl`; install with `JULIA_NUM_PRECOMPILE_TASKS=1` under `heavy -m 12G`, parallel precompilation of its Lustre toolchain gets OOM-killed at 8G). Clocked sublanguage: `component X@[input clk extends Discrete]`, `RealInput@[clk]()`, `y@clk = y@(clk-1) + u@clk`, `initial y@(clk-1) = 0`; `PeriodicClock(dt)`, `ZeroCrossingClock` (event clock on a signal's zero crossings), `Sampler`, `ZeroOrderHold`, `DiscreteBinaryController` (hysteresis relay). Read clocked signals with `sol[m.path.to.clocked_output]` (one entry per tick); the `ZeroOrderHold` continuous output reads back as a constant through symbolic indexing. See `Machinery.EventPeakSampler`.
- **IncompressibleFlowComponents circuits:** declare `path medium::IncompressibleFlowComponents.Media.AbstractMedium = IncompressibleFlowComponents.Media.ConstantIncompressible(density = 1025.0)` in the component that owns the circuit and `continuity(medium, some.port.medium)`. Keep the circuit flat: a `FluidPort` passed through a wrapper component's own port makes MTK's stream expansion turn `h_outflow` into a runaway state (`Unstable` within a step). Put the circuit in a `partial component` that the test `extends` (multiple `extends` are fine). The library `Pump` is a flow source with a speed floor and chatters at its `max(0, …)` zero-flow kink when the head across an idle branch changes sign; `Ship6DOF.CentrifugalPump` (affinity law `H = s² H0 - k Q|Q| - k_lin Q`, check valve via `k_back`, laminar term so the Jacobian is non-singular at rest) is robust. Enum-variant fields (`LinearFlowCoefficient.K(value = …)`) only accept literals; size valves with `m_flow_nominal` / `dp_nominal`. Never reverse a running transfer instantly: latch the direction in the clocked relay and ramp valve openings with `SlewRateLimiter`s.
- **Clocked models: discrete outputs must not feed algebraic equations directly, and FBDF is unusable.** The DiscreteComponents callbacks (every step and every tick) use `CheckInit`, which re-checks the algebraic residual without re-solving it, so a clocked signal that jumps into an algebraic loop (pump speed from a latched direction) aborts with "DAE initialization failed … normresid > abstol" at an integer-second tick even though initialization succeeded; route such signals through a `SlewRateLimiter`/`FirstOrder` state. `alg = ODEAlg.FBDF()` stalls at the initial step size on any clocked model because the per-step callback resets the BDF history. Diagnose by building the problem with `DyadInterface.get_simplified_model`/`setup_prob`, `init(prob, Rodas5P())` and stepping manually until it throws, then printing the algebraic residuals.
- **Prismatic-slide devices (tuned mass damper):** `MultibodyComponents.Prismatic(n = [0, 1, 0], s(initial = 0), v(initial = 0))` on a `FixedTranslation`, a `Body` on `frame_b`, and `TranslationalComponents.Components.SpringDamper` between `slide.support` and `slide.axis` (see `RollTunedMassDamper`).
- **Time-varying masses on the hull:** `Ship6DOF.VariableMass` is `Body` without an inertia tensor and with `m` and the mass-centre height as inputs (`f = m (resolve2(R, a_0 - g_0) + z×r_cm + w×(w×r_cm))`, `tau = r_cm × f`); hang it from a `FixedTranslation` on `ShipBody.frame_a`. Used by `BallastTank`; the same pattern serves fuel, cargo or crane loads that change mass.
- **Hysteresis / latches / counters on clocks:** write a small clocked component (`component X@[input clk extends Discrete]`, `RealInput@[clk]()`, `x@clk = ifelse(x@(clk-1) > 0.5, …)`, `initial x@(clk-1) = …`) and wrap it in a continuous component as `Sampler → X → ZeroOrderHold` with a `PeriodicClock(dt)` connected to the sampler output (`connect(sampler.y, x.u, clock.y)`). Inside the clocked component only sampled inputs and previous values are available (no `time`); gate by time outside. Examples: `ZigZagRelay`, `AntiHeelingRelay`, `CableBreakLatch`, `WaypointIndex`. Cost: one solver step per tick, so keep `dt` as large as the physics allows. The old continuous relay trick (`der(s) = (target - s)/T` with `target` an `ifelse` on `s`) is no longer needed.
- **Autopilots:** an `ifelse`-frozen integrator with a `clamp` on the output chatters at the saturation boundary (320k solver steps for 2400 s); `BlockComponents.Continuous.LimPID` with back-calculation anti-windup (`u_s` = wrapped course error, `u_m` = 0) takes 5k steps.
- **Sign-check every hydrodynamic port with a physical experiment,** not only with a compile: the yaw-rate persistence bug was a rudder inflow-angle sign, caught by the rudder pulse-and-return analysis, and the planar wind loads had both lateral signs flipped relative to upstream.

### Frame conventions (caught during code review)

- **`PlanarMechanics.Body.frame_a.fx`/`fy` are world-frame forces** (`m·der(v) = f`
  where `v = der(r)` is the world-frame velocity). When a force-producing component
  computes a body-frame `Force_X`/`Force_Y`, it MUST rotate by `frame_a.phi` before
  writing to `frame_a.fx`/`fy` — otherwise the force is interpreted as world-frame and
  the result is silently wrong (compiles, runs, drifts in the wrong direction).
  The correct pattern (see `Rudder`, `ShipWind`, `WingSail`, `Propeller1Q`):
  ```dyad
  frame_a.fx + cos(frame_a.phi) * Force_X_body - sin(frame_a.phi) * Force_Y_body = 0
  frame_a.fy + sin(frame_a.phi) * Force_X_body + cos(frame_a.phi) * Force_Y_body = 0
  frame_a.tau + Moment = 0   # tau is rotation-invariant in 2D, no rotation needed
  ```
  Same convention as `MultibodyComponents.PlanarMechanics.WorldForceTorque(resolve_in_frame
  = FrameB())` — read that component for the canonical implementation.
- **`World` is only required if a `Body` (or any component referencing
  `gravity_acceleration_2d()`) is present in the model.** Tests that only use `Fixed`
  + `frame_a` (no `Body`) should NOT include a `World`, otherwise MTK errors with
  `Expected an Initial parameter to exist for variable g`. Counter-intuitive — the
  natural assumption is "always include a World". For sub-component analyses with no
  body, drop it.
