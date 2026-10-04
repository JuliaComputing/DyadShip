#!/usr/bin/env julia
# Focused checks of the synthetic SimpleDieselEngine parameter presets: the files and
# their sidecars, that each wrapper is the one engine with the file's values (including a
# table length other than the default), that both presets solve, and that the metered fuel
# is what the file's table gives at every accepted step.
#
# Runs three engine analyses only; never the generated test suite. The files that must be
# refused are checked by test/invalid_engine_presets/check.jl.
#
# Usage:  julia --project=<environment with DyadShip> scripts/validate_engine_presets.jl

using DyadShip, Printf, Test, TOML
using DyadInterface: symbolic_container
using ModelingToolkit: assertions, equations, get_systems, nameof

const PR = DyadShip.Propulsion
const PRESETS = joinpath(pkgdir(DyadShip), "assets", "presets", "Synthetic")

successful_retcode(sol) = string(sol.retcode) == "Success"
trapz(ts, ys) = sum((ts[i + 1] - ts[i]) * (ys[i] + ys[i + 1]) / 2 for i in 1:length(ts) - 1)

# Independent piecewise-linear lookup (no held ends needed: callers stay inside the table).
function lookup(xs, ys, x)
    i = clamp(searchsortedlast(xs, x), 1, length(xs) - 1)
    return ys[i] + (ys[i + 1] - ys[i]) * (x - xs[i]) / (xs[i + 1] - xs[i])
end

core_of(sys) = only(filter(s -> nameof(s) == :core, get_systems(sys)))

presets = [
    (label = "A", file = "diesel_engine_a.toml", wrapper = PR.SyntheticDieselEngineA, analysis = PR.SyntheticDieselEngineATransient, n = 6),
    (label = "B", file = "diesel_engine_b.toml", wrapper = PR.SyntheticDieselEngineB, analysis = PR.SyntheticDieselEngineBTransient, n = 3),
]

@testset "SimpleDieselEngine parameter presets" begin

# The canonical run both presets are compared with (same load, default table).
base = PR.SimpleDieselEngineInRangeTransient()
be = symbolic_container(base).engine
@testset "default component unchanged" begin
    @test successful_retcode(base.sol)
    @test base.sol[be.Fuel][end] ≈ 1.21982503 rtol = 1e-7
    @test base.sol.ps[be.SFOC_P] == [605, 907.5, 1028.5, 1210]
    @test base.sol.ps[be.SFOC_g] == [185, 179, 178, 182]
    @test base.sol.ps[be.m_dot_idle] == 0
end

@testset "table length is checked at construction" begin
    # A table whose length disagrees with n_sfoc is refused before a problem is built.
    @test_throws ArgumentError PR.SimpleDieselEngine(; name = :engine, SFOC_P = [300.0, 605.0, 907.5, 1028.5, 1210.0])
    @test_throws ArgumentError PR.SimpleDieselEngine(; name = :engine, SFOC_g = [185.0, 179.0, 178.0])
    @test_throws ArgumentError PR.SimpleDieselEngine(; name = :engine, n_sfoc = 5)
    @test_throws ArgumentError PR.checked_table_length(4, [1.0, 2.0, 3.0], [1.0, 2.0, 3.0, 4.0])
    @test PR.checked_table_length(3, [1.0, 2.0, 3.0], [1.0, 2.0, 3.0]) == 3
    # Consistent lengths are accepted, default and non-default.
    @test PR.SimpleDieselEngine(; name = :engine) !== nothing
    @test PR.SimpleDieselEngine(; name = :engine, n_sfoc = 5, SFOC_P = [300.0, 605.0, 907.5, 1028.5, 1210.0],
        SFOC_g = [185.0, 185.0, 179.0, 178.0, 182.0]) !== nothing
end

fuel = Dict{String, Float64}()

for p in presets
@testset "preset $(p.label)" begin
    asset = TOML.parsefile(joinpath(PRESETS, p.file))
    side = TOML.parsefile(joinpath(PRESETS, replace(p.file, ".toml" => ".provenance.toml")))

    @testset "file and sidecar" begin
        # Exactly the table and the idle rate: nothing else of the engine is overridden.
        @test sort(collect(keys(asset))) == ["SFOC_P", "SFOC_g", "m_dot_idle", "n_sfoc"]
        @test asset["n_sfoc"] == p.n != 4
        @test length(asset["SFOC_P"]) == p.n && length(asset["SFOC_g"]) == p.n
        @test PR.table_knots_valid(Float64.(asset["SFOC_P"])) == 1
        @test PR.table_values_valid(Float64.(asset["SFOC_g"])) == 1
        @test PR.rate_valid(asset["m_dot_idle"]) == 1
        @test side["asset"] == p.file
        @test side["calibration_status"] == "synthetic"
        @test side["manufacturer"] == "none" && side["model"] == "none"
        @test occursin("SYNTHETIC", read(joinpath(PRESETS, p.file), String))
    end

    @testset "wrapper is the one engine" begin
        w = p.wrapper(; name = :engine)
        core = core_of(w)
        ref = PR.SimpleDieselEngine(; name = :core, n_sfoc = p.n,
            SFOC_P = Float64.(asset["SFOC_P"]), SFOC_g = Float64.(asset["SFOC_g"]), m_dot_idle = asset["m_dot_idle"])
        # Same equations and assertions as the component built directly; the wrapper
        # itself contributes connections only.
        @test string.(equations(core)) == string.(equations(ref))
        @test length(assertions(core)) == 3
        @test sort(nameof.(get_systems(w))) == [:core, :flange]
    end

    res = p.analysis(); sol = res.sol; m = symbolic_container(res)
    e = m.engine; c = m.engine.core
    @testset "solve with the file's table" begin
        @test successful_retcode(sol)
        # The structural length and the values delivered by the file reached the solve.
        @test length(sol.ps[c.SFOC_P]) == p.n
        @test sol.ps[c.SFOC_P] == asset["SFOC_P"]
        @test sol.ps[c.SFOC_g] == asset["SFOC_g"]
        @test sol.ps[c.m_dot_idle] == asset["m_dot_idle"]
        # Only the fuel changes: speed and power are those of the default-table run.
        @test sol[c.rpm][end] ≈ 1800 rtol = 1e-6
        @test sol[c.ShaftPower][end] ≈ base.sol[be.ShaftPower][end] rtol = 1e-6
        @test sol[c.KWh][end] ≈ base.sol[be.KWh][end] rtol = 1e-5
        # Wrapper ports carry the core's signals.
        for name in (:KWh, :Fuel, :Fuel_extrapolated, :Inst_Fuel, :ShaftPower, :SFOC_valid)
            @test sol[getproperty(e, name)] == sol[getproperty(c, name)]
        end
    end

    @testset "metering at accepted steps" begin
        ts = sol.t; P = sol[e.ShaftPower] ./ 1000; mdot = sol[e.Inst_Fuel]; F = sol[e.Fuel]
        xs = Float64.(asset["SFOC_P"]); ys = Float64.(asset["SFOC_g"]); idle = asset["m_dot_idle"]
        # Inside this preset's table for the whole run.
        @test all(xs[1] .<= P .<= xs[end])
        @test all(sol[e.SFOC_valid] .== 1)
        @test all(sol[e.Fuel_extrapolated] .== 0)
        # Rate at every accepted step from an independent lookup of the file's table.
        expected = [lookup(xs, ys, P[i]) * P[i] / 3.6e6 + idle for i in eachindex(P)]
        @test maximum(abs.(mdot .- expected) ./ expected) <= 1e-12
        @test all(diff(F) .>= 0) && F[1] == 0
        # Counter against a tight-tolerance run and against the rate samples (trapezoid
        # gap within twice the trapezoid-versus-Simpson estimate plus the tolerance).
        ref = p.analysis(abstol = 1e-10, reltol = 1e-10)
        Fref = ref.sol[symbolic_container(ref).engine.Fuel][end]
        @test F[end] ≈ Fref rtol = 1e-6
        T = trapz(ts, mdot)
        S = sum((ts[i + 1] - ts[i]) * (mdot[i] + 4 * sol((ts[i] + ts[i + 1]) / 2; idxs = e.Inst_Fuel) + mdot[i + 1]) / 6
                for i in 1:length(ts) - 1)
        @test abs(T - F[end]) <= 2 * abs(T - S) + 1e-6 * F[end]
        # The idle share is exactly the idle rate times the elapsed time.
        @test F[end] - idle * ts[end] > 0
        @test idle == 0 || mdot[end] > lookup(xs, ys, P[end]) * P[end] / 3.6e6
        fuel[p.label] = F[end]
        @printf("  preset %s: %d knots, idle %.4f kg/s, %d accepted steps, power %.1f .. %.1f kW, Fuel(30) = %.8f kg (tight-tolerance run %.8f)\n",
            p.label, p.n, idle, length(ts), minimum(P), maximum(P), F[end], Fref)
    end
end
end

@testset "presets are distinct parameter sets" begin
    @test !isapprox(fuel["A"], fuel["B"]; rtol = 1e-3)
    @test !isapprox(fuel["A"], base.sol[be.Fuel][end]; rtol = 1e-3)
    @test !isapprox(fuel["B"], base.sol[be.Fuel][end]; rtol = 1e-3)
end

end
