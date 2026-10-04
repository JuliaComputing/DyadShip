#!/usr/bin/env julia
# Parameter files that must be refused. Three kinds:
#  - files the loader cannot apply (unknown key, wrong type, a provenance sidecar, a
#    missing file): the compiler reports each as an `urn:instantiate:apply-…` diagnostic
#    while still exiting 0, and the component cannot be constructed;
#  - arrays whose length disagrees with `n_sfoc`: the compiler does not report it (the
#    length is structural), and the engine refuses it with an `ArgumentError` when the
#    component is constructed, before any problem is built;
#  - files that load but cannot be metered (unordered, duplicate or non-finite knots,
#    non-positive or non-finite SFOC, negative idle rate): the component is constructed
#    and the engine's assertions end the run without a successful return code.
#
# Usage:
#   julia --project=<environment with DyadShip and this package> check.jl [compile-log]
# The optional argument is the log of `dyad compile` run on this package.

using DyadShipInvalidEnginePresets, DyadShip, Printf, Test, TOML

const FX = DyadShipInvalidEnginePresets
const PR = DyadShip.Propulsion
const FILES = joinpath(pkgdir(FX), "assets", "presets", "Invalid")

successful_retcode(sol) = string(sol.retcode) == "Success"

loader_cases = [
    ("unknown_key", FX.InvalidUnknownKey),
    ("wrong_type", FX.InvalidWrongType),
    ("sidecar.provenance", FX.InvalidSidecar),
    ("does_not_exist", FX.InvalidMissingFile),
]
length_cases = [
    ("length_mismatch", FX.InvalidLengthMismatch, FX.InvalidLengthMismatchTransient),
    ("longer_without_count", FX.InvalidLongerWithoutCount, FX.InvalidLongerWithoutCountTransient),
]
value_cases = [
    ("decreasing_knots", FX.InvalidDecreasingKnots, FX.InvalidDecreasingKnotsTransient, :knots),
    ("duplicate_knot", FX.InvalidDuplicateKnot, FX.InvalidDuplicateKnotTransient, :knots),
    ("inf_knot", FX.InvalidInfKnot, FX.InvalidInfKnotTransient, :knots),
    ("nan_sfoc", FX.InvalidNanSfoc, FX.InvalidNanSfocTransient, :values),
    ("zero_sfoc", FX.InvalidZeroSfoc, FX.InvalidZeroSfocTransient, :values),
    ("negative_idle", FX.InvalidNegativeIdle, FX.InvalidNegativeIdleTransient, :idle),
]

@testset "invalid engine parameter files" begin

@testset "refused by the loader" begin
    for (file, ctor) in loader_cases
        err = try
            ctor(; name = :engine); nothing
        catch e
            e
        end
        @printf("  %-22s construction: %s\n", file, err === nothing ? "ACCEPTED" : first(replace(sprint(showerror, err), '\n' => ' '), 110))
        @test err !== nothing
    end
end

run_outcome(analysis) = try
    successful_retcode(analysis().sol) ? :accepted : :failed
catch
    :threw
end

@testset "length that disagrees with n_sfoc" begin
    for (file, ctor, analysis) in length_cases
        asset = TOML.parsefile(joinpath(FILES, file * ".toml"))
        @test !(length(asset["SFOC_P"]) == get(asset, "n_sfoc", 4) == length(asset["SFOC_g"]))
        # Refused when the component is constructed, with the three lengths named.
        err = try
            ctor(; name = :engine); nothing
        catch e
            e
        end
        @printf("  %-22s construction: %s\n", file, err === nothing ? "ACCEPTED" : first(replace(sprint(showerror, err), '\n' => ' '), 110))
        @test err isa ArgumentError
        @test err isa ArgumentError && occursin("SFOC table length mismatch", err.msg)
        # And so the analysis never reaches a solve.
        @test run_outcome(analysis) == :threw
    end
end

@testset "refused by the engine's checks" begin
    for (file, ctor, analysis, kind) in value_cases
        asset = TOML.parsefile(joinpath(FILES, file * ".toml"))
        # The file is well-formed for the loader and the plain check names the defect.
        @test length(asset["SFOC_P"]) == asset["n_sfoc"] == length(asset["SFOC_g"])
        ok = (knots = PR.table_knots_valid(Float64.(asset["SFOC_P"])),
              values = PR.table_values_valid(Float64.(asset["SFOC_g"])),
              idle = PR.rate_valid(asset["m_dot_idle"]))
        @test ok[kind] == 0
        @test sum(values(ok)) == 2
        @test ctor(; name = :engine) !== nothing
        outcome = run_outcome(analysis)
        @printf("  %-22s run: %s\n", file, outcome)
        @test outcome != :accepted
    end
end

if !isempty(ARGS)
    @testset "compiler diagnostics" begin
        log = read(ARGS[1], String)
        diags = filter(l -> occursin("urn:instantiate:apply-", l), split(log, '\n'))
        @printf("  %d apply diagnostics in %s\n", length(diags), ARGS[1])
        source = split(read(joinpath(pkgdir(FX), "dyad", "InvalidEnginePresets.dyad"), String), '\n')
        reported(file) = (line = findfirst(l -> occursin("presets/Invalid/$(file).toml", l), source);
                          any(d -> occursin("InvalidEnginePresets.dyad:$(line):", d), diags))
        for (file, _) in loader_cases
            @test reported(file)
        end
        for (file, _, _) in length_cases
            @test !reported(file)
        end
        # The files that load are not reported.
        for (file, _, _, _) in value_cases
            @test !reported(file)
        end
    end
end

end
