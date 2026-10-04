#!/usr/bin/env julia
# Focused checks of the propeller and hotel-load parameter presets: the files and their
# sidecars, that each wrapper is the one component with the file's values, that the
# propeller presets give the thrust and torque coefficients of the B-series regression at
# their geometry, and that each hotel-load preset draws the cited power.
#
# Runs one small analysis only; never the generated test suite.
#
# Usage:  julia --project=<environment with DyadShip> scripts/validate_component_presets.jl

using DyadShip, Test, TOML
using DyadInterface: symbolic_container
using ModelingToolkit: ModelingToolkit, equations, get_systems, nameof, initial_conditions

const PRESETS = DyadShip.Presets
const S6 = DyadShip.Ship6DOF
const ASSETS = joinpath(pkgdir(DyadShip), "assets", "presets")
const SymbolicUtils = ModelingToolkit.Symbolics.SymbolicUtils

core_of(sys) = only(filter(s -> nameof(s) == :core, get_systems(sys)))
# Value a parameter of `core` was constructed with, as a plain number.
function loaded(sys, name)
    values = Dict(string(k) => SymbolicUtils.unwrap_const(v) for (k, v) in initial_conditions(sys))
    Float64(values["core₊" * name])
end
# A wrapper adds no equation of its own: everything it states is a connection, and its
# `core` has exactly the equations of the bare component.
only_connections(sys) = !isempty(ModelingToolkit.get_eqs(sys)) &&
    all(eq -> occursin("connect(", string(eq)), ModelingToolkit.get_eqs(sys))
same_equations(core, bare) = string.(equations(core)) == string.(equations(bare))
function applied_files(dir)
    files = filter(f -> endswith(f, ".toml") && !endswith(f, ".provenance.toml"), readdir(dir))
    @test sort(files) == sort(filter(f -> isfile(joinpath(dir, replace(f, ".toml" => ".provenance.toml"))), files))
    files
end
function sidecar_of(dir, file, keys_)
    sidecar = TOML.parsefile(joinpath(dir, replace(file, ".toml" => ".provenance.toml")))
    @test sidecar["schema"] == "dyadship.provenance/1"
    @test sidecar["asset"] == file
    @test sort(collect(keys(sidecar["parameters"]))) == keys_
    @test all(p -> haskey(p, "unit") && haskey(p, "origin"), values(sidecar["parameters"]))
    @test !isempty(sidecar["source"]) && all(s -> !isempty(s["title"]) && !isempty(s["role"]), sidecar["source"])
    sidecar
end

@testset "propeller and hotel-load parameter presets" begin

@testset "Wageningen B-series propeller presets" begin
    dir = joinpath(ASSETS, "Propeller")
    presets = [
        (file = "wageningen_b4_55.toml", wrapper = PRESETS.WageningenB4_55Propeller, Z = 4.0, Ae_Ao = 0.55, P_D = 1.0, Diameter = 4.0),
        (file = "wageningen_b4_70.toml", wrapper = PRESETS.WageningenB4_70Propeller, Z = 4.0, Ae_Ao = 0.70, P_D = 0.8, Diameter = 6.0),
        (file = "wageningen_b5_75.toml", wrapper = PRESETS.WageningenB5_75Propeller, Z = 5.0, Ae_Ao = 0.75, P_D = 1.0, Diameter = 7.0),
    ]
    @test sort(applied_files(dir)) == sort([p.file for p in presets])
    # The regression the component evaluates has the published number of terms.
    @test size(S6.WAGENINGEN_KT) == (39, 5) && size(S6.WAGENINGEN_KQ) == (47, 5)
    bare = S6.Propeller1Q(; name = :core)
    for p in presets
        data = TOML.parsefile(joinpath(dir, p.file))
        @test sort(collect(keys(data))) == ["Ae_Ao", "Diameter", "P_D", "Z"]
        @test (data["Z"], data["Ae_Ao"], data["P_D"], data["Diameter"]) == (p.Z, p.Ae_Ao, p.P_D, p.Diameter)
        # Inside the range the component states for the regression.
        @test 2 <= data["Z"] <= 7 && 0.30 <= data["Ae_Ao"] <= 1.05 && 0.5 <= data["P_D"] <= 1.4
        @test data["Diameter"] > 0
        sidecar = sidecar_of(dir, p.file, ["Ae_Ao", "Diameter", "P_D", "Z"])
        @test sidecar["component"] == "DyadShip.Ship6DOF.Propeller1Q"
        @test occursin("Oosterveld", sidecar["source"][1]["title"])
        @test occursin("not from any document", sidecar["parameters"]["Diameter"]["origin"])

        # The wrapper is the one Propeller1Q with the file's values and no equation of its own.
        sys = p.wrapper(; name = :prop)
        core = core_of(sys)
        @test same_equations(core, bare)
        @test Set(nameof.(get_systems(sys))) == Set([:core, :flange, :frame_a])
        @test only_connections(sys)
        for name in ("Z", "Ae_Ao", "P_D", "Diameter")
            @test loaded(sys, name) == data[name]
        end

        # KT and KQ of the preset geometry at a few advance ratios: the regression evaluated
        # at the values the wrapper loaded, against the polynomial summed term by term here.
        Z, A, P = loaded(sys, "Z"), loaded(sys, "Ae_Ao"), loaded(sys, "P_D")
        polynomial(T, J) = sum(T[i, 1] * J^T[i, 2] * P^T[i, 3] * A^T[i, 4] * Z^T[i, 5] for i in axes(T, 1))
        for J in (0.0, 0.3, 0.6)
            kt, kq = S6.wageningen_kt(J, P, A, Z), S6.wageningen_kq(J, P, A, Z)
            @test kt ≈ polynomial(S6.WAGENINGEN_KT, J) rtol = 1e-12
            @test kq ≈ polynomial(S6.WAGENINGEN_KQ, J) rtol = 1e-12
            # Open-water sanity inside the working range: positive thrust and torque,
            # both falling with advance ratio, efficiency below one.
            @test kt > 0 && kq > 0
            J > 0 && @test J / (2π) * kt / kq < 1
        end
        @test S6.wageningen_kt(0.0, P, A, Z) > S6.wageningen_kt(0.3, P, A, Z) > S6.wageningen_kt(0.6, P, A, Z)
        @test S6.wageningen_kq(0.0, P, A, Z) > S6.wageningen_kq(0.3, P, A, Z) > S6.wageningen_kq(0.6, P, A, Z)
    end
    # The first preset is the propeller the sample ship already uses, value for value.
    @test (presets[1].Diameter, presets[1].P_D, presets[1].Ae_Ao, presets[1].Z) == (4.0, 1.0, 0.55, 4.0)
    # Presets differ where their files differ: more pitch and area give more thrust at J = 0.3.
    kt(p) = S6.wageningen_kt(0.3, p.P_D, p.Ae_Ao, p.Z)
    @test kt(presets[1]) > kt(presets[2])       # P/D 1.0 against 0.8
end

@testset "hotel-load presets" begin
    dir = joinpath(ASSETS, "HotelLoad")
    # Auxiliary engine power output in kW as printed in the cited table (at berth, anchored,
    # manoeuvring, at sea), typed here independently of the files.
    table = Dict("bulk_carrier_60k_100k_dwt" => (240, 400, 1100, 410),
                 "container_ship_8k_12k_teu" => (1150, 1600, 2900, 1800))
    presets = [
        (file = "hotel_load_bulk_carrier_60k_100k_dwt_at_sea.toml", wrapper = PRESETS.BulkCarrierHotelLoadAtSea, signal = :bulk_at_sea, kw = table["bulk_carrier_60k_100k_dwt"][4]),
        (file = "hotel_load_bulk_carrier_60k_100k_dwt_at_berth.toml", wrapper = PRESETS.BulkCarrierHotelLoadAtBerth, signal = :bulk_at_berth, kw = table["bulk_carrier_60k_100k_dwt"][1]),
        (file = "hotel_load_container_ship_8k_12k_teu_at_sea.toml", wrapper = PRESETS.ContainerShipHotelLoadAtSea, signal = :container_at_sea, kw = table["container_ship_8k_12k_teu"][4]),
        (file = "hotel_load_container_ship_8k_12k_teu_at_berth.toml", wrapper = PRESETS.ContainerShipHotelLoadAtBerth, signal = :container_at_berth, kw = table["container_ship_8k_12k_teu"][1]),
    ]
    @test sort(applied_files(dir)) == sort([p.file for p in presets])
    result = PRESETS.HotelLoadPresetTransient()
    @test string(result.sol.retcode) == "Success"
    rig = symbolic_container(result)
    for p in presets
        data = TOML.parsefile(joinpath(dir, p.file))
        @test sort(collect(keys(data))) == ["CycleAmplitude", "Kr", "NominalPower"]
        @test data["NominalPower"] == 1000.0 * p.kw && data["Kr"] == 1.0 && data["CycleAmplitude"] == 0.0
        sidecar = sidecar_of(dir, p.file, ["CycleAmplitude", "Kr", "NominalPower"])
        @test sidecar["component"] == "DyadShip.Machinery.OnOffConsumer"
        @test occursin("not a ship", sidecar["calibration_status"])
        source = only(sidecar["source"])
        @test occursin("Fourth IMO GHG Study 2020", source["title"]) && source["page"] == "68"
        @test occursin("Table 17", source["table"]) && length(source["sha256"]) == 64
        ship = first(k for k in keys(table) if occursin(k, p.file))
        row = sidecar["table_row"]
        @test (row["at_berth"], row["anchored"], row["manoeuvring"], row["at_sea"]) == table[ship]
        @test occursin("$(p.kw) kW", sidecar["parameters"]["NominalPower"]["origin"])

        sys = p.wrapper(; name = :load)
        @test nameof.(get_systems(sys)) == [:core]
        @test only_connections(sys)
        @test same_equations(core_of(sys), DyadShip.Machinery.OnOffConsumer(; name = :core))
        @test loaded(sys, "NominalPower") == data["NominalPower"]
        @test loaded(sys, "Kr") == 1.0 && loaded(sys, "CycleAmplitude") == 0.0
        # Switched on, the consumer draws the cited power once its start-up ramp is over.
        load = getproperty(rig, p.signal)
        @test result.sol[load.y][end] ≈ 1000.0 * p.kw rtol = 1e-6
        @test result.sol[load.y][1] == 0
        @test maximum(result.sol[load.y]) <= 1000.0 * p.kw * (1 + 1e-9)
    end
end

end
