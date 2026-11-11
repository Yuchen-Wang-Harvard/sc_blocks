# Directed validation of the integration driver; deliberately below paper cutoffs.
include("run_cross_channel.jl")
include("compare_cross_channel.jl")
using Test

"""Validate series evaluation, all geometries and isolated-worker reduction on a small grid.

The resulting integrals are smoke tests only. Their partial block cutoffs and
two small quadrature orders are never promoted to the Appendix C release check.
"""
function run_cross_smoke(directory::String)
    ispath(directory) && throw(ArgumentError("refusing to overwrite smoke-test evidence"))
    mkpath(directory)
    root=normpath(joinpath(@__DIR__,".."))
    report=Dict{String,Any}("passed"=>false,"scope"=>"driver_smoke_not_paper_cutoffs",
        "timestamp_utc"=>string(now(UTC)),"julia"=>string(VERSION),"cpu"=>Sys.CPU_NAME,
        "source_sha256"=>source_hashes(root))
    started=time_ns()
    try
        unit=@timed include(joinpath(root,"test","test_cross_driver.jl"))
        report["unit_seconds"]=unit.time
        config=TOML.parsefile(joinpath(root,"examples","cross_channel.toml"))
        config["source_total_level"]=2
        config["target_total_level"]=3
        config["block_precision_bits"]=96
        config["source_orders"]=[1,2]
        config["target_orders"]=[2,3]
        config["constant_precision_bits"]=96
        config["upsilon_order"]=256
        config["workers"]=1
        serial_input=joinpath(directory,"serial_input.toml")
        serial_output=joinpath(directory,"serial.toml")
        save_cross_report(serial_input,config)
        serial=@timed run_cross_channel(serial_input,serial_output)
        report["serial_seconds"]=serial.time
        @test serial.value["completed"]
        @test !serial.value["release_validated"]
        @test length(serial.value["geometry"])==10
        @test length(last(serial.value["runs"])["ratios"])==20
        @test last(serial.value["runs"])["source"]["nodes"]==8
        @test last(serial.value["runs"])["target"]["nodes"]==27
        config["workers"]=2
        parallel_input=joinpath(directory,"parallel_input.toml")
        parallel_output=joinpath(directory,"parallel.toml")
        save_cross_report(parallel_input,config)
        parallel=@timed run_cross_channel(parallel_input,parallel_output)
        report["parallel_seconds"]=parallel.time
        comparison=compare_cross_channel(serial_output,parallel_output,joinpath(directory,"serial_parallel.toml");tolerance=1e-12)
        @test comparison["passed"]
        # A separate target-method comparison at one momentum, without reading
        # either production or reference coefficients from a saved result.
        config["surface"]=config["surface"][1:1]
        config["workers"]=1
        config["block_precision_bits"]=136
        surfaces,_=prepare_cross_surfaces(config)
        native=@timed integrate_cross_channel(config,surfaces,false,1)
        config["target_method"]="double_virasoro"
        doubled=@timed integrate_cross_channel(config,surfaces,false,1)
        relative=maximum(abs,native.value.values-doubled.value.values)/maximum(abs,native.value.values)
        @test relative<1e-12
        report["target_method_relative_difference"]=relative
        report["target_native_seconds"]=native.time
        report["target_double_virasoro_seconds"]=doubled.time
        config["reuse_vertex_label_products"]=true
        cached=@timed integrate_cross_channel(config,surfaces,true,1)
        config["reuse_vertex_label_products"]=false
        uncached=@timed integrate_cross_channel(config,surfaces,true,1)
        reuse_error=maximum(abs,cached.value.values-uncached.value.values)/maximum(abs,uncached.value.values)
        @test reuse_error<1e-12
        @test cached.value.report["product_cache_hits"]>0
        @test cached.value.report["ccy_transitions"]<uncached.value.report["ccy_transitions"]
        report["source_reuse_maximum_scaled_difference"]=reuse_error
        report["source_cached_seconds"]=cached.time
        report["source_uncached_seconds"]=uncached.time
        report["source_product_cache_hits"]=cached.value.report["product_cache_hits"]
        report["source_cached_transitions"]=cached.value.report["ccy_transitions"]
        report["source_uncached_transitions"]=uncached.value.report["ccy_transitions"]
        report["passed"]=true
    catch error
        report["failure"]=sprint(showerror,error)
        rethrow()
    finally
        report["elapsed_seconds"]=(time_ns()-started)/1e9
        save_cross_report(joinpath(directory,"summary.toml"),report)
    end
    report
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==1 || error("Usage: run_cross_smoke.jl fresh_output_directory")
    run_cross_smoke(only(ARGS))
end
