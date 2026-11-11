# Fresh numerical plans with cold and compiled-code timings. No coefficient
# cache is carried between runs; only Julia's compiled methods are reused.
include(joinpath(@__DIR__,"..","bin","blocks.jl"))

"""Time fresh block computations while separating JIT warmup from numerical-cache reuse."""
function benchmark_blocks(input::String,prefix::String;repetitions::Int=2)
    repetitions>=1 || throw(ArgumentError("benchmark repetitions must be positive"))
    for iteration in 1:repetitions
        GC.gc(true)
        label=iteration==1 ? "cold" : "warm$(iteration-1)"
        output=prefix*"_"*label*".toml"
        isfile(output) && throw(ArgumentError("refusing to overwrite benchmark $output"))
        measurement=@timed run_blocks(input,output)
        report=measurement.value
        get(report,"validated_recovery",true) || throw(ErrorException("diagnostic recovery output is not a validated benchmark"))
        report["benchmark_iteration"]=iteration
        report["fresh_numerical_plan"]=true
        report["complete_driver_seconds"]=measurement.time
        report["harness_sha256"]=bytes2hex(sha256(read(@__FILE__)))
        report["timing_scope"]="fresh numerical caches; iteration one includes cold JIT; later iterations reuse compiled methods only"
        temporary=output*".tmp"
        open(temporary,"w") do io
            TOML.print(io,report;sorted=true)
        end
        mv(temporary,output;force=true)
        println("Benchmark ",label,": numerical_s=",report["seconds"],
            " complete_driver_s=",measurement.time," compilation_s=",get(report,"compile_time",0))
        flush(stdout)
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    2<=length(ARGS)<=3 || error("Usage: benchmark_blocks.jl input.toml output_prefix [repetitions]")
    benchmark_blocks(ARGS[1],ARGS[2];repetitions=length(ARGS)==3 ? parse(Int,ARGS[3]) : 2)
end
