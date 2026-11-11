#= Run the component checks and save reproducible timing/provenance.

Full-cutoff graph comparisons and cross-channel integrals have separate
drivers in this directory; see README.md and docs/CODE_MAP.md.
=#
using Dates
using SHA
using TOML

"""Hash every Julia source, CLI, check and test file relative to the portable package root."""
function source_hashes(root)
    hashes=Dict{String,String}()
    for folder in ("src","bin","checks","test"), (directory,_,files) in walkdir(joinpath(root,folder)), file in sort(files)
        endswith(file,".jl") || continue
        path=joinpath(directory,file)
        hashes[relpath(path,root)]=bytes2hex(sha256(read(path)))
    end
    hashes
end

"""Run one fresh-process suite and write elapsed/allocation/compilation measurements."""
function run_foundations(output)
    ispath(output) && throw(ArgumentError("refusing to overwrite foundation-check evidence"))
    root=normpath(joinpath(@__DIR__,".."))
    report=Dict{String,Any}(
        "scope"=>"foundations_only_not_full_production_validation",
        "timestamp_utc"=>string(now(UTC)),"julia"=>string(VERSION),
        "machine"=>Sys.MACHINE,"cpu"=>Sys.CPU_NAME,"threads"=>Threads.nthreads(),
        "precision_bits"=>192,
        "precision_note"=>"main algebra/PBW comparisons use 192 bits; directed precision tests use 96/136/192 bits and analytic tests also use Float64, as specified in each test file",
        "source_sha256"=>source_hashes(root),
        "passed"=>false)
    started=time_ns()
    try
        measurement=@timed include(joinpath(root,"test","runtests.jl"))
        report["passed"]=true
        report["seconds"]=measurement.time
        report["allocated_bytes"]=measurement.bytes
        report["gc_seconds"]=measurement.gctime
        for field in (:compile_time,:recompile_time)
            hasproperty(measurement,field) && (report[string(field)]=getproperty(measurement,field))
        end
    catch error
        report["seconds"]=(time_ns()-started)/1e9
        report["failure"]=sprint(showerror,error)
        rethrow()
    finally
        mkpath(dirname(abspath(output)))
        open(output,"w") do io
            TOML.print(io,report;sorted=true)
        end
    end
    println("Saved foundation check record: ",output)
    report
end

if abspath(PROGRAM_FILE)==@__FILE__
    output=isempty(ARGS) ? joinpath(@__DIR__,"..","results","foundations.toml") : only(ARGS)
    run_foundations(output)
end
