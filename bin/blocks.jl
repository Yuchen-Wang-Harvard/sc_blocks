# Standalone coefficient driver for closed all-NS and mixed NSNSNS/NSRR graphs.
using SuperconformalBlocks
using TOML
using Dates
using SHA

"""Parse an integer, decimal or rational string directly at the working precision."""
function parse_real(value)
    value isa Integer && return BigFloat(value)
    value isa AbstractString || throw(ArgumentError("noninteger numeric inputs must be quoted strings"))
    pieces=split(value,'/')
    length(pieces)==1 && return parse(BigFloat,only(pieces))
    length(pieces)==2 || throw(ArgumentError("invalid real parameter: $value"))
    parse(BigFloat,pieces[1])/parse(BigFloat,pieces[2])
end

"""Parse a real scalar or [real,imaginary] string pair without Float64 rounding."""
function parse_complex(value)
    value isa AbstractVector && return length(value)==2 ?
        Complex{BigFloat}(parse_real(value[1]),parse_real(value[2])) :
        throw(ArgumentError("complex parameter needs [real,imaginary]"))
    Complex{BigFloat}(parse_real(value))
end

"""Hash the production Julia sources, independently of git or the parent repository."""
function production_hashes(root)
    hashes=Dict{String,String}()
    for folder in ("src","bin"),(directory,_,files) in walkdir(joinpath(root,folder)),file in sort(files)
        endswith(file,".jl") || continue
        path=joinpath(directory,file)
        hashes[relpath(path,root)]=bytes2hex(sha256(read(path)))
    end
    hashes
end

"""Compute one graph configuration and atomically save coefficients and measured provenance."""
function run_blocks(input,output)
    any(ispath,(output,output*".failed.toml")) && throw(ArgumentError("refusing to overwrite block-computation evidence"))
    config=TOML.parsefile(input)
    sector=get(config,"sector","NS")
    sector in ("NS","NSRR") || throw(ArgumentError("sector must be NS or NSRR"))
    bits=Int(get(config,"precision_bits",136))
    bits>=64 || throw(ArgumentError("precision_bits must be at least 64"))
    root=normpath(joinpath(@__DIR__,".."))
    hashes=production_hashes(root)
    graphname=config["graph"]
    graph=graphname=="theta" ? theta_graph() : graphname=="glasses" ? glasses_graph() :
          graphname=="tetrahedron" ? tetrahedron_graph() : throw(ArgumentError("unknown graph"))
    ne=length(graph.ends)
    setprecision(BigFloat,bits) do
        momenta=Tuple(parse_complex(value) for value in config["momenta"])
        parities=Tuple(Int.(get(config,"primary_parities",zeros(Int,ne))))
        forms=Int.(config["forms"])
        plan=if sector=="NS"
            NSBlockPlan(BigFloat,graph,parse_complex(config["b"]),momenta,forms;primary_parities=parities)
        else
            RamondBlockPlan(BigFloat,graph,parse_complex(config["b"]),momenta,Tuple(Bool.(config["ramond"])),
                forms,Int.(config["eta"]);primary_parities=parities,
                recovery_tolerance=haskey(config,"recovery_tolerance") ? parse_real(config["recovery_tolerance"]) : inv(BigFloat(10)^8),
                recovery_policy=Symbol(get(config,"recovery_policy","error")))
        end
        cutoff=config["cutoff"]
        kind=cutoff["kind"]
        domain=kind=="total" ? total_domain(ne,2Int(cutoff["level"])) : kind=="per_edge" ?
            box_domain(ntuple(_->2Int(cutoff["level"]),ne)) : throw(ArgumentError("unknown cutoff kind"))
        started=time_ns()
        measurement=try
            @timed sector=="NS" ? ns_block(plan,domain) : ramond_block(plan,domain)
        catch error
            failure=Dict{String,Any}("completed"=>false,"validated_recovery"=>false,
                "convention"=>CONVENTION,"input"=>config,"precision_bits"=>bits,
                "seconds_before_failure"=>(time_ns()-started)/1e9,"failure"=>sprint(showerror,error),
                "source_sha256"=>hashes,"timestamp_utc"=>string(now(UTC)))
            if sector=="NSRR"
                failure["stage_seconds"]=plan.timings
                failure["recovery_diagnostics"]=plan.recovery_diagnostics
            end
            mkpath(dirname(abspath(output)))
            open(output*".failed.toml","w") do io
                TOML.print(io,failure;sorted=true)
            end
            rethrow()
        end
        answer=measurement.value
        all(value->value isa Number ? isfinite(value) : all(isfinite,value),values(answer)) ||
            throw(ErrorException("nonfinite block coefficient; use generic input parameters or an explicit limiting prescription"))
        records=if sector=="NS"
            [Dict("level2"=>collect(key),"real"=>string(real(get(answer,key,0))),
                  "imaginary"=>string(imag(get(answer,key,0)))) for key in domain.indices]
        else
            [Dict("level2"=>collect(key),"parity_real"=>string.(real.(get(answer,key,zeros(Complex{BigFloat},2^ne)))),
                  "parity_imaginary"=>string.(imag.(get(answer,key,zeros(Complex{BigFloat},2^ne))))) for key in domain.indices]
        end
        report=Dict{String,Any}("convention"=>CONVENTION,"implementation_status"=>"partial_release",
            "input"=>config,"timestamp_utc"=>string(now(UTC)),"julia"=>string(VERSION),
            "machine"=>Sys.MACHINE,"cpu"=>Sys.CPU_NAME,"threads"=>Threads.nthreads(),
            "precision_bits"=>bits,"seconds"=>measurement.time,"allocated_bytes"=>measurement.bytes,
            "gc_seconds"=>measurement.gctime,"branches"=>plan.branches,"ccy_transitions"=>plan.ccy_transitions,
            "source_sha256"=>hashes,"coefficients"=>records)
        if sector=="NSRR"
            report["validated_recovery"]=plan.recovery_policy==:error && plan.recovery_diagnostics["passed"]
            report["recovery_diagnostics"]=plan.recovery_diagnostics
            report["inserted_edges"]=plan.marked
            report["forward_residual"]=string(plan.forward_residual)
            report["stage_seconds"]=plan.timings
            report["parity_order"]="index=1+mask; edge e uses bit e-1; includes intrinsic NS parity"
        end
        for field in (:compile_time,:recompile_time)
            hasproperty(measurement,field) && (report[string(field)]=getproperty(measurement,field))
        end
        mkpath(dirname(abspath(output)))
        temporary=output*".tmp"
        open(temporary,"w") do io
            TOML.print(io,report;sorted=true)
        end
        mv(temporary,output;force=true)
        println("Saved ",length(records)," coefficients to ",output,
                "; seconds=",measurement.time,"; branches=",plan.branches)
        report
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==2 || error("Usage: julia --project=. bin/blocks.jl input.toml output.toml")
    run_blocks(ARGS[1],ARGS[2])
end
