# Local-process integration worker. Communication uses atomic TOML files, so
# the portable parallel path needs no sockets or external service.
include("run_cross_channel.jl")

"""Decode freshly prepared geometry/spin data from a coordinator's explicit job file."""
function worker_cross_surfaces(records)
    [(id=String(row["id"]),source_q=Tuple(ComplexF64(z...) for z in row["source_q"]),
        target_q=Tuple(ComplexF64(z...) for z in row["target_q"]),
        source_spins=Tuple(Tuple(Int.(s)) for s in row["source_spins"]),
        target_spins=Tuple(Tuple(Int.(s)) for s in row["target_spins"]),
        free_ratio=Tuple(Float64.(row["free_ratio"]))) for row in records]
end

"""Compute one disjoint grid chunk and retain partial values or a complete failure record."""
function run_cross_worker(input::String,output::String)
    isfile(output) && throw(ArgumentError("refusing to overwrite a worker result"))
    payload=TOML.parsefile(input)
    surfaces=worker_cross_surfaces(payload["surfaces"])
    indices=[Tuple(Int.(point)) for point in payload["point_indices"]]
    started=time_ns()
    try
        worker_progress=state->save_cross_report(output,state)
        result=integrate_cross_channel(payload["config"],surfaces,Bool(payload["ramond"]),Int(payload["order"]);
            progress=worker_progress,point_indices=indices)
        result.report["completed"]=true
        result.report["nodes_done"]=result.report["nodes"]
        result.report["process_id"]=getpid()
        result.report["input_sha256"]=bytes2hex(sha256(read(input)))
        save_cross_report(output,result.report)
    catch error
        failure=isfile(output) ? TOML.parsefile(output) : Dict{String,Any}()
        failure["completed"]=false
        failure["failure"]=sprint(showerror,error)
        failure["elapsed_seconds"]=(time_ns()-started)/1e9
        save_cross_report(output,failure)
        rethrow()
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==2 || error("Usage: cross_worker.jl coordinator_input.toml result.toml")
    run_cross_worker(ARGS[1],ARGS[2])
end
