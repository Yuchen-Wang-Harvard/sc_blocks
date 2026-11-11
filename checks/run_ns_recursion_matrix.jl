# Full all-NS double-Virasoro versus independent SCA c-recursion release gate.
using SuperconformalBlocks
using TOML
using Dates
using SHA
include("run_foundations.jl")

"""Write a partial or completed NS comparison record atomically."""
function save_ns_matrix(output,report)
    mkpath(dirname(abspath(output)))
    temporary=output*".tmp"
    open(temporary,"w") do io
        TOML.print(io,report;sorted=true)
    end
    mv(temporary,output;force=true)
end

"""Compare every allowed vertex form at the explicit theta/tetrahedron cutoff.

The independent reference computes all vertex-form sectors in one pass. Each
grade has one NS parity mask, so the signed coefficient comparison certifies
every fixed-spin evaluation; the report distinguishes that fact from executed
spin sums. No stored block coefficients are read.
"""
function run_ns_recursion_matrix(family::String,level::Int,output::String;bits::Int=192,primarymask::Int=0,
                                 momentum_mode::Symbol=:real)
    family in ("theta","tetrahedron") || throw(ArgumentError("family must be theta or tetrahedron"))
    isfile(output) && throw(ArgumentError("refusing to overwrite an existing comparison"))
    level>=0 && bits>=64 || throw(ArgumentError("invalid cutoff or precision"))
    momentum_mode in (:real,:imaginary) || throw(ArgumentError("momentum_mode must be real or imaginary"))
    root=normpath(joinpath(@__DIR__,".."))
    report=Dict{String,Any}("family"=>family,"level"=>level,"precision_bits"=>bits,
        "cutoff_kind"=>(family=="theta" ? "per_edge" : "total"),"primary_mask"=>primarymask,
        "completed"=>false,"passed"=>false,"convention"=>CONVENTION,"momentum_mode"=>string(momentum_mode),
        "timestamp_utc"=>string(now(UTC)),"julia"=>string(VERSION),"cpu"=>Sys.CPU_NAME,
        "process_id"=>getpid(),
        "source_sha256"=>source_hashes(root),"cases"=>Dict{String,Any}[])
    started=time_ns()
    try
        setprecision(BigFloat,bits) do
            graph=family=="theta" ? theta_graph() : tetrahedron_graph()
            count=length(graph.ends)
            primarymask in 0:2^count-1 || throw(ArgumentError("invalid primary parity mask"))
            parities=ntuple(e->(primarymask>>(e-1))&1,count)
            momenta=map(Complex{BigFloat},family=="theta" ? (11//23,13//29,17//31) :
                (11//23,13//29,17//31,19//37,23//41,29//43))
            momentum_mode==:imaginary && (momenta=map(p->im*p,momenta))
            b=Complex{BigFloat}(7//5)
            Q=b+inv(b)
            weights=map(p->Q*Q/8-p*p/2,momenta)
            central=Complex{BigFloat}(3//2)+3Q*Q
            domain=family=="theta" ? box_domain(ntuple(_->2level,count)) : total_domain(count,2level)
            reference_plan=NSRecursionPlan(BigFloat,graph,central,weights;primary_parities=parities)
            report["reference_compact_theta"]=reference_plan.compact_theta
            report["stage"]="independent_ns_recursion"
            report["b"]=[string(real(b)),string(imag(b))]
            report["momenta"]=[[string(real(p)),string(imag(p))] for p in momenta]
            report["multigrades"]=length(domain.indices)
            save_ns_matrix(output,report)
            reference_progress=state->begin
                report["reference_progress"]=Dict(string(key)=>value for (key,value) in pairs(state))
                report["elapsed_seconds"]=(time_ns()-started)/1e9
                save_ns_matrix(output,report)
                if state.stage=="completed_parity"
                    println("NS reference mask=",state.output_mask," terms=",state.seed_terms,
                        " cached_vertices=",state.cached_vertices," elapsed_s=",report["elapsed_seconds"])
                    flush(stdout)
                end
            end
            reference=@timed ns_c_block(reference_plan,domain;progress=reference_progress)
            report["reference_seconds"]=reference.time
            report["reference_allocated_bytes"]=reference.bytes
            report["reference_compile_seconds"]=hasproperty(reference,:compile_time) ? reference.compile_time : 0.0
            report["reference_transitions"]=reference_plan.transitions
            report["reference_seed_terms"]=reference_plan.seed_terms
            report["reference_cached_inverse_norms"]=length(reference_plan.inverse_norms)
            tolerance=BigFloat(10)^(-18)
            report["tolerance"]=string(tolerance)
            report["spin_policy"]="each grade has one supported NS parity; equality of that coefficient certifies all 2^edges spin lifts"
            plan=NSBlockPlan(BigFloat,graph,b,momenta,zeros(Int,length(graph.slots));primary_parities=parities)
            forms_=unique([Tuple(mod(sum(key[e] for e in slots),2) for slots in graph.slots) for key in domain.indices])
            sort!(forms_)
            for forms in forms_
                plan.forms.=forms
                foreach(empty!,plan.vertices)
                report["stage"]="double_virasoro"
                report["current_forms"]=collect(forms)
                save_ns_matrix(output,report)
                production=@timed ns_block(plan,domain)
                maximum_error=zero(BigFloat)
                worst=first(domain.indices)
                supported=0
                for key in domain.indices
                    allowed=all(v->mod(sum(key[e] for e in graph.slots[v]),2)==forms[v],eachindex(forms))
                    expected=allowed ? get(reference.value,key,zero(Complex{BigFloat})) : zero(Complex{BigFloat})
                    actual=get(production.value,key,zero(Complex{BigFloat}))
                    error=abs(actual-expected)/max(big(1),abs(actual),abs(expected))
                    isfinite(error) || throw(ErrorException("nonfinite NS comparison at $key"))
                    if error>maximum_error
                        maximum_error,worst=error,key
                    end
                    supported+=allowed
                end
                case=Dict{String,Any}("forms"=>collect(forms),"maximum_scaled_error"=>string(maximum_error),
                    "worst_level2"=>collect(worst),"supported_multigrades"=>supported,
                    "production_seconds"=>production.time,"production_allocated_bytes"=>production.bytes,
                    "production_compile_seconds"=>(hasproperty(production,:compile_time) ? production.compile_time : 0.0),
                    "passed"=>maximum_error<=tolerance)
                push!(report["cases"],case)
                report["elapsed_seconds"]=(time_ns()-started)/1e9
                save_ns_matrix(output,report)
                maximum_error<=tolerance || throw(ErrorException("NS recursion comparison failed: $maximum_error at $worst"))
                println("NS comparison forms=",forms," production_s=",production.time," error=",maximum_error)
                flush(stdout)
            end
            report["completed"]=true
            report["passed"]=true
            delete!(report,"current_forms")
            report["stage"]="completed"
        end
    catch error
        report["failure"]=sprint(showerror,error)
        rethrow()
    finally
        report["elapsed_seconds"]=(time_ns()-started)/1e9
        save_ns_matrix(output,report)
    end
    report
end

if abspath(PROGRAM_FILE)==@__FILE__
    3<=length(ARGS)<=6 || error("Usage: run_ns_recursion_matrix.jl theta|tetrahedron level output.toml [bits] [primary_mask] [real|imaginary]")
    run_ns_recursion_matrix(ARGS[1],parse(Int,ARGS[2]),ARGS[3];bits=length(ARGS)>=4 ? parse(Int,ARGS[4]) : 192,
        primarymask=length(ARGS)>=5 ? parse(Int,ARGS[5]) : 0,momentum_mode=length(ARGS)>=6 ? Symbol(ARGS[6]) : :real)
end
