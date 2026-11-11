# Independent integral checks and explicit arithmetic/quadrature escalation of
# the Appendix C constants, including momenta in the integration-grid tail.
using SuperconformalBlocks
using TOML
using Dates
using SHA
include("run_foundations.jl")

"""Compare the fast constants with a higher-order, higher-precision calculation.

Independent tanh-sinh Upsilon integrals separately test the Gaussian evaluator
at oscillatory arguments. Both constant evaluators use the same displayed
Appendix C normalization; this is a numerical precision check, not a new
normalization convention or an archived-value comparison.
"""
function run_liouville_precision(output::String)
    isfile(output) && throw(ArgumentError("refusing to overwrite constant-precision evidence"))
    root=normpath(joinpath(@__DIR__,".."))
    report=Dict{String,Any}("passed"=>false,"convention"=>CONVENTION,
        "timestamp_utc"=>string(now(UTC)),"julia"=>string(VERSION),"cpu"=>Sys.CPU_NAME,
        "source_sha256"=>source_hashes(root),"cases"=>Dict{String,Any}[],
        "fast_precision_bits"=>53,"reference_precision_bits"=>192,
        "fast_order"=>512,"reference_order"=>1024,"fast_tmax"=>120,"reference_tmax"=>180,
        "b_input"=>1.4,"parameter_policy"=>"the Float64 input b and momenta are represented exactly in BigFloat to isolate numerical evaluation errors",
        "constant_relative_tolerance"=>1e-10,"upsilon_absolute_tolerance"=>1e-10)
    started=time_ns()
    try
        setprecision(BigFloat,192) do
            fast=LiouvillePlan(Float64,7/5;order=512,tmax=120)
            high=LiouvillePlan(BigFloat,BigFloat(fast.b);order=1024,tmax=180)
            maximum_error=big(0)
            for p in ((0.05,0.15,0.3),(0.3,0.5,0.7),(1.2,1.7,2.1),(2.8,3.1,3.4))
                low_time=@timed (liouville_ns(fast,p),liouville_nsrr(fast,p))
                high_time=@timed (liouville_ns(high,p),liouville_nsrr(high,p))
                error=maximum(abs(Complex{BigFloat}(a)-b)/max(abs(b),BigFloat(10)^(-100))
                    for (low,ref) in zip(low_time.value,high_time.value) for (a,b) in zip(low,ref))
                maximum_error=max(maximum_error,error)
                push!(report["cases"],Dict("momenta"=>collect(p),"maximum_relative_error"=>string(error),
                    "fast_seconds"=>low_time.time,"reference_seconds"=>high_time.time,
                    "passed"=>error<=big"1e-10"))
                println("Liouville momenta=",p," relative error=",error)
                flush(stdout)
                error<=big"1e-10" || throw(ErrorException("constant precision comparison failed at $p: $error"))
            end
            report["maximum_relative_constant_error"]=string(maximum_error)
            integral_checks=Dict{String,Any}[]
            for x in (Complex{BigFloat}(high.Q/4,BigFloat(93)/20),
                      Complex{BigFloat}(inv(high.b)/2,BigFloat(17)/5),
                      Complex{BigFloat}(inv(high.b),BigFloat(17)/5))
                independent=@timed reference_log_upsilon(high.b,x;bits=120)
                error=abs(log_upsilon(fast,ComplexF64(x))-independent.value.value)
                push!(integral_checks,Dict("argument"=>[string(real(x)),string(imag(x))],
                    "absolute_error"=>string(error),"reference_seconds"=>independent.time,
                    "reference_evaluations"=>independent.value.evaluations,"passed"=>error<=big"1e-10"))
                error<=big"1e-10" || throw(ErrorException("independent oscillatory Upsilon integral failed: $error"))
            end
            report["independent_integrals"]=integral_checks
        end
        report["passed"]=true
    catch error
        report["failure"]=sprint(showerror,error)
        rethrow()
    finally
        report["elapsed_seconds"]=(time_ns()-started)/1e9
        mkpath(dirname(abspath(output)))
        open(output,"w") do io
            TOML.print(io,report;sorted=true)
        end
    end
    report
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==1 || error("Usage: run_liouville_precision.jl output.toml")
    run_liouville_precision(only(ARGS))
end
