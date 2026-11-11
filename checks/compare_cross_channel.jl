# Compare freshly generated integration runs; never treat this as an independent
# coefficient oracle or as evidence beyond the explicit input cutoffs.
using TOML
using SHA
using Dates

"""Read keyed partition integrals or ratios from one completed quadrature entry."""
function cross_comparison_values(entry,channel)
    records=channel=="ratios" ? entry["ratios"] : entry[channel]["integrals"]
    Dict((string(row["id"]),Int(get(row,"level2",-1)),Int(row["spin_index"]))=>
        ComplexF64((channel=="ratios" ? row["ratio"] : row["value"])...) for row in records)
end

"""Compare precision, Gaussian reference, target method or quadrature variants.

Geometry seeds, spin labels, b and block cutoffs must agree. Both input runs
must have completed; their ratio pass flag is reported separately. The last
quadrature order in each is compared, allowing an order-convergence check.
"""
function compare_cross_channel(left_path::String,right_path::String,output::String;tolerance::Float64=1e-6)
    isfile(output) && throw(ArgumentError("refusing to overwrite a comparison"))
    isfinite(tolerance) && tolerance>0 || throw(ArgumentError("positive finite tolerance required"))
    left,right=TOML.parsefile(left_path),TOML.parsefile(right_path)
    left["completed"] && right["completed"] || throw(ArgumentError("both integrations must have completed"))
    left["convention"]==right["convention"] || throw(ArgumentError("different block conventions"))
    for key in ("b","source_total_level","target_total_level","surface")
        left["input"][key]==right["input"][key] || throw(ArgumentError("incompatible input $key"))
    end
    report=Dict{String,Any}("timestamp_utc"=>string(now(UTC)),"tolerance"=>tolerance,
        "left_sha256"=>bytes2hex(sha256(read(left_path))),"right_sha256"=>bytes2hex(sha256(read(right_path))),
        "left_input"=>left["input"],"right_input"=>right["input"],
        "left_ratio_checks_passed"=>left["passed"],"right_ratio_checks_passed"=>right["passed"],
        "passed"=>true,"comparisons"=>Dict{String,Any}[])
    for channel in ("source","target","ratios")
        a=cross_comparison_values(last(left["runs"]),channel)
        b=cross_comparison_values(last(right["runs"]),channel)
        Set(keys(a))==Set(keys(b)) || throw(ArgumentError("different integral/spin domains"))
        maximum_error=0.0
        worst=first(keys(a))
        for key in keys(a)
            scale=channel=="ratios" ? 1.0 : max(abs(a[key]),abs(b[key]),floatmin(Float64))
            error=abs(a[key]-b[key])/scale
            isfinite(error) || throw(ErrorException("nonfinite comparison"))
            if error>maximum_error
                maximum_error,worst=error,key
            end
        end
        push!(report["comparisons"],Dict("channel"=>channel,"maximum_error"=>maximum_error,
            "worst_surface"=>worst[1],"worst_level2"=>worst[2],"worst_spin_index"=>worst[3],
            "passed"=>maximum_error<=tolerance))
        report["passed"] &= maximum_error<=tolerance
    end
    mkpath(dirname(abspath(output)))
    open(output,"w") do io
        TOML.print(io,report;sorted=true)
    end
    report["passed"] || throw(ErrorException("cross-channel convergence failed; see $output"))
    report
end

if abspath(PROGRAM_FILE)==@__FILE__
    3<=length(ARGS)<=4 || error("Usage: compare_cross_channel.jl left.toml right.toml output.toml [tolerance]")
    compare_cross_channel(ARGS[1],ARGS[2],ARGS[3];tolerance=length(ARGS)==4 ? parse(Float64,ARGS[4]) : 1e-6)
end
