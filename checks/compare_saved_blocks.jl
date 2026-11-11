# Regression comparison for two standalone CLI outputs in the same convention.
# This supplements independent PBW checks; it is not a substitute for them.
using TOML
using SHA
using Dates

"""Read saved parity coefficients at higher precision and compare every grade/component."""
function compare_saved_blocks(reference::String,candidate::String,output::String;
                              tolerance::String="1e-28",allow_precision_difference::Bool=false)
    a,b=TOML.parsefile(reference),TOML.parsefile(candidate)
    input_a,input_b=deepcopy(a["input"]),deepcopy(b["input"])
    if allow_precision_difference
        pop!(input_a,"precision_bits",nothing)
        pop!(input_b,"precision_bits",nothing)
    end
    input_a==input_b || throw(ArgumentError("different block inputs"))
    a["convention"]==b["convention"] || throw(ArgumentError("different conventions"))
    report=Dict{String,Any}("scope"=>"coefficient_regression_not_independent_PBW",
        "reference"=>reference,"candidate"=>candidate,"tolerance"=>tolerance,
        "reference_precision_bits"=>a["precision_bits"],"candidate_precision_bits"=>b["precision_bits"],
        "reference_sha256"=>bytes2hex(sha256(read(reference))),
        "candidate_sha256"=>bytes2hex(sha256(read(candidate))),
        "checker_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "timestamp_utc"=>string(now(UTC)),"passed"=>false)
    setprecision(BigFloat,256) do
        rows=Dict(Tuple(row["level2"])=>row for row in a["coefficients"])
        Set(keys(rows))==Set(Tuple(row["level2"]) for row in b["coefficients"]) ||
            throw(ArgumentError("different coefficient domains"))
        largest=zero(BigFloat)
        worst=Int[]
        entries=0
        for row in b["coefficients"]
            original=rows[Tuple(row["level2"])]
            x=complex.(parse.(BigFloat,original["parity_real"]),parse.(BigFloat,original["parity_imaginary"]))
            y=complex.(parse.(BigFloat,row["parity_real"]),parse.(BigFloat,row["parity_imaginary"]))
            error=maximum(abs,x-y)/max(big(1),maximum(abs,x))
            isfinite(error) || throw(ErrorException("nonfinite comparison"))
            if error>largest
                largest,worst=error,row["level2"]
            end
            entries+=length(x)
        end
        report["maximum_scaled_error"]=string(largest)
        report["worst_level2"]=worst
        report["compared_entries"]=entries
        report["passed"]=largest<parse(BigFloat,tolerance)
    end
    open(output,"w") do io
        TOML.print(io,report;sorted=true)
    end
    report["passed"] || error("saved-coefficient regression: "*report["maximum_scaled_error"])
    println("Compared ",report["compared_entries"]," entries; maximum scaled error=",report["maximum_scaled_error"])
    report
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==3 || error("Usage: compare_saved_blocks.jl reference.toml candidate.toml report.toml")
    compare_saved_blocks(ARGS...)
end
