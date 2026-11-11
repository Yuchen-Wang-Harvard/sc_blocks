# Independent full-cutoff production/PBW comparison driver. Each family/cutoff
# is explicit, and every case is checkpointed with separate measured timings.
using SuperconformalBlocks
using Dates
using SHA
using TOML
include("run_foundations.jl")

"""Timing metadata without retaining a numerical return value in the report."""
function measurement_fields(measurement)
    result=Dict{String,Any}("seconds"=>measurement.time,"allocated_bytes"=>measurement.bytes,"gc_seconds"=>measurement.gctime)
    for field in (:compile_time,:recompile_time)
        hasproperty(measurement,field) && (result[string(field)]=getproperty(measurement,field))
    end
    result
end

"""Atomic checkpoint after each completed comparison, retaining partial and failed runs."""
function save_matrix_report(path,report)
    mkpath(dirname(abspath(path)))
    temporary=path*".tmp"
    open(temporary,"w") do io
        TOML.print(io,report;sorted=true)
    end
    mv(temporary,path;force=true)
end

"""Construct production and independent PBW modules from the same input parameters."""
function pbw_matrix_setup(family::String,primarymask::Int)
    tetra=family=="tetrahedron"
    graph=tetra ? tetrahedron_graph() : occursin("glasses",family) ? glasses_graph() : theta_graph()
    allns=family in ("ns_theta","ns_glasses")
    family in ("ns_theta","ns_glasses","nsrr_theta","mixed_glasses","double_r_glasses","tetrahedron") ||
        throw(ArgumentError("unknown matrix family $family"))
    ramond=allns ? (false,false,false) : family=="mixed_glasses" ? (false,false,true) :
        tetra ? (false,false,false,true,true,true) : (false,true,true)
    count=length(ramond)
    primarymask>=0 && primarymask<2^count || throw(ArgumentError("invalid intrinsic-parity mask"))
    parities=ntuple(e->(primarymask>>(e-1))&1,count)
    all(e->!ramond[e] || parities[e]==0,1:count) || throw(ArgumentError("intrinsic parity only applies to NS edges"))
    b=Complex{BigFloat}(7//5)
    momenta=map(Complex{BigFloat},tetra ? (11//23,13//29,17//31,19//37,23//41,29//43) : (11//23,13//29,17//31))
    nv=length(graph.slots)
    plan=allns ? NSBlockPlan(BigFloat,graph,b,momenta,zeros(Int,nv);primary_parities=parities) :
        RamondBlockPlan(BigFloat,graph,b,momenta,ramond,zeros(Int,nv),ones(Int,nv);primary_parities=parities)
    q=b+inv(b)
    central=Complex{BigFloat}(3//2)+3q*q
    modules=ntuple(e->ramond[e] ? SCAModule(BigFloat,:R,central;beta=momenta[e]/sqrt(big(2))) :
        SCAModule(BigFloat,:NS,central;h=q*q/8-momenta[e]^2/2,parity=parities[e]),count)
    vertices=[ramond[slots[2]] ? nsrr_vertex(ntuple(j->modules[slots[j]],3)) :
        ns_vertex(ntuple(j->modules[slots[j]],3)) for slots in graph.slots]
    rvertices=findall(slots->ramond[slots[2]],graph.slots)
    cases=Tuple{Vector{Int},Vector{Int}}[]
    for fmask in 0:2^nv-1,etamask in 0:2^length(rvertices)-1
        iseven(count_ones(fmask)) || continue
        forms=[(fmask>>(j-1))&1 for j in 1:nv]
        eta=ones(Int,nv)
        for (j,v) in enumerate(rvertices)
            eta[v]=isodd(etamask>>(j-1)) ? -1 : 1
        end
        push!(cases,(forms,eta))
    end
    plan,vertices,cases,ramond,parities
end

"""Independently sew every PBW multigrade and compare parity plus all fixed-spin coefficients."""
function compare_pbw_matrix(plan,vertices,answer,domain,forms,eta,tolerance)
    maximum_error=zero(BigFloat)
    worst=first(domain.indices)
    evaluations=0
    for level in domain.indices
        reference=pbw_coefficient(plan.graph,vertices,level,forms,eta)
        if plan isa NSBlockPlan
            actual=zeros(Complex{BigFloat},length(reference))
            mask=sum(mod(level[e]+plan.primary_parities[e],2)<<(e-1) for e in eachindex(level))
            actual[mask+1]=get(answer,level,zero(Complex{BigFloat}))
        else
            actual=get(answer,level,zeros(Complex{BigFloat},length(reference)))
        end
        for (a,b) in ((actual,reference),(walsh!(copy(actual)),walsh!(copy(reference))))
            discrepancy=maximum(abs,a-b)/max(big(1),maximum(abs,b))
            isfinite(discrepancy) || throw(ErrorException("nonfinite PBW comparison at $level"))
            if discrepancy>maximum_error
                maximum_error,worst=discrepancy,level
            end
            evaluations+=length(a)
        end
    end
    maximum_error<tolerance || throw(ErrorException("PBW mismatch $maximum_error at $worst (tolerance $tolerance)"))
    Dict{String,Any}("maximum_scaled_error"=>string(maximum_error),"worst_level2"=>collect(worst),
        "multigrades"=>length(domain.indices),"parity_and_spin_entries"=>evaluations)
end

"""Run one complete graph-family matrix at its explicitly requested cutoff and precision."""
function run_pbw_matrix(family::String,level::Int,output::String;primarymask::Int=0,bits::Int=192,
                        cache_products::Bool=false)
    level>=0 && bits>=64 || throw(ArgumentError("invalid cutoff/precision"))
    isfile(output) && throw(ArgumentError("refusing to overwrite PBW comparison evidence"))
    root=normpath(joinpath(@__DIR__,".."))
    report=Dict{String,Any}("family"=>family,"level"=>level,
        "cutoff_kind"=>family=="tetrahedron" ? "total" : "per_edge",
        "intrinsic_primary_mask"=>primarymask,"precision_bits"=>bits,"convention"=>CONVENTION,
        "timestamp_utc"=>string(now(UTC)),"julia"=>string(VERSION),"machine"=>Sys.MACHINE,
        "cpu"=>Sys.CPU_NAME,"threads"=>Threads.nthreads(),"source_sha256"=>source_hashes(root),
        "process_id"=>getpid(),
        "reuse_vertex_label_products"=>cache_products,
        "completed"=>false,"passed"=>false,"cases"=>Dict{String,Any}[])
    started=time_ns()
    try
        setprecision(BigFloat,bits) do
            plan,vertices,cases,ramond,parities=pbw_matrix_setup(family,primarymask)
            cache_products && !(plan isa RamondBlockPlan) && throw(ArgumentError("product reuse applies to the Ramond plan"))
            plan isa RamondBlockPlan && (plan.cache_products=cache_products)
            report["ramond_edges"]=collect(ramond)
            report["primary_parities"]=collect(parities)
            report["b"]=[string(real(plan.b)),string(imag(plan.b))]
            report["momenta"]=[[string(real(p)),string(imag(p))] for p in plan.momenta]
            report["expected_cases"]=length(cases)
            domain=family=="tetrahedron" ? total_domain(length(ramond),2level) : box_domain(ntuple(_->2level,length(ramond)))
            tolerance=BigFloat(2)^(-floor(Int,55bits/100))
            report["tolerance"]=string(tolerance)
            save_matrix_report(output,report)
            for (number,(forms,eta)) in enumerate(cases)
                if plan isa NSBlockPlan
                    plan.forms.=forms
                    foreach(empty!,plan.vertices)
                else
                    set_vertex_labels!(plan,forms,eta)
                end
                report["current_case"]=Dict("number"=>number,"forms"=>forms,"eta"=>eta,"stage"=>"production")
                save_matrix_report(output,report)
                production=@timed plan isa NSBlockPlan ? ns_block(plan,domain) : ramond_block(plan,domain)
                report["current_case"]["stage"]="independent_pbw"
                report["current_production"]=measurement_fields(production)
                save_matrix_report(output,report)
                reference=@timed compare_pbw_matrix(plan,vertices,production.value,domain,forms,eta,tolerance)
                record=Dict{String,Any}("forms"=>forms,"eta"=>eta,"passed"=>true,
                    "production"=>measurement_fields(production),"reference"=>measurement_fields(reference),
                    "comparison"=>reference.value)
                push!(report["cases"],record)
                if plan isa RamondBlockPlan
                    report["product_cache_hits"]=plan.product_cache_hits
                    report["cached_products"]=length(plan.products)
                    report["ccy_transitions"]=plan.ccy_transitions
                end
                delete!(report,"current_production")
                report["elapsed_seconds"]=(time_ns()-started)/1e9
                save_matrix_report(output,report)
                println("Case ",number,"/",length(cases)," forms=",forms," eta=",eta,
                    " production_s=",production.time," reference_s=",reference.time,
                    " error=",reference.value["maximum_scaled_error"])
                flush(stdout)
            end
        end
        report["completed"]=true
        report["passed"]=true
        delete!(report,"current_case")
    catch error
        report["failure"]=sprint(showerror,error)
        rethrow()
    finally
        report["elapsed_seconds"]=(time_ns()-started)/1e9
        save_matrix_report(output,report)
    end
    report
end

if abspath(PROGRAM_FILE)==@__FILE__
    3<=length(ARGS)<=6 || error("Usage: run_pbw_matrix.jl family level output.toml [intrinsic_primary_mask] [bits] [cache_products]")
    run_pbw_matrix(ARGS[1],parse(Int,ARGS[2]),ARGS[3];primarymask=length(ARGS)>=4 ? parse(Int,ARGS[4]) : 0,
        bits=length(ARGS)>=5 ? parse(Int,ARGS[5]) : 192,
        cache_products=length(ARGS)>=6 ? parse(Bool,ARGS[6]) : false)
end
