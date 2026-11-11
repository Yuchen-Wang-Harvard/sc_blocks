# Fresh Appendix C integrals. No archived blocks, periods or fitted scales enter.
include(joinpath(@__DIR__,"..","bin","blocks.jl"))
include("run_foundations.jl")
using LinearAlgebra

"""Encode complex numerical values as real/imaginary pairs for portable TOML."""
cross_encode(z::Number)=[Float64(real(z)),Float64(imag(z))]
cross_encode(values::AbstractArray)=[cross_encode(value) for value in values]
cross_encode(values::Tuple)=[cross_encode(value) for value in values]

"""Atomically retain progress, including partial integrals and numerical failures."""
function save_cross_report(output,report)
    mkpath(dirname(abspath(output)))
    temporary=output*".tmp"
    open(temporary,"w") do io
        TOML.print(io,report;sorted=true)
    end
    mv(temporary,output;force=true)
end

"""Recompute every period, refine printed seeds, and verify both exact spin transports."""
function prepare_cross_surfaces(config)
    plan=ThetaGeometryPlan(Float64,Int(config["geometry_mode_cutoff"]))
    comparison=ThetaGeometryPlan(Float64,plan.cutoff+4)
    cutoff=Int(config["theta_cutoff"])
    surfaces=NamedTuple[]
    metadata=Dict{String,Any}[]
    for row in config["surface"]
        source=Tuple(ComplexF64(parse_complex(z)) for z in row["source_q"])
        seed=Tuple(ComplexF64(parse_complex(z)) for z in row["target_q"])
        refined=refine_theta_surface(plan,source,seed)
        source_spins=((1,-1,1),(-1,-1,1))
        target_spins=Tuple(Tuple(Int.(s)) for s in row["target_spins"])
        length(target_spins)==2 || throw(ArgumentError("two spin structures required per surface"))
        for (s,t) in zip(source_spins,target_spins)
            alpha,beta=theta_characteristic(true,s)
            transport_characteristic(refined.matrix,alpha,beta)==theta_characteristic(false,t) ||
                throw(ErrorException("spin transport failed for surface $(row["id"]), spin $s"))
        end
        higher_source=theta_geometry(comparison,source)
        higher_target=theta_geometry(comparison,refined.target_q)
        geometry_error=max(maximum(abs,higher_source.omega-refined.source.omega),
            maximum(abs,higher_target.omega-refined.target.omega),
            abs(higher_source.log_abs_P-refined.source.log_abs_P),
            abs(higher_target.log_abs_P-refined.target.log_abs_P))
        geometry_error<1e-11 || throw(ErrorException("geometry cutoff is insufficient: $geometry_error"))
        frame=ntuple(j->free_partition(refined.target,false,target_spins[j];theta_cutoff=cutoff)/
            free_partition(refined.source,true,source_spins[j];theta_cutoff=cutoff),2)
        all(x->isfinite(x)&&x>0,frame) || throw(ErrorException("invalid free-theory frame"))
        push!(surfaces,(id=string(row["id"]),source_q=source,target_q=refined.target_q,
            source_spins=source_spins,target_spins=target_spins,free_ratio=frame))
        push!(metadata,Dict{String,Any}("id"=>string(row["id"]),"source_q"=>cross_encode(source),
            "target_seed"=>cross_encode(seed),"refined_target_q"=>cross_encode(refined.target_q),
            "symplectic_matrix"=>[collect(refined.matrix[j,:]) for j in 1:4],
            "initial_period_residual"=>refined.initial_period_residual,
            "final_period_residual"=>refined.final_period_residual,"iterations"=>refined.iterations,
            "mode_cutoff"=>plan.cutoff,"geometry_cutoff_change"=>geometry_error,
            "source_spins"=>[collect(s) for s in source_spins],"target_spins"=>[collect(s) for s in target_spins],
            "free_target_over_source"=>collect(frame),"spin_transport_passed"=>true))
    end
    allunique(surface.id for surface in surfaces) || throw(ArgumentError("surface identifiers must be unique"))
    surfaces,metadata
end

"""Precompute analytic plumbing powers and spin characters at the working block precision."""
function cross_evaluation_plan(surfaces,ramond::Bool,domain,levels)
    keys=domain.indices
    n=length(surfaces)
    powers=zeros(Complex{BigFloat},length(keys),n)
    signs=zeros(Int,ramond ? 8 : length(keys),n,2)
    for (s,surface) in enumerate(surfaces)
        q=ramond ? surface.source_q : surface.target_q
        spins=ramond ? surface.source_spins : surface.target_spins
        logarithms=map(x->log(Complex{BigFloat}(x)),q)
        for (k,key) in enumerate(keys)
            powers[k,s]=exp(sum(logarithms[e]*key[e]/2 for e in 1:3))
            if !ramond
                for spin in 1:2
                    signs[k,s,spin]=prod(spins[spin][e]^mod(key[e],2) for e in 1:3)
                end
            end
        end
        if ramond
            for mask in 0:7,spin in 1:2
                signs[mask+1,s,spin]=prod(spins[spin][e]^((mask>>(e-1))&1) for e in 1:3)
            end
        end
    end
    (keys=keys,powers=powers,signs=signs,levels=levels,ramond=ramond)
end

"""Evaluate a fresh coefficient series at all surfaces/spins and stated total cutoffs.

An all-NS c-recursion series contains both a sectors; the form filter separates
them before taking absolute squares. Ramond parity coefficients are summed
against the spin character before taking absolute squares.
"""
function evaluate_cross_series(series,plan;form::Int=0)
    result=zeros(Complex{BigFloat},size(plan.powers,2),length(plan.levels),2)
    for (k,key) in enumerate(plan.keys)
        grade=sum(key)
        !plan.ramond && mod(grade,2)!=form && continue
        haskey(series,key) || continue
        value=series[key]
        for s in axes(result,1),spin in 1:2
            coefficient=if plan.ramond
                sum(value[j]*plan.signs[j,s,spin] for j in 1:8)
            else
                value*plan.signs[k,s,spin]
            end
            term=coefficient*plan.powers[k,s]
            for l in eachindex(plan.levels)
                grade<=plan.levels[l] && (result[s,l,spin]+=term)
            end
        end
    end
    Float64.(abs2.(result))
end

"""Reference Gaussian nodes and exact measure compensation for dp_1 dp_2 dp_3/pi^3."""
function cross_momentum_grid(surfaces,ramond::Bool,order::Int,mode::String;
                             reference_q=nothing,alpha_scheme::String="sector")
    order>=1 || throw(ArgumentError("positive quadrature order required"))
    rates=[-log(abs((ramond ? surface.source_q : surface.target_q)[e])) for surface in surfaces,e in 1:3]
    reference=vec(minimum(rates;dims=1))
    if mode=="fixed"
        isnothing(reference_q) && throw(ArgumentError("fixed Gaussian requires three reference q magnitudes"))
        length(reference_q)==3 && all(q->isfinite(q)&&0<q<1,reference_q) ||
            throw(ArgumentError("reference q magnitudes must lie strictly between zero and one"))
        reference=-log.(Float64.(reference_q))
    elseif mode=="geometric"
        reference=sqrt.(reference.*vec(maximum(rates;dims=1)))
    elseif mode!="minimum"
        throw(ArgumentError("reference_gaussian must be fixed, minimum or geometric"))
    end
    alpha_scheme in ("sector","uniform_minus_half") || throw(ArgumentError("invalid Laguerre alpha scheme"))
    nodes,weights=Vector{Float64}[],Vector{Float64}[]
    alphas=Rational{Int}[]
    for e in 1:3
        alpha=alpha_scheme=="uniform_minus_half" || (ramond && e>1) ? -1//2 : 1//2
        push!(alphas,alpha)
        x,w=gauss_laguerre(Float64,order,alpha)
        p=sqrt.(x./reference[e])
        push!(nodes,p)
        push!(weights,w./(2reference[e]^(alpha+1))./p.^(2alpha+1)./pi)
    end
    (rates=rates,reference=reference,nodes=nodes,weights=weights,alphas=alphas)
end

"""Validate independent source/target grid orders, retaining the legacy common-order input."""
function cross_order_schedule(config)
    if haskey(config,"source_orders") || haskey(config,"target_orders")
        haskey(config,"orders") && throw(ArgumentError("specify common orders or separate channel orders, not both"))
        haskey(config,"source_orders") && haskey(config,"target_orders") ||
            throw(ArgumentError("both source_orders and target_orders are required"))
        source,target=Int.(config["source_orders"]),Int.(config["target_orders"])
    else
        source=Int.(config["orders"])
        target=copy(source)
    end
    length(source)==length(target) && !isempty(source) || throw(ArgumentError("channel order lists must have equal positive lengths"))
    all(orders->all(>(0),orders)&&issorted(orders)&&allunique(orders),(source,target)) ||
        throw(ArgumentError("quadrature orders must be positive and strictly increasing"))
    collect(zip(source,target))
end

"""Count successive stable changes; no preceding grid or a failed change resets the count."""
function cross_stability_steps(change::Real,tolerance::Real,previous_steps::Int)
    isfinite(tolerance) && tolerance>0 && previous_steps>=0 || throw(ArgumentError("invalid convergence controls"))
    isfinite(change) && 0<=change<=tolerance ? previous_steps+1 : 0
end

"""Pack channel integrals by surface, cutoff, and spin; no array-axis convention is implicit."""
function cross_integral_records(surfaces,levels,values,control)
    [Dict{String,Any}("id"=>surfaces[s].id,"level2"=>levels[l],
        "spin_index"=>spin,"value"=>cross_encode(values[s,l,spin]),
        "incorrect_same_label_control"=>cross_encode(control[s,l,spin]))
        for s in eachindex(surfaces) for l in eachindex(levels) for spin in 1:2]
end

"""Integrate one channel with fresh momentum plans and shared surfaces/spin evaluations.

Only combinatorial and momentum-local data are reused. The optional NS target
reference is independent SCA recursion; source blocks always use production
double Virasoro and the split-edge algorithm where prescribed.
"""
function integrate_cross_channel(config,surfaces,ramond::Bool,order::Int;
                                 progress=(state)->nothing,point_indices=nothing)
    BLAS.set_num_threads(1)
    bits=Int(config["block_precision_bits"])
    cutoff=2Int(config[ramond ? "source_total_level" : "target_total_level"])
    cutoff>=0 && bits>=64 || throw(ArgumentError("invalid cutoff/precision"))
    levels=sort!(unique([0,max(0,cutoff-2),cutoff]))
    grid=cross_momentum_grid(surfaces,ramond,order,String(config["reference_gaussian"]);
        reference_q=get(config,ramond ? "source_reference_q" : "target_reference_q",nothing),
        alpha_scheme=String(get(config,"laguerre_alpha_scheme","sector")))
    points=isnothing(point_indices) ? vec(collect(Iterators.product(1:order,1:order,1:order))) : collect(point_indices)
    !isempty(points) && allunique(points) && all(key->all(i->1<=i<=order,key),points) ||
        throw(ArgumentError("invalid or repeated momentum-grid indices"))
    values=zeros(ComplexF64,length(surfaces),length(levels),2)
    control=zeros(ComplexF64,size(values))
    source=ramond ? "NSRR" : "NSNSNS"
    started=time_ns()
    largest_recovery=0.0
    block_seconds=0.0
    block_compile_seconds=0.0
    constants_seconds=0.0
    product_cache_hits=0
    ccy_transitions=0
    setprecision(BigFloat,bits) do
        b=parse_real(config["b"])
        Q=b+inv(b)
        graph=theta_graph()
        domain=total_domain(3,cutoff)
        evaluation=cross_evaluation_plan(surfaces,ramond,domain,levels)
        method=String(config["target_method"])
        method in ("ns_reference","double_virasoro") || throw(ArgumentError("invalid target_method"))
        vacuum=!ramond && method=="ns_reference" ? ns_schottky(graph,domain) : nothing
        constant_bits=Int(get(config,"constant_precision_bits",53))
        constant_bits==53 || constant_bits>=64 || throw(ArgumentError("constant precision must be 53 or at least 64 bits"))
        constants=if constant_bits==53
            LiouvillePlan(Float64,Float64(b);order=Int(config["upsilon_order"]),tmax=config["upsilon_tmax"])
        else
            setprecision(BigFloat,constant_bits) do
                LiouvillePlan(BigFloat,parse_real(config["b"]);order=Int(config["upsilon_order"]),tmax=config["upsilon_tmax"])
            end
        end
        primary=[exp(sum(log(abs((ramond ? surface.source_q : surface.target_q)[e]))*
            (Float64(Q)^2/4+(ramond && e>1 ? 1/8 : 0)) for e in 1:3)) for surface in surfaces]
        for (number,indices) in enumerate(points)
            p=ntuple(e->grid.nodes[e][indices[e]],3)
            measure=prod(grid.weights[e][indices[e]] for e in 1:3)
            surface_weights=[measure*primary[s]*exp(-sum((grid.rates[s,e]-grid.reference[e])*p[e]^2 for e in 1:3))
                for s in eachindex(surfaces)]
            momenta=map(x->Complex{BigFloat}(0,BigFloat(x)),p)
            coupling_time=@timed setprecision(BigFloat,constant_bits==53 ? bits : constant_bits) do
                ramond ? liouville_nsrr(constants,p) : liouville_ns(constants,p)
            end
            couplings=coupling_time.value
            constants_seconds+=coupling_time.time
            if ramond
                plan=RamondBlockPlan(BigFloat,graph,b,momenta,(false,true,true),[0,0],[1,1];
                    cache_products=Bool(get(config,"reuse_vertex_label_products",true)))
                for f in 0:1,eta1 in 1:2,eta2 in 1:2
                    set_vertex_labels!(plan,[f,f],[eta1==1 ? 1 : -1,eta2==1 ? 1 : -1])
                    computation=@timed ramond_block(plan,domain)
                    block_seconds+=computation.time
                    block_compile_seconds+=hasproperty(computation,:compile_time) ? computation.compile_time : 0.0
                    plan.recovery_diagnostics["passed"] || throw(ErrorException("unvalidated source recovery"))
                    largest_recovery=max(largest_recovery,Float64(parse(BigFloat,plan.recovery_diagnostics["maximum_sector_residual"])))
                    squared=evaluate_cross_series(computation.value,evaluation)
                    coupling=(-1)^f*couplings[f+1,eta1]*couplings[f+1,eta2]
                    incorrect=couplings[1,eta1]*couplings[1,eta2]
                    for s in eachindex(surfaces),l in eachindex(levels),spin in 1:2
                        term=surface_weights[s]*squared[s,l,spin]
                        values[s,l,spin]+=coupling*term
                        control[s,l,spin]+=incorrect*term
                    end
                end
                product_cache_hits+=plan.product_cache_hits
                ccy_transitions+=plan.ccy_transitions
            else
                reference=if method=="ns_reference"
                    weights=ntuple(e->Q*Q/8+BigFloat(p[e])^2/2,3)
                    plan=NSRecursionPlan(BigFloat,graph,BigFloat(3)/2+3Q*Q,weights)
                    computation=@timed ns_c_block(plan,domain;vacuum=vacuum)
                    block_seconds+=computation.time
                    block_compile_seconds+=hasproperty(computation,:compile_time) ? computation.compile_time : 0.0
                    computation.value
                else
                    nothing
                end
                plan=method=="double_virasoro" ? NSBlockPlan(BigFloat,graph,b,momenta,[0,0]) : nothing
                for a in 0:1
                    series=if isnothing(reference)
                        plan.forms.=(a,a)
                        foreach(empty!,plan.vertices)
                        computation=@timed ns_block(plan,domain)
                        block_seconds+=computation.time
                        block_compile_seconds+=hasproperty(computation,:compile_time) ? computation.compile_time : 0.0
                        computation.value
                    else
                        reference
                    end
                    squared=evaluate_cross_series(series,evaluation;form=a)
                    coupling=(-1)^a*couplings[a+1]^2
                    for s in eachindex(surfaces),l in eachindex(levels),spin in 1:2
                        values[s,l,spin]+=surface_weights[s]*coupling*squared[s,l,spin]
                    end
                end
            end
            all(isfinite,values) || throw(ErrorException("nonfinite partition integral"))
            if number%max(1,div(length(points),10))==0 || number==length(points)
                state=Dict{String,Any}("channel"=>source,"order"=>order,"nodes_done"=>number,
                    "nodes_total"=>length(points),"elapsed_seconds"=>(time_ns()-started)/1e9,
                    "integrals"=>cross_integral_records(surfaces,levels,values,control))
                progress(state)
                println("Cross channel=",source," order=",order," nodes=",number,"/",length(points),
                    " seconds=",state["elapsed_seconds"])
                flush(stdout)
            end
        end
    end
    result=Dict{String,Any}("channel"=>source,"order"=>order,"nodes"=>length(points),
        "level2_cutoffs"=>levels,"reference_rates"=>grid.reference,
        "laguerre_alphas"=>Float64.(grid.alphas),
        "reference_gaussian"=>config["reference_gaussian"],"block_seconds"=>block_seconds,
        "block_compile_seconds"=>block_compile_seconds,"constant_seconds"=>constants_seconds,
        "elapsed_seconds"=>(time_ns()-started)/1e9,"maximum_recovery_residual"=>largest_recovery,
        "product_cache_hits"=>product_cache_hits,"ccy_transitions"=>ccy_transitions,
        "integrals"=>cross_integral_records(surfaces,levels,values,control))
    (report=result,values=values,control=control,levels=levels)
end

"""Reconstruct one freshly completed worker result, checking every labelled array entry."""
function read_cross_piece(path,surfaces)
    report=TOML.parsefile(path)
    get(report,"completed",false) || throw(ErrorException("incomplete momentum worker: $path"))
    levels=Int.(report["level2_cutoffs"])
    indices=Dict(surface.id=>i for (i,surface) in enumerate(surfaces))
    level_indices=Dict(level=>i for (i,level) in enumerate(levels))
    values=zeros(ComplexF64,length(surfaces),length(levels),2)
    control=zeros(ComplexF64,size(values))
    filled=falses(size(values))
    for row in report["integrals"]
        s,l,j=indices[row["id"]],level_indices[row["level2"]],Int(row["spin_index"])
        filled[s,l,j] && throw(ErrorException("duplicate worker integral entry"))
        values[s,l,j]=ComplexF64(row["value"]...)
        control[s,l,j]=ComplexF64(row["incorrect_same_label_control"]...)
        filled[s,l,j]=true
    end
    all(filled) && all(isfinite,values) && all(isfinite,control) || throw(ErrorException("incomplete or nonfinite worker integrals"))
    (report=report,values=values,control=control,levels=levels)
end

"""Distribute disjoint momentum nodes to isolated processes, retaining every spin at each node.

The reduction order is deterministic. Worker-local progress files preserve
partial integrals; their caches and MPFR precision contexts are never shared.
Serial execution calls the same integration function on the complete grid.
"""
function dispatch_cross_channel(config,surfaces,ramond::Bool,order::Int,workers::Int,output;
                                progress=(state)->nothing)
    workers==1 && return integrate_cross_channel(config,surfaces,ramond,order;progress=progress)
    points=vec(collect(Iterators.product(1:order,1:order,1:order)))
    count=min(workers,length(points))
    chunks=[points[j:count:end] for j in 1:count]
    prefix=output*"."*(ramond ? "NSRR" : "NSNSNS")*".order$(order)"
    started=time_ns()
    root=normpath(joinpath(@__DIR__,".."))
    jobs=NamedTuple[]
    pieces=NamedTuple[]
    try
        for (chunk,indices) in enumerate(chunks)
            input=abspath(prefix*".chunk$(chunk).input.toml")
            result=abspath(prefix*".chunk$(chunk).toml")
            logfile=abspath(prefix*".chunk$(chunk).log")
            any(isfile,(input,result,logfile)) && throw(ArgumentError("refusing to overwrite a worker checkpoint"))
            payload=Dict{String,Any}("config"=>config,"ramond"=>ramond,"order"=>order,
                "point_indices"=>[collect(key) for key in indices],
                "surfaces"=>[Dict("id"=>s.id,"source_q"=>cross_encode(s.source_q),"target_q"=>cross_encode(s.target_q),
                    "source_spins"=>[collect(x) for x in s.source_spins],"target_spins"=>[collect(x) for x in s.target_spins],
                    "free_ratio"=>collect(s.free_ratio)) for s in surfaces])
            save_cross_report(input,payload)
            io=open(logfile,"w")
            command=`$(Base.julia_cmd()) --project=$root --threads=1 $(joinpath(@__DIR__,"cross_worker.jl")) $input $result`
            process=try
                run(pipeline(command;stdout=io,stderr=io);wait=false)
            catch
                close(io)
                rethrow()
            end
            push!(jobs,(process=process,io=io,result=result,logfile=logfile))
        end
        last_progress=time()
        while any(job->process_running(job.process),jobs)
            for job in jobs
                !process_running(job.process) && !success(job.process) &&
                    throw(ErrorException("momentum worker failed; see $(job.logfile)"))
            end
            if time()-last_progress>=5
                done=sum(isfile(job.result) ? get(TOML.parsefile(job.result),"nodes_done",0) : 0 for job in jobs)
                progress(Dict{String,Any}("channel"=>(ramond ? "NSRR" : "NSNSNS"),"order"=>order,
                    "nodes_done"=>done,"nodes_total"=>order^3,"elapsed_seconds"=>(time_ns()-started)/1e9))
                last_progress=time()
            end
            sleep(0.25)
        end
        for job in jobs
            success(job.process) || throw(ErrorException("momentum worker failed; see $(job.logfile)"))
            push!(pieces,read_cross_piece(job.result,surfaces))
        end
    finally
        for job in jobs
            process_running(job.process) && kill(job.process)
            wait(job.process)
            close(job.io)
        end
    end
    values=reduce(+,piece.values for piece in pieces)
    control=reduce(+,piece.control for piece in pieces)
    levels=first(pieces).levels
    all(piece->piece.levels==levels,pieces) || throw(ErrorException("worker cutoff mismatch"))
    report=copy(first(pieces).report)
    report["nodes"]=sum(piece.report["nodes"] for piece in pieces)
    report["nodes"]==order^3 || throw(ErrorException("incomplete momentum grid"))
    report["workers"]=length(pieces)
    report["elapsed_seconds"]=(time_ns()-started)/1e9
    report["worker_seconds"]=[piece.report["elapsed_seconds"] for piece in pieces]
    for key in ("block_seconds","block_compile_seconds","constant_seconds","product_cache_hits","ccy_transitions")
        report[key]=sum(piece.report[key] for piece in pieces)
    end
    report["component_timing_scope"]="sum over isolated worker processes; elapsed_seconds is wall time"
    report["maximum_recovery_residual"]=maximum(piece.report["maximum_recovery_residual"] for piece in pieces)
    report["integrals"]=cross_integral_records(surfaces,levels,values,control)
    progress(Dict{String,Any}("channel"=>report["channel"],"order"=>order,
        "nodes_done"=>report["nodes"],"nodes_total"=>order^3,"elapsed_seconds"=>report["elapsed_seconds"]))
    (report=report,values=values,control=control,levels=levels)
end

"""Run all supplied surfaces/orders, recording convergence without fitting a normalization.

A completed run is not automatically a release pass. The explicit ratio and
quadrature tolerances must pass, and separate Gaussian/precision comparisons
are required by the release checklist even when these local criteria pass.
"""
function run_cross_channel(input::String,output::String)
    isfile(output) && throw(ArgumentError("refusing to overwrite cross-channel evidence"))
    config=TOML.parsefile(input)
    root=normpath(joinpath(@__DIR__,".."))
    report=Dict{String,Any}("completed"=>false,"passed"=>false,"release_validated"=>false,
        "convention"=>CONVENTION,"input"=>config,"source_sha256"=>source_hashes(root),
        "timestamp_utc"=>string(now(UTC)),"julia"=>string(VERSION),"cpu"=>Sys.CPU_NAME,
        "threads"=>Threads.nthreads(),"momentum_measure"=>"product dp_e/pi",
        "arithmetic"=>"BigFloat blocks and series summation; explicit constant precision; Float64 geometry and final quadrature",
        "constant_precision_bits"=>Int(get(config,"constant_precision_bits",53)),
        "reuse_vertex_label_products"=>Bool(get(config,"reuse_vertex_label_products",true)),
        "independence"=>"fresh production source, independent NS target when selected, fresh constants and geometry; no archived numerical input",
        "runs"=>Dict{String,Any}[])
    started=time_ns()
    try
        surfaces,metadata=prepare_cross_surfaces(config)
        report["geometry"]=metadata
        b=Float64(parse_real(config["b"]))
        exponent=1+2(b+inv(b))^2
        frame=[surface.free_ratio[spin]^exponent for surface in surfaces,spin in 1:2]
        previous=nothing
        previous_source=nothing
        previous_target=nothing
        ratio_tolerance=Float64(get(config,"ratio_tolerance",1e-4))
        quadrature_tolerance=Float64(get(config,"quadrature_tolerance",1e-6))
        schedule=cross_order_schedule(config)
        required_steps=Int(get(config,"required_successive_steps",1))
        required_steps>=1 || throw(ArgumentError("at least one successive convergence step is required"))
        all(t->isfinite(t)&&t>0,(ratio_tolerance,quadrature_tolerance)) || throw(ArgumentError("positive finite tolerances required"))
        stable_steps=0
        report["ratio_tolerance"]=ratio_tolerance
        report["quadrature_tolerance"]=quadrature_tolerance
        report["required_successive_steps"]=required_steps
        workers=Int(get(config,"workers",1))
        workers>=1 || throw(ArgumentError("workers must be positive"))
        report["workers"]=workers
        for (source_order,target_order) in schedule
            entry=Dict{String,Any}("source_order"=>source_order,"target_order"=>target_order)
            source_order==target_order && (entry["order"]=source_order)
            report["current_run"]=entry
            save_cross_report(output,report)
            progress=state->begin
                entry["progress"]=state
                save_cross_report(output,report)
            end
            source=dispatch_cross_channel(config,surfaces,true,source_order,workers,output;progress=progress)
            entry["source"]=source.report
            save_cross_report(output,report)
            target=dispatch_cross_channel(config,surfaces,false,target_order,workers,output;progress=progress)
            entry["target"]=target.report
            delete!(entry,"progress")
            source_final=source.values[:,end,:]
            target_final=target.values[:,end,:]
            all(z->real(z)>0 && abs(imag(z))<1e-10*real(z),source_final) || throw(ErrorException("source partition function is not positive real"))
            all(z->real(z)>0 && abs(imag(z))<1e-10*real(z),target_final) || throw(ErrorException("target partition function is not positive real"))
            ratios=source_final./target_final.*frame
            controls=source.control[:,end,:]./target_final.*frame
            entry["ratios"]=[Dict("id"=>surfaces[s].id,"spin_index"=>spin,
                "ratio"=>cross_encode(ratios[s,spin]),"distance_from_one"=>abs(ratios[s,spin]-1),
                "incorrect_same_label_control"=>cross_encode(controls[s,spin])) for s in eachindex(surfaces) for spin in 1:2]
            entry["maximum_ratio_distance_from_one"]=maximum(abs,ratios.-1)
            entry["maximum_source_cutoff_change"]=maximum(abs,source_final./source.values[:,max(1,end-1),:].-1)
            entry["maximum_target_cutoff_change"]=maximum(abs,target_final./target.values[:,max(1,end-1),:].-1)
            convergence=Inf
            if !isnothing(previous)
                entry["maximum_ratio_change_from_previous_order"]=maximum(abs,ratios-previous)
                entry["maximum_source_relative_change"]=maximum(abs,source_final./previous_source.-1)
                entry["maximum_target_relative_change"]=maximum(abs,target_final./previous_target.-1)
                convergence=max(entry["maximum_ratio_change_from_previous_order"],
                    entry["maximum_source_relative_change"],entry["maximum_target_relative_change"])
            end
            stable_steps=cross_stability_steps(convergence,quadrature_tolerance,stable_steps)
            entry["successive_stable_steps"]=stable_steps
            entry["quadrature_converged"]=stable_steps>=required_steps
            entry["ratio_tolerance_passed"]=entry["maximum_ratio_distance_from_one"]<=ratio_tolerance
            push!(report["runs"],entry)
            delete!(report,"current_run")
            save_cross_report(output,report)
            previous,previous_source,previous_target=ratios,source_final,target_final
            println("Cross source order=",source_order," target order=",target_order,
                " maximum |Q-1|=",entry["maximum_ratio_distance_from_one"],
                " quadrature change=",convergence," successive stable steps=",stable_steps)
            flush(stdout)
        end
        final=last(report["runs"])
        report["passed"]=final["quadrature_converged"] && final["ratio_tolerance_passed"]
        report["completed"]=true
    catch error
        report["failure"]=sprint(showerror,error)
        rethrow()
    finally
        report["elapsed_seconds"]=(time_ns()-started)/1e9
        save_cross_report(output,report)
    end
    report
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==2 || error("Usage: run_cross_channel.jl input.toml output.toml")
    run_cross_channel(ARGS[1],ARGS[2])
end
