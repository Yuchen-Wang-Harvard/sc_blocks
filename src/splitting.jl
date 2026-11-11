"""Theta psi_r on an auxiliary Ramond state, with Theta outside the psi insertion."""
function theta_psi_action(module_::FermionModule{T},r2::Int,state::FermionState) where T
    module_.ramond || throw(ArgumentError("Theta psi requires a Ramond module"))
    out=Dict{FermionState,Complex{T}}()
    for (intermediate,value) in act_psi(module_,r2,state)
        target=FermionState(intermediate.modes,!intermediate.ground)
        out[target]=get(out,target,zero(Complex{T}))+value*
            parity_sign(length(intermediate.modes)+Int(intermediate.ground)+1)
    end
    out
end

"""Even Theta psi operator between possibly unequal oscillator grades."""
function theta_psi_apply(module_::FreeField{T},r2::Int,vector::Dict{UInt32,Complex{T}}) where T
    out=Dict{UInt32,Complex{T}}()
    for (id,value) in vector
        state=module_.states[id]
        for (auxiliary,coefficient) in theta_psi_action(module_.auxiliary,r2,state.auxiliary)
            ff_add!(module_,out,OscillatorState(state.bosons,state.fermions,auxiliary),value*coefficient)
        end
    end
    prune_zeros!(out)
end

"""The two Virasoro weights of the external psi insertion."""
function psi_weights(b::Complex{T}) where T
    (-(1+2b*b)/(2(1-b*b)),(b*b+2)/(2(1-b*b)))
end

"""Physical L_1 splitting recurrence with directly evaluated ground anchors."""
mutable struct SplittingPlan{T<:AbstractFloat}
    actions::BranchActions{T}
    values::Dict{NTuple{3,Int},Complex{T}}
    factors::Dict{Tuple{Int,Int,Int,Int,Tuple{Vararg{Int16}}},Complex{T}}
    alternative_pivots::Int
    recursive_values::Int
end
function SplittingPlan(actions::BranchActions{T}) where T
    actions.positive.ramond || throw(ArgumentError("splitting needs a Ramond action plan"))
    SplittingPlan{T}(actions,Dict(),Dict(),0,0)
end

"""Direct unnormalized primary splitting matrix element in the current full tensor metric."""
function direct_splitting(plan::SplittingPlan,out::Int,input::Int,alpha::Int)
    module_=plan.actions.positive
    bra,ket=ff_primary(module_,out,alpha),ff_primary(module_,input,alpha)
    r2=branch_level2(true,input)-branch_level2(true,out)
    ff_pairing(module_,bra,theta_psi_apply(module_,r2,ket))
end

"""Cached one-leg descendant reduction at the insertion vertex (out,psi,input)."""
function splitting_factor(plan::SplittingPlan,out::Int,input::Int,slot::Int,term::BranchActionTerm)
    module_=plan.actions.positive
    value=term.coefficient
    external=psi_weights(module_.b)
    for copy in 1:2
        word=copy==1 ? term.first : term.second
        isempty(word) && continue
        value*=get!(plan.factors,(out,input,slot,copy,word)) do
            weights=(branch_weight(module_.b,module_.p,out,copy),external[copy],
                     branch_weight(module_.b,module_.p,input,copy))
            single_leg_vertex(slot,word,weights)
        end
    end
    value
end

"""Alternative physical L_0 Ward identity when the first splitting pivot is unresolved."""
function splitting_zero_pivot(plan::SplittingPlan{T},out::Int,input::Int,alpha::Int) where T
    denominator,numerator,bound=zero(Complex{T}),zero(Complex{T}),zero(T)
    for slot in (1,3)
        source=slot==1 ? out : input
        for term in branch_action(plan.actions,source,alpha,0)
            changed_out,changed_in=slot==1 ? (term.label,input) : (out,term.label)
            abs(changed_out-changed_in)==2 || continue
            coefficient=(slot==1 ? 1 : -1)*splitting_factor(plan,changed_out,changed_in,slot,term)
            if term.label==source
                denominator+=coefficient
                bound+=abs(coefficient)
            else
                numerator-=coefficient*splitting_raw(plan,changed_out,changed_in,alpha)
            end
        end
    end
    denominator,numerator,bound
end

"""Unnormalized splitting coefficient from L_±1 support and the physical Ward recurrence."""
function splitting_raw(plan::SplittingPlan{T},out::Int,input::Int,alpha::Int) where T
    isodd(out) && isodd(input) && alpha in (0,1) || throw(ArgumentError("invalid splitting labels"))
    abs(out-input)==2 || return zero(Complex{T})
    get!(plan.values,(out,input,alpha)) do
        abs(out)==abs(input)==1 && return direct_splitting(plan,out,input,alpha)
        denominator,numerator,bound=zero(Complex{T}),zero(Complex{T}),zero(T)
        if abs(out)>abs(input)
            for term in branch_action(plan.actions,input,alpha,-1)
                if term.label==input
                    factor=splitting_factor(plan,out,input,3,term)
                    denominator+=factor
                    bound+=abs(factor)
                elseif abs(out-term.label)==2
                    numerator-=splitting_factor(plan,out,term.label,3,term)*splitting_raw(plan,out,term.label,alpha)
                end
            end
            for term in branch_action(plan.actions,out,alpha,1)
                numerator+=splitting_factor(plan,term.label,input,1,term)*splitting_raw(plan,term.label,input,alpha)
            end
        else
            for term in branch_action(plan.actions,out,alpha,-1)
                if term.label==out
                    factor=splitting_factor(plan,out,input,1,term)
                    denominator+=factor
                    bound+=abs(factor)
                elseif abs(term.label-input)==2
                    numerator-=splitting_factor(plan,term.label,input,1,term)*splitting_raw(plan,term.label,input,alpha)
                end
            end
            for term in branch_action(plan.actions,input,alpha,1)
                numerator+=splitting_factor(plan,out,term.label,3,term)*splitting_raw(plan,out,term.label,alpha)
            end
        end
        if abs(denominator)<=ward_pivot_tolerance(T)*bound
            denominator,numerator,bound=splitting_zero_pivot(plan,out,input,alpha)
            plan.alternative_pivots+=1
        end
        abs(denominator)>ward_pivot_tolerance(T)*bound ||
            throw(ErrorException("both splitting pivots failed for out=$out, input=$input, alpha=$alpha"))
        plan.recursive_values+=1
        numerator/denominator
    end
end
