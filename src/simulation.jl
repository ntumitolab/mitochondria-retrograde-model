struct RTGoutput{T}
    IsNucAccum::Bool
    Conc_Cyt::T
    Conc_Nuc::T
    threshold::T
end

IsNucAccum(res::RTGoutput) = res.IsNucAccum

"""
Validate responses
"""
function isValid(m::RTGmodel; kwags...)
    valid, _ = try_conditions(m; kwags...)
    return valid
end

function getOutput(m::RTGmodel, gfp; kwags...)
    return getOutput(m.u, gfp, m.protein_lookup; kwags...)
end

"""
Compare the cytosolic concentration of either 'rtg1' or 'rtg3' with their nucleus concentrations.
Output:
1 : nucleus concentration is higher than the cytosolic one
0 : otherwise
"""
function getOutput(sol, gfp, protein_lookup; threshold=TRANS_THRESHOLD)
    gfp = lowercase(gfp)
    if gfp == "rtg1"
        cyt_index = protein_lookup[:Rtg1_c]
        nuc_index = protein_lookup[:Rtg1_n]
    elseif gfp == "rtg3"
        cyt_index = protein_lookup[:Rtg3_c]
        nuc_index = protein_lookup[:Rtg3_n]
    else
        throw(MathodError)
    end

    total_conc_cyt = sum(sol[cyt_index])
    total_conc_nuc = sum(sol[nuc_index])

    NucAccum = total_conc_nuc > threshold * total_conc_cyt ? true : false

    return RTGoutput(NucAccum, total_conc_cyt, total_conc_nuc, threshold)
end

"""
Knockout specified protein with name [`prName`](@ref), and set as low concentration equal to [`DEL_CONC`](@ref)
"""
function knockout(m::RTGmodel, prNames; del_conc=DEL_CONC, WildExpLevels::NamedTuple=getExpLevels(; condition=DefaultCondition), kwags...)

    newConcs = fill(del_conc, length(prNames))
    knockout_exp = setTuples(WildExpLevels, prNames, newConcs)
    u = init_u(m; expLevels=knockout_exp, kwags...)
    return u
end



"""
Solve steady state of [`RTGmodel`](@ref) based on given initial variables (`u`). Noted that fieldname `u` is ignored.
"""
function getSteadySol(m::RTGmodel, u; ssmethod=SSMETHOD)
    #todo
    prob = DEsteady(func=m.model, u0=u, p=m.p, method=ssmethod)
    sol = solve(prob)
    return sol
end

"""
Returm [`RTGmodel`](@ref) with `u` in steady state.
"""
function getSteady(m::RTGmodel; warning=true, kwags...)
    # consider
    sol_ss = getSteadySol(m, m.u; kwags...)
    model = get_ctor(m)
    m_ss = model(m, sol_ss.u, m.p, m.protein_lookup)

    # Warning
    warning && sol_ss.retcode == :Success ? @warn("Steady state not found with resid=$(sol_ss.resid)") : nothing

    return m_ss
end


"""
Try conditions of the model with given initial variables and parameters.
"""
function try_conditions(m::RTGmodel; expLevels=getExpLevels(; condition=DefaultCondition), S_SPAN=S_SPAN, del_conc=DEL_CONC, TRANS_THRESHOLD=TRANS_THRESHOLD, SSMETHOD=SSMETHOD, conditions=getConditions(; df=DataTables.BoolCond), Break=true, Cond_random=true, input_name="s")

    s_idx = findfirst(x -> x == input_name, catalyst_name(m.model))
    valid = true
    num_cond = size(conditions)[1]

    # Knockout Candidates
    iExs = findall(x -> x in Symbol.(names(conditions)), keys(expLevels))
    prNames = collect(keys(expLevels)[iExs])

    # Shuffling conditions
    trial_cond_order = Cond_random ? shuffle(1:num_cond) : 1:num_cond
    test_log = Vector{Union{Missing,Int16}}(missing, num_cond)

    # Get steady state
    df_pr = conditions[!, prNames]

    for tn in trial_cond_order
        conPr = df_pr[tn, :]
        cond = conditions[tn, :]

        # Get knockout protein names
        Pr_Exts = values(df_pr[tn, prNames])
        i_dels = findall(x -> x == 0, Pr_Exts)
        prName_dels = prNames[i_dels]

        # Knockout
        u_k = knockout(m, prName_dels; del_conc=del_conc)
        u_k[s_idx] = S_SPAN[cond[:s]+1]

        # Steady state
        u_ss = getSteadySol(m, u_k; ssmethod=SSMETHOD)

        # Status of RTG switch
        output = getOutput(u_ss, cond[:gfp], m.protein_lookup; threshold=TRANS_THRESHOLD)

        test_log[tn] = IsNucAccum(output)
        if IsNucAccum(output) == cond[:Trans2Nuc]
            continue
        else
            valid = false
            if Break
                break
            end
        end
    end
    return valid, test_log
end
