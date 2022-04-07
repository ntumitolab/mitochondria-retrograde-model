# Parameter searching to meet the conditions
using DataFrames
using CSV
using Tables
using DifferentialEquations

"""
Scan for parameters that meet the boolean conditions in the retrograde (RTG) signalling model.
"""
function param_scan(
    Model=rtgMTK;
    trajectories=1000,
    batch_size=trajectories,
    proteinlevels=STRESSED,
    knockoutlevel=1e-6,
    datafile=joinpath(@__DIR__, "data", "boolean_table_RTG13.csv"),
    rollparams=() -> exp10(rand(-3:0.1:3)),
    rollhill=() -> rand(0.5:0.5:5.0),
    nuclearRatioThreshold=10,
    steadyStateSolver=DynamicSS(Rodas5()),
    ensembleSolver=EnsembleThreads(),
    ntarget=100
)

    conds = Tables.rowtable(dropmissing(CSV.read(datafile, DataFrame), disallowmissing=true))

    @named sys = Model(ONE_SIGNAL; proteinlevels=proteinlevels)

    param2idx = Dict(k => i for (i, k) in enumerate(parameters(sys)))
    idxnRuns = param2idx[nRuns]

    prob = SteadyStateProblem(sys, resting_u0(sys))

    # Initalize parameters for the first batch
    batchParams = exp10.(rand(-3:0.1:3, batch_size, length(parameters(sys))))
    batchParams[:, param2idx[n_S]] .= rand(0.5:0.5:5.0, batch_size, 1)


    # select random param for a fresh problem
    # repeat for each conditions
    function prob_func(prob, i, iter)
        params = prob.p
        params .= batchParams[i%batch_size+1, :]
        params[idxnRuns] = iter

        # Adjust params according to conditions
        cond = conds[iter]
        params[param2idx[ΣRtg1]] = ifelse(cond[:Rtg1] == 0, knockoutlevel, proteinlevels[ΣRtg1])
        params[param2idx[ΣRtg2]] = ifelse(cond[:Rtg2] == 0, knockoutlevel, proteinlevels[ΣRtg2])
        params[param2idx[ΣRtg3]] = ifelse(cond[:Rtg3] == 0, knockoutlevel, proteinlevels[ΣRtg3])
        params[param2idx[ΣMks]] = ifelse(cond[:Mks] == 0, knockoutlevel, proteinlevels[ΣMks])
        params[param2idx[mul_S]] = cond[:s]

        remake(prob, p=params)
    end

    function pass_cond(idx, conds, sol)
        cond = conds[idx]
        if cond[:gfp] == "rtg3"
            if cond[:Trans2Nuc] == 1
                passed = rtg3_nucleus(sol) > nuclearRatioThreshold * rtg3_cytosol(sol)
            else
                passed = nuclearRatioThreshold * rtg3_nucleus(sol) < rtg3_cytosol(sol)
            end
        elseif cond[:gfp] == "rtg1"
            if cond[:Trans2Nuc] == 1
                passed = rtg1_nucleus(sol) > nuclearRatioThreshold * rtg1_cytosol(sol)
            else
                passed = nuclearRatioThreshold * rtg1_nucleus(sol) < rtg1_cytosol(sol)
            end
        end
        return passed
    end

    # rerun when passed a condition and number for runs < number of conditions
    # otherwise, do not rerun
    function output_func(sol, i)
        idx = Int(sol.prob.p[param2idx[nRuns]])
        rerun = idx < length(conds) && pass_cond(idx, conds, sol)
        return (sol, rerun)
    end

    # Only add valid, all tests apassed results
    function reduction(u, batch, I)

        for sol in batch
            idx = Int(sol.prob.p[param2idx[nRuns]])
            if idx == length(conds) && pass_cond(idx, conds, sol)
                u = push!(u, sol)
            end
        end

        # Repopulate new parameters for the batch
        batchParams = exp10.(rand(-3:0.1:3, batch_size, length(parameters(sys))))
        batchParams[:, param2idx[n_S]] .= rand(0.5:0.5:5.0, batch_size, 1)

        return (u, length(u) >= ntarget)
    end

    ensprob = EnsembleProblem(prob; output_func, prob_func, reduction)

    sim = solve(ensprob, steadyStateSolver, ensembleSolver; trajectories, batch_size)
end
