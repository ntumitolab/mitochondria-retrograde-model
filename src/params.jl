# Parameter searching to meet the conditions
using DataFrames
using CSV
using Tables
using DifferentialEquations

"""
Scan for parameters that meet the boolean conditions
"""
function param_scan(Model=rtgMTK;
    trajectories=1000,
    batch_size=trajectories,
    proteinlevels=STRESSED,
    knockoutlevel=1e-6,
    datafile=joinpath(@__DIR__, "data", "boolean_table_RTG13.csv"),
    rollparams=() -> exp10(rand(-3:0.1:3)),
    rollhill=() -> rand(0.5:0.5:5.0),
    nuclearRatioThreshold=10,
    steadyStateSolver=DynamicSS(Rodas5()),
    ensembleSolver=EnsembleThreads()
)

    conds = Tables.rowtable(dropmissing(CSV.read(datafile, DataFrame), disallowmissing=true))

    @named sys = Model(ONE_SIGNAL; proteinlevels=proteinlevels)

    param2idx = Dict(k => i for (i, k) in enumerate(parameters(sys)))

    # select random param for a fresh problem
    # rerun for each conditions
    function prob_func(prob, i, repeat)
        cond_idx = repeat + 1
        params = prob.p
        params[param2idx[nRuns]] = cond_idx

        # Initalize parameters for the first run
        if cond_idx == 1
            for k in (ksV, ksD, k2I, k2M, kn2M, kBM, knBM, k13I, k13IV, k13ID, k3A_c, k3I_c, k3I_n, k13_c, kn13_c, k13_n, kn13_n, k1in, k1out, k3inA, k3outA, k3inI, k3outI)
                params[param2idx[k]] = rollparams()
            end
            params[param2idx[n_S]] = rollhill()
        end

        # Adjust params according to conditions
        cond = conds[cond_idx]

        # Knockout Rtg1
        if cond[:Rtg1] == 0
            params[param2idx[ΣRtg1]] = knockoutlevel
        end

        # Knockout Rtg2
        if cond[:Rtg2] == 0
            params[param2idx[ΣRtg2]] = knockoutlevel
        end

        # Knockout Rtg3
        if cond[:Rtg3] == 0
            params[param2idx[ΣRtg3]] = knockoutlevel
        end

        # Knockout Mks
        if cond[:Mks] == 0
            params[param2idx[ΣMks]] = knockoutlevel
        end

        # Turn off signal
        if cond[:s] == 0
            params[param2idx[mul_S]] = 0
        end

        remake(prob, p=params)
    end

    function pass_cond(idx, conds, sol)
        cond = conds[idx]
        if cond[:gfp] == "rtg3"
            if cond[:Trans2Nuc] == 1
                rerun = rtg3_nucleus(sol) > nuclearRatioThreshold * rtg3_cytosol(sol)
            else
                rerun = nuclearRatioThreshold * rtg3_nucleus(sol) < rtg3_cytosol(sol)
            end
        elseif cond[:gfp] == "rtg1"
            if cond[:Trans2Nuc] == 1
                rerun = rtg1_nucleus(sol) > nuclearRatioThreshold * rtg1_cytosol(sol)
            else
                rerun = nuclearRatioThreshold * rtg1_nucleus(sol) < rtg1_cytosol(sol)
            end
        end
        return rerun
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
                push!(u, sol)
            end
        end

        return (u, false)
    end

    prob = EnsembleProblem(SteadyStateProblem(sys, resting_u0(sys), []); output_func, prob_func, reduction)

    sim = solve(prob, steadyStateSolver, ensembleSolver; trajectories, batch_size)
end
