# Parameter searching to meet the conditions
using DifferentialEquations
using ModelingToolkit
using Optim

"""
Find parameters that satisfy the boolean conditions in the retrograde (RTG) signalling model using Optimization methods.
"""
function optim_params(
    Model=RtgMTK;
    datafile=joinpath(@__DIR__, "data", "boolean_table_RTG13.csv"),
    knockoutlevel=1e-4,
    proteinlevels=STRESSED,
    desolver=DynamicSS(Rodas5()),
    ensemble_method=EnsembleThreads(),
    lowerbound=1e-3,
    upperbound=1e3,
    hilllowerbound=2.0,
    hillupperbound=7.0,
    xinit=1.0,
    optimsolver=Optim.SAMIN(),
    optimoptions=Optim.Options(iterations=10^5, show_trace=true, show_every=1000),
    targetratio=2)

    # Boolean conditions
    conds = load_data(datafile)

    @named sys = Model(ONE_SIGNAL; proteinlevels)
    prob = SteadyStateProblem(sys, resting_u0(sys))

    # Indices to parameters in the system
    param2idx = Dict(k => i for (i, k) in enumerate(parameters(sys)))
    iΣRtg1 = param2idx[ΣRtg1]
    iΣRtg2 = param2idx[ΣRtg2]
    iΣRtg3 = param2idx[ΣRtg3]
    iΣMks = param2idx[ΣMks]
    imul_S = param2idx[mul_S]

    # Parameters to be optimized
    params_optim = [k for k in parameters(sys) if !any(isequal(k), (ΣRtg1, ΣRtg2, ΣRtg3, ΣMks, ΣBmh, mul_S))]
    # Mapping indices of the x vector in Optim to params in the system
    xi2params = [param2idx[k] for k in params_optim]
    xidxnS = findfirst(isequal(n_S), params_optim)

    # The cost is bounded by target protein concentration ratio
    scorecap = -log10(targetratio)

    # Cost function
    function cost(x)

        function output_func(sol, i)
            cond = conds[i]
            if cond[:gfp] == "rtg3"
                score = -log10(rtg3_nucleus(sol) / rtg3_cytosol(sol)) * ifelse(cond[:Trans2Nuc] == 1, 1, -1)
            elseif cond[:gfp] == "rtg1"
                score = -log10(rtg1_nucleus(sol) / rtg1_cytosol(sol)) * ifelse(cond[:Trans2Nuc] == 1, 1, -1)
            else
                score = 0.0
            end
            return (max(score - scorecap, 0.0), false)
        end

        function prob_func(prob, i, repeat)
            cond = conds[i]
            p = prob.p
            # Adjust params according to conditions
            p[iΣRtg1] = cond[:Rtg1] == 0 ? knockoutlevel : p[iΣRtg1]
            p[iΣRtg2] = cond[:Rtg2] == 0 ? knockoutlevel : p[iΣRtg2]
            p[iΣRtg3] = cond[:Rtg3] == 0 ? knockoutlevel : p[iΣRtg3]
            p[iΣMks] = cond[:Mks] == 0 ? knockoutlevel : p[iΣMks]
            p[imul_S] = cond[:s]

            # Assign optim vector to ODE parameters
            for i in 1:length(x)
                p[xi2params[i]] = x[i]
            end
            return remake(prob, p=p)
        end

        eprob = EnsembleProblem(prob; output_func, prob_func)
        esol = solve(eprob, desolver, ensemble_method, trajectories=length(conds))
        return sum(esol) / length(conds)
    end

    # Initial conditions and lower / upper bounds for Optim
    x0 = ones(length(xi2params)) .* xinit
    lb = ones(length(xi2params)) .* lowerbound
    ub = ones(length(xi2params)) .* upperbound
    ub[xidxnS] = hillupperbound
    lb[xidxnS] = hilllowerbound
    x0 .= clamp.(x0, lb, ub)

    res = Optim.optimize(cost, lb, ub, x0, optimsolver, optimoptions)
    parammap = Dict(params_optim .=> Optim.minimizer(res))

    return (res=res, parammap=parammap)
end
