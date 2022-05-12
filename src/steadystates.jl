using ModelingToolkit
using OrdinaryDiffEq
using SteadyStateDiffEq
using SciMLBase


"""
Find steady states of the ODE system `sys` by randomly
"""
function find_steady_states(;
    S=ZERO_SIGNAL, proteins=STRESSED, params=Dict(),
    trajectories=10, batch_size=trajectories, ntarget=trajectories,
    solver=SSRootfind(), ensemble_method=EnsembleThreads()
)

    @named sys = RtgMTK(S; proteinlevels=proteins)

    statesmap = Dict(k => i for (i, k) in enumerate(states(sys)))

    # Create a set of initial conditions respecting conservation relationships for the ensemble
    function make_u0()
        u0 = zeros(length(statesmap))

        # Randomly split total into `n` fractions
        function _frac(n::Int, total=1.0)
            total .* diff([zero(total); sort([rand(typeof(total)) for _ in 1:(n-1)]); one(total)])
        end

        # Split ΣMks into 3
        frac = _frac(3, proteins[ΣMks])
        # u0[statesmap[Mks]] = frac[1]
        u0[statesmap[BmhMks]] = frac[2]
        u0[statesmap[Rtg2Mks_c]] = frac[3]

        # Remainder of Bmh
        # u0[statesmap[Bmh]] = proteins[ΣBmh] - u0[statesmap[BmhMks]]

        # Split ΣRtg2 into 2
        frac = _frac(2, proteins[ΣRtg2] - u0[statesmap[Rtg2Mks_c]])
        # u0[statesmap[Rtg2I_c]] = frac[1]
        u0[statesmap[Rtg2A_c]] = frac[2]

        # Split ΣRtg1 into 6
        frac = _frac(6, proteins[ΣRtg1])
        u0[statesmap[Rtg13A_c]] = frac[1]
        u0[statesmap[Rtg13I_c]] = frac[2]
        # u0[statesmap[Rtg1_c]] = frac[3]
        u0[statesmap[Rtg1_n]] = frac[4]
        u0[statesmap[Rtg13A_n]] = frac[5]
        u0[statesmap[Rtg13I_n]] = frac[6]

        # Split ΣRtg3 not in Rtg13 complex into 4
        frac = _frac(4, proteins[ΣRtg3] - u0[statesmap[Rtg13A_c]] - u0[statesmap[Rtg13I_c]] - u0[statesmap[Rtg13A_n]] - u0[statesmap[Rtg13I_n]])

        # u0[statesmap[Rtg3I_c]] = frac[1]
        u0[statesmap[Rtg3A_c]] = frac[2]
        u0[statesmap[Rtg3A_n]] = frac[3]
        u0[statesmap[Rtg3I_n]] = frac[4]

        return u0
    end

    """Rerun when the solution is invalid (NaNs and negatives)"""
    function output_func(sol, i; negtol=-1e-6)
        rerun = any(isnan.(sol.u)) || any(sol.u .< negtol)
        clamp!(sol.u, 0.0, Inf)
        return (sol, rerun)
    end

    function prob_func(prob, i, repeat)
        # Roll a random initial condition
        remake(prob, u0=make_u0())
    end

    """Only save unique solutions"""
    function reduction(us, batch, I; norm_atol=1e-4)
        for res in batch
            isUnique = true
            for existing in us
                if isapprox(existing.u, res.u, atol=norm_atol)
                    isUnique = false
                    break
                end
            end

            if isUnique
                us = push!(us, res)
            end
        end
        return (us, length(us) >= ntarget)
    end

    prob = SteadyStateProblem(sys, resting_u0(sys), params)
    ensprob = EnsembleProblem(prob; output_func, prob_func, reduction)
    return solve(ensprob, solver, ensemble_method; trajectories, batch_size)
end
