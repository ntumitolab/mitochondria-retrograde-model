import RetroSignalModel as rs
using DifferentialEquations

# Construct the model
rtgM4 = rs.rtgM4(1)

# Reaction system
rn = rtgM4.model

# Reference initial conditions
u0 = rtgM4.u
statemap = Dict(k => i for (i,) in enumerate(keys(u0)))

# System parameters
params = rtgM4.p

# Conservation relationships
Σbmh = u0.Bmh + u0.BmhMks
ΣMks = u0.Mks + u0.BmhMks + u0.Rtg2Mks_c
ΣRtg1 = u0.Rtg13_a_c + u0.Rtg13_i_c + u0.Rtg1_c + u0.Rtg1_n + u0.Rtg13_a_n + u0.Rtg13_i_n
ΣRtg2 = u0.Rtg2Mks_c + u0.Rtg2_ina_c + u0.Rtg2_act_c
ΣRtg3 = u0.Rtg13_a_c + u0.Rtg13_i_c + u0.Rtg13_a_n + u0.Rtg13_i_n + u0.Rtg3_i_c + u0.Rtg3_a_c + u0.Rtg3_a_n + u0.Rtg3_i_n

"""
Create a set of initial conditions respecting conservation relationships for the ensemble

    `rand_func`: a function that returns a valeu between 0 and 1
"""
function make_u0!(u0, Σbmh, ΣMks, ΣRtg1, ΣRtg2, ΣRtg3; rand_func=rand)

    # Randomly make n fractions of total
    function _frac(n::Int, total=1)
        total .* diff([0.0; sort([rand_func() for _ in 1:(n-1)]); 1.0])
    end

    # Split ΣMks into 3
    frac = _frac(3, ΣMks)
    u0[statemap[:Mks]] = frac[1]
    u0[statemap[:BmhMks]] = frac[2]
    u0[statemap[:Rtg2Mks_c]] = frac[3]
    u0[statemap[:Bmh]] = Σbmh - u0[statemap[:BmhMks]]

    remainder = ΣRtg2 - u0[statemap[:Rtg2Mks_c]]
    u0[statemap[:Rtg2_ina_c]] = remainder * rand_func()
    u0[statemap[:Rtg2_act_c]] = remainder - u0[statemap[:Rtg2_ina_c]]

    # Split ΣRtg1 into 6
    frac = _frac(6, ΣRtg1)

    u0[statemap[:Rtg13_a_c]] = frac[1]
    u0[statemap[:Rtg13_i_c]] = frac[2]
    u0[statemap[:Rtg1_c]] = frac[3]
    u0[statemap[:Rtg1_n]] = frac[4]
    u0[statemap[:Rtg13_a_n]] = frac[5]
    u0[statemap[:Rtg13_i_n]] = frac[6]

    # Split ΣRtg3 into 4
    remainder = ΣRtg3 - (u0[statemap[:Rtg13_a_c]] + u0[statemap[:Rtg13_i_c]] + u0[statemap[:Rtg13_a_n]] + u0[statemap[:Rtg13_i_n]])
    frac = _frac(4, remainder)

    u0[statemap[:Rtg3_i_c]] = frac[1]
    u0[statemap[:Rtg3_a_c]] = frac[2]
    u0[statemap[:Rtg3_a_n]] = frac[3]
    u0[statemap[:Rtg3_i_n]] = frac[4]

    return u0
end

function prob_func(prob, i, repeat)
    make_u0!(prob.u0, Σbmh, ΣMks, ΣRtg1, ΣRtg2, ΣRtg3)
    prob
end

"""Reject invalid (NaNs and negatives) and duplicate results"""
function reduction(u, batch, I; negtol=-1e-6)
    for result in batch
        # Skip invalid (NaNs and negatives) results
        if any(isnan.(result)) || any(result .< negtol)
            continue
        end

        result .= max.(0, result)
        isUnique = true

        # Only save unique solutions
        for existing in u
            if all(isapprox.(existing, result, rtol=1e-4))
                isUnique = false
                break
            end
        end

        if isUnique
            u = push!(u, result)
        end
    end
    return (u, false)
end

## Ensemble simulation
oprob = ODEProblem(rtgM4.model, rtgM4.u, (0.0, 1000.0), rtgM4.p)
ssprob = SteadyStateProblem(oprob)
method = SSRootfind()
ensprob = EnsembleProblem(ssprob, prob_func=prob_func, reduction=reduction)
trajectories = 120000
batch_size = 300

@time sim = solve(ensprob, method, EnsembleThreads(); trajectories=trajectories, batch_size=batch_size)

#=
**TODO**

Steady states at different levels of input signal `S`. (Use a loop).
=#
