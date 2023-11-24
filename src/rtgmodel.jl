using Catalyst
using DifferentialEquations

"""
```julia
make_rtg_rn(S = DEFAULT_SIGNAL; name, proteinlevels=STRESSED)
```

Make a RTG signalling ModelingToolkit (mtk) system with a time-dependent damage signal function `S` and `proteinlevels` for Rtg1, Rtg2, Rtg3, Bmh, and Mks proteins.
"""
function make_rtg_rn(s)
    rn = @reaction_network begin
        (hill($s, ksV, ksD, 7), k2I), Rtg2I_c <--> Rtg2A_c
        (k2M, kn2M), Rtg2A_c + Mks <--> Rtg2Mks_c
        # (kBM, knBM), Bmh + Mks <--> BmhMks
        # k13I + mm(BmhMks, k13IV, k13ID), Rtg13A_c --> Rtg13I_c
        # (k3A_c, k3I_c), Rtg3I_c <--> Rtg3A_c
        # k3I_n, Rtg3A_n --> Rtg3I_n
        # (k13_c, kn13_c), Rtg1_c + Rtg3A_c <--> Rtg13A_c
        # (k13_c, kn13_c), Rtg1_c + Rtg3I_c <--> Rtg13I_c
        # (k13_n, kn13_n), Rtg1_n + Rtg3A_n <--> Rtg13A_n
        # (k13_n, kn13_n), Rtg1_n + Rtg3I_n <--> Rtg13I_n
        # (k1in, k1out), Rtg1_c <--> Rtg1_n
        # (k3inA, k3outA), Rtg3A_c <--> Rtg3A_n
        # (k3inI, k3outI), Rtg3I_c <--> Rtg3I_n
    end

    return rn
end

# test

import HomotopyContinuation
using RetroSignalModel
import RetroSignalModel as rs

params = load_parameters("solution_rtgM4.csv")[1]

#==
u0 = [
    :Rtg1_c=>rs.STRESSED[rs.ΣRtg1],
    :Rtg1_n=>0.0, :Rtg13A_c=>0.0, :Rtg13I_c=>0.0, :Rtg13A_n=>0.0, :Rtg13I_n=>0.0,
    :Rtg2I_c=>rs.STRESSED[rs.ΣRtg2], :Rtg2A_c=>0.0, :Rtg2Mks_c=>0.0,
    :Rtg3I_c=>rs.STRESSED[rs.ΣRtg3], :Rtg3A_c=>0.0, :Rtg3A_n=>0.0, :Rtg3I_n=>0.0,
    :Mks=>rs.STRESSED[rs.ΣMks], :BmhMks=>0.0,
    :Bmh=>rs.STRESSED[rs.ΣBmh]
]

==#

rn = make_rtg_rn(0)

u0 = [
    :Rtg2I_c=>rs.STRESSED[rs.ΣRtg2], :Rtg2A_c=>0.0,
    :Mks=>rs.STRESSED[rs.ΣMks], :Rtg2Mks_c=>0.0
]


prob = SteadyStateProblem(rn, u0, params)

sol = solve(prob, DynamicSS(Rodas5()))


hc_steady_states(rn, params; u0)

wilhelm_2009_model = @reaction_network begin
    k1, Y --> 2X
    k2, 2X --> X + Y
    k3, X + Y --> Y
    k4, X --> 0
end

ps = [:k1 => 8.0, :k2 => 2.0, :k3 => 1.0, :k4 => 1.5]

hc_steady_states(wilhelm_2009_model, ps)
