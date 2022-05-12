# RTG model(s)
using ModelingToolkit

"""
```julia
hil(x, k, n=1)
```

Hill function

```math
hil(x, k, n) := \frac{x^n}{x^n + k^n}
```
"""
hil(x, k) = x / (x + k)
hil(x, k, n) = hil(x^n, k^n)

"""Constant damage signal"""
CONSTANT_SIGNAL(s) = t -> s

"""Damage signal = 0"""
ZERO_SIGNAL = CONSTANT_SIGNAL(0)

"""Damage signal = 1"""
ONE_SIGNAL = CONSTANT_SIGNAL(1)

# Stressed protein levels
const STRESSED = Dict(
    ΣRtg1 => 17.2068,
    ΣRtg2 => 221.9025,
    ΣRtg3 => 25.7976,
    ΣMks => 47.8059,
    ΣBmh => 1611.5946
)

# Unstressed protein levels
const UNSTRESSED = Dict(
    ΣRtg1 => 8.0878,
    ΣRtg2 => 80.7769,
    ΣRtg3 => 50.7607,
    ΣMks => 33.158,
    ΣBmh => 1239.5666
)

"""
```julia
rtgMTK(S = DEFAULT_SIGNAL; name, proteinlevels=STRESSED)
```

Make a RTG signalling ModelingToolkit (mtk) system with a time-dependent damage signal function `S` and `proteinlevels` for Rtg1, Rtg2, Rtg3, Bmh, and Mks proteins.
"""
function RtgMTK(S=ZERO_SIGNAL; name,
    proteinlevels=STRESSED, simplify=true)

    D = Differential(t)
    @variables v[1:13](t)

    # MTK does not support elimination in overdetermined systems
    # ODEs are eliminated manually
    eqs = [
        s ~ S(t),
        # Rtg2I_c <=> Rtg2A_c
        v[1] ~ ksV * hil(s * mul_S, ksD, n_S) * Rtg2I_c - k2I * Rtg2A_c,
        # Rtg2A_c + Mks <=> Rtg2Mks_c
        v[2] ~ k2M * Rtg2A_c * Mks - kn2M * Rtg2Mks_c,
        # Bmh + Mks <=> BmhMks
        v[3] ~ kBM * Bmh * Mks - knBM * BmhMks,
        # Rtg13A_c => Rtg13I_c
        v[4] ~ (k13I + k13IV * hil(BmhMks, k13ID)) * Rtg13A_c,
        # Rtg3I_c <=> Rtg3A_c
        v[5] ~ k3A_c * Rtg3I_c - k3I_c * Rtg3A_c,
        # Rtg3A_n => Rtg3I_n
        v[6] ~ k3I_n * Rtg3A_n,
        # Rtg1_c + Rtg3A_c <=> Rtg13A_c
        v[7] ~ k13_c * Rtg1_c * Rtg3A_c - kn13_c * Rtg13A_c,
        # Rtg1_c + Rtg3I_c <=> Rtg13I_c
        v[8] ~ k13_c * Rtg1_c * Rtg3I_c - kn13_c * Rtg13I_c,
        # Rtg1_n + Rtg3A_n <=> Rtg13A_n
        v[9] ~ k13_n * Rtg1_n * Rtg3A_n - kn13_n * Rtg13A_n,
        # Rtg1_n + Rtg3I_n <=> Rtg13I_n
        v[10] ~ k13_n * Rtg1_n * Rtg3I_n - kn13_n * Rtg13I_n,
        # Rtg1_c <=> Rtg1_n
        v[11] ~ k1in * Rtg1_c - k1out * Rtg1_n,
        # Rtg3A_c <=> Rtg3A_n
        v[12] ~ k3inA * Rtg3A_c - k3outA * Rtg3A_n,
        # Rtg3I_c <=> Rtg3I_n
        v[13] ~ k3inI * Rtg3I_c - k3outI * Rtg3I_n,

        # D(Rtg2I_c) ~ -v[1],
        D(Rtg2A_c) ~ v[1] - v[2],
        # D(Mks) ~ -v[2] - v[3],
        D(Rtg2Mks_c) ~ v[2],
        # D(Bmh) ~ -v[3],
        D(BmhMks) ~ v[3],
        D(Rtg13A_c) ~ -v[4] + v[7],
        D(Rtg13I_c) ~ v[4] + v[8],
        # D(Rtg3I_c) ~ -v[5] - v[8] - v[13],
        D(Rtg3A_c) ~ v[5] - v[7] - v[12],
        D(Rtg3A_n) ~ -v[6] - v[9] + v[12],
        D(Rtg3I_n) ~ v[6] - v[10] + v[13],
        # D(Rtg1_c) ~ -v[7] - v[8] - v[11],
        D(Rtg1_n) ~ -v[9] - v[10] + v[11],
        D(Rtg13A_n) ~ v[9],
        D(Rtg13I_n) ~ v[10],

        # Conservation relationships
        ΣRtg1 ~ Rtg1_c + Rtg1_n + Rtg13A_c + Rtg13I_c + Rtg13A_n + Rtg13I_n,
        ΣRtg2 ~ Rtg2I_c + Rtg2A_c + Rtg2Mks_c,
        ΣRtg3 ~ Rtg13A_c + Rtg13I_c + Rtg13A_n + Rtg13I_n + Rtg3I_c + Rtg3A_c + Rtg3A_n + Rtg3I_n,
        ΣMks ~ Mks + Rtg2Mks_c + BmhMks,
        ΣBmh ~ Bmh + BmhMks,
    ]

    sys = ODESystem(eqs, t; name, defaults=proteinlevels)

    if simplify
        sys = structural_simplify(sys)
    end
    return sys
end

"""
Make resting initial conditions fot the system. Default to all zeros.
"""
resting_u0(sys, val=0.0) = Dict(states(sys) .=> val)
