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

"""Damage signal = 0"""
ZERO_SIGNAL(t) = 0
"""Damage signal = 1"""
ONE_SIGNAL(t) = 1


const STRESSED = Dict(
    ΣRtg1 => 17.2068,
    ΣRtg2 => 221.9025,
    ΣRtg3 => 25.7976,
    ΣMks => 47.8059,
    ΣBmh => 1611.5946
)

const UNSTRESSED = Dict(
    ΣRtg1 => 8.0878,
    ΣRtg2 => 80.7769,
    ΣRtg3 => 50.7607,
    ΣMks => 33.158,
    ΣBmh => 1239.5666
)

# Distributions of RTG1 and RTG3 proteins

rtg13_nucleus(sol) = sol[Rtg13I_n] + sol[Rtg13A_n]
rtg13_cytosol(sol) = sol[Rtg13I_c] + sol[Rtg13A_c]
rtg3_nucleus(sol) = rtg13_nucleus(sol) + sol[Rtg3A_n] + sol[Rtg3I_n]
rtg3_cytosol(sol) = rtg13_cytosol(sol) + sol[Rtg3A_c] + sol[Rtg3I_c]
rtg1_nucleus(sol) = rtg13_nucleus(sol) + sol[Rtg1_n]
rtg1_cytosol(sol) = rtg13_cytosol(sol) + sol[Rtg1_c]

"""
```julia
rtgMTK(S = DEFAULT_SIGNAL; name, proteinlevels=STRESSED)
```

Make a RTG signalling model with a time-dependent damage signal function `S`.
"""
function rtgMTK(S=ZERO_SIGNAL; name, proteinlevels=STRESSED, simplify=true)

    D = Differential(t)

    # Modulation Layer
    @variables v_01(t) v_02(t) v_03(t)

    # Nuclear Localization Signal (NLS) activation
    @variables v_04(t) v_05(t) v_06(t)

    # Rtg1Binding: Rtg13 formation
    @variables v_07(t) v_08(t) v_09(t) v_10(t)

    # Translocation
    @variables v_11(t) v_12(t) v_13(t)

    # MTK does not support elimination in overdetermined systems
    # ODEs are eliminated manually
    eqs = [
        s ~ S(t),
        # Rtg2I_c <=> Rtg2A_c
        v_01 ~ ksV * hil(s * mul_S, ksD, n_S) * Rtg2I_c - k2I * Rtg2A_c,
        # Rtg2A_c + Mks <=> Rtg2Mks_c
        v_02 ~ k2M * Rtg2A_c * Mks - kn2M * Rtg2Mks_c,
        # Bmh + Mks <=> BmhMks
        v_03 ~ kBM * Bmh * Mks - knBM * BmhMks,
        # Rtg13A_c => Rtg13I_c
        v_04 ~ (k13I + k13IV * hil(BmhMks, k13ID)) * Rtg13A_c,
        # Rtg3I_c <=> Rtg3A_c
        v_05 ~ k3A_c * Rtg3I_c - k3I_c * Rtg3A_c,
        # Rtg3A_n => Rtg3I_n
        v_06 ~ k3I_n * Rtg3A_n,
        # Rtg1_c + Rtg3A_c <=> Rtg13A_c
        v_07 ~ k13_c * Rtg1_c * Rtg3A_c - kn13_c * Rtg13A_c,
        # Rtg1_c + Rtg3I_c <=> Rtg13I_c
        v_08 ~ k13_c * Rtg1_c * Rtg3I_c - kn13_c * Rtg13I_c,
        # Rtg1_n + Rtg3A_n <=> Rtg13A_n
        v_09 ~ k13_n * Rtg1_n * Rtg3A_n - kn13_n * Rtg13A_n,
        # Rtg1_n + Rtg3I_n <=> Rtg13I_n
        v_10 ~ k13_n * Rtg1_n * Rtg3I_n - kn13_n * Rtg13I_n,
        # Rtg1_c <=> Rtg1_n
        v_11 ~ k1in * Rtg1_c - k1out * Rtg1_n,
        # Rtg3A_c <=> Rtg3A_n
        v_12 ~ k3inA * Rtg3A_c - k3outA * Rtg3A_n,
        # Rtg3I_c <=> Rtg3I_n
        v_13 ~ k3inI * Rtg3I_c - k3outI * Rtg3I_n,

        # D(Rtg2I_c) ~ -v_01,
        D(Rtg2A_c) ~ v_01 - v_02,
        # D(Mks) ~ -v_02 - v_03,
        D(Rtg2Mks_c) ~ v_02,
        # D(Bmh) ~ -v_03,
        D(BmhMks) ~ v_03,
        D(Rtg13A_c) ~ -v_04 + v_07,
        D(Rtg13I_c) ~ v_04 + v_08,
        # D(Rtg3I_c) ~ -v_05 - v_08 - v_13,
        D(Rtg3A_c) ~ v_05 - v_07 - v_12,
        D(Rtg3A_n) ~ -v_06 - v_09 + v_12,
        D(Rtg3I_n) ~ v_06 - v_10 + v_13,
        # D(Rtg1_c) ~ -v_07 - v_08 - v_11,
        D(Rtg1_n) ~ -v_09 - v_10 + v_11,
        D(Rtg13A_n) ~ v_09,
        D(Rtg13I_n) ~ v_10,

        # Conservation relationships
        ΣRtg1 ~ Rtg1_c + Rtg1_n + Rtg13A_c + Rtg13I_c + Rtg13A_n + Rtg13I_n,
        ΣRtg2 ~ Rtg2I_c + Rtg2A_c + Rtg2Mks_c,
        ΣRtg3 ~ Rtg13A_c + Rtg13I_c + Rtg13A_n + Rtg13I_n + Rtg3I_c + Rtg3A_c + Rtg3A_n + Rtg3I_n,
        ΣMks ~ Mks + Rtg2Mks_c + BmhMks,
        ΣBmh ~ Bmh + BmhMks,
        nRuns ~ nRuns
    ]

    sys = ODESystem(eqs, t; name, defaults=proteinlevels)

    if simplify
        sys = structural_simplify(sys)
    end
    return sys
end

"""
Resting initial conditions = all zeroes.
"""
resting_u0(sys) = zeros(length(states(sys)))
