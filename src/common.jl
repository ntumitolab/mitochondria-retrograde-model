using ModelingToolkit

@variables t
@variables Rtg13A_c(t) Rtg13I_c(t) Rtg13A_n(t) Rtg13I_n(t) Rtg3I_c(t) Rtg3A_c(t) Rtg3A_n(t) Rtg3I_n(t) Rtg1_c(t) Rtg1_n(t)
@variables s(t) Rtg2I_c(t) Rtg2A_c(t)
@variables Mks(t) Rtg2Mks_c(t) Bmh(t) BmhMks(t)

@parameters ΣRtg1 ΣRtg2 ΣRtg3 ΣMks ΣBmh
@parameters ksV = 1.0 ksD = 1.0 n_S = 1.0 k2I = 1.0 mul_S = 1.0 nRuns = 1.0
@parameters k2M = 1.0 kn2M = 1.0 kBM = 1.0 knBM = 1.0
@parameters k13I = 1.0 k13IV = 1.0 k13ID = 1.0 k3A_c = 1.0 k3I_c = 1.0 k3I_n = 1.0
@parameters k13_c = 1.0 kn13_c = 1.0 k13_n = 1.0 kn13_n = 1.0
@parameters k1in = 1.0 k1out = 1.0 k3inA = 1.0 k3outA = 1.0 k3inI = 1.0 k3outI = 1.0
