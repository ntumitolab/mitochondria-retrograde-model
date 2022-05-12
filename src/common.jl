using ModelingToolkit
using DataFrames
using CSV
using Tables

# State variable and parameters
@variables t
@variables Rtg13A_c(t) Rtg13I_c(t) Rtg13A_n(t) Rtg13I_n(t) Rtg3I_c(t) Rtg3A_c(t) Rtg3A_n(t) Rtg3I_n(t) Rtg1_c(t) Rtg1_n(t)
@variables s(t) Rtg2I_c(t) Rtg2A_c(t)
@variables Mks(t) Rtg2Mks_c(t) Bmh(t) BmhMks(t)
@parameters ΣRtg1 ΣRtg2 ΣRtg3 ΣMks ΣBmh
@parameters ksV = 1.0 ksD = 1.0 n_S = 1.0 k2I = 1.0 mul_S = 1.0
@parameters k2M = 1.0 kn2M = 1.0 kBM = 1.0 knBM = 1.0
@parameters k13I = 1.0 k13IV = 1.0 k13ID = 1.0 k3A_c = 1.0 k3I_c = 1.0 k3I_n = 1.0
@parameters k13_c = 1.0 kn13_c = 1.0 k13_n = 1.0 kn13_n = 1.0
@parameters k1in = 1.0 k1out = 1.0 k3inA = 1.0 k3outA = 1.0 k3inI = 1.0 k3outI = 1.0

# Distributions of RTG1 and RTG3 proteins
rtg13_nucleus(u) = u[Rtg13I_n] + u[Rtg13A_n]
rtg13_cytosol(u) = u[Rtg13I_c] + u[Rtg13A_c]
rtg3_nucleus(u) = rtg13_nucleus(u) + u[Rtg3A_n] + u[Rtg3I_n]
rtg3_cytosol(u) = rtg13_cytosol(u) + u[Rtg3A_c] + u[Rtg3I_c]
rtg1_nucleus(u) = rtg13_nucleus(u) + u[Rtg1_n]
rtg1_cytosol(u) = rtg13_cytosol(u) + u[Rtg1_c]

# Mapping from symbols to MTK parameters
const symbol2mtk = Dict(Symbol(k) => k for k in (ΣRtg1, ΣRtg2, ΣRtg3, ΣMks, ΣBmh, ksV, ksD, n_S, k2I, mul_S, k2M, kn2M, kBM, knBM, k13I, k13IV, k13ID, k3A_c, k3I_c, k3I_n, k13_c, kn13_c, k13_n, kn13_n, k1in, k1out, k3inA, k3outA, k3inI, k3outI))

# Utility functions

"""Loads a CSV file and returns a RowTable (An tuple of NamedTuple for each row)"""
load_data(filename) = Tables.rowtable(dropmissing(CSV.read(filename, DataFrame), disallowmissing=true))

"""Load boolean conditions for model validation."""
load_conditions(filename="boolean_table_RTG13.csv") = load_data(joinpath(@__DIR__, "data", filename))

"""Load model parameters"""
function load_parameters(filename="solution_rtgMTK_optim.csv")
    paramsets = load_data(joinpath(@__DIR__, "data", filename))
    map(paramsets) do p
        Dict(symbol2mtk[k] => i for (k, i) in pairs(p))
    end
end

