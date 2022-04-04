using LinearAlgebra
using DataFrames
using CSV

"""Read a CSV file and return a DataFrame"""
read_csv(filename; kwargs...) = DataFrame(CSV.File(filename); kwargs...)

"""Merge csv files under a directory (first layer)."""
function merge_csv(dirname; format=".csv")

    fns_ = readdir(dirname)

    fns = filter_files(fns_, format)

    length(fns) == 0 ? @warn("No file is found.") : nothing

    fns = joinpath.(dirname, fns)

    df = CSV.read(fns[1], DataFrame)
    for fn in fns[2:end]
        df_ = CSV.read(fn, DataFrame)
        df = vcat(df, df_)
    end

    return df
end


function catalyst_name(m::ReactionSystem; fcall=Catalyst.species, remove_t=true)
    spec_names_t = string.(fcall(m))
    spec_names = remove_t ? replace.(spec_names_t, "(t)" => "") : spec_names_t
    return spec_names
end

function get_tables(csv_fns::NamedTuple)
    dfs = []
    for fn in csv_fns
        push!(dfs, read_csv(fn))
    end
    return NamedTuple{keys(csv_fns)}(tuple(dfs...))
end

"""
Get constructor from a data type
"""
get_ctor(datatype) = typeof(datatype).name.wrapper


function find_id(sNames, features; convert_lower=true)
    sNames = convert_lower ? lowercase.(sNames) : sNames
    features = convert_lower ? lowercase.(features) : features
    idxs = Vector{Int}()
    for (ind, i) in enumerate(sNames)
        isEle = true
        for f in features
            if !occursin(f, string(i))
                isEle = false
                break
            end
        end
        if isEle
            push!(idxs, ind)
        end
    end
    return idxs
end

"""
Scale the maximum element of vector `v` to `capVal`.
"""
cap_vec(v::AbstractVector, capVal) = v .* (capVal / norm(v, Inf))

# Default Settings

const COND_DATA_PATH = joinpath(@__DIR__, "data", "boolean_table_RTG13.csv")
const DEL_CONC = 1e-4

"""
Steady state solver.

A Stiff solver is chosen `TRBDF2` (stand for Tranpezoidal Backward Differeital Formula) from **DifferentialEquations.jl** for robust simulation of stability.

`TRBF2` is a one-step method based on trapezoidal rule and the backward differentiaion formula of order 2, and it is strongly L-stable that means `TRBFS` is good at integrating stiff equations [1].

References
----------
1. Hosea, M. E., & Shampine, L. F. (1996). Analysis and implementation of TR-BDF2. Applied Numerical Mathematics, 20(1-2), 21-37. URL: https://www.sciencedirect.com/science/article/pii/0168927495001158?via%3Dihub
2. DifferentialEqiations.jl - Stiff Problems. [[link](https://diffeq.sciml.ai/stable/solvers/ode_solve/#Stiff-Problems)]
"""
const DEFAULT_SSMETHOD = DynamicSS(TRBDF2())

const DataFiles = (;
    RNAseq=joinpath(@__DIR__, "data", "RNAseq_RTG_expression.csv"),
    BoolCond=joinpath(@__DIR__, "data", "boolean_table_RTG13.csv"),
    solution_rtgM4=joinpath(@__DIR__, "data", "solution_rtgM4.csv")
)

const DEFAULT_SOL_IDX = 1

"""
Parameter searching
"""
const K_dist = Exponential(1000)
const K_N_dist = Geometric(0.3)
const S_SPAN = (1e-6, 1.0)
const Max_expr = 10000.0
const TRANS_THRESHOLD = 1.5
const DataTables = get_tables(DataFiles)
const DefaultCondition = "stressed"

function get_expression_levels(; df=DataTables.RNAseq, condition="unstressed", Max_expr=Max_expr, idxName=1)
    prNames = df[!, idxName]
    ex = df[!, condition]
    ex = cap_vec(ex, Max_expr)
    return NamedTuple{tuple(Symbol.(prNames)...)}(tuple(ex...))
end


function get_conditions(; df=DataTables.BoolCond, drop_missing=true)
    if drop_missing
        return dropmissing(df)
    else
        return df
    end
end
