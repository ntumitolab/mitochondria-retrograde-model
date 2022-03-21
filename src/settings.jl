const COND_DATA_PATH = joinpath(@__DIR__,"data","boolean_table_RTG13.csv")
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
const SSMETHOD = DynamicSS(TRBDF2())

const DataFiles = (;
    RNAseq= joinpath(@__DIR__,"data","RNAseq_RTG_expression.csv"),
    BoolCond = joinpath(@__DIR__,"data","boolean_table_RTG13.csv"),
    solution_rtgM4 = joinpath(@__DIR__,"data","solution_rtgM4.csv")
)

const default_sol_i = 1

"""
Parameter searching
"""
const K_dist = Exponential(1000)
const K_N_dist = Geometric(0.3) 
const S_SPAN = (1e-6,1.0)
const Max_expr = 10000.
const TRANS_THRESHOLD = 1.5
const DataTables = getTables(DataFiles)
const DefaultCondition = "stressed"