using ModelingToolkit
using RetroSignalModel
using DifferentialEquations
using DataFrames
using Tables
using CSV

import RetroSignalModel as rs
import ModelingToolkit as mtk

sim = param_scan(ensembleSolver=EnsembleThreads(), trajectories=2^30, batch_size=2^20, ntarget=100)
