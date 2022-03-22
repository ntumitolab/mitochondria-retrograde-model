
#' Environment
using Pkg
Pkg.activate(@__DIR__)
cd(@__DIR__)
Pkg.status()
Pkg.instantiate()

#' Load library to main worker
@time using Distributed, RetroSignalModel, CSV
import RetroSignalModel as rs
#' Add workers
addprocs(exeflags="--project=$(Base.active_project())");
@show nprocs();

#' Load library to workers
@time @everywhere using RetroSignalModel;

#' Parameter searching
@time df = rs.paramSearching(rs.rtgM4(); num_sim=1e9, save_iter=1e6, distributed=true, saveall=false)

#'Output
CSV.write("../../src/data/solution_rtgM4.csv", df)
