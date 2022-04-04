#' Load library to main worker
using CSV
import RetroSignalModel as rs

#' Parameter searching
@time df = rs.search_params(rs.rtgM4(); num_sim=1e9, save_iter=1e6, distributed=true, saveall=false)

#'Output
CSV.write("../../src/data/solution_rtgM4.csv", df)
