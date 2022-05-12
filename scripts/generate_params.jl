using RetroSignalModel
using Optim
using CSV

optim_params(targetratio=5)

dfparams = optim_params_threads(120; targetratio=5)

CSV.write("solution_rtgMTK_optim_5x.csv", dfparams)
