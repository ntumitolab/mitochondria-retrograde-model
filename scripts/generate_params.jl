using RetroSignalModel
using Optim
using CSV

versioninfo()

optim_params(targetratio=2)

dfparams = optim_params_threads(200; targetratio=2)

CSV.write("solution_rtgMTK_optim_2x.csv", dfparams)
