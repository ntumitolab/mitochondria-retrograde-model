using RetroSignalModel
using ModelingToolkit
using DataFrames
using Optim
using CSV

versioninfo()

optim_params(targetratio=5)

res = optim_params_threads(120; targetratio=5);

paramaps = map(res) do x
    Dict(Symbol(k) => v for (k, v) in x.parammap)
end

CSV.write("solution_rtgMTK_optim_5x.csv", DataFrame(paramaps))
