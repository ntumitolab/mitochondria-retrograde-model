# Parameter searching to meet the conditions
using DataFrames
using CSV
using Tables

function param_scan(Model=rtgMTK;
    proteinlevels=STRESSED,
    knockoutlevel=1e-6,
    datafile=joinpath(@__DIR__, "data", "boolean_table_RTG13.csv"))

    conds = Tables.rowtable(dropmissing(CSV.read(datafile, DataFrame), disallowmissing=true))

    @named sys = Model(ONE_SIGNAL; proteinlevels=proteinlevels)

    # prob_func()
    # select random param for a fresh problem
    # rerun for each conditions

    # output_func()
    # passed a condition and number for runs < number of conditions : rerun
    # failed a condition: write negative numbers as output and quit
    # passed all conditions: write params (or prob) as output and quit


    # reduction()
    # Only add results with valid params
end
