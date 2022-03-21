using JuMP
using GLPK
model = Model(GLPK.Optimizer)
n = 5
JuMP.@variables(model, begin
    y[1:n] >= 0
end)
rtg1_idx = [1,2,4]
rtg2_idx = [3]
rtg3_idx = [5]




a = [@constraint(model,sum(y[rtg1_idx])==12),
     @constraint(model,sum(y[rtg2_idx])==12),
     @constraint(model,sum(y[rtg3_idx])==12)
]

set_name(a[1], "efefef")
set_name(a[2], "efef")
set_name(a[3], "efefe")


optimize!(model)

model.obj_dict
value.(y)