import RetroSignalModel as rs

models = [
    rs.actM1(),
    rs.rtgM1(),
    rs.rtgM2(),
    rs.rtgM3(),
    rs.rtgM4(),
]

m4 = rs.rtgM4()
m4rand = rs.rtgM4(m4; u=rs.init_u(m4))
rs.get_steadysol(m4rand)

@show rs.get_protein_lookup(models[end].model)

rs.get_steadysol.(models)
[rs.knockout(m, [:Rtg1]) for m in models]

@test rs.is_valid(rs.rtgM4()) == true

valid, _ = rs.try_conditions(models[end])
@test valid == true
