# Validation

using CairoMakie
import RetroSignalModel as rs

rtgM4 = rs.rtgM4(1)
pr_table = rs.get_protein_lookup(rtgM4)

# Set expression level
explvls = rs.get_expression_levels()
# Knock out Bmh (should knockout Total Bmh)

expΔBmh = (exp..., Bmh=1e-4)

# New initial conditions
u_low = rs.init_u(rtgM4.model, pr_table; expLevels=exp_dBmh, init_s=0.2)

# Steady state solutions

sol = rs.get_steadysol(rtgM4, u_low)

# Show outputs

# Translocation

function measure_trans(m, p, dmg_sig, explvl, gfp)
    pr_table = rs.get_protein_lookup(m)
    u = rs.init_u(m.model, pr_table; expLevels=explvl, init_s=dmg_sig)
    prob = DEsteady(func=m.model, u0=u, p=p, method=rs.SSMETHOD)
    sol = solve(prob)
    return rs.check_nu_accumulation(sol, gfp, pr_table)
end

measure_trans(rtgM4, rtgM4.p, 0.4, exp, "Rtg3")

measure_trans(rtgM4, rtgM4.p, 0.1, exp, "Rtg3")
measure_trans(rtgM4, rtgM4.p, 0.9, (exp..., Bmh=1e-4), "Rtg3")

# Grouped table for plotting
tbl = (x=[1, 1, 1, 1, 2, 2, 2, 2],
    height=[
        measure_trans(rtgM4, rtgM4.p, 0.1, exp, "Rtg3").Conc_Cyt,
        measure_trans(rtgM4, rtgM4.p, 0.1, exp, "Rtg3").Conc_Nuc,
        measure_trans(rtgM4, rtgM4.p, 0.9, exp, "Rtg3").Conc_Cyt,
        measure_trans(rtgM4, rtgM4.p, 0.9, exp, "Rtg3").Conc_Nuc,
        measure_trans(rtgM4, rtgM4.p, 0.1, (exp..., Bmh=1), "Rtg3").Conc_Cyt,
        measure_trans(rtgM4, rtgM4.p, 0.1, (exp..., Bmh=1), "Rtg3").Conc_Nuc,
        measure_trans(rtgM4, rtgM4.p, 0.9, (exp..., Bmh=1), "Rtg3").Conc_Cyt,
        measure_trans(rtgM4, rtgM4.p, 0.9, (exp..., Bmh=1), "Rtg3").Conc_Nuc
    ],
    grp=[1, 1, 2, 2, 3, 3, 4, 4], # Genotype
    grp1=[1, 2, 1, 2, 1, 2, 1, 2] # Mitochondrial condition
)

# Makie plots
# Could be replaced by PyPlot?

fig = Figure()
colors = Makie.wong_colors()
ax = Axis(fig[1, 1], xticks=(1:2, ["Control (WT/mtDamage)", "Bmh-del (WT/mtDamage)"]),
    title="Rtg3 translocation")

barplot!(ax, tbl.x, tbl.height,
    dodge=tbl.grp,
    stack=tbl.grp1,
    color=colors[tbl.grp1]
)

# Legend
labels = ["Rtg3-GFP (Cytosol)", "Rtg3-GFP (Nucleus)"]
elements = [PolyElement(polycolor=colors[i]) for i in 1:length(labels)]
Legend(fig[1, 2], elements, labels, "Rtg3 Concentration")
# save("validation.pdf", fig)
fig
