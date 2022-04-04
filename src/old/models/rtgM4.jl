@with_kw struct rtgM4{A,B,C} <: RTGmodel
    model::A = rtgM4_model()
    protein_lookup = get_protein_lookup(model)
    u::B = get_u(DataFiles.solution_rtgM4, model, DEFAULT_SOL_IDX)
    p::C = get_p(DataFiles.solution_rtgM4, model, DEFAULT_SOL_IDX)
end

function rtgM4(sol_id::Integer; sol_csv=DataFiles.solution_rtgM4)
    model = rtgM4_model()
    u = get_u(sol_csv, model, sol_id)
    p = get_p(sol_csv, model, sol_id)
    return rtgM4(u=u, p=p)
end

"""
Modified version from rtgM3
1. Change input to Hill kinetics
"""
function rtgM4_model()
    model = @reaction_network begin
        #=
        Input
        =#
        (1), s → s

        #=
        Modulation Layer
        =#
        (hill(s, ksV, ksD, n_s), k2I), Rtg2_ina_c ↔ Rtg2_act_c #[1]
        (k2M, kn2M), Rtg2_act_c + Mks ↔ Rtg2Mks_c
        (kBM, knBM), Bmh + Mks ↔ BmhMks

        #=
        Activation Model:
        =#
        # NLSact: Nuclear Localization Signal (NLS) activation
        (k13I + mm(BmhMks, k13IV, k13ID)), Rtg13_a_c → Rtg13_i_c
        (k3A_c, k3I_c), Rtg3_i_c ↔ Rtg3_a_c
        (k3I_n), Rtg3_a_n → Rtg3_i_n

        # Rtg1Binding: Rtg13 formation
        (k13_c, kn13_c), Rtg1_c + Rtg3_a_c ↔ Rtg13_a_c
        (k13_c, kn13_c), Rtg1_c + Rtg3_i_c ↔ Rtg13_i_c

        (k13_n, kn13_n), Rtg1_n + Rtg3_a_n ↔ Rtg13_a_n #[3]
        (k13_n, kn13_n), Rtg1_n + Rtg3_i_n ↔ Rtg13_i_n #[3]

        # Translocation
        (k1in, k1out), Rtg1_c ↔ Rtg1_n #[2]
        (k3inA, k3outA), Rtg3_a_c ↔ Rtg3_a_n #[2]
        (k3inI, k3outI), Rtg3_i_c ↔ Rtg3_i_n

        #(k4, kn4), Rtg13_a_c ↔ Rtg13_a_n # Deleted [2]
    end n_s ksV ksD k2I k2M kn2M kBM knBM k13I k13IV k13ID k3A_c k3I_c k3I_n k13_c kn13_c k13_n kn13_n k1in k1out k3inA k3outA k3inI k3outI
    return model
end
