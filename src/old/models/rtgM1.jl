@with_kw struct rtgM1{A,B,C} <: RTGmodel
    model::A = rtgM1_model()
    protein_lookup = get_protein_lookup(model)
    p::C = init_p(model; idx_hill_coefs=[])
    u::B = init_u(model, protein_lookup)
end



"""
Rtg Model: version 1
"""
function rtgM1_model()
    rtgM1 = @reaction_network  begin
        #=
        Input
        =#
        (1), s→s
        #=
        Modulation Layer
        =#
        (k1md*s, kn1md), Rtg2_ina_c ↔ Rtg2_act_c
        (k2md, kn2md), Rtg2_act_c + Mks ↔ Rtg2Mks_c
        (k3md, kn3md), Bmh + Mks ↔ BmhMks
        #=
        Activation Model:
        =#
        # NLSact: Nuclear Localization Signal (NLS) activation
        (k2n + k1*(BmhMks)/(k1m+BmhMks)), Rtg13_a_c → Rtg13_i_c
        (k2, k2n), Rtg3_i_c ↔ Rtg3_a_c
        # Rtg1Binding: Rtg13 formation
        (k3, kn3), Rtg1_c + Rtg3_a_c ↔ Rtg13_a_c
        (k3, kn3), Rtg1_c + Rtg3_i_c ↔ Rtg13_i_c
        # Translocation
        (k4, kn4), Rtg13_a_c ↔ Rtg13_a_n
        (k4, kn4), Rtg3_a_c ↔ Rtg3_a_n

    end k1md kn1md k2md kn2md k3md kn3md k1 k1m k2 k2n k3 kn3 k4 kn4
    return rtgM1
end