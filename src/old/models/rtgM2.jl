@with_kw struct rtgM2{A,B,C} <: RTGmodel
    model::A = rtgM2_model()
    protein_lookup = get_protein_lookup(model)
    p::C = init_p(model; idx_hill_coefs=[])
    u::B = init_u(model, protein_lookup)
end


"""
This model is modified version of rtgM1.
Major changes:
1. Change input to be michaelis menten
2. Change the translocation to be Michaelis menten
3. Rtg1p and Rtg3p enter nucleus separately
"""
function rtgM2_model()
    rtgM2 = @reaction_network  begin
        #=
        This model is modified version of rtgM1.
        Major changes:
        1. Change input to be michaelis menten
        2. Change the translocation to be Michaelis menten
        3. Rtg1p and Rtg3p enter nucleus separately
        =#
        #=
        Input
        =#
        (1), s→s
        #=
        Modulation Layer
        =#
        (k1md*s*(k1mk+s), kn1md), Rtg2_ina_c ↔ Rtg2_act_c #[1]
        (k2md, kn2md), Rtg2_act_c + Mks ↔ Rtg2Mks_c
        (k3md, kn3md), Bmh + Mks ↔ BmhMks
        #=
        Activation Model:
        =#
        # NLSact: Nuclear Localization Signal (NLS) activation
        (k2n + mm(BmhMks, k1, k1m)), Rtg13_a_c → Rtg13_i_c
        (k2, k2n), Rtg3_i_c ↔ Rtg3_a_c

        # Rtg1Binding: Rtg13 formation
        (k3, kn3), Rtg1_c + Rtg3_a_c ↔ Rtg13_a_c
        (k3, kn3), Rtg1_c + Rtg3_i_c ↔ Rtg13_i_c

        (k3nu, kn3nu), Rtg1_n + Rtg3_a_n ↔ Rtg13_a_n #[3]

        # Translocation
        ( mm(Rtg1_c, kt1, ktd1)  , mm(Rtg1_c, kt1nu, ktd1nu)), Rtg1_c ↔ Rtg1_n #[2]
        ( mm(Rtg3_a_c, kt2, ktd2) , mm(Rtg3_a_n, kt2nu, ktd2nu)), Rtg3_a_c ↔ Rtg3_a_n #[2]

        #(k4, kn4), Rtg13_a_c ↔ Rtg13_a_n # Deleted [2]
    end k1mk k1md k1mk kn1md k2md kn2md k3md kn3md k1 k1m k2 k2n k3 kn3 k4 kn4 kt1 ktd1 kt1nu ktd1nu kt2 ktd2 kt2nu ktd2nu k3nu kn3nu

    return rtgM2
end