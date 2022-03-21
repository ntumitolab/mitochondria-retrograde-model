@with_kw struct rtgM3{A,B,C} <: RTGmodel
    model::A = rtgM3_model()
    protein_lookup = get_protein_lookup(model)
    p::C = init_p(model; idx_hill_coefs=[])
    u::B = init_u(model, protein_lookup)
end


"""
This model is modified version of rtgM2.
Major changes:
1. Change input to be michaelis menten
2. Change the translocation to be Michaelis menten: (which is not ture). See [supplemental material](https://www.pnas.org/content/pnas/suppl/2006/06/30/0604085103.DC1/04085SuppText.pdf)
3. Rtg1p and Rtg3p enter nucleus separately

References 
-----------

1. Hao, N., Budnik, B. A., Gunawardena, J., & O'Shea, E. K. (2013). Tunable signal processing through modular control of transcription factor translocation. Science, 339(6118), 460-464. [[Article](https://science.sciencemag.org/content/339/6118/460), [SOM](https://www.pnas.org/content/pnas/suppl/2006/06/30/0604085103.DC1/04085SuppText.pdf)]

2. Mitochondrial Signaling: The Retrograde Response. 2004. (Reaction of Bmh-Mks on Rtg13 dimer) [[Article](https://www.cell.com/molecular-cell/pdf/S1097-2765(04)00179-0.pdf)]
"""
function rtgM3_model()
    rtgM3= @reaction_network  begin
        #=
        Input
        =#
        (1), s→s

        #=
        Modulation Layer
        =#
        (mm(s, ksV, ksD), k2I), Rtg2_ina_c ↔ Rtg2_act_c #[1]
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
        (k1in , k1out), Rtg1_c ↔ Rtg1_n #[2]
        (k3inA, k3outA), Rtg3_a_c ↔ Rtg3_a_n #[2]
        (k3inI, k3outI), Rtg3_i_c ↔ Rtg3_i_n

        #(k4, kn4), Rtg13_a_c ↔ Rtg13_a_n # Deleted [2]
    end ksV ksD k2I k2M kn2M kBM knBM k13I k13IV k13ID k3A_c k3I_c k3I_n k13_c kn13_c k13_n kn13_n k1in k1out k3inA k3outA k3inI k3outI

    ##
    return rtgM3
end