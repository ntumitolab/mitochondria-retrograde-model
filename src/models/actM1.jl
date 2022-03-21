@with_kw struct actM1{A,B,C} <: RTGmodel
    model::A = actM1_model()
    protein_lookup = get_protein_lookup(model)
    p::C = init_p(model; idx_hill_coefs=[])
    u::B = init_u(model, protein_lookup)
end


"""
Activation Model:
"""
function actM1_model()
    actM1 = @reaction_network begin
        # NLSact: Nuclear Localization Signal (NLS) activation
        (1), BmhMks → BmhMks
        (k2n + k1*(BmhMks)/(k1m+BmhMks)), Rtg13_a_c → Rtg13_i_c
        (k2, k2n), Rtg3_i_c ↔ Rtg3_a_c
        # Rtg1Binding: Rtg13 formation
        (k3, kn3), Rtg1_c + Rtg3_a_c ↔ Rtg13_a_c
        (k3, kn3), Rtg1_c + Rtg3_i_c ↔ Rtg13_i_c
        # Translocation
        (k4, kn4), Rtg13_a_c ↔ Rtg13_a_n
        (k4, kn4), Rtg3_a_c ↔ Rtg3_a_n
    end k1 k1m k2 k2n k3 kn3 k4 kn4
    return actM1
end