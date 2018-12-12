import sympy

funcs = {
    # Strain growth rate
    'mu_#N#': '( mu_max_#N# * S_#S# / ( K_#N# + S_#S# ) )',

    # Bacteriocin induced death
    'omega_#B#': '- (omega_max_#B# * B_#B# ^ n_omega_#B# ) / ( K_omega_#B#^n_omega_#B#  + B_#B#^n_omega_#B# )',

    # Induction of bacteriocin expression by AHL
    'k_b_ind_#B#': '( A_#A# ^ nB_#B# / ( KB_#B# ^ nB_#B# + A_#A# ^ nB_#B# ) )',

    # Repression of bacteriocin expression by AHL
    'k_b_repr_#B#': '( KB_#B# ^ nB_#B# /  \
                    ( KB_#B# ^ nB_#B# + A_#A# ^ nB_#B# ) )'
}

# Base eqs contain the description of species before interactions with other species
base_eqs = {
    'N_#N#': '( - D * N_#N# )',
    'S_#S#': '( - D * ( S0_#S# - S_#S# ) )',
    'B_#B#': '( - D * B_#B# )',
    'A_#A#': '( - D * A_#A# )'
}


def gen_diff_strain_growth(strain_list):
    pass

def gen_diff_eq_substrate(substrate_id, strain_list):
    dS_dt = base_eqs['S_#S#'].replace('#S#', substrate_id)

    # Term defining consumption of substrate by a strain
    strain_growth_rate = funcs['mu_#N#']
    # strain_consumption = ' '.join([strain_growth_rate, ' * N_#N# ', ' / g_#N# '])
    strain_consumption = strain_growth_rate + ' * N_#N# / g_#N#'

    for strain in strain_list:
        substrates = strain.substrate_dependence
        if substrate_id in substrates:
            dS_dt = dS_dt + ' - ' + strain_consumption
            dS_dt = dS_dt.replace('#N#', strain.id)



    dS_dt = dS_dt.replace('#S#', substrate_id)

    return {'S_#S#'.replace('#S#', substrate_id): dS_dt}

def gen_AHL_diff_eq(AHL_ids, strain_list):
