import sympy
import numpy as np

funcs = {
    # Strain growth rate
    'mu_#N#': '( mu_max_#N# * S_#S# / ( K_#N# + S_#S# ) )',


    # Induction of bacteriocin expression by AHL
    'k_b_ind_#B#': '( A_#A# ^ nB_#B# / ( KB_#B# ^ nB_#B# + A_#A# ^ nB_#B# ) )',

    # Repression of bacteriocin expression by AHL
    'k_b_repr_#B#': '( KB_#B# ^ nB_#B# / ( KB_#B# ^ nB_#B# + A_#A# ^ nB_#B# ) )',

    # Production of an AHL species
    'A_production': 'kA_#A# * N_#N# * C',

    # Function defining sensitivity to microcin
    'omega': '( omega_max_#B# * B_#B#^n_omega_#B# / ( K_omega_#B# ^ n_omega_#B# + B_#B# ^ n_omega_#B# ) )'
}

# Base eqs contain the description of species before interactions with other species.
# This usually consists of the dilution term which we can apply production terms to.
base_eqs = {
    'N_#N#': '( - D * N_#N# )',
    'S_#S#': '( D * ( S0_#S# - S_#S# ) )',
    'B_#B#': '( - D * B_#B# )',
    'A_#A#': '( - D * A_#A# )'
}

def gen_strain_growth_diff(strain_id, strain_list):

    dN_dt = base_eqs['N_#N#']

    for strain in strain_list:
        if strain.id is strain_id:
            dN_dt = dN_dt + ' + N_#N# '
            for s in strain.substrate_dependences:
                dN_dt = dN_dt + ' * ' + funcs['mu_#N#'].replace('#S#', s.id)

            for m in strain.sensitivities:
                dN_dt = dN_dt + ' - ' + funcs['omega'] + ' * N_#N#'
                dN_dt = dN_dt.replace('#B#', m)

    dN_dt = dN_dt.replace('#N#', strain_id)
    N_key = 'N_#N#'.replace('#N#', strain_id)

    return {N_key: dN_dt}


def gen_diff_eq_substrate(substrate_id, strain_list):
    dS_dt = base_eqs['S_#S#'].replace('#S#', substrate_id)

    # Term defining consumption of substrate by a strain
    strain_growth_rate = funcs['mu_#N#']
    strain_consumption = strain_growth_rate + ' * N_#N# * C / g_#N#'
    strain_production = ' N_#N# * p_#S#'

    # Sum of all consumption by strains
    for strain in strain_list:
        substrates = strain.substrate_dependences
        for consume_s in substrates:
            if consume_s.id is substrate_id:
                dS_dt = dS_dt + ' - ' + strain_consumption
                dS_dt = dS_dt.replace('#N#', strain.id)

    for strain in strain_list:
        for produce_s in strain.substrate_production:
            if produce_s.id is substrate_id:
                dS_dt = dS_dt + ' + ' + strain_production
                dS_dt = dS_dt.replace('#N#', strain.id)

    dS_dt = dS_dt.replace('#S#', substrate_id)
    S_key = 'S_#S#'.replace('#S#', substrate_id)
    return {S_key: dS_dt}


def gen_AHL_diff_eq(AHL_id, strain_list):
    dA_dt = base_eqs['A_#A#']
    production_term = funcs['A_production']

    for strain in strain_list:
        for a in strain.AHLs:
            if a.id is AHL_id:
                dA_dt = dA_dt + ' + ' + production_term.replace('#N#', strain.id)

    dA_dt = dA_dt.replace('#A#', AHL_id)
    AHL_key = 'A_#A#'.replace('#A#', AHL_id)

    return {AHL_key: dA_dt}


def gen_microcin_diff_eq(microcin_id, strain_list):
    dB_dt = base_eqs['B_#B#']

    for strain in strain_list:

        for b in strain.microcins:
            if b.id is microcin_id:
                dB_dt = dB_dt + ' + ' + ' kBmax_#B# '

                if b.AHL_inducers is not np.nan:
                    # Induction terms
                    for a in b.AHL_inducers:
                        dB_dt = dB_dt + ' * ' + funcs['k_b_ind_#B#'].replace('#A#', a.id)

                if b.AHL_repressors is not np.nan:
                    for a in b.AHL_repressors:
                        dB_dt = dB_dt + ' * ' + funcs['k_b_repr_#B#'].replace('#A#', a.id)

                dB_dt = dB_dt + ' * N_#N# * C'
                dB_dt = dB_dt.replace('#N#', strain.id)

    dB_dt = dB_dt.replace('#B#', microcin_id)
    B_key = 'B_#B#'.replace('#B#', microcin_id)

    return {B_key: dB_dt}
