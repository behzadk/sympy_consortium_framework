import sympy
import numpy as np

funcs = {
    # Strain growth rate
    'mu_#N#': '( S_#S# / ( K_mu_#S# + S_#S# ) )',

    # Induction of bacteriocin expression by AHL
    'k_b_ind_#B#': '( A_#A# ^ nB_#B# / ( KB_#B# ^ nB_#B# + A_#A# ^ nB_#B# ) )',

    # Repression of bacteriocin expression by AHL
    'k_b_repr_#B#': '( KB_#B# ^ nB_#B# / ( KB_#B# ^ nB_#B# + A_#A# ^ nB_#B# ) )',

    # Production of an AHL species
    'A_production': 'kA_#A# * N_#N# * C',

    # Function defining sensitivity to microcin
    'omega': '( omega_max_#B# * B_#B# ^ n_omega_#B# / ( K_omega_#B# ^ n_omega_#B# + B_#B# ^ n_omega_#B# ) )',

    # Function defining production by antitoxin V
    'V_antitoxin': '( K_V_#V# / ( #V# + K_V_#V# ) )',

    # Induction of antitoxin expression by AHL
    'k_v_ind_#V#': '( A_#A# ^ nV_#V# / ( kV_#V# ^ nV_#V#  + A_#A# ^ nV_#V#) )',

    # Repression of antitoxin expression by AHL
    'k_v_repr_#V#': '( KB_#V# ^ nV_#V#  / ( kV_#V# ^ nV_#V#  + A_#A# ^ nV_#V# ) )'

}

# Base eqs contain the description of species before interactions with other species.
# This usually consists of the dilution term which we can apply production terms to.
base_eqs = {
    'N_#N#': '( - D * N_#N# )',
    'S_#S#': '( D * ( S0_#S# - S_#S# ) )',
    'B_#B#': '( - D * B_#B# )',
    'A_#A#': '( - D * A_#A# )',
    'V_#V#': '( - D * V_#V# )'
}

def gen_strain_growth_diff(strain_id, strain_list):
    dN_dt = base_eqs['N_#N#']

    # Get antitoxins
    antitoxin_list = []
    for strain in strain_list:
        antitoxin_list += strain.antitoxins

    antitoxin_list = list(set([x.id for x in antitoxin_list]))

    for strain in strain_list:
        if strain.id is strain_id:
            dN_dt = dN_dt + ' + N_#N# * mu_max_#N# '
            for s in strain.substrate_dependences:
                dN_dt = dN_dt + ' * ' + funcs['mu_#N#'].replace('#S#', s.id)

            # Strain sensitive to microcin, protected by 'antitoxins'
            for m in strain.sensitivities:
                dN_dt = dN_dt + ' * ' + funcs['omega']
                dN_dt = dN_dt.replace('#B#', m)

                # Apply antitoxin function if cognate antitoxin is present.
                for v in antitoxin_list:
                    if v.split('_')[-1] == m:
                        dN_dt = dN_dt + ' * ' + funcs['V_antitoxin']
                        dN_dt = dN_dt.replace('#V#', v)
                        break

                dN_dt = dN_dt + ' * N_#N#'


    dN_dt = dN_dt.replace('#N#', strain_id)
    N_key = 'N_#N#'.replace('#N#', strain_id)

    return {N_key: dN_dt}

def gen_diff_eq_antitoxin(antitoxin_id, strain_list):
    dV_dt = base_eqs['V_#V#']

    for strain in strain_list:

        for v in strain.antitoxins:
            if v.id is antitoxin_id:
                dV_dt = dV_dt + ' + ' + ' kV_max_#V# '

                if v.AHL_inducers is not np.nan:
                    # Induction terms
                    for a in v.AHL_inducers:
                        dV_dt = dV_dt + ' * ' + funcs['k_v_ind_#V#'].replace('#A#', a.id)

                if v.AHL_repressors is not np.nan:
                    for a in v.AHL_repressors:
                        dV_dt = dV_dt + ' * ' + funcs['k_v_repr_#V#'].replace('#A#', a.id)

                dV_dt = dV_dt + ' * N_#N# * C'
                dV_dt = dV_dt.replace('#N#', strain.id)

    dV_dt = dV_dt.replace('#V#', antitoxin_id)
    V_key = 'V_#V#'.replace('#V#', antitoxin_id)

    return {V_key: dV_dt}



def gen_diff_eq_substrate(substrate_id, strain_list):
    dS_dt = base_eqs['S_#S#'].replace('#S#', substrate_id)

    # Term defining consumption of substrate by a strain
    strain_growth_rate = funcs['mu_#N#']
    strain_consumption = strain_growth_rate + ' * N_#N# * C / g_#N#'
    strain_production = ' C * N_#N# * p_#S#'

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
