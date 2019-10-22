from species import Microcin
from species import AHL
from species import Strain
from species import Substrate
from species import Toxin
from species import Antitoxin

from model import Model
from cpp_output import Cpp_source_output
from cpp_output import Cpp_header_output

import model_space_generator


def known_two_strain_systems():
    AHL_1 = AHL('1')
    AHL_2 = AHL('2')
    S_glu = Substrate('glu')

    # Model a
    B_1 = Microcin(0, '1', [AHL_1], [])
    B_2 = Microcin(1, '2', [AHL_2], [])
    nx = Strain('x', [B_1], [AHL_1], [S_glu], [B_1.id])
    nc = Strain('c', [B_2], [AHL_2], [S_glu],  [B_2.id])
    model_a = Model(0, [nx, nc])

    # Model b
    B_1 = Microcin(0, '1', [], [AHL_2])
    B_2 = Microcin(1, '2', [AHL_1], [])

    nx = Strain('x', [B_1], [AHL_1], [S_glu], [B_1.id])
    nc = Strain('c', [B_2], [AHL_2], [S_glu], [B_2.id])
    model_b = Model(1, [nx, nc])

    # Model c
    B_1 = Microcin(0, '1', [], [AHL_1])

    nx = Strain('x',  [B_1], [AHL_1], [S_glu], [])
    nc = Strain('c', [], [], [S_glu],  [B_1.id])
    model_c = Model(2, [nx, nc])

    # Model D
    B_1 = Microcin(0, '1', [AHL_1], [])

    nx = Strain('x', [B_1], [], [S_glu],  [])
    nc = Strain('c', [], [AHL_1], [S_glu], [B_1.id])
    model_d = Model(3, [nx, nc])

    # Model E
    B_1 = Microcin(0, '1', [AHL_2], [AHL_1])

    nx = Strain('x',  [B_1], [AHL_1], [S_glu], [])
    nc = Strain('c', [], [AHL_2],  [S_glu], [B_1.id])
    model_e = Model(4, [nx, nc])

    # Model F
    B_1 = Microcin(0, '1', [AHL_1], [])
    B_2 = Microcin(1, '2', [], [AHL_1])

    nx = Strain('x',  [B_1, B_2], [AHL_1, AHL_2], [S_glu], [B_1.id])
    nc = Strain('c', [], [],  [S_glu], [B_2.id])
    model_f = Model(5, [nx, nc])

    # Mode g
    B_1 = Microcin(0, '1', [AHL_1], [])
    B_2 = Microcin(1, '2', [AHL_2], [AHL_1])

    nx = Strain('x',  [B_1, B_2], [AHL_1], [S_glu], [B_1.id])
    nc = Strain('c', [], [AHL_2],  [S_glu], [B_2.id])
    model_g = Model(6, [nx, nc])

    model_list = [model_a, model_b, model_c, model_d, model_e, model_f, model_g]

    for idx, m in enumerate(model_list):
        m.build_equations()
        m.build_symbolic_equations()
        m.build_jacobian()
        m.extract_species()
        m.extract_params()

        default_params_path = './default_params/default_params.csv'
        default_init_species_path = './default_params/default_init_species.csv'

        m.write_prior_parameter_dict(default_params_path, './output/known_2_species/params_' + str(idx) + '.csv')
        m.write_init_species_dict(default_init_species_path, './output/known_2_species/species_' + str(idx) + '.csv')
        print("Model ", idx, "Is legal?: ", m.is_legal())

    print("writing source and header files")
    cpp_out = Cpp_source_output(model_list)
    cpp_out.write_source_file("./output/known_2_species/model.cpp")
    header = Cpp_header_output(model_list)
    header.write_header_file("./output/known_2_species/model.h")


def one_strain_tests():
    pass


def balagadde(model_idx, adj_mat_out_dir):

    # Set species IDs
    substrate_ids = ['glu']
    S_glu = Substrate(substrate_ids[0])
    # S_trp = Substrate(substrate_ids[1])

    substrate_objects = [S_glu]

    AHL_ids = ['1', '2']
    AHL_1 = AHL(AHL_ids[0])
    AHL_2 = AHL(AHL_ids[1])

    AHL_objects = [AHL_1, AHL_2]

    toxin_ids = ['1T', '2T']
    antitoxin_ids = ['1T', '2T']  # Name of antitoxins must match name of toxin.(B_#ID# format)
    strain_ids = ['1', '2']

    A_1 = AHL(AHL_ids[0])
    A_2 = AHL(AHL_ids[1])

    T_1 = Toxin(config_idx=0, toxin_id=toxin_ids[0], AHL_inducer_list=[], AHL_repressor_list=[], constitutive_expression=True)

    # AHL from strain 1 induces toxin in strain 2
    T_2 = Toxin(config_idx=1, toxin_id=toxin_ids[1], AHL_inducer_list=[A_1], AHL_repressor_list=[], constitutive_expression=False)

    # Immunity in strain 1 is induced by AHL 2
    V_1 = Antitoxin(config_idx=0, antitoxin_id=antitoxin_ids[0], AHL_inducer_list=[A_2], AHL_repressor_list=[], constitutive_expression=False)

    N_1 = Strain(strain_id=strain_ids[0],
                 microcin_expression=[],
                 AHL_expression=[A_1],
                 substrate_dependences=[S_glu],
                 microcin_sensitivities=[],
                 substrate_production=[],
                 antitoxins=[V_1],
                 immunity_expression=[],
                 toxin_expression=[T_1])

    N_2 = Strain(strain_id=strain_ids[1],
                 microcin_expression=[],
                 AHL_expression=[A_2],
                 substrate_dependences=[S_glu],
                 microcin_sensitivities=[],
                 substrate_production=[],
                 antitoxins=[],
                 immunity_expression=[],
                 toxin_expression=[T_2])

    bala_model = Model(model_idx=model_idx, strain_list=[N_1, N_2])
    bala_model.generate_adjacency_matrix(len(substrate_ids), len(AHL_ids), 0, len(strain_ids), len(antitoxin_ids), 0, len(toxin_ids))
    bala_model.write_adj_matrix(adj_mat_out_dir, [], AHL_ids, strain_ids, substrate_ids, antitoxin_ids, [], toxin_ids)

    return bala_model

def scott(model_idx, adj_mat_out_dir):

    # Set species IDs
    substrate_ids = ['glu']
    S_glu = Substrate(substrate_ids[0])
    # S_trp = Substrate(substrate_ids[1])

    substrate_objects = [S_glu]

    AHL_ids = ['1', '2']
    AHL_1 = AHL(AHL_ids[0])
    AHL_2 = AHL(AHL_ids[1])

    AHL_objects = [AHL_1, AHL_2]

    toxin_ids = ['1T', '2T']
    strain_ids = ['1', '2']

    A_1 = AHL(AHL_ids[0])
    A_2 = AHL(AHL_ids[1])

    T_1 = Toxin(config_idx=0, toxin_id=toxin_ids[0], AHL_inducer_list=[A_1], AHL_repressor_list=[], constitutive_expression=True)

    # AHL from strain 1 induces toxin in strain 2
    T_2 = Toxin(config_idx=1, toxin_id=toxin_ids[1], AHL_inducer_list=[A_2], AHL_repressor_list=[], constitutive_expression=False)


    N_1 = Strain(strain_id=strain_ids[0],
                 microcin_expression=[],
                 AHL_expression=[A_1],
                 substrate_dependences=[S_glu],
                 microcin_sensitivities=[],
                 substrate_production=[],
                 antitoxins=[],
                 immunity_expression=[],
                 toxin_expression=[T_1])

    N_2 = Strain(strain_id=strain_ids[1],
                 microcin_expression=[],
                 AHL_expression=[A_2],
                 substrate_dependences=[S_glu],
                 microcin_sensitivities=[],
                 substrate_production=[],
                 antitoxins=[],
                 immunity_expression=[],
                 toxin_expression=[T_2])

    scott_model = Model(model_idx=model_idx, strain_list=[N_1, N_2])
    scott_model.generate_adjacency_matrix(len(substrate_ids), len(AHL_ids), 0, len(strain_ids), 0, len(toxin_ids))
    scott_model.write_adj_matrix(adj_mat_out_dir, [], AHL_ids, strain_ids, substrate_ids, [], [], toxin_ids)

    return scott_model


def mccardell(model_idx, adj_mat_out_dir):

    # Set species IDs
    substrate_ids = ['glu']
    S_glu = Substrate(substrate_ids[0])
    # S_trp = Substrate(substrate_ids[1])

    substrate_objects = [S_glu]

    AHL_ids = ['1', '2']
    AHL_1 = AHL(AHL_ids[0])
    AHL_2 = AHL(AHL_ids[1])

    AHL_objects = [AHL_1, AHL_2]

    toxin_ids = ['1T', '2T']
    strain_ids = ['1', '2']

    A_1 = AHL(AHL_ids[0])
    A_2 = AHL(AHL_ids[1])

    T_1 = Toxin(config_idx=0, toxin_id=toxin_ids[0], AHL_inducer_list=[A_1], AHL_repressor_list=[], constitutive_expression=True)

    # AHL from strain 1 induces toxin in strain 2
    T_2 = Toxin(config_idx=1, toxin_id=toxin_ids[1], AHL_inducer_list=[A_2], AHL_repressor_list=[A_1], constitutive_expression=False)


    N_1 = Strain(strain_id=strain_ids[0],
                 microcin_expression=[],
                 AHL_expression=[A_1],
                 substrate_dependences=[S_glu],
                 microcin_sensitivities=[],
                 substrate_production=[],
                 antitoxins=[],
                 immunity_expression=[],
                 toxin_expression=[T_1])

    N_2 = Strain(strain_id=strain_ids[1],
                 microcin_expression=[],
                 AHL_expression=[A_2],
                 substrate_dependences=[S_glu],
                 microcin_sensitivities=[],
                 substrate_production=[],
                 antitoxins=[],
                 immunity_expression=[],
                 toxin_expression=[T_2])


    mccardell = Model(model_idx=model_idx, strain_list=[N_1, N_2])

    mccardell.generate_adjacency_matrix(len(substrate_ids), len(AHL_ids), 0, len(strain_ids), 0, 0, len(toxin_ids))
    mccardell.write_adj_matrix(adj_mat_out_dir, [], AHL_ids, strain_ids, substrate_ids, [], [], toxin_ids)

    return mccardell


if __name__ == "__main__":
    balagadde()
    # known_two_strain_systems()
