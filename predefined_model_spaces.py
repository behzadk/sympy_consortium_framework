from species import Microcin
from species import  AHL
from species import Strain
from species import Substrate

from model import Model
from cpp_output import Cpp_source_output
from cpp_output import Cpp_header_output


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


if __name__== "__main__":
    known_two_strain_systems()