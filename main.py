from species import Substrate
from species import AHL
from species import Microcin
from species import Strain

from model import Model
from cpp_output import Cpp_source_output
from cpp_output import Cpp_header_output

import model_space_generator

def two_species_no_symm():
    output_dir = "./output/two_species_no_symm/input_files/"
    adj_matrix_out_dir = "./output/two_species_no_symm/adj_matricies/"
    S_glu = Substrate('glu')
    AHL_ids = ['1', '2']
    AHL_1 = AHL(AHL_ids[0])
    AHL_2 = AHL(AHL_ids[1])

    AHL_objects = [AHL_1, AHL_2]
    substrate_objects = [S_glu]
    microcin_ids = ['mccV', 'mccB']
    strain_ids = ['x', 'c']

    microcin_objects, microcin_configs_df = model_space_generator.generate_microcin_combinations(microcin_ids,
                                                                                                 AHL_objects,
                                                                                                 microcin_induced=True,
                                                                                                 microcin_repressed=False)
    max_substrate_parts = len(substrate_objects)
    max_microcin_parts = len(AHL_objects)
    max_AHL_parts = len(microcin_ids)
    max_strains_parts = len(strain_ids)


    model_space = model_space_generator.model_space(strain_ids, microcin_objects,
                                                    AHL_objects, substrate_objects,
                                                    max_microcin_parts, max_AHL_parts,
                                                    max_substrate_parts, max_microcin_sensitivities=2)

    part_combos = model_space.generate_part_combinations(1, 2, 1, 2)

    print("Number of part combinations: ", len(part_combos))

    model_list = model_space.generate_models()

    models_ref_df = model_space.generate_model_reference_table(max_microcin_parts, max_AHL_parts,
                                                               max_substrate_parts, max_microcin_sensitivities=2)

    microcin_configs_df.to_csv(output_dir + 'microcin_config.csv')
    models_ref_df.to_csv(output_dir + 'model_ref.csv')
    print("Number of legal models: ", len(model_list))
    print("building equations and writing input files")

    for idx, m in enumerate(model_list):
        m.build_equations()
        m.build_symbolic_equations()
        m.build_jacobian()
        m.extract_species()
        m.extract_params()

        default_params_path = './default_params/default_params.csv'
        default_init_species_path = './default_params/default_init_species.csv'

        adj_mat_path = adj_matrix_out_dir + 'model_' + str(idx) + '_adj_mat.csv'
        m.write_adj_matrix(adj_mat_path, microcin_ids, AHL_ids, strain_ids)
        m.write_prior_parameter_dict(default_params_path, output_dir + 'params_' + str(idx) + '.csv')
        m.write_init_species_dict(default_init_species_path, output_dir + 'species_' + str(idx) + '.csv')

    print("writing source and header files")
    cpp_out = Cpp_source_output(model_list)
    cpp_out.write_source_file(output_dir + 'model.cpp')
    header = Cpp_header_output(model_list)
    header.write_header_file(output_dir + 'model.h')

def three_species_no_symm():
    output_dir = "./output/3_species_space/"
    S_glu = Substrate('glu')
    AHL_1 = AHL('1')
    AHL_2 = AHL('2')

    AHL_objects = [AHL_1, AHL_2]
    substrate_objects = [S_glu]
    microcin_ids = ['mccV', 'mccB']
    strain_ids = ['x', 'c', 'y']
    microcin_objects, microcin_configs_df = model_space_generator.generate_microcin_combinations(microcin_ids,
                                                                                                 AHL_objects,
                                                                                                 microcin_induced=True,
                                                                                                 microcin_repressed=False)

    model_space = model_space_generator.model_space(strain_ids, microcin_objects, AHL_objects, substrate_objects)
    part_combos = model_space.generate_part_combinations(1, 1, 1, max_microcin_sensitivities=2)
    model_list = model_space.generate_models()

    models_ref_df = model_space.generate_model_reference_table(1, 1, 1, max_microcin_sensitivities=2)

    microcin_configs_df.to_csv(output_dir + 'microcin_config.csv')
    models_ref_df.to_csv(output_dir + 'model_ref.csv')

    print("Number of part combinations: ", len(part_combos))
    print("Number of legal models: ", len(model_list))
    print("building equations and writing input files")

    for idx, m in enumerate(model_list):
        m.build_equations()
        m.build_symbolic_equations()
        m.build_jacobian()
        m.extract_species()
        m.extract_params()

        default_params_path = './default_params/default_params.csv'
        default_init_species_path = './default_params/default_init_species.csv'

        m.write_prior_parameter_dict(default_params_path, output_dir + 'params_' + str(idx) + '.csv')
        m.write_init_species_dict(default_init_species_path, output_dir + 'species_' + str(idx) + '.csv')

    print("writing source and header files")
    cpp_out = Cpp_source_output(model_list)
    cpp_out.write_source_file(output_dir + 'model.cpp')
    header = Cpp_header_output(model_list)
    header.write_header_file(output_dir + 'model.h')


def main():
    two_species_no_symm()
    # three_species_no_symm()

if __name__ == "__main__":
    main()
