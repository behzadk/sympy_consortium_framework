from species import Substrate
from species import AHL
from species import Microcin
from species import Strain

from model import Model
from cpp_output import Cpp_source_output
from cpp_output import Cpp_header_output

import model_space_generator
import csv
import ast

import utils
from tqdm import tqdm

class PowForDoubleStar(ast.NodeTransformer):
    def visit_BinOp(self, node):
        node.left = self.visit(node.left)
        node.right = self.visit(node.right)
        pow_func = ast.parse("math.pow", mode="eval").body
        if isinstance(node.op, ast.Pow):
            node = ast.copy_location(
                       ast.Call(func=pow_func,
                                args=[node.left, node.right],
                                keywords=[]
                               ),
                       node
                   )

        return node


def generate_simulation_files(model_list, params_path, init_species_path, output_dir):
    utils.make_folder(output_dir)
    print("Number of legal models: ", len(model_list))
    print("building equations and writing input files")

    print("Building cpp files.. ")
    for idx, m in enumerate(tqdm(model_list)):
        m.build_equations()
        m.build_symbolic_equations()
        m.build_jacobian()
        m.extract_species()
        m.extract_params()
        m.write_python_equations(output_dir)

        m.write_prior_parameter_dict(params_path, output_dir)
        m.write_init_species_dict(init_species_path, output_dir)

    print("Writing source and header files")
    cpp_out = Cpp_source_output(model_list)
    cpp_out.write_source_file(output_dir)
    header = Cpp_header_output(model_list)
    header.write_header_file(output_dir)

def generate_adjacency_matricies(model_list, substrate_ids, microcin_ids, AHL_ids, strain_ids, antitoxin_ids, immunity_ids, toxin_ids, output_dir):
    utils.make_folder(output_dir)

    for idx, m in enumerate(model_list):
        m.write_adj_matrix(output_dir, microcin_ids, AHL_ids, strain_ids, substrate_ids, antitoxin_ids, immunity_ids, toxin_ids)


def spock_manu_no_symm():
    output_dir = "./output/input_files_two_species_spock_manu_1/"

    default_params_path = './default_params/spock_manu/default_params_spock_manu.csv'
    default_init_species_path = './default_params/spock_manu/default_init_species_spock_manu.csv'

    # Set species IDs
    substrate_ids = ['glu']
    S_glu = Substrate(substrate_ids[0])
    # S_trp = Substrate(substrate_ids[1])

    substrate_objects = [S_glu]

    AHL_ids = ['1']
    AHL_1 = AHL(AHL_ids[0])
    # AHL_2 = AHL(AHL_ids[1])

    AHL_objects = [AHL_1]


    microcin_ids = ['1']
    toxin_ids = ['1T']
    antitoxin_ids = ['1T'] # Name of antitoxins must match name of toxin.(B_#ID# format)
    immunity_ids = ['1']
    strain_ids = ['1', '2']

    # Numerical maximum number of parts
    max_substrate_parts = len(substrate_objects)
    max_microcin_parts = len(microcin_ids)
    max_AHL_parts = len(AHL_objects)
    max_strains_parts = len(strain_ids)
    max_toxin_parts = len(toxin_ids)
    max_antitoxins = len(antitoxin_ids)
    max_immunity_parts = len(immunity_ids)

    # Generate microcin expression objects from AHLs and microcins
    microcin_objects, microcin_configs_df = model_space_generator.generate_microcin_combinations(microcin_ids,
                                                                                                 AHL_objects,
                                                                                                 microcin_induced=True,
                                                                                                 microcin_repressed=True, microcin_constitutive=True)
    toxin_objects, toxin_configs_df = model_space_generator.generate_toxin_combinations(toxin_ids,
                                                                                                 AHL_objects,
                                                                                                 toxin_induced=False,
                                                                                                 toxin_repressed=False, toxin_constitutive=False)


    antitoxin_objects, antitoxin_configs_df = model_space_generator.generate_antitoxin_combinations(antitoxin_ids,
                                                                                                 AHL_objects,
                                                                                                 antitoxin_induced=False,
                                                                                                 antitoxin_repressed=False, antitoxin_constitutive=False)

    immunity_objects, immunity_configs_df = model_space_generator.generate_immunity_combinations(immunity_ids,
                                                                                                 AHL_objects,
                                                                                                 immunity_induced=True,
                                                                                                 immunity_repressed=True, immunity_constitutive=True)


    model_space = model_space_generator.model_space(strain_ids, microcin_objects,
                                                    AHL_objects, substrate_objects, antitoxin_objects, immunity_objects, toxin_objects,
                                                    max_microcin_parts, max_AHL_parts,
                                                    max_substrate_parts, max_antitoxins, max_immunity_parts, max_toxin_parts, max_microcin_sensitivities=2)

    part_combos = model_space.generate_part_combinations(
        strain_max_microcin=1, strain_max_AHL=1, strain_max_sub_dependencies=1, 
        strain_max_microcin_sens=1, strain_max_sub_production=0, strain_max_antitoxin=1, 
        strain_max_immunity=1, strain_max_toxin=1
        )


    model_space.generate_models()

    # clean_models = []
    # for m in model_space.models_list:
    #     clean = True
    #     for strain in m.strains:
    #         if len(strain.sensitivities) == 0:
    #             clean = False

    #     if clean:
    #         clean_models.append(m)

    # model_space.models_list = clean_models

    model_space.spock_manu_model_filter()

    model_space.remove_symmetries()
    model_space.reset_model_indexes()

    model_list = model_space.models_list

    generate_adjacency_matricies(model_list, substrate_ids, microcin_ids, AHL_ids, strain_ids, antitoxin_ids, immunity_ids, toxin_ids, output_dir)
    generate_simulation_files(model_list, default_params_path, default_init_species_path, output_dir)

def single_strain_test():
    output_dir = "./output/input_files_one_species_0/"

    default_params_path = './default_params/spock_manu/default_params_spock_manu.csv'
    default_init_species_path = './default_params/spock_manu/default_init_species_spock_manu.csv'

    # Set species IDs
    substrate_ids = ['glu']
    S_glu = Substrate(substrate_ids[0])

    substrate_objects = [S_glu]

    AHL_ids = ['1']
    AHL_1 = AHL(AHL_ids[0])

    AHL_objects = [AHL_1]


    microcin_ids = ['1']
    antitoxin_ids = ['1'] # Name of antitoxins must match name of microcins.(B_#ID# format)
    strain_ids = ['1']

    # Numerical maximum number of parts
    max_substrate_parts = len(substrate_objects)
    max_microcin_parts = len(microcin_ids)
    max_AHL_parts = len(AHL_objects)
    max_strains_parts = len(strain_ids)
    max_antitoxins = len(antitoxin_ids)

    # Generate microcin expression objects from AHLs and microcins
    microcin_objects, microcin_configs_df = model_space_generator.generate_microcin_combinations(microcin_ids,
                                                                                                 AHL_objects,
                                                                                                 microcin_induced=True,
                                                                                                 microcin_repressed=True, microcin_constitutive=False)

    antitoxin_objects, antitoxin_configs_df = model_space_generator.generate_antitoxin_combinations(antitoxin_ids,
                                                                                                 AHL_objects,
                                                                                                 antitoxin_induced=True,
                                                                                                 antitoxin_repressed=True, antitoxin_constitutive=False)


    model_space = model_space_generator.model_space(strain_ids, microcin_objects,
                                                    AHL_objects, substrate_objects, antitoxin_objects,
                                                    max_microcin_parts, max_AHL_parts,
                                                    max_substrate_parts, max_antitoxins, max_microcin_sensitivities=2)

    part_combos = model_space.generate_part_combinations(strain_max_microcin=1, strain_max_AHL=1, strain_max_sub_dependencies=1, strain_max_microcin_sens=1, strain_max_sub_production=1, strain_max_antitoxin=1)

    print("Number of part combinations: ", len(part_combos))

    model_space.generate_models()
    model_space.one_species_filter()

    model_space.remove_symmetries()

    model_space.reset_model_indexes()

    model_list = model_space.models_list

    generate_adjacency_matricies(model_list, substrate_ids, microcin_ids, AHL_ids, strain_ids, antitoxin_ids, output_dir)
    generate_simulation_files(model_list, default_params_path, default_init_species_path, output_dir)


def three_species_no_symm():
    output_dir = "./output/input_files_three_species_0/"

    out_inputs_dir = output_dir + "input_files/"
    adj_matrix_out_dir = output_dir + "adj_matricies/"

    S_glu = Substrate('glu')
    AHL_ids = ['1', '2']

    AHL_1 = AHL(AHL_ids[0])
    AHL_2 = AHL(AHL_ids[1])

    AHL_objects = [AHL_1, AHL_2]
    substrate_objects = [S_glu]
    microcin_ids = ['1', '2']
    strain_ids = ['1', '2', '3']

    print("generating microcin combinations")
    microcin_objects, microcin_configs_df = model_space_generator.generate_microcin_combinations(microcin_ids,
                                                                                                 AHL_objects,
                                                                                                 microcin_induced=True,
                                                                                                 microcin_repressed=True,
                                                                                                 microcin_constitutive=True)

    max_substrate_parts = len(substrate_objects)
    max_microcin_parts = len(microcin_ids)
    max_AHL_parts = len(AHL_objects)
    max_strains_parts = len(strain_ids)

    print("generating model space")
    model_space = model_space_generator.model_space(strain_ids, microcin_objects, AHL_objects, substrate_objects,
                                                    max_microcin_parts=max_microcin_parts, max_AHL_parts=max_AHL_parts,
                                                    max_substrate_dependencies=max_substrate_parts,
                                                    max_microcin_sensitivities=max_microcin_parts)

    print("generating_part_combinations")
    part_combos = model_space.generate_part_combinations(1, 2, 1, 2)

    print("Generating model list")
    model_list = model_space.generate_models()
    print("generating reference table")

    # models_ref_df = model_space.generate_model_reference_table(max_microcin_parts, max_AHL_parts,
    #                                                            max_substrate_parts, max_microcin_parts)

    print("generating simulation files")
    generate_simulation_files(model_list, output_dir)

    print("generating generate_adjacency_matricies")

    generate_adjacency_matricies(model_list, substrate_ids, microcin_ids, AHL_ids, strain_ids, antitoxin_ids, output_dir)


def main():
    # single_strain_test()
    spock_manu_no_symm()
    # three_species_no_symm()

if __name__ == "__main__":
    main()
