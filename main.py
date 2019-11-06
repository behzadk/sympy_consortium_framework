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

import predefined_model_spaces

# class PowForDoubleStar(ast.NodeTransformer):
#     def visit_BinOp(self, node):
#         node.left = self.visit(node.left)
#         node.right = self.visit(node.right)
#         pow_func = ast.parse("math.pow", mode="eval").body
#         if isinstance(node.op, ast.Pow):
#             node = ast.copy_location(
#                        ast.Call(func=pow_func,
#                                 args=[node.left, node.right],
#                                 keywords=[]
#                                ),
#                        node
#                    )

#         return node


def spock_manu_no_symm():
    output_dir = "./output/input_files_two_species_spock_manu_1/"

    default_params_path = './default_params/spock_manu/default_params_spock_manu.csv'
    default_init_species_path = './default_params/spock_manu/default_init_species_spock_manu.csv'

    model_space = predefined_model_spaces.spock_manu_no_symm(default_params_path, default_init_species_path, output_dir)

def single_strain_test():
    output_dir = "./output/input_files_one_species_0/"
    default_params_path = './default_params/default/default_params.csv'
    default_init_species_path = './default_params/default/default_init_species.csv'

    model_space = predefined_model_spaces.three_species_no_symm(default_params_path, default_init_species_path, output_dir)

def two_species_no_symm():
    output_dir = "./output/input_files_two_species_0/"
    default_params_path = './default_params/default/default_params.csv'
    default_init_species_path = './default_params/default/default_init_species.csv'

    model_space = predefined_model_spaces.two_species_no_symm(default_params_path, default_init_species_path, output_dir)

def two_species_no_symm_auxos():
    output_dir = "./output/input_files_two_species_auxos_0/"
    default_params_path = './default_params/default/default_params.csv'
    default_init_species_path = './default_params/default/default_init_species.csv'

    model_space = predefined_model_spaces.two_species_no_symm_auxos(default_params_path, default_init_species_path, output_dir)

def three_species_no_symm_auxos():
    output_dir = "./output/input_files_three_species_auxos_0/"
    default_params_path = './default_params/default/default_params.csv'
    default_init_species_path = './default_params/default/default_init_species.csv'

    model_space = predefined_model_spaces.three_species_no_symm_auxos(default_params_path, default_init_species_path, output_dir)


def three_species_no_symm():
    output_dir = "./output/input_files_three_species_0/"
    default_params_path = './default_params/default/default_params.csv'
    default_init_species_path = './default_params/default/default_init_species.csv'

    model_space = predefined_model_spaces.three_species_no_symm(default_params_path, default_init_species_path, output_dir)

def main():
    # single_strain_test()
    # spock_manu_no_symm()
    # three_species_no_symm()
    two_species_no_symm()
    # two_species_no_symm_auxos()
    # three_species_no_symm_auxos()

if __name__ == "__main__":
    main()
