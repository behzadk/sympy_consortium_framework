import sympy as sp
from collections import OrderedDict

from cpp_output import Cpp_source_output
from cpp_output import Cpp_header_output

def modified_LV_equations():
    species_list = ['G_1', 'H_1', 'H_2', 'P_1']
    diff_eqs = OrderedDict()

    dG_1 = 'a_g  - (H_1 * G_1 * a_1)  - (H_2 * G_1 * a_2)'
    dH_1 = 'H_1 * G_1 * a_1 - (H_1^0.5) * (alpha_1 * P_1)'
    dH_2 = 'H_2 * G_1 * a_2 - (H_2^0.5) * (alpha_2 * P_1)'
    dP_1 = 'P_1 * ( (beta_1 * H_1^(0.5)) + (beta_2 * H_2^(0.5))  - b)'

    zeros_list = [0 for i in range(len(species_list))]
    symbolic_equations = sp.Matrix(zeros_list)
   
    for idx, eq in enumerate([dG_1, dH_1, dH_2, dP_1]):
        species_name = species_list[idx]
        symbolic_equations[idx] = sp.sympify(eq, locals=locals())
        diff_eqs[species_name] = symbolic_equations[idx]

    return symbolic_equations, species_list, diff_eqs

def LV_equations():
    species_list = ['H_1', 'H_2', 'P_1']
    diff_eqs = OrderedDict()

    dH_1 = 'H_1 * (a_1 - (alpha_1 * P_1))'
    dH_2 = 'H_2 * (a_2 - (alpha_2 * P_1))'
    dP_1 = 'P_1 * ( (beta_1 * H_1) + (beta_2 * H_2)  - b)'

    zeros_list = [0 for i in range(len(species_list))]
    symbolic_equations = sp.Matrix(zeros_list)
   
    for idx, eq in enumerate([dH_1, dH_2, dP_1]):
        species_name = species_list[idx]
        symbolic_equations[idx] = sp.sympify(eq, locals=locals())
        diff_eqs[species_name] = symbolic_equations[idx]

    return symbolic_equations, species_list, diff_eqs


def van_der_pol():
    pass

class Model:
    def __init__(self, model_idx, equations_func):
        self.idx = model_idx
        self.symbolic_equations, self.species_list, self.diff_eqs = equations_func()
        
        self.jac = self.make_jac()
        self.params_list = self.extract_params()

    def make_diff_eqs(self):
        dH_1 = 'H_1 * (a_1 - (alpha_1 * P_1))'
        dH_2 = 'H_2 * (a_2 - (alpha_2 * P_1))'
        dP_1 = 'P_1 * ( (beta_1 * H_1) + (beta_2 * H_2)  - b)'

        zeros_list = [0 for i in range(len(self.species_list))]
        symbolic_equations = sp.Matrix(zeros_list)
       
        for idx, eq in enumerate([dH_1, dH_2, dP_1]):
            species_name = self.species_list[idx]
            symbolic_equations[idx] = sp.sympify(eq, locals=locals())
            self.diff_eqs[species_name] = symbolic_equations[idx]

        return symbolic_equations


    def make_jac(self):
        J = self.symbolic_equations.jacobian(self.species_list)
        # for idx_i in range(len(self.species_order)):
        #     for idx_j in range(len(self.species_order)):
        #         print(idx_i, idx_j)
        #         print(J[idx_i, idx_j])
        #         print("")
        return J

    def extract_params(self):
        all_params = []

        for eq in self.symbolic_equations:
            free_symbols = eq.free_symbols
            for symbol in free_symbols:
                if str(symbol) not in self.species_list:
                    all_params.append(str(symbol))

        all_params = sorted(list(set(all_params)), key=str.lower) # Alphabetical order!

        return all_params


if __name__ == "__main__":
    x = Model(model_idx=0, equations_func=LV_equations)

    model_list = [x]
    output_dir = './lv_output/'
    cpp_out = Cpp_source_output(model_list)
    cpp_out.write_source_file(output_dir)
    header = Cpp_header_output(model_list)
    header.write_header_file(output_dir)
