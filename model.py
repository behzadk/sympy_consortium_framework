import equation_builder
import sympy
from sympy.printing.cxxcode import cxxcode
import numpy as np
import pandas as pd
import csv
import species
import os
import utils


class Model:
    def __init__(self, model_idx, strain_list):
        self.idx = model_idx
        self.strains = strain_list

        self.strain_ids = self.get_strain_species()
        self.AHL_ids = self.get_AHL_species()
        self.microcin_ids = self.get_microcin_species()

        self.substrate_ids = self.get_substrate_species()

        self.all_ids = self.microcin_ids + self.AHL_ids + self.substrate_ids + self.strain_ids

        self.diff_eqs = {}
        self.symbolic_equations = None
        self.jac = None

        self.params_list = []
        self.species_list = []

        # Model at this stage can be split into terms
        # Substrate rhs
        # Cell number rhs
        # AHL rhs

    # edges show interactions between species from j (column) to i (row)
    def generate_adjacency_matrix(self, max_sub, max_AHL, max_mic, max_strains):
        # Cell1, Cell2, M1, M2, AHL1, AHL2,
        total_species = max_AHL + max_mic + max_strains + max_sub
        # total_species = len(self.AHL_ids) + len(self.microcin_ids) + len(self.strain_ids)
        adjacency_matrix = np.zeros([total_species, total_species])

        strain_init_idx = 0
        substrate_init_idx = strain_init_idx + max_strains
        microcin_init_idx = substrate_init_idx + max_sub
        AHL_init_idx = microcin_init_idx + max_mic


        all_microcin_objects = []
        all_microcin_ids = []
        all_AHLs = []
        all_substrates = []

        # Collect species objects
        for strain in self.strains:
            for microcin in strain.microcins:
                if microcin not in all_microcin_objects:
                    all_microcin_objects.append(microcin)

                if microcin.id not in all_microcin_ids:
                    all_microcin_ids.append(microcin.id)

            for AHL in strain.AHLs:
                if AHL not in all_AHLs:
                    all_AHLs.append(AHL)

            for s_dependence in strain.substrate_dependences:
                if s_dependence not in all_substrates:
                    all_substrates.append(s_dependence)

            for s_production in strain.substrate_production:
                if s_production not in all_substrates:
                    all_substrates.append(s_production)


        # Fill strain sensitivities to microcin sensitivity is a negative interaction from
        # microcin to strain. (i strain row, j mic col )
        for idx_strain, strain in enumerate(self.strains):
            for sens in strain.sensitivities:
                for idx_mic, mic_id in enumerate(all_microcin_ids):
                    if sens == mic_id:
                        from_node = idx_mic + microcin_init_idx
                        to_node = idx_strain + strain_init_idx

                        adjacency_matrix[to_node, from_node] = -1
        

        # Fill strain substrate dependency and production
        for idx_strain, strain in enumerate(self.strains):
            # Dependency
            for sub_dependence in strain.substrate_dependences:
                sub_indexs = [all_substrates.index(x) for i, x in enumerate(all_substrates) if x == sub_dependence]
                for s_idx in sub_indexs:
                    from_node = s_idx + substrate_init_idx
                    to_node = idx_strain + strain_init_idx

                    adjacency_matrix[to_node, from_node] = 1


            # Production
            for sub_production in strain.substrate_production:
                sub_indexs = [all_substrates.index(x) for i, x in enumerate(all_substrates) if x == sub_production]
                for s_idx in sub_indexs:
                    from_node = idx_strain + strain_init_idx
                    to_node = s_idx + substrate_init_idx

                    adjacency_matrix[to_node, from_node] = 1


        # Fill strain production of microcin, AHL and substrate
        for idx_strain, strain in enumerate(self.strains):
            # Microcin production
            for idx_strain_mic, strain_mic in enumerate(strain.microcins):
                mic_produced_idx = [all_microcin_ids.index(x.id) for i, x in enumerate(all_microcin_objects) if x == strain_mic]

                for idx_mic in mic_produced_idx:
                    from_node = idx_strain + strain_init_idx
                    to_node = idx_mic + microcin_init_idx

                    adjacency_matrix[to_node, from_node] = 1


            # AHL production
            for idx_strain_AHL, strain_AHL in enumerate(strain.AHLs):
                AHL_produced_idx = [i for i, x in enumerate(all_AHLs) if x == strain_AHL]
                for idx_AHL in AHL_produced_idx:
                    from_node = idx_strain + strain_init_idx
                    to_node = idx_AHL + AHL_init_idx
                    adjacency_matrix[to_node, from_node] = 1

            # AHL mic interactions. from AHL to mic
            for mic in all_microcin_objects:
                idx_mic = all_microcin_ids.index(mic.id)

                # Repressors
                try:#
                    if mic.AHL_repressors is np.nan:
                        continue

                    repressor = mic.AHL_repressors[0]
                    AHL_repressor_idx = [i for i, x in enumerate(all_AHLs) if x == repressor]

                    for idx_AHL in AHL_repressor_idx:
                        from_node = idx_AHL + AHL_init_idx
                        to_node = idx_mic + microcin_init_idx

                        adjacency_matrix[to_node, from_node] = -1

                except(IndexError):
                    pass

                # Inducers
                try:
                    if mic.AHL_repressors is np.nan:
                        continue

                    inducer = mic.AHL_inducers[0]
                    AHL_inducer_idx = [i for i, x in enumerate(all_AHLs) if x == inducer]
                    for idx_AHL in AHL_inducer_idx:
                        from_node = idx_AHL + AHL_init_idx
                        to_node = idx_mic + microcin_init_idx

                        adjacency_matrix[to_node, from_node] = 1

                except(IndexError):
                    pass

        self.adjacency_matrix = adjacency_matrix

    def get_strain_species(self):
        strain_id_list = []
        for strain in self.strains:
            strain_id_list.append(strain.id)

        return list(set(strain_id_list))

    def get_microcin_species(self):
        microcin_id_list = []
        for strain in self.strains:
            strain_microcins = strain.microcins
            for m in strain_microcins:
                microcin_id_list.append(m.id)

        return list(set(microcin_id_list))

    def get_AHL_species(self):
        AHL_id_list = []
        for strain in self.strains:
            strain_AHLs = strain.AHLs
            for a in strain_AHLs:
                AHL_id_list.append(a.id)

        return list(set(AHL_id_list))

    def get_substrate_species(self):
        substrate_id_list = []
        for strain in self.strains:
            strain_susbtrates = strain.substrate_dependences
            for s in strain_susbtrates:
                substrate_id_list.append(s.id)

        for strain in self.strains:
            strain_susbtrates = strain.substrate_production
            for s in strain_susbtrates:
                substrate_id_list.append(s.id)

        return list(set(substrate_id_list))

    def is_legal(self):
        required_microcin = []
        required_AHL = []
        required_sub = []

        # Legal requirements
        for s in self.strains:
            required_sub += s.substrate_dependences

            for m in s.microcins:
                if m.AHL_inducers is not np.nan:
                    required_AHL += m.AHL_inducers
                if m.AHL_inducers is not np.nan:
                    required_AHL += m.AHL_repressors

            for m_sens in s.sensitivities:
                required_microcin += [m_sens]

        for a in required_AHL:
            if a.id not in self.AHL_ids:
                return False

        for m in required_microcin:
            if m not in self.microcin_ids:
                return False

        # Remove redundancies
        # Remove models where AHL has no action
        required_AHL_ids = [a.id for a in required_AHL]
        for a_expressed in self.AHL_ids:
            if a_expressed not in required_AHL_ids:
                return False

        # Expression of microcin which no strain is sensitive to
        system_sensitivities = []
        for s in self.strains:
            system_sensitivities = system_sensitivities + s.sensitivities

        for m in self.microcin_ids:
            if m not in system_sensitivities:
                return False

        # Remove models where a produced substrate is not consumed
        for strain in self.strains:
            for sub in strain.substrate_production:
                if sub not in required_sub:
                    return False

        return True

    def build_equations(self):

        # For each strain
        for n in self.strains:
            dN_dt = equation_builder.gen_strain_growth_diff(n.id, self.strains)
            self.diff_eqs.update(dN_dt)

        # For each substrate
        for s in self.substrate_ids:
            dS_dt = equation_builder.gen_diff_eq_substrate(s, self.strains)
            self.diff_eqs.update(dS_dt)

        # For each microcin
        for b in self.microcin_ids:
            dB_dt = equation_builder.gen_microcin_diff_eq(b, self.strains)
            self.diff_eqs.update(dB_dt)

        # For each AHL
        for a in self.AHL_ids:
            dA_dt = equation_builder.gen_AHL_diff_eq(a, self.strains)
            self.diff_eqs.update(dA_dt)

    def build_jacobian(self):
        species_names = list(self.diff_eqs.keys())
        order = sympy.symbols(species_names)
        J = self.symbolic_equations.jacobian(order)

        self.jac = J

    def build_symbolic_equations(self):
        species_names = list(self.diff_eqs.keys())
        order = sympy.symbols(species_names)

        zeros_list = [0 for i in range(len(order))]
        symbolic_equations = sympy.Matrix(zeros_list)

        for idx, eq_key in enumerate(species_names):
            symbolic_equations[idx] = sympy.sympify(self.diff_eqs[eq_key], locals=locals())

        self.symbolic_equations = symbolic_equations


    def extract_species(self):
        self.species_list = list(self.diff_eqs.keys())


    def extract_params(self):
        all_params = []

        for eq in self.symbolic_equations:
            free_symbols = eq.free_symbols
            for symbol in free_symbols:
                if str(symbol) not in self.species_list:
                    all_params.append(str(symbol))

        all_params = sorted(list(set(all_params))) # Alphabetical order!

        self.params_list = all_params

    def config_data(self):
        # N, S, B, A
        for N in self.strains:
            print(N.id)
            for S in N.substrate_dependences:
                print(S.id)

            for B in N.microcins:
                print(B.config_idx)

            for A in N.AHLs:
                print(A.id)

    def write_adj_matrix(self, output_dir, mic_ids, AHL_ids, strain_ids, substrate_ids):
        adj_mat_dir = output_dir + "adj_matricies/"
        utils.make_folder(adj_mat_dir)

        adj_mat_path = adj_mat_dir + 'model_' + str(self.idx) + '_adj_mat.csv'

        new_mic_ids = []
        for idx, i in enumerate(mic_ids):
            new_mic_ids.append("B_" + i)

        new_AHL_ids = []
        for idx, i in enumerate(AHL_ids):
            new_AHL_ids.append("A_" + i)

        new_strain_ids = []
        for idx, i in enumerate(strain_ids):
            new_strain_ids.append("N_" + i)

        new_substrate_ids = []
        for idx, i in enumerate(substrate_ids):
            new_strain_ids.append("S_" + i)


        adj_matrix = self.adjacency_matrix

        with open(adj_mat_path, 'w') as f:
            w = csv.writer(f)
            adj_mat_species = new_strain_ids + new_substrate_ids + new_mic_ids + new_AHL_ids
            header = [None] + adj_mat_species
            w.writerow(header)

            for row_idx in range(len(adj_matrix)):
                w.writerow([adj_mat_species[row_idx]] + adj_matrix[row_idx].tolist())


    ##
    # Writes upper and lower boundaries for uniform priors to a csv, separately for
    # parameters and species using a csv containing default parameters.
    ##
    def write_prior_parameter_dict(self, default_params_path, output_dir):
        sim_inputs_dir = output_dir + "input_files/"
        utils.make_folder(sim_inputs_dir)

        sim_params_path = sim_inputs_dir + 'params_' + str(self.idx) + ".csv"

        model_parameters = self.params_list
        default_params = pd.read_csv(default_params_path)

        prior_dict = {}
        for idx, row in default_params.iterrows():
            p = row['parameter']

            # Adds param to dict if not linked to a particular species
            if p in model_parameters:
                prior_dict[p] = [row['lower_bound'], row['upper_bound']]
                continue

            # Adds param to dict if it is linked to a particulr species, identified by the id tag
            p = row['parameter'] + '_#ID#'
            for id in self.all_ids:
                param_species = p.replace('#ID#', id)
                if param_species in model_parameters:
                    prior_dict[param_species] = [row['lower_bound'], row['upper_bound']]
                    continue


        if len(model_parameters) != len(prior_dict):
            print('mismatch in length of prior dict and model parameters.', len(prior_dict), len(model_parameters))


        with open(sim_params_path, 'w') as f:  # Just use 'w' mode in 3.x
            w = csv.writer(f)
            for key, value in prior_dict.items():
                w.writerow([key, value[0], value[1]])


    def write_init_species_dict(self, default_init_species_path, output_dir):
        sim_inputs_dir = output_dir + "input_files/"
        utils.make_folder(sim_inputs_dir)

        sim_species_path = sim_inputs_dir + 'species_' + str(self.idx) + ".csv"

        model_species = self.species_list
        default_species = pd.read_csv(default_init_species_path)

        prior_dict = {}
        for idx, row in default_species.iterrows():
            p = row['species']

            # Adds param to dict if not linked to a particular species
            if p in model_species:
                prior_dict[p] = [row['lower_bound'], row['upper_bound']]
                continue

            # Adds param to dict if it is linked to a particulr species, identified by the id tag
            p = row['species'] + '_#ID#'
            for id in self.all_ids:
                param_species = p.replace('#ID#', id)

                if param_species in model_species:
                    prior_dict[param_species] = [row['lower_bound'], row['upper_bound']]
                    continue

        if len(model_species) != len(prior_dict):
            print('mismatch in length of prior dict and model parameters.', len(prior_dict), len(model_species))
            print(prior_dict)
            print("")
            print(model_species)
            exit()



        with open(sim_species_path, 'w') as f:  # Just use 'w' mode in 3.x
            w = csv.writer(f)
            for key, value in prior_dict.items():
                w.writerow([key, value[0], value[1]])


    def write_python_equations(self, output_path):
        py_eqs_dir = output_path + "py_eqs_txt_files/"
        utils.make_folder(py_eqs_dir)

        model_py_eqs_path = py_eqs_dir + "model_" + str(self.idx) + "_eqs.txt"

        with open(model_py_eqs_path, 'w') as txt_file:  # Just use 'w' mode in 3.x
            for eq in self.diff_eqs:
                res = eq + " = " + str(self.diff_eqs[eq]) + "\n\n"
                res = res.replace("^", "**")
                txt_file.write(res)
