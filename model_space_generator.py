import species
import pandas as pd
import itertools
import model
import numpy as np

##
# Generates microcin objects for all combinations, given a list of AHL objects and information on
# if microcin can be induced, repressed, or constitutively expressed. Currently assumes only one AHl can mediate
# expression at any one time.
#
# @param microcin_ids - List of possible ids for microcin, a proxy for the number of independent microcin
# @param AHL_objects - A list containing AHL objects
# @param microcin_induced - Determines whether microcin can be induced by AHL
# @param microcin_repressed - Determines whether microcin can be repressed by AHL
# @param microcin_constitutive - Determines whether microcin can be constitutively expressed
##
def generate_microcin_combinations(microcin_ids, AHL_objects, microcin_induced=False,
                                   microcin_repressed=False, microcin_constitutive=False):
    microcin_config_idx = 0
    microcin_objects = []
    microcin_config_data = []  # Microcin reference contains - config index, microcin_id, inducer_id, repressor_id, constitutive

    if microcin_induced is True:
        for AHL in AHL_objects:
            for m_id in microcin_ids:
                microcin_objects.append(species.Microcin(microcin_config_idx, m_id, [AHL], []))
                config_data = [microcin_config_idx, m_id, AHL.id, np.nan]
                microcin_config_data.append(config_data)
                microcin_config_idx += 1

    if microcin_repressed is True:
        for AHL in AHL_objects:
            for m_id in microcin_ids:
                microcin_objects.append(species.Microcin(microcin_config_idx, m_id, [], [AHL]))

                config_data = [microcin_config_idx, m_id, np.nan, AHL.id]
                microcin_config_data.append(config_data)

                microcin_config_idx += 1

    if microcin_constitutive is True:
        for m_id in microcin_ids:
            microcin_objects.append(species.Microcin(microcin_config_idx, m_id, np.nan, np.nan, constitutive_expression=True))
            microcin_config_idx += 1

    microcin_config_df = pd.DataFrame(columns=["microcin_idx", "microcin_id", "inducer_id", "repressor_id"],
                                      data=microcin_config_data)

    return microcin_objects, microcin_config_df


##
# Generates strain objects for all combinations, given a list of microcin objects and substrate objects
##
class model_space():
    def __init__(self, strain_ids, microcin_objects, AHL_objects, substrate_objects,
                 max_microcin_parts, max_AHL_parts, max_substrate_dependencies, max_microcin_sensitivities=1):
        self.strain_ids = strain_ids
        self.strain_objects = []
        self.microcin_objects = microcin_objects
        self.AHL_objects = AHL_objects
        self.substrate_objects = substrate_objects

        self.microcin_ids = list(set([m.id for m in microcin_objects]))

        self.part_combinations = []
        self.models_list = []

        # Maximum parts for each strain
        self.max_microcin_parts = max_microcin_parts
        self.max_AHL_parts = max_AHL_parts
        self.max_substrate_dependencies = max_substrate_dependencies
        self.max_microcin_sensitivities = max_microcin_sensitivities

    def generate_part_combinations(self, strain_max_microcin, strain_max_AHL, strain_max_sub, strain_max_microcin_sens):
        # Construct possible combinations for each part. [None] added to AHL production, microcin production and
        # microcin sensitivity represent empty part. The purpose of this is so we generate combinations with one or more
        # of each part.

        # Consider making this into a lambda function
        microcin_production_lists = [list(i for i in m if i != None) for m in
                                     itertools.combinations(self.microcin_objects + [None], strain_max_microcin)]

        AHL_production_lists = [list(i for i in a if i != None) for a in
                                itertools.combinations(self.AHL_objects + [None], strain_max_AHL)]

        substrate_dependencies_list = [list(i for i in s if i != None) for s in
                                       itertools.combinations(self.substrate_objects, strain_max_sub)]

        microcin_sensitivities_list = [list(i for i in m_id if i != None) for m_id in
                                       itertools.combinations(self.microcin_ids + [None], strain_max_microcin_sens)]

        # Append empty list representing no production or sensitivity, only necessarry if less than two max parts
        if strain_max_AHL > 1: AHL_production_lists.append([])
        if strain_max_microcin > 1: microcin_production_lists.append([])
        if strain_max_microcin_sens > 1: microcin_sensitivities_list.append([])


        # Generate all different combinations of parts
        for idx_m, m in enumerate(microcin_production_lists):
            for idx_a, a in enumerate(AHL_production_lists):
                for idx_s, s in enumerate(substrate_dependencies_list):
                    for idx_sensi, sensi in enumerate(microcin_sensitivities_list):
                        self.part_combinations.append([m, a, s, sensi])


        return self.part_combinations

    def remove_symmetries(self):
        all_adj_mats = []
        keep_idx_0 = []

        # Check for direct matches
        for idx, model in enumerate(self.models_list):
            if any(np.array_equal(x, model.adjacency_matrix) for x in all_adj_mats):
                continue


            else:
                keep_idx_0.append(idx)
                all_adj_mats.append(model.adjacency_matrix)

        strain_init_idx = 0
        microcin_init_idx = strain_init_idx + len(self.strain_ids)
        AHL_init_idx = microcin_init_idx + self.max_microcin_parts

        strains_index_range = range(strain_init_idx, len(self.strain_ids))
        mic_index_range = range(microcin_init_idx, microcin_init_idx + self.max_microcin_parts)
        AHL_index_range = range(AHL_init_idx, AHL_init_idx + self.max_AHL_parts)

        # Flip strain columns and remove symmetrical strains
        keep_idx_1 = []
        for idx in keep_idx_0:
            permute_strains = list(itertools.permutations(strains_index_range))
            original_config = permute_strains[0]
            model_adj = self.models_list[idx].adjacency_matrix

            match = False
            for perm in permute_strains[1:]:
                new_adj = model_adj.copy()
                new_adj.T[[original_config]] = new_adj.T[[perm]]
                new_adj[[original_config]] = new_adj[[perm]]

                if any(np.array_equal(x, new_adj) for x in all_adj_mats):
                    match = True

            if match is False:
                keep_idx_1.append(idx)



        return keep_idx_1


    def generate_models(self):
        model_idx = 0

        system_combinations = itertools.combinations(self.part_combinations, len(self.strain_ids))
        total_sys = 0

        for sys in system_combinations:
            model_strains = []
            for idx, N_id in enumerate(self.strain_ids):
                new_strain = species.Strain(N_id, *sys[idx])
                model_strains.append(new_strain)

            new_model = model.Model(model_idx, model_strains)


            new_model.generate_adjacency_matrix(self.max_AHL_parts, self.max_microcin_parts, len(self.strain_ids))
            total_sys +=1

            if new_model.is_legal():
                self.models_list.append(new_model)
                model_idx += 1
        print("Number of systems: ", total_sys)

        keep_idx = self.remove_symmetries()
        self.models_list = [self.models_list[i] for i in keep_idx]

        # reset model indexes
        model_idx = 0
        for m in self.models_list:
            m.idx = model_idx
            model_idx +=1

        return self.models_list

    def generate_model_reference_table(self, max_microcin_parts, max_AHL_parts,
                                   max_substrate_dependencies, max_microcin_sensitivities):
        # Make column headers

        cell_prefix = 'cell_IDX_'
        microcin_prefix = 'M_IDX'
        AHL_prefix = 'AHL_IDX'
        substrate_prefix = 'S_IDX'
        sensitivity_prefix = 'Sens_IDX'

        models_datasheet_cols = []

        for cell_idx, i in enumerate(self.strain_ids):
            for m_idx, m in enumerate(range(max_microcin_parts)):
                cell = cell_prefix.replace('IDX', str(cell_idx))
                j = microcin_prefix.replace('IDX', str(m_idx))
                models_datasheet_cols.append(cell + j)

            for a_idx, a in enumerate(range(max_AHL_parts)):
                cell = cell_prefix.replace('IDX', str(cell_idx))
                j = AHL_prefix.replace('IDX', str(a_idx))
                models_datasheet_cols.append(cell + j)

            for s_idx, s in enumerate(range(max_substrate_dependencies)):
                cell = cell_prefix.replace('IDX', str(cell_idx))
                j = substrate_prefix.replace('IDX', str(s_idx))
                models_datasheet_cols.append(cell + j)

            for sens_idx, sens in enumerate(range(max_microcin_sensitivities)):
                cell = cell_prefix.replace('IDX', str(cell_idx))
                j = sensitivity_prefix.replace('IDX', str(sens_idx))
                models_datasheet_cols.append(cell + j)

        models_datasheet = pd.DataFrame(columns=models_datasheet_cols)

        # NaN's used to fill an empty cell, where no action takes place
        for model_idx, model in enumerate(self.models_list):
            model_data = []
            for strain_idx, strain in enumerate(model.strains):

                for idx_m, m in enumerate(range(max_microcin_parts)):
                    try:
                        model_data.append(strain.microcins[idx_m].config_idx)
                    except(IndexError):
                        model_data.append(np.nan)

                for idx_a, a in enumerate(range(max_AHL_parts)):
                    try:
                        model_data.append(strain.AHLs[idx_a].id)
                    except(IndexError):
                        model_data.append(np.nan)

                for idx_s, s in enumerate(range(max_substrate_dependencies)):
                    try:
                        model_data.append(strain.substrate_dependences[idx_s].id)
                    except(IndexError):
                        model_data.append(np.nan)

                for idx_sens, sens in enumerate(range(max_microcin_sensitivities)):
                    try:
                        model_data.append(strain.sensitivities[idx_sens])

                    except(IndexError):
                        model_data.append(np.nan)

            models_datasheet.loc[model.idx] = model_data

        return models_datasheet


