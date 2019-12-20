from species import Microcin
from species import  AHL
from species import Strain
from species import Substrate

from model import Model
from cpp_output import Cpp_source_output
from cpp_output import Cpp_header_output
from species import Substrate
from species import AHL
from species import Microcin
from species import Strain
from species import Immunity

from model import Model
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


def single_strain_test(default_params_path, default_init_species_path, output_dir):
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
    strain_ids = ['1']

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
                                                                                                 microcin_induced=False,
                                                                                                 microcin_repressed=False, microcin_constitutive=True)
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
                                                                                                 immunity_induced=False,
                                                                                                 immunity_repressed=False, immunity_constitutive=False)


    model_space = model_space_generator.model_space(strain_ids, microcin_objects,
                                                    AHL_objects, substrate_objects, antitoxin_objects, immunity_objects, toxin_objects,
                                                    max_microcin_parts, max_AHL_parts,
                                                    max_substrate_parts, max_antitoxins, max_immunity_parts, max_toxin_parts, max_microcin_sensitivities=1)

    part_combos = model_space.generate_part_combinations(
        strain_max_microcin=1, strain_max_AHL=1, strain_max_sub_dependencies=1, 
        strain_max_microcin_sens=1, strain_max_sub_production=0, strain_max_antitoxin=1, 
        strain_max_immunity=1, strain_max_toxin=1
        )

    print("Generating model list")
    model_list = model_space.generate_models()
    print("generating reference table")

    # models_ref_df = model_space.generate_model_reference_table(max_microcin_parts, max_AHL_parts,
    #                                                            max_substrate_parts, max_microcin_parts)

    model_space.remove_symmetries()
    model_space.reset_model_indexes()

    model_list = model_space.models_list

    generate_adjacency_matricies(model_list, substrate_ids, microcin_ids, AHL_ids, strain_ids, antitoxin_ids, immunity_ids, toxin_ids, output_dir)
    generate_simulation_files(model_list, default_params_path, default_init_species_path, output_dir)


def spock_manu_no_symm(default_params_path, default_init_species_path, output_dir):
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


def two_species_no_symm(default_params_path, default_init_species_path, output_dir):
    # Set species IDs
    substrate_ids = ['glu']
    S_glu = Substrate(substrate_ids[0])
    substrate_objects = [S_glu]

    AHL_ids = ['1', '2']
    AHL_1 = AHL(AHL_ids[0])
    AHL_2 = AHL(AHL_ids[1])

    AHL_objects = [AHL_1, AHL_2]

    microcin_ids = ['1', '2']
    strain_ids = ['1', '2']

    toxin_ids = []
    antitoxin_ids = [] # Name of antitoxins must match name of toxin.(B_#ID# format)
    immunity_ids = []

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
                                                                                                 immunity_induced=False,
                                                                                                 immunity_repressed=False, immunity_constitutive=False)


    model_space = model_space_generator.model_space(strain_ids, microcin_objects,
                                                    AHL_objects, substrate_objects, antitoxin_objects, immunity_objects, toxin_objects,
                                                    max_microcin_parts, max_AHL_parts,
                                                    max_substrate_parts, max_antitoxins, max_immunity_parts, max_toxin_parts, max_microcin_sensitivities=2)

    part_combos = model_space.generate_part_combinations(
        strain_max_microcin=1, strain_max_AHL=1, strain_max_sub_dependencies=1, 
        strain_max_microcin_sens=2, strain_max_sub_production=0, strain_max_antitoxin=1, 
        strain_max_immunity=0, strain_max_toxin=0
        )

    print("Generating model list")
    model_list = model_space.generate_models()
    print("generating reference table")
    model_space.remove_symmetries()
    print(len(model_space.models_list))

    model_space.reset_model_indexes()

    model_list = model_space.models_list

    model_space_generator.generate_adjacency_matricies(model_list, substrate_ids, microcin_ids, AHL_ids, strain_ids, antitoxin_ids, immunity_ids, toxin_ids, output_dir)
    model_space_generator.generate_simulation_files(model_list, default_params_path, default_init_species_path, output_dir)

    return model_space


def two_species_no_symm_auxos(default_params_path, default_init_species_path, output_dir):
    # Set species IDs
    substrate_ids = ['glu', 'aa1', 'aa2']
    S_glu = Substrate(substrate_ids[0])
    S_trp = Substrate(substrate_ids[1])
    S_val = Substrate(substrate_ids[2])

    substrate_objects = [S_glu, S_trp, S_val]

    AHL_ids = ['1', '2']
    AHL_1 = AHL(AHL_ids[0])
    AHL_2 = AHL(AHL_ids[1])

    AHL_objects = [AHL_1, AHL_2]

    microcin_ids = ['1', '2']
    strain_ids = ['1', '2']

    toxin_ids = []
    antitoxin_ids = [] # Name of antitoxins must match name of toxin.(B_#ID# format)
    immunity_ids = []

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
                                                                                                 microcin_induced=False,
                                                                                                 microcin_repressed=False, microcin_constitutive=True)
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
                                                                                                 immunity_induced=False,
                                                                                                 immunity_repressed=False, immunity_constitutive=False)


    model_space = model_space_generator.model_space(strain_ids, microcin_objects,
                                                    AHL_objects, substrate_objects, antitoxin_objects, immunity_objects, toxin_objects,
                                                    max_microcin_parts, max_AHL_parts,
                                                    max_substrate_parts, max_antitoxins, max_immunity_parts, max_toxin_parts, max_microcin_sensitivities=2)

    part_combos = model_space.generate_part_combinations(
        strain_max_microcin=1, strain_max_AHL=0, strain_max_sub_dependencies=3, 
        strain_max_microcin_sens=2, strain_max_sub_production=3, strain_max_antitoxin=0, 
        strain_max_immunity=0, strain_max_toxin=0
        )

    clean_part_combos = []
    for combo in part_combos:
        keep = True

        # Does not produce glucose
        for sub in combo[4]:
            if sub.id == "glu":
                keep = False

        # Must be dependent on glucose
        if "glu" not in [sub.id for sub in combo[2]]:
            keep = False

        if keep:
            clean_part_combos.append(combo)

        print("")

    model_space.part_combinations = clean_part_combos

    print("Generating model list")
    model_list = model_space.generate_models()
    model_space.aux_filter()

    print("generating reference table")
    print(len(model_space.models_list))

    model_space.remove_symmetries()
    print(len(model_space.models_list))

    model_space.reset_model_indexes()

    model_list = model_space.models_list

    model_space_generator.generate_adjacency_matricies(model_list, substrate_ids, microcin_ids, AHL_ids, strain_ids, antitoxin_ids, immunity_ids, toxin_ids, output_dir)
    model_space_generator.generate_simulation_files(model_list, default_params_path, default_init_species_path, output_dir)

    return model_space

def three_species_no_symm_auxos(default_params_path, default_init_species_path, output_dir):
    # Set species IDs
    substrate_ids = ['glu', 'aa1', 'aa2', 'aa3']
    S_glu = Substrate(substrate_ids[0])
    S_trp = Substrate(substrate_ids[1])
    S_val = Substrate(substrate_ids[2])

    substrate_objects = [S_glu, S_trp, S_val]

    AHL_ids = ['1', '2']
    AHL_1 = AHL(AHL_ids[0])
    AHL_2 = AHL(AHL_ids[1])

    AHL_objects = [AHL_1, AHL_2]

    microcin_ids = ['1', '2', '3']
    strain_ids = ['1', '2', '3']

    toxin_ids = []
    antitoxin_ids = [] # Name of antitoxins must match name of toxin.(B_#ID# format)
    immunity_ids = []

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
                                                                                                 microcin_induced=False,
                                                                                                 microcin_repressed=False, microcin_constitutive=True)
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
                                                                                                 immunity_induced=False,
                                                                                                 immunity_repressed=False, immunity_constitutive=False)


    model_space = model_space_generator.model_space(strain_ids, microcin_objects,
                                                    AHL_objects, substrate_objects, antitoxin_objects, immunity_objects, toxin_objects,
                                                    max_microcin_parts, max_AHL_parts,
                                                    max_substrate_parts, max_antitoxins, max_immunity_parts, max_toxin_parts, max_microcin_sensitivities=3)

    part_combos = model_space.generate_part_combinations(
        strain_max_microcin=1, strain_max_AHL=0, strain_max_sub_dependencies=4, 
        strain_max_microcin_sens=3, strain_max_sub_production=1, strain_max_antitoxin=0, 
        strain_max_immunity=0, strain_max_toxin=0
        )

    clean_part_combos = []
    for combo in part_combos:
        keep = True

        # Does not produce glucose
        for sub in combo[4]:
            if sub.id == "glu":
                keep = False

        # Must be dependent on glucose
        if "glu" not in [sub.id for sub in combo[2]]:
            keep = False

        if keep:
            clean_part_combos.append(combo)

        print("")

    model_space.part_combinations = clean_part_combos

    print("Generating model list")
    model_list = model_space.generate_models()
    model_space.aux_filter()

    print("generating reference table")
    print(len(model_space.models_list))

    model_space.remove_symmetries()
    print(len(model_space.models_list))


    model_space.reset_model_indexes()

    model_list = model_space.models_list

    model_space_generator.generate_adjacency_matricies(model_list, substrate_ids, microcin_ids, AHL_ids, strain_ids, antitoxin_ids, immunity_ids, toxin_ids, output_dir)
    model_space_generator.generate_simulation_files(model_list, default_params_path, default_init_species_path, output_dir)

    return model_space


def three_species_no_symm(default_params_path, default_init_species_path, output_dir):
    # Set species IDs
    substrate_ids = ['glu']
    S_glu = Substrate(substrate_ids[0])
    substrate_objects = [S_glu]

    AHL_ids = ['1', '2']
    AHL_1 = AHL(AHL_ids[0])
    AHL_2 = AHL(AHL_ids[1])

    AHL_objects = [AHL_1, AHL_2]

    microcin_ids = ['1', '2', '3']
    strain_ids = ['1', '2', '3']

    toxin_ids = []
    antitoxin_ids = [] # Name of antitoxins must match name of toxin.(B_#ID# format)
    immunity_ids = []

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
                                                                                                 immunity_induced=False,
                                                                                                 immunity_repressed=False, immunity_constitutive=False)


    model_space = model_space_generator.model_space(strain_ids, microcin_objects,
                                                    AHL_objects, substrate_objects, antitoxin_objects, immunity_objects, toxin_objects,
                                                    max_microcin_parts, max_AHL_parts,
                                                    max_substrate_parts, max_antitoxins, max_immunity_parts, max_toxin_parts, max_microcin_sensitivities=3)

    part_combos = model_space.generate_part_combinations(
        strain_max_microcin=1, strain_max_AHL=1, strain_max_sub_dependencies=1, 
        strain_max_microcin_sens=1, strain_max_sub_production=0, strain_max_antitoxin=1, 
        strain_max_immunity=1, strain_max_toxin=1
        )
    print(len(part_combos))

    print("Generating model list")
    model_list = model_space.generate_models()
    print(len(model_space.models_list))
    model_space.max_immunity_filter(1)
    print(len(model_space.models_list))

    model_space.remove_symmetries()
    model_space.reset_model_indexes()

    model_list = model_space.models_list
    print(len(model_space.models_list))

    model_space_generator.generate_adjacency_matricies(model_list, substrate_ids, microcin_ids, AHL_ids, strain_ids, antitoxin_ids, immunity_ids, toxin_ids, output_dir)
    model_space_generator.generate_simulation_files(model_list, default_params_path, default_init_species_path, output_dir)

    return model_space

def three_species_one_pred_two_prey(default_params_path, default_init_species_path, output_dir):
    # Set species IDs
    substrate_ids = ['glu']
    S_glu = Substrate(substrate_ids[0])
    substrate_objects = [S_glu]

    AHL_ids = ['1']
    AHL_1 = AHL(AHL_ids[0])

    AHL_objects = [AHL_1]

    microcin_ids = ['1B', '2B']
    strain_ids = ['1', '2', '3']
    
    immunity_ids = ['1B']
    imm_species = Immunity(config_idx=1, immunity_id=immunity_ids[0], AHL_inducer_list=[AHL_objects[0]], AHL_repressor_list=[], constitutive_expression=False)
    m_species = Microcin(config_idx=1, microcin_id=microcin_ids[0], AHL_inducer_list=[AHL_objects[0]], AHL_repressor_list=[], constitutive_expression=False)

    N_1 = Strain(strain_id='1', microcin_expression=[m_species], AHL_expression=[], 
        substrate_dependences=[substrate_objects[0]], microcin_sensitivities=[microcin_ids[0]], substrate_production=[], 
        antitoxins=[], immunity_expression=[imm_species], toxin_expression=[])

    N_2 = Strain(strain_id='2', microcin_expression=[], AHL_expression=[AHL_1], 
        substrate_dependences=[substrate_objects[0]], microcin_sensitivities=[microcin_ids[0]], substrate_production=[], 
        antitoxins=[], immunity_expression=[imm_species], toxin_expression=[])

    N_3 = Strain(strain_id='3', microcin_expression=[], AHL_expression=[AHL_1], 
        substrate_dependences=[substrate_objects[0]], microcin_sensitivities=[microcin_ids[0]], substrate_production=[], 
        antitoxins=[], immunity_expression=[imm_species], toxin_expression=[])


    model_1 = Model(1, [N_1, N_2, N_3])


    model_list = [model_1]

    for idx, m in enumerate(model_list):
        m.build_equations()
        m.build_symbolic_equations()
        m.build_jacobian()
        m.extract_species()
        m.extract_params()
        m.generate_adjacency_matrix(1, 1, 2, 3, 0, 1, 0)

        print("Model ", idx, "Is legal?: ", m.is_legal())

    model_space_generator.generate_adjacency_matricies(model_list, substrate_ids, microcin_ids, AHL_ids, strain_ids, [], immunity_ids, [], output_dir)
    model_space_generator.generate_simulation_files(model_list, default_params_path, default_init_species_path, output_dir)


    # return model_space
