import sympy
import equation_builder

"""
Following flags for a species type
#N# - strain
#S# - substrate
#B# - bacteriocin
#A# - AHL

The flags are replaced with a index
"""



class Strain:
    def __init__(self, strain_id, substrate_dependence, microcin_expression, AHL_expression):
        self.id = strain_id
        self.substrate_dependence = substrate_dependence
        self.microcins = microcin_expression
        self.AHLs = AHL_expression
        self.diff_eqs = {}


class Microcin:
    def __init__(self, microcin_id, inducer_id_list, repressor_id_list):
        self.id = microcin_id
        self.inducer_id_list = inducer_id_list
        self.repressor_id_list = repressor_id_list


class AHL:
    def __init__(self, AHL_id):
        self.id = AHL_id


class Model:
    def __init__(self, strain_list):
        self.strains = strain_list
        self.microcin_ids = self.get_microcin_species()
        self.AHL_ids = self.get_AHL_species()
        self.substrate_ids = self.get_substrate_species()

        self.diff_eqs = {}

        # Model at this stage can be split into terms
        # Substrate rhs
        # Cell number rhs
        # AHL rhs


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
            strain_susbtrates = strain.substrate_dependence
            for s in strain_susbtrates:
                substrate_id_list.append(s)

        return list(set(substrate_id_list))

    def build_equations(self):

        # For each species
        for s in self.substrate_ids:
            dS_dt = equation_builder.gen_diff_eq_substrate(s, self.strains)
            self.diff_eqs.update(dS_dt)



        print(self.diff_eqs)


def main():
    nc = Strain('c', ['glu'], [], [])
    nx = Strain('x', ['glu'], [], [])
    m = Model([nx, nc])
    m.build_equations()

    print("hello world")
    # kb_ind = funcs['k_b_ind_#B#']
    # kb_ind = kb_ind.replace('#B#', '0')
    # kb_ind = kb_ind.replace('#A#', '1')
    #
    # print(kb_ind)
    # i = sympy.sympify(kb_ind)
    # print(i)


if __name__ == "__main__":
    main()
