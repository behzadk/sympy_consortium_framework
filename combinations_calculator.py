import math
import numpy as np
import itertools

def spock_combs():
	cst = 0
	ind = 1
	rpr = 2
	no_exp = 3

	B = [0, 1, 2, 3]
	I = [0, 1, 2, 3]
	T = [0, 1, 2, 3]
	V = [1, 2, 3]
	Q = [0, 3]


	numerical_sol = 4**4 * 2
	print("unique combs", numerical_sol)

	Q_not_exp = 2**4

	print("Q not expressed", Q_not_exp)
	numerical_sol = numerical_sol - Q_not_exp
	print("num_sol", numerical_sol)
	print("")

	I_exp = 3**2 * 4**2 * 2
	print("I expressed", I_exp)
	numerical_sol = numerical_sol - I_exp
	print("num_sol", numerical_sol)

	combs = []
	for b in B:
		for i in I:
			for t in T:
				for v in V:
					for q in Q:
						combs.append([b, i, t, v, q])

	print(len(combs))
	C_1 = []
	C_2 = []
	C_3 = []
	C_4 = []
	C_5 = []

	invalid = []

	for comb in combs:
		failed = False
		is_expressed = [0, 1, 2]
		is_regulated = [1, 2]

		# Inducer not expressed means other species can only be constitutive or unexpressed
		no_inducer = [0, 3]
		if comb[4] == no_exp:
			C_1_criteria = [True if i in no_inducer else False for i in comb[0:4] ]
			if all(C_1_criteria):
				C_1.append(comb)

			else:
				failed = True

		# Immunity expression requires bacteriocin expression
		if comb[1] in is_expressed:
			C_2_criteria = [True if comb[0] in is_expressed else False]
			if all(C_2_criteria):
				C_2.append(comb)

			else:
				failed = True

		# Antitoxin expression requires toxin expression
		if comb[3] in is_expressed:
			C_3_criteria = [True if comb[2] in is_expressed else False]
			if all(C_3_criteria):
				C_3.append(comb)
			else:
				failed = True

		# Quorum expression requires any induced
		if comb[4] in is_expressed:
			C_4_criteria = [True if i in is_regulated else False for i in comb[0:4] ]
			if any(C_4_criteria):
				C_4.append(comb)
			else:
				failed = True

		if comb[4] != no_exp and comb[1] not in is_expressed and comb[3] not in is_expressed and comb[4] not in is_expressed:
			C_5.append(comb)
			print(C_5)

		if failed:
			invalid.append(comb)
			
	print(len(C_1))
	print(len(C_2))
	print(len(invalid))
	print(len(combs) - len(invalid))

if __name__ == "__main__":
	spock_combs()