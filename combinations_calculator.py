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
	V = [0, 1, 2, 3]
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

	for comb in combs:

	print("Part combinations: ", len(combs))



if __name__ == "__main__":
	spock_combs()