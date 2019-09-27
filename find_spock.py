import pandas as pd
import numpy as np
import glob

def main():
	adj_mat_dir = "./output/input_files_two_species_spock_manu_1/adj_matricies/"
	adj_mat_files = glob.glob(adj_mat_dir + '*.csv')

	for f in adj_mat_files:
		df = pd.read_csv(f)
		header = list(df)[1:]
		mat = df.values
		mat = mat[:, 1:]

		# Is producing AHL
		A_1_idx = header.index("A_1")
		B_1_idx = header.index("B_1")
		N_1_idx = header.index("N_1")
		N_2_idx = header.index("N_2")
		I_1_idx =  header.index("I_1")
		V_1_idx =  header.index("V_1T")
		T_1_idx = header.index("T_1T")

		# mat[To][from]
		is_spock = True
		# No immunity
		if not mat[I_1_idx][N_1_idx] == 1:
			is_spock = False

		# Constitutive expression
		if not mat[I_1_idx][A_1_idx] == 0:
			is_spock = False

		# Produces bacteriocin
		if not mat[B_1_idx][N_1_idx] == 1:
			is_spock = False

		# Kills N1
		if not mat[N_1_idx][B_1_idx] == -1:
			is_spock = False

		# Kills N2
		if not mat[N_2_idx][B_1_idx] == -1:
			is_spock = False

		# Produces AHL
		if not mat[A_1_idx][N_1_idx] == 1:
			is_spock = False

		# Represses bacteriocin
		if not mat[B_1_idx][A_1_idx] == -1:
			is_spock = False

		if sum(mat[:, N_1_idx]) >  3:
			is_spock = False

		if is_spock:
			print(header)
			print(mat)
			print(f)
			# exit()


if __name__ == "__main__":
	main()