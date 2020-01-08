import pandas as pd
import numpy as np
import glob

def find_spock():
	adj_mat_dir = "./output/input_files_two_species_spock_manu_0/adj_matricies/"
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


		# mat[To][from]
		is_spock = True
		# No immunity
		if not mat[I_1_idx][N_1_idx] == 1:
			is_spock = False

		if not mat[I_1_idx][A_1_idx] == 0:
			is_spock = False

		if not mat[B_1_idx][N_1_idx] == 1:
			is_spock = False

		if not mat[N_1_idx][B_1_idx] == -1:
			is_spock = False

		if not mat[N_2_idx][B_1_idx] == -1:
			is_spock = False

		if not mat[A_1_idx][N_1_idx] == 1:
			is_spock = False

		if not mat[B_1_idx][A_1_idx] == -1:
			is_spock = False

		if is_spock:
			print(header)
			print(mat)
			print(f)
			exit()


def find_this_thing():
	adj_mat_dir = "./output/input_files_three_species_1/adj_matricies/"
	adj_mat_files = glob.glob(adj_mat_dir + '*.csv')

	for f in adj_mat_files:
		df = pd.read_csv(f)
		header = list(df)[1:]
		mat = df.values
		mat = mat[:, 1:]

		# Is producing AHL
		A_1_idx = header.index("A_1")
		A_2_idx = header.index("A_2")

		B_1_idx = header.index("B_1")
		B_2_idx = header.index("B_2")
		B_3_idx = header.index("B_3")

		N_1_idx = header.index("N_1")
		N_2_idx = header.index("N_2")
		N_3_idx = header.index("N_3")

		is_thing = True

		if not sum(mat[:, B_1_idx]) < 0:
			is_thing = False

		if not sum(mat[:, B_2_idx]) < 0:
			is_thing = False

		if not sum(mat[:, B_3_idx]) < 0:
			is_thing = False

		if is_thing:
			print(header)
			print(mat)
			print(f)
			exit()







if __name__ == "__main__":
	find_this_thing()