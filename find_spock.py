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


def find_another_thing():
	adj_mat_dir = "./output/input_files_three_species_3/adj_matricies/"
	adj_mat_files = glob.glob(adj_mat_dir + '*.csv')

	models = []

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

		N_1_idx = header.index("N_1")
		N_2_idx = header.index("N_2")
		N_3_idx = header.index("N_3")

		is_thing = True

		# N1 produces B1
		if not mat[B_1_idx][N_1_idx] == 1:
			is_thing = False


		# # A1 produced by N3
		# if not mat[A_1_idx][N_3_idx] == 1:
		# 	is_thing = False

		# N2 produces B2
		if not mat[B_2_idx][N_2_idx] == 1:
			is_thing = False

		# N3 produces B2
		if not mat[B_2_idx][N_2_idx] == 1:
			is_thing = False

		# # B3 not produced
		# if not sum(mat[B_3_idx]) == 0:
		# 	is_thing = False

		# N1 sensitive to B2
		if not mat[N_1_idx][B_2_idx] == -1:
			is_thing = False

		# N2 sensitive to B1
		if not mat[N_2_idx][B_1_idx] == -1:
			is_thing = False

		# N3 sensitive to B1
		if not mat[N_3_idx][B_1_idx] == -1:
			is_thing = False

		# # N3 sensitive to B1
		if not mat[N_3_idx][B_1_idx] == -1:
			is_thing = False

		# B1 induced by A1
		if not mat[B_1_idx][A_1_idx] == 1:
			is_thing = False

		# B2 induced by A2
		if not mat[B_2_idx][A_2_idx] == 1:
			is_thing = False

		# A2 produced by N2
		if not mat[A_2_idx][N_2_idx] == 1:
			is_thing = False

		# A2 produced by N3
		if not mat[A_2_idx][N_3_idx] == 1:
			is_thing = False

		# A1 produced by N1
		if not mat[A_1_idx][N_1_idx] == 1:
			is_thing = False


		# # A1 not produced
		# if not sum(mat[A_1_idx]) == 0:
		# 	is_thing = False

		# # A2 not produced
		# if not sum(mat[A_2_idx]) == 0:
		# 	is_thing = False


		if is_thing:
			print(header)
			print(mat)
			print(mat[B_1_idx])
			print(mat[B_1_idx][N_1_idx])
			print(f)
			print(df)
			models.append(header)
			# exit()
	print(np.shape(models))


def find_rings():
	adj_mat_dir = "./output/input_files_three_species_3/adj_matricies/"
	adj_mat_files = glob.glob(adj_mat_dir + '*.csv')

	models = []

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

		# N1 produces B1
		if not mat[B_1_idx][N_1_idx] == 1:
			is_thing = False

		# N2 produces B2
		if not mat[B_2_idx][N_2_idx] == 1:
			is_thing = False

		# N3 produces B3
		if not mat[B_3_idx][N_3_idx] == 1:
			is_thing = False

		# N1 sensitive to B3
		if not mat[N_1_idx][B_3_idx] == -1:
			is_thing = False

		# N2 sensitive to B1
		if not mat[N_2_idx][B_1_idx] == -1:
			is_thing = False

		# N3 sensitive to B1
		if not mat[N_3_idx][B_2_idx] == -1:
			is_thing = False


		# B1 produced by one strain
		if not sum([ mat[B_1_idx][N_1_idx] , mat[B_1_idx][N_2_idx] , mat[B_1_idx][N_3_idx] ]) == 1:
			is_thing = False

		# B2 produced by one strain
		if not sum([ mat[B_2_idx][N_1_idx] , mat[B_2_idx][N_2_idx] , mat[B_2_idx][N_3_idx] ]) == 1:
			is_thing = False

		# B3 produced by one strain
		if not sum([ mat[B_3_idx][N_1_idx] , mat[B_3_idx][N_2_idx] , mat[B_3_idx][N_3_idx] ]) == 1:
			is_thing = False

		## No AHL
		if not sum([ mat[A_1_idx][N_1_idx] , mat[A_1_idx][N_2_idx] , mat[A_1_idx][N_3_idx] ]) == 0:
			is_thing = False

		if not sum([ mat[A_2_idx][N_1_idx] , mat[A_2_idx][N_2_idx] , mat[A_2_idx][N_3_idx] ]) == 0:
			is_thing = False


		## Pure constitutive ring
		if not sum(mat[:, B_1_idx]) == -1:
			is_thing = False

		## Pure constitutive ring
		if not sum(mat[:, B_2_idx]) == -1:
			is_thing = False

		## Pure constitutive ring
		if not sum(mat[:, B_3_idx]) == -1:
			is_thing = False


		if is_thing:
			# print(header)
			# print(df)
			models.append(f.split('_')[-3])
			print(f)
			# models.append(header)
			# exit()
	print(np.shape(models))
	print(models)


def find_AHL_ring():
	adj_mat_dir = "./output/input_files_three_species_3/adj_matricies/"
	adj_mat_files = glob.glob(adj_mat_dir + '*.csv')

	models = []

	for f in adj_mat_files:
		df = pd.read_csv(f)
		header = list(df)[1:]
		mat = df.values
		mat = mat[:, 1:]

		# Is producing AHL
		A_1_idx = header.index("A_1")
		A_2_idx = header.index("A_2")
		A_3_idx = header.index("A_3")

		B_1_idx = header.index("B_1")
		B_2_idx = header.index("B_2")
		B_3_idx = header.index("B_3")

		N_1_idx = header.index("N_1")
		N_2_idx = header.index("N_2")
		N_3_idx = header.index("N_3")

		is_thing = True

		# N1 produces B1
		if not mat[B_1_idx][N_1_idx] == 1:
			is_thing = False

		# N2 produces B2
		if not mat[B_2_idx][N_2_idx] == 1:
			is_thing = False

		# N3 produces B3
		if not mat[B_3_idx][N_3_idx] == 1:
			is_thing = False

		# N1 sensitive to B3
		if not mat[N_1_idx][B_3_idx] == -1:
			is_thing = False

		# N2 sensitive to B1
		if not mat[N_2_idx][B_1_idx] == -1:
			is_thing = False

		# N3 sensitive to B1
		if not mat[N_3_idx][B_2_idx] == -1:
			is_thing = False

		# B1 produced by one strain
		if not sum([ mat[B_1_idx][N_1_idx] , mat[B_1_idx][N_2_idx] , mat[B_1_idx][N_3_idx] ]) == 1:
			is_thing = False

		# B2 produced by one strain
		if not sum([ mat[B_2_idx][N_1_idx] , mat[B_2_idx][N_2_idx] , mat[B_2_idx][N_3_idx] ]) == 1:
			is_thing = False

		# B3 produced by one strain
		if not sum([ mat[B_3_idx][N_1_idx] , mat[B_3_idx][N_2_idx] , mat[B_3_idx][N_3_idx] ]) == 1:
			is_thing = False

		# AHL production
		if not mat[A_1_idx][N_1_idx] == 1:
			is_thing = False
		if not mat[A_2_idx][N_2_idx] == 1:
			is_thing = False
		if not mat[A_3_idx][N_3_idx] == 1:
			is_thing = False

		# AHL repression
		if not mat[B_1_idx][A_1_idx] == -1:
			is_thing = False
		if not mat[B_2_idx][A_2_idx] == -1:
			is_thing = False
		if not mat[B_3_idx][A_3_idx] == -1:
			is_thing = False


		if is_thing:
			# print(header)
			# print(df)
			models.append(f.split('_')[-3])
			print(f)
			# models.append(header)
			# exit()


	print(np.shape(models))
	print(models)


if __name__ == "__main__":
	find_AHL_ring()