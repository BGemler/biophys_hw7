import csv
import re


def gen_all_possible_paranthesis_seqlenonly(seq_len):
	"""
	disregard the sequence itself, just generate all 
	possible parantheses representations "(", ")", or "." at every position
	"""
	parantheses_structures = set()
	possible_states = [")", "(", "."]

	# initialize the first state of parantheses_structures:
	for state in possible_states:
		parantheses_structures.add(state)

	print("seq len:", seq_len, ", N^3:", 3 ** seq_len)

	done_making_structs = False
	while done_making_structs == False:
		done_making_structs = True
		# need to initialize set for each loop step to avoid set size changing
		loop_paranthesis_structures = set(parantheses_structures)

		for seq in loop_paranthesis_structures:
			# check if seq is of required len
			if len(seq) < seq_len:
				done_making_structs = False

				# add possible variations
				for state in possible_states:
					parantheses_structures.add(seq + state)

				# remove seq from dict
				parantheses_structures.remove(seq)

	print("# of structures:", len(parantheses_structures))

	return parantheses_structures


def find_bps_from_seq(sequence, structure):
	"""
	for a structure to be passed here, equal number of ( and )
	need to check if "." between ( and )
	"""
	print(structure)
	bp_poses = []

	for i in range(len(sequence)):
		state = structure[i]

		if state == "(":
			# initialize a new bp
			bp_poses.append([i, ""])
		elif state == ")":
			# find most recent unpaired bp
			placed_bp = False
			for j in range(len(bp_poses)):
				bp = bp_poses[len(bp_poses) - j - 1]
				if bp[1] == "" and placed_bp == False:
					placed_bp = True
					bp[1] = i
					bp_poses[len(bp_poses) - j - 1] = bp


			# if can't place bp, exit with infinite energy bp (seq rejected)
			if placed_bp == False:
				return [["N", "N"]]

	bps = []
	for bp_pos in bp_poses:
		bps.append([sequence[bp_pos[0]], sequence[bp_pos[1]]])

	return bps

	"""

	for_indexes = 

	# find indexes of ( and ) in sequence
	match_for = re.finditer("(", structure)
	for_indexes = [match.start() for match in match_for]

	match_rev = re.finditer(")", structure)
	rev_indexes = [match.start() for match in match_rev]

	# reverse the rev_indexes
	rev_indexes.reverse()

	# the first index in the ( links to the last index in the ), etc...
	for i in range(len(for_indexes)):
		for_index = for_indexes[i]
		rev_index = rev_indexes[i]

		bps.append([sequence[for_index], sequence[rev_index]])
	"""


def evaluate_structures(all_parantheses_structures, \
													sequence, bp_energies, \
													min_allowed_hairpin_size):
	"""
	"""
	passing_parantheses_structures = {}
	for structure in all_parantheses_structures:
		# find number of "(" and number of ")"
		num_for = structure.count("(")
		num_rev = structure.count(")")

		# find minimum hairpin size in the structure
		min_hairpin_size = ""
		seq = structure.split(")")
		for seq in seq:
			se = seq.split("(")
			for s in se:
				if len(s) > 0:
					if min_hairpin_size == "":
						min_hairpin_size = len(s)
					if len(s) < min_hairpin_size:
						min_hairpin_size = len(s)

		# ( can't be right next to )
		no_unallowed_bp_loops = True
		for i in range(len(structure) - 1):
			first_item = structure[i]
			second_item = structure[i + 1]

			if first_item == ")" and second_item == "(" or \
					first_item == "(" and second_item == ")":
					no_unallowed_bp_loops = False

		# find base paired energies
		inf_energy_bp = False
		tot_energy = 0
		bps_out = []
		
		if num_for == num_rev and \
				min_hairpin_size != "" and \
				min_hairpin_size >= min_allowed_hairpin_size and \
				no_unallowed_bp_loops == True:

			bps = find_bps_from_seq(sequence, structure)
			for bp in bps:
				bp_string = bp[0] + bp[1]
				bps_out.append(bp_string)

				if bp_string in bp_energies:
					tot_energy += bp_energies[bp_string]
				else:
					inf_energy_bp = True


		# check conditions
		if num_for == num_rev and \
					min_hairpin_size != "" and \
					min_hairpin_size >= min_allowed_hairpin_size and \
					inf_energy_bp == False and \
					no_unallowed_bp_loops == True:

			print(structure, tot_energy, inf_energy_bp)
			passing_parantheses_structures[structure] = [tot_energy, bps_out]

	return passing_parantheses_structures


def write_out_good_structures(passing_parantheses_structures, sequence):
	"""
	only writing out structures that pass all conditions
	"""
	with open("results/passing-structures.csv", "w") as f:
		out = csv.writer(f)
		out.writerow(["Paranthesis Representation for Sequence:\n" + sequence, \
										"Base Pairs", "Total Energy"])

		for structure in passing_parantheses_structures:
			tot_energy, bps_out = passing_parantheses_structures[structure]

			out.writerow([structure, bps_out, tot_energy])
	f.close()

	return


def main(sequence, bp_energies, min_allowed_hairpin_size):
	"""
	find all possible parantheses structures
	Evaluate parantheses structures by:
	- number of "(" and ")" must be equal
	- a hairpin must exist between "(" and ")"
	- runs of "."s must be >= min_allowed_hairpin_size
	- no base pairs allowed outside of those in bp_energies
	"""
	# need to walk along the sequence, try every possible combination 
	# note: no restraints applied here
	all_parantheses_structures = gen_all_possible_paranthesis_seqlenonly(len(sequence))

	# apply restraints, only get passing seqs
	passing_parantheses_structures = evaluate_structures(all_parantheses_structures, \
																													sequence, bp_energies, \
																													min_allowed_hairpin_size)

	# write passing structures
	write_out_good_structures(passing_parantheses_structures, sequence)

	return


sequence = "UGCGCAAAGUAUCA"
bp_energies = {
	"GU": -1,
	"UG": -1,
	"AU": -2,
	"UA": -2,
	"GC": -3,
	"CG": -3
}
min_allowed_hairpin_size = 3


main(sequence, bp_energies, min_allowed_hairpin_size)