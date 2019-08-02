# -*- coding: utf-8 -*-

"""
`` **module description**:
This module allows to simulate biological sequence evolution throught sequence and splicing structure evolution events.
.. moduleauthor:: Esaie Kuitche - CoBIUS Team - Univeristy of Sherbrooke.
Juillet 2019

"""

import copy
import numpy as np
import random
from random import shuffle
import itertools
from ete3 import Tree
import argparse
import os

def build_arg_parser():
	parser = argparse.ArgumentParser(description="SimSpliceEvol program parameters")
	parser.add_argument('-i', '--input_tree_file', help = "input guide tree (default = Example/input/large.nw)", default = "Example/input/large.nw")    
	parser.add_argument('-n', '--number_of_simulation', default= 10, help = "number of simulations") 	

	parser.add_argument('-k_indel', '--k_indel', help="multiplicative constant for codon indel rate (default =0.5)", default = 0.5) 
	parser.add_argument('-k_ei_c', '--k_ei_c',  help="multiplicative constant for exon-intron change (eic) rate  (default =5)",  default = 20) 
	parser.add_argument('-k_t_c', '--k_t_c', default =5, help="multiplicative constant for transcript change (tc) rate (default =5)") 
	parser.add_argument('-k_intron', '--k_intron', default = 1.5,  help = "multiplicative constant for substitution rate in intron (default =1.5)") 
	parser.add_argument('-k_nb_exons', '--k_nb_exons', default = 1.5, help="multiplicative constant for number of exons in gene (default =1.5)")  


	parser.add_argument('-eic_l', '--exon_i_change_loss',  default =0.4, help="relative frequence of exon loss in eic (default =0.4)") 
	parser.add_argument('-eic_g', '--exon_i_change_gain',  default = 0.5, help="relative frequence of exon gain in eic (default =0.5)") 
	parser.add_argument('-eic_d', '--exon_i_change_dup',  default =0.1, help="relative frequence of exon duplication in eic (default =0.1)")  

	parser.add_argument('-tc_rs', '--random_selection',  default = 0.5, help="relative frequence of random selection in tc (default =0.5)")
	parser.add_argument('-tc_a5', '--alternative_five_prime', default =  0.1, help="relative frequence of alternative five prime in tc (default =0.1)")
	parser.add_argument('-tc_a3', '--alternative_three_prime', default = 0.1, help="relative frequence of alternative three prime in tc (default =0.1)")
	parser.add_argument('-tc_ek', '--exon_skipping', default = 0.2, help="relative frequence of exon skipping in tc (default =0.1)")
	parser.add_argument('-tc_me', '--mutually_exclusive', default = 0.1, help="relative frequence of mutually exclusive in tc (default =0.1)")
	parser.add_argument('-tc_ir', '--intron_retention', default = 0.05, help="relative frequence of intron retention in tc (default =0.05)")
	parser.add_argument('-tc_tl', '--transcript_loss', default = 0.4, help="relative frequence of transcript loss in tc (default =0.4)")

	
	
	

	return parser



def initCodonMatrix(codonMatrixFile):
	with open(codonMatrixFile, 'r') as f:
		lines = f.readlines()
		size = len(lines)
		matrix = np.zeros((size, size))
		for i, line in enumerate(lines):
			for j, elem in enumerate(line.split()):
				matrix[int(i), int(j)] = matrix[int(j), int(i)] = elem
		for i in range(size):
			for j in range(size):
				matrix[i, j] = 10 * round(matrix[i, j])
		for i in range(size):
			s = 0
			for j in range(size):
				s += matrix[i, j]
			matrix[i, i] = 1000 - s
		return matrix


def init(one_nt_dist, two_nt_dist, three_nt_dist):
	one_nt_dist = {'A': 24, 'C': 32, 'T': 19, 'G': 25}

	two_nt_dist = {
		'A': {'A': 21, 'C': 28, 'T': 21, 'G': 30},
		'C': {'A': 27, 'C': 38, 'T': 21, 'G': 14},
		'T': {'A': 12, 'C': 35, 'T': 20, 'G': 33},
		'G': {'A': 29, 'C': 30, 'T': 19, 'G': 22}
	}

	three_nt_dist = {
		'AA': {'A': 16, 'C': 35, 'T': 13, 'G': 36},
		'AC': {'A': 22, 'C': 49, 'T': 14, 'G': 15},
		'GT': {'A': 8, 'C': 30, 'T': 16, 'G': 46},
		'AG': {'A': 12, 'C': 55, 'T': 20, 'G': 13},
		'CC': {'A': 26, 'C': 35, 'T': 28, 'G': 10},
		'CA': {'A': 9, 'C': 22, 'T': 12, 'G': 57},
		'CG': {'A': 17, 'C': 37, 'T': 20, 'G': 26},
		'TT': {'A': 10, 'C': 49, 'T': 22, 'G': 20},
		'GG': {'A': 19, 'C': 39, 'T': 13, 'G': 29},
		'GC': {'A': 15, 'C': 55, 'T': 18, 'G': 12},
		'AT': {'A': 2, 'C': 40, 'T': 14, 'G': 43},
		'GA': {'A': 20, 'C': 25, 'T': 18, 'G': 36},
		'TG': {'A': 5, 'C': 41, 'T': 30, 'G': 23},
		'TA': {'A': 1, 'C': 65, 'T': 34, 'G': 1},
		'TC': {'A': 20, 'C': 48, 'T': 20, 'G': 12},
		'CT': {'A': 7, 'C': 26, 'T': 12, 'G': 56}
	}
	nt_inton_dist = {'A': 29, 'C': 20, 'T': 30, 'G': 21}
	mean_exon_gene = 0
	sd_exon_gene = 0
	mean_exon_cds = 9
	sd_exon_cds = 8
	mean_exon_len = 170
	sd_exon_len = 258
	mean_len_introns = 3760
	sd_len_intons = 20126
	return mean_len_introns, sd_len_intons, mean_exon_cds, sd_exon_cds, mean_exon_gene, sd_exon_gene, mean_exon_len, sd_exon_len, one_nt_dist, two_nt_dist, three_nt_dist, nt_inton_dist


def generate_nt_for_random(one_nt_dist, two_nt_dist, three_nt_dist, nt_inton_dist):
	one_nt_dist_nt = []
	one_nt_dist_nt = list("".join([key * value for key, value in one_nt_dist.items()]))
	shuffle(one_nt_dist_nt)

	nt_inton_dist_nt = []
	nt_inton_dist_nt = list("".join([key * value for key, value in nt_inton_dist.items()]))
	shuffle(nt_inton_dist_nt)

	two_nt_dist_nt = {}
	for k, v in two_nt_dist.items():
		tmp_list = []
		tmp_list = list("".join([key * value for key, value in v.items()]))
		shuffle(tmp_list)
		two_nt_dist_nt[k] = tmp_list

	three_nt_dist_nt = {}
	for k, v in three_nt_dist.items():
		tmp_list = []
		tmp_list = list("".join([key * value for key, value in v.items()]))
		shuffle(tmp_list)
		three_nt_dist_nt[k] = tmp_list
	return one_nt_dist_nt, two_nt_dist_nt, three_nt_dist_nt, nt_inton_dist_nt


def generate_number_by_normal_law(mu, sigma):
	val = int(np.random.normal(mu, sigma, 1)[0])
	if val > 0:
		return val
	else:
		return generate_number_by_normal_law(mu, sigma)

def generate_number_max_exon_transcript(mu, sigma):
	val = int(np.random.normal(mu, sigma, 1)[0])
	if val > 0:
		return val
	else:
		return generate_number_max_exon_transcript(mu, sigma)

def generate_number_max_exon_gene(number_max_exon_transcript, k):   
	if k >= 1:
		return int(round(k*number_max_exon_transcript))
	else:
		return number_max_exon_transcript



def generate_exon(len_exon, one_nt_dist_nt, two_nt_dist_nt, three_nt_dist_nt, end_exon):
	i = 0
	one = ""
	two = ""
	three = ""
	#exon = "ATG"
	exon = ""
	stopCodon = ['TAG', 'TAA', 'TGA']

	if len(end_exon) == 1:
		while True:
			two = random.choice(two_nt_dist_nt[end_exon])
			exon = two + random.choice(three_nt_dist_nt[end_exon + two])
			if (two + exon in stopCodon):
				pass
			else:
				break

	elif len(end_exon) == 2:
		while True:
			exon = random.choice(three_nt_dist_nt[end_exon])
			if (end_exon + exon in stopCodon):
				pass
			else:
				break

	while i < len_exon:
		if i % 3 == 0:
			one = random.choice(one_nt_dist_nt)
		elif i % 3 == 1:
			two = one + random.choice(two_nt_dist_nt[one])

		else:
			while True:
				three = two + random.choice(three_nt_dist_nt[two])
				if three in stopCodon:
					pass
				else:
					exon = exon + three
					break
		i = i + 1

	if i % 3 == 0:
		return exon, ""
	elif i % 3 == 1:
		return exon + one, one
	else:
		return exon + two, two


def generate_intron(intron_lenght, nt_inton_dist_nt):
	i = 0
	splicesites = [["GT-AG"]*90,  ["GC-AG"]*7, [""]*3]
	splicesites = list(itertools.chain.from_iterable(splicesites))
	splicesite = random.choice(splicesites)

	if len(splicesite)>0:
		splicesitetmp = splicesite.split("-")
		donnor = splicesitetmp[0]
		acceptor = splicesitetmp[1]
	else:
		donnor = random.choice(["A", "G", "C", "T"]) + "" + random.choice(["A", "G", "C", "T"])
		acceptor = random.choice(["A", "G", "C", "T"]) + "" + random.choice(["A", "G", "C", "T"])

	intron = donnor

	while i < (intron_lenght - 4):
		intron += random.choice(nt_inton_dist_nt)
		i += 1

	intron += acceptor   
	return intron


def write_init(seq_dict, name):
	keys = seq_dict.keys()
	listKeyOrder = []
	n = len(keys)
	#file = open(name + "_set_at_root.fasta", "w")

	for i in range(n):
		if name + "_" + str(i) in keys:
			listKeyOrder.append(name + "_" + str(i))
			#file.write(">" + name + "_" + str(i) + "\n")
			#file.write(seq_dict[name + "_" + str(i)] + "\n")
	return listKeyOrder


def adjust_exons_extremities(exons_dict, end_exon):
	stopCodon = ['TAG', 'TAA', 'TGA']
	seq = exons_dict["exon_0"]
	#le premier exon commence par un ATG
	if seq[:3] != "ATG":
		seq = "ATG" + seq
		exons_dict["exon_0"] = seq

	#le dernier exon finir par un codon stop
	id_last = len(exons_dict.keys()) - 1
	seq = exons_dict["exon_" + str(id_last)]
	if len(end_exon) == 0:
		seq_tmp = seq[:-3] + random.choice(stopCodon)
	else:
		seq_tmp = seq[:-len(end_exon)] + random.choice(stopCodon)
	exons_dict["exon_" + str(id_last)] = seq_tmp
	
	if len(exons_dict.keys()) > 2:
		i = 1 
		while i < (len(exons_dict.keys())-1):
			if random.choice([True, False]):
				exons_dict["exon_" + str(i)] = "ATG" + exons_dict["exon_" + str(i)]
			i += 1

def make_structure_evol(tree, exon_loss_dup_gain, prob_loss, prob_dup, prob_gain, exons_dict, exon_to_gain, listKeyOrder, k2):
	gain = []
	dup = []
	loss = []
	freq_gain = 0
	freq_dup = 0
	freq_loss = 0
	i = 0
	for node in tree.traverse("preorder"):
		if node.is_root():
			node.name = [copy.deepcopy(exon_loss_dup_gain), copy.deepcopy(listKeyOrder), "", ""]

		else:
			cs = node.dist

			leaf_name = node.name

			node.name = copy.deepcopy(node.get_ancestors()[0].name)

			node.name[3] = leaf_name

			nb_exon = 0
			exon_may_be_lost = []
			for e in node.name[0].keys():
				if node.name[0][e][0] == 0:
					nb_exon +=1
					exon_may_be_lost.append(e)

				if node.name[2] == "":
					node.name[2] = {}
					node.name[2][e] = exons_dict[e]
				else:
					node.name[2][e] = exons_dict[e]


			nb_loss = int(nb_exon * prob_loss * k2 * cs) 
			exon_to_delete = random.choices(exon_may_be_lost, k=nb_loss)
			for e in exon_to_delete:
				node.name[0][e][0] = 1
			

			nb_dup = int((nb_exon-nb_loss) * prob_loss * k2 * cs) 
			exon_may_be_dup = [e for e in exon_may_be_lost if e not in exon_to_delete]
			exon_to_duplicate = random.choices(exon_may_be_dup, k=nb_dup)
			for e in exon_to_duplicate:
				dup_state = [0, e]
				exons_dict[e + "_dup"] = exons_dict[e]
				exon_loss_dup_gain[e + "_dup"] = [0, [0, e], 0]
				node.name[1].insert(node.name[1].index(e) + 1, e + "_dup")
				node.name[0][e + "_dup"] = [0, [1, e], 0]
				node.name[0][e][1] = [dup_state]
				if node.name[2] == "":
					node.name[2] = {}
					node.name[2][e + "_dup"] = exons_dict[e]
				else:
					node.name[2][e + "_dup"] = exons_dict[e]



			nb_gain = int((nb_exon-nb_loss+nb_dup) * prob_loss * k2 * cs) 
			for cmpt in range(nb_gain):
				exons_dict["exon_new_" + str(i)] = exon_to_gain.pop()
				if node.name[2] == "":
					node.name[2] = {}
					node.name[2]["exon_new_" + str(i)] = exons_dict["exon_new_" + str(i)]
				else:
					node.name[2]["exon_new_" + str(i)] = exons_dict["exon_new_" + str(i)]
				
				exon_loss_dup_gain["exon_new_" + str(i)] = [0, [0, ""], 1]
				node.name[1].insert(random.randint(0, len(node.name[1])), "exon_new_" + str(i))
				i = i + 1			


def makeDistribution(codon_matrix, codonslist):
	codon_dict = {}
	for i in range(len(codonslist)):
		l = []
		for j in range(len(codonslist)):
			l += [codonslist[j]] * int(codon_matrix[i, j])
		codon_dict[codonslist[i]] = l
	codon_dict['TAG'] = ['TAG']
	codon_dict['TAA'] = ['TAA']
	codon_dict['TGA'] = ['TGA']
	return codon_dict


def make_substition_evol(tree, exons_dict, codon_distribution, listKeyOrder, k1):
	gene_exon_state = {}
	for node in tree.traverse("preorder"):
		if node.is_root():
			node.name[2] = copy.deepcopy(exons_dict)
		else:					
			parent = node.get_ancestors()[0]
			parent_exon_id = parent.name[1]

			parentSeq = "".join([exons_dict[x] for x in parent_exon_id])
			
			#l = int(len(parentSeq)/3.0)
			l = int(len(parentSeq)/3.0)

			nodeSeq = ""

			cs = node.dist

			nb_subt = int(l * cs * k1)
			all_data = range(l)
			subt_positions = choice_distint(all_data, nb_subt)
			subt_positions = [3*x for x in subt_positions]
			#print(all_data)
			#print ("===========", nb_subt, len(parentSeq), sorted(subt_positions))
			start = 0
			end = 0
			for exon_id in parent_exon_id:
				nb = 0
				start = end 

				end = start + len(parent.name[2][exon_id])
				exon_tmp = parent.name[2][exon_id]
				#print("exons_len", len(exon_tmp))
				for position in sorted(subt_positions):

					if position >= start and position+3 <= end:						
				
						position = position - start

						if False : #exon_tmp[position:position+3] == "ATG":
							codon_subt = "ATG"
						else:	
							codon_subt = random.choice(codon_distribution[exon_tmp[position:position+3]])
							nb +=1
						exon_tmp = exon_tmp[0:position] + codon_subt + exon_tmp[position+3:]
				
				if node.name[3] in gene_exon_state.keys():
					gene_exon_state[node.name[3]][exon_id] = exon_tmp
				else:
					gene_exon_state[node.name[3]] =  {}
					gene_exon_state[node.name[3]][exon_id] = exon_tmp

				
				node.name[2][exon_id] = exon_tmp
				#print("-------------------------",nb)

	
def make_indel_evol(tree, exons_dict, codon_distribution, listKeyOrder, mean_nb_indel, mean_len_indel, sd_len_indel, exon_to_gain,
					dict_indel_exons, k1, prob_insertion, prob_deletion):
	indel_exon_state = {}

	for node in tree.traverse("preorder"):		
		if node.is_root():
			pass
		elif node.is_leaf:
			
			cs = node.dist
			
			seqLen = 0
			nb_exon = 0
			for e, seq in node.name[0].items():
				if node.name[0][e][0] == 0:
					seqLen += len(seq)
					nb_exon += 1

			seqLen = int(seqLen/3)

			nb_deletion = int(seqLen*prob_deletion*k1*cs)

			nb_insertion = int(seqLen*prob_insertion*k1*cs)

			deletion_positions = choice_distint(range(seqLen), nb_deletion)
			insertion_positions = choice_distint(range(seqLen), nb_insertion)


			start = 0
			end = 0
			for exon_id in node.name[0].keys():			
				start = end 
				end = start + len(node.name[2][exon_id])
				exon_tmp = node.name[2][exon_id]

				for position in sorted(insertion_positions):
					
					if position >= start and position < end:						
						len_indel = generate_number_by_normal_law(mean_len_indel, sd_len_indel)
						begin = int(random.choice(range(len(exon_to_gain)))/2)

						exon_gained = exon_to_gain[begin: begin+len_indel]

						if node.name[3] in indel_exon_state.keys():
							if exon_id in indel_exon_state[node.name[3]].keys():
								indel_exon_state[node.name[3]][exon_id].append(["insertion", position,exon_gained])
							else:
								indel_exon_state[node.name[3]][exon_id] = [["insertion", position,exon_gained]]
						else:
							indel_exon_state[node.name[3]] =  {}
							indel_exon_state[node.name[3]][exon_id] = [["insertion", position,exon_gained]]

			start = 0
			end = 0
			for exon_id in node.name[0].keys():
				start = end 
				end = start + len(node.name[2][exon_id])
				exon_tmp = node.name[2][exon_id]

				for position in sorted(deletion_positions):
					#print (position, start, end, len(exon_tmp))
					len_indel = generate_number_by_normal_law(mean_len_indel, sd_len_indel)
					if position >= start and position +len_indel < end:												

						if node.name[3] in indel_exon_state.keys():
							if exon_id in indel_exon_state[node.name[3]].keys():
								indel_exon_state[node.name[3]][exon_id].append(["deletion", position, len_indel])
							else:
								indel_exon_state[node.name[3]][exon_id] = [["deletion", position, len_indel]]

						else:
							indel_exon_state[node.name[3]] =  {}
							indel_exon_state[node.name[3]][exon_id] = [["deletion", position, len_indel]]

	
	return indel_exon_state

		




def make_list(liste, n):
	"""
	Cette fonction calcule l'ensemble des plages des séquences utilisé dans un exon
	:param liste: liste[0] == position d'insertion et liste[1] == liste des positions des deletes dans l'exon
	:param n: longueur total de la séquences de l'exon
	:return: liste des positions triés en ordre croissant
	"""
	exclude_values = []
	position = -1
	# print liste
	for e in liste[0]:
		exclude_values = exclude_values + range(e[0], e[0] + e[1])
	for e in liste[1]:
		exclude_values.append(e[0])

	while (True):
		position = random.choice(range(n + 1))
		if not (position in exclude_values):
			break

	return exclude_values, position


def generate_bloc_insertion(exon, exon_use_position, seq_node, exon_to_gain, len_indel):
	position = 0
	seq = ""
	if len(exon_use_position[0]) == 0 and len(exon_use_position[1]) == 0:
		position = random.randint(0, len(seq_node))
		seq = exon_to_gain[0:len_indel]
		exon_to_gain = exon_to_gain[len_indel:]
	else:
		liste = []
		n = len(seq_node[exon])
		if n == 0:
			exit("len exon == 0")
		liste, position = make_list(copy.deepcopy(exon_use_position), n)
		seq = exon_to_gain[0:len_indel]
		exon_to_gain = exon_to_gain[len_indel:]
	return position, seq, exon_to_gain


def generate_block_deletion(exon, exon_use_position, len_indel, seq_node):
	position = 0
	extention = 0
	range_block = 0
	liste = []
	if len(exon_use_position[0]) == 0 and len(exon_use_position[1]) == 0:		
		position = random.randint(0, len(seq_node[exon]))
	else:
		liste = []
		liste, position = make_list(copy.deepcopy(exon_use_position), len(seq_node[exon]))
	return position, extention


def multiple_alignment(tree, dict_indel_exons, gene_datas):
	# parcours des feuilles de l'arbre tree
	for n in tree.traverse("preorder"):
		if n.is_leaf():
			exons_include_in_node = []
			structure_event = n.name[0]
			for exon_id, events in structure_event.items():
				if events[0] == 0:
					exons_include_in_node.append(exon_id)
			# Prendre tous les noeuds des protéines

			for exon_id in exons_include_in_node:
				to_delete = []
				to_insert = []
				if n.name[0][exon_id][0] == 1:
					continue
				else:
					to_delete = n.name[4][exon_id][0]
					to_insert = n.name[4][exon_id][1]
					exon_seq = n.name[2][exon_id]
					lon = len(exon_seq)
					insert_block_list = {}
					for insert in to_insert:
						insert_block_list[insert[0]] = insert[1]  # cle = position et val = seq

					if len(to_delete) > 0:
						for block in to_delete:
							if block[0] in insert_block_list:
								exit()
								exon_seq = exon_seq[:block[0]] + insert_block_list[block[0]] + exon_seq[
																							   block[0] + block[1]:]
							else:
								exon_seq = exon_seq[:block[0]] + "-" * block[1] + exon_seq[block[0] + block[1]:]

					if lon != len(exon_seq):
						print ("EROOOOR")
						exit()
					position_seq = {}
					list_position_insertion = []
					for insert_block in dict_indel_exons[exon_id][1]:
						if insert_block[0] in position_seq.keys():
							print ("dup postion to insert")
							exit()
						position_seq[insert_block[0]] = insert_block[1]
						list_position_insertion.append(insert_block[0])

					if len(list_position_insertion) != len(list(set(list_position_insertion))):
						print ("Error")
						exit(-1)
					list_position_insertion.sort(reverse=True)
					len_seq = len(exon_seq)

					for pos in list_position_insertion:
						seq = position_seq[pos]
						lon2 = len(exon_seq)
						if pos < len_seq - 1:
							if pos in insert_block_list:
								exon_seq = exon_seq[0:pos] + seq + exon_seq[pos:]
							else:
								exon_seq = exon_seq[0:pos] + "-" * len(seq) + exon_seq[pos:]
						else:
							return -1
						if lon2 + len(seq) != len(exon_seq):
							print ("error insert len")
							exit()

					n.name[2][exon_id] = exon_seq
					# print exon_seq

	return 1


def write_gene_in_file(tree, indel_exon_state, intron_dict, src, src2, name=""):
	genes_seq = {}
	gene = open(src + "genes/" + src2 + name + "gene.fasta", "w")

	exon_intron_of_genes = {}	
	resulting_gene_exon = {}

	for node in tree.traverse("preorder"):
		if node.is_leaf():

			resulting_gene_exon[node.name[3]] = {}

			l_datas = {}
			gene_ids = node.name[3]
			l = []
			
			for e, k in node.name[0].items():
				l.append(e)					

			exon_ids = sorted(l)

			seq = ""			

			exon_intron = []

			intron_id_list = list(intron_dict.keys())

			for exon_id in exon_ids:
				seq_exon = node.name[2][exon_id]
				if exon_id in indel_exon_state.keys():
					exon_events = indel_exon_state[node.name[3]][exon_id]
					
					for e in exon_events:
					
						if e[0] == "insertion":
							seq_exon = seq_exon[0:e[1]] + e[2] + seq_exon[e[1]:]
					
					for e in exon_events:
						print(e)
						if e[0] == "deletion":
							seq_exon = seq_exon[0:e[1]] + seq_exon[e[1]+e[2]:]
				
				resulting_gene_exon[node.name[3]][exon_id] = seq_exon

				exon_intron.append(exon_id)
				intron = random.choice(intron_id_list)
				intron_id_list.remove(intron)

				exon_intron.append(intron)				

				seq += seq_exon + intron_dict[intron]

			gene.write(">" + node.name[3] + "\n")
			gene.write(seq + "\n")
			exon_intron_of_genes[node.name[3]] = exon_intron
			genes_seq[node.name[3]] = seq
	
	return exon_intron_of_genes, resulting_gene_exon, genes_seq




def choice_distint(all_data, k):
	choices = []
	if k > len(all_data):
		k = len(all_data)
	while len(choices) < k:
		selection = random.choice(all_data)
		if selection not in choices:
			choices.append(selection)
	return choices

def make_protein_evol(tree, exon_intron_of_genes, resulting_gene_exon, tc5, tc3, tcsk, tcme, tc6, k3, intron_dict, path, iteration, genes_seq):
	
	leaves = {}

	cds_list = {}
	exon_events = []
	for node in tree.traverse("preorder"):
		cds_id =1
		l = 1.5
		if node.is_root():
			pass
		elif node.is_leaf():
			used_cds =  []
			gene_id = node.name[3]

			cds_list[node.name[3]] = {}
			cds_list_gene = {}
			exon_intron = exon_intron_of_genes[node.name[3]]
			exon_sequences = resulting_gene_exon[node.name[3]]	

			exon_ids = []

			for e in exon_intron:
				if e.startswith("exon"):
					exon_ids.append(e)
			cds_0 = choice_distint(exon_ids, int(len(exon_ids)-len(exon_ids)/5))
			
			cds_0 = [e for e in cds_0 if node.name[0][e][0] == 0]
			used_cds.append(cds_0)
			#for e in exon
			cds_list[node.name[3]]["cds_0"] = [cds_0, {}, {}]
			cds_list_gene["cds_0"] = [cds_0, {}, {}]

			cs = node.dist

			
			nb_random = int(round(l * tc6 * k3 * cs *15))
			#print(nb_random)
			nb_random_retain = 0
			for i in range(nb_random):
				exons = choice_distint(exon_ids, int(len(exon_ids)-len(exon_ids)/5))
				exons_not_delete = [e for e in cds_0 if node.name[0][e][0] == 0]
				if len(exons) != len(exons_not_delete) or exons_not_delete in used_cds:
					pass
				else:					
					cds_list[node.name[3]]["cds_" + str(cds_id)] = [exons, {}, {}]
					cds_list_gene["cds_" + str(cds_id)] = [exons, {}, {}]
					cds_id +=1
					nb_random_retain += 1
					used_cds.append(exons)					
			
			l+= nb_random_retain

			nb_5 = int(round(l * tc5 * k3 * cs ))
			for i in range(nb_5):				
				exon =random.choices(cds_0, k=1)
				len_5 = random.choice(range(len(exon_sequences[exon[0]])))	
				#len_5 = int(len_5/2)
				cds_list[node.name[3]]["cds_" + str(cds_id)] = [cds_0, {exon[0]:len_5}, {}]				
				cds_list_gene["cds_" + str(cds_id)] = [cds_0, {exon[0]:len_5}, {}]				
				
				cds_id += 1				

			l += nb_5
			"""
			nb_3 = int(round(l * tc3 * k3 * cs ))			
			nb_5 = int(round(l * tc5 * k3 * cs ))
			
			for i in range(nb_5):
				exon =random.choices(cds_0, k=1)
				len_5 = random.choice(range(len(exon_sequences[exon[0]])))	
				cds_list[node.name[3]]["cds_" + str(cds_id)] = [cds_0, {exon[0]:len_5}, {}]				
				cds_list_gene["cds_" + str(cds_id)] = [cds_0, {exon[0]:len_5}, {}]				
				cds_id += 1

			l += nb_5
			"""

			nb_3 = 0#int(round(l * tc3 * k3 * cs ))
			
			for i in range(nb_3):
				cds_seleted = choice_distint(list(cds_list_gene.keys()), nb_3)
				for cds in cds_seleted:
					exon =random.choices(cds_list_gene[cds][0], k=1)
					
					if cds in cds_list_gene.keys() and len(cds_list_gene[cds][1].keys())>0 and exon in list(cds_list_gene[cds][1].keys()) :						

						len_3 = random.choice(range(len(exon_sequences[exon[0]])-cds_list_gene[cds][1][exon[0]]))					
						cds_list[node.name[3]]["cds_" + str(cds_id)] = [cds_list_gene[cds][0], cds_list_gene[cds][1], {exon[0]:len_3}]
						cds_list_gene["cds_" + str(cds_id)] = [cds_list_gene[cds][0], cds_list_gene[cds][1], {exon[0]:len_3}]
					else:
						len_3 = random.choice(range(len(exon_sequences[exon[0]])))					
						cds_list[node.name[3]]["cds_" + str(cds_id)] = [cds_list_gene[cds][0], {}, {exon[0]:len_3}]
						cds_list_gene["cds_" + str(cds_id)] = [cds_list_gene[cds][0], {}, {exon[0]:len_3}]
					
					cds_id += 1

			l += nb_3
			nb_skipping = int(round(l * tcsk * k3 * cs ))
			cds_seleted = choice_distint(list(cds_list_gene.keys()), nb_skipping)
			for cds in cds_seleted:
				exon =random.choices(cds_list_gene[cds][0], k=1)
				cds_list[node.name[3]]["cds_" + str(cds_id)] = [[x for x in cds_list_gene[cds][0] if x!= exon[0]], {}, {}]			
				cds_list_gene["cds_" + str(cds_id)] = [[x for x in cds_list_gene[cds][0] if x!= exon[0]], {}, {}]			
				cds_id +=1

			l += nb_skipping
			nb_mutually_exclusive = int(round(l * tcme * k3 * cs ))
			cds_seleted = choice_distint(list(cds_list_gene.keys()), nb_mutually_exclusive)
			for cds in cds_seleted:
				if len(cds_list_gene[cds][0])>2:
					exons =choice_distint(cds_list_gene[cds][0], 2)
					cds_list[node.name[3]]["cds_" + str(cds_id)] = [[x for x in cds_list_gene[cds][0] if x != exons[0]], {}, {}]			
					cds_list_gene["cds_" + str(cds_id)] = [[x for x in cds_list_gene[cds][0] if x != exons[0]], {}, {}]			
					cds_id +=1
					cds_list[node.name[3]]["cds_" + str(cds_id)] = [[x for x in cds_list_gene[cds][0] if x != exons[1]], {}, {}]			
					cds_list_gene["cds_" + str(cds_id)] = [[x for x in cds_list_gene[cds][0] if x != exons[1]], {}, {}]			
					cds_id +=1

			leaves[gene_id] = node.name


	write_datas(leaves, cds_list, exon_intron_of_genes, resulting_gene_exon,intron_dict, path, iteration, genes_seq)
	
	
		
def write_datas(leaves, cds_list, exon_intron_of_genes, resulting_gene_exon,intron_dict, path, iteration, genes_seq):
	cds_file = open(path + "cds/" + iteration + "cds.fasta", "w")
	exon_positions_file = open(path + "positions/" + iteration + "exon_positions.fasta", "w")
	cds_to_gene_file = open(path + "cds_gene/" + iteration + "cds_gene.txt", "w")	
	alignement_pair = open(path + "pairwise_alignment/" + iteration + "alignement.fasta", "w")
	cluster = open(path + "cluster/" + iteration + "cluster.fasta", "w")
	
	alignement_multiple = open(path + "multiple_alignment/" + iteration + "multiple.fasta", "w")

	list_cds = {}
	all_exons = []

	all_exons_seq = {}
	for gene_id, node in leaves.items():		
		exon_intron = exon_intron_of_genes[gene_id]
		for e in exon_intron:
			if e.startswith("exon"):
				all_exons.append(e)
		for id, seq in node[2].items():
			all_exons_seq[id] = seq

	all_exons = sorted(list(set(all_exons)))

	all_cds_datas = {}
	for gene_id, node in leaves.items():		
		exon_intron = exon_intron_of_genes[gene_id]
		exon_sequences = resulting_gene_exon[gene_id]	

		alignement_pair.write(">" + gene_id + "\n")
		alignement_pair.write(genes_seq[gene_id] + "\n")

		cds_of_this_gene = cds_list[gene_id]
		for cds_id, cds_infos in cds_of_this_gene.items():
			cds_datas = {}
			all_cds_datas[cds_id + "_" +gene_id] = {}
			cds_seq = ""
			cds_position = []			
			cds_exons = sorted(cds_infos[0])

			cds_5prime = cds_infos[1]
			cds_3prime = cds_infos[2]
			
			exon_start_cds = 0
			exon_end_cds = 0
			cds_seq_aligned_to_is_gene = ""
			
			last_position = -1
			for exon_id in cds_exons:
				
				exon_start_gene = 0
				exon_end_gene = 0				

				exon_seq = node[2][exon_id]
				len_init_exon_seq = len(exon_seq)
				position_start_3prime = 0
				position_start_5prime = 0
				if exon_id in cds_3prime.keys():
					position_start_3prime = cds_3prime[exon_id]
					exon_seq = exon_seq[0:len(exon_seq)-position_start_3prime]

				if exon_id in cds_5prime.keys():
					position_start_5prime = cds_5prime[exon_id]
					exon_seq = exon_seq[position_start_5prime:]

				exon_start_cds = exon_end_cds
				exon_end_cds += len(exon_seq)

				cds_datas[exon_id] = [exon_seq, exon_start_cds, exon_end_cds, position_start_5prime, position_start_3prime, len_init_exon_seq]				
				all_cds_datas[cds_id + "_" +gene_id][exon_id] = [exon_seq, exon_start_cds, exon_end_cds, position_start_5prime, position_start_3prime, len_init_exon_seq]				
				cds_seq += exon_seq

			exon_start_gene = 0
			exon_end_gene = 0
			cds_seq_aligned_to_is_gene = ""

			for e_i in exon_intron:
				if e_i in cds_datas.keys():
					exon_start_cds = cds_datas[e_i][1]
					exon_end_cds = cds_datas[e_i][2]
					exon_seq = cds_datas[e_i][0]
					position_start_5prime = cds_datas[e_i][3]
					position_start_3prime = cds_datas[e_i][4]
					len_init_exon_seq = cds_datas[e_i][5]

					exon_start_gene = exon_end_gene + position_start_5prime
					exon_end_gene = exon_start_gene + len(exon_seq) - position_start_3prime
					cds_position.append([exon_start_cds, exon_end_cds, exon_start_gene, exon_end_gene])

					cds_seq_aligned_to_is_gene += position_start_5prime*"-"
					cds_seq_aligned_to_is_gene += exon_seq
					cds_seq_aligned_to_is_gene += position_start_3prime*"-"						
				else:
					exon_start_gene = exon_end_gene
					if e_i.startswith("exon"):
						cds_seq_aligned_to_is_gene += len(node[2][e_i])*"-"
						exon_end_gene += len(node[2][e_i])
					else:
						exon_end_gene += len(intron_dict[e_i])					
						cds_seq_aligned_to_is_gene += len(intron_dict[e_i])*"-"


			alignement_pair.write(">" + cds_id + "_" +gene_id + "\n")
			alignement_pair.write(cds_seq_aligned_to_is_gene + "\n")
			cds_file.write(">" + cds_id + "_" +gene_id +"\n")
			cds_file.write(cds_seq +"\n")

			exon_positions_file.write(">" + cds_id + "_" +gene_id +"\n")
			for  position in cds_position:
				exon_positions_file.write(str(position[0])+ " " + str(position[1]) + " " + str(position[2]) + " " +str(position[3]) + "\n")

			cds_to_gene_file.write( cds_id + "_" +gene_id + " " + gene_id + "\n")
			list_cds[cds_id + "_" +gene_id ] = cds_exons


			cds_seq_multi_aligned = ""

			for e_i in all_exons:
				if e_i in cds_datas.keys():
					exon_seq = cds_datas[e_i][0]
					position_start_5prime = cds_datas[e_i][3]
					position_start_3prime = cds_datas[e_i][4]

					cds_seq_multi_aligned += position_start_5prime*"-"
					cds_seq_multi_aligned += exon_seq
					cds_seq_multi_aligned += position_start_3prime*"-"						
					exon_i = position_start_5prime + len(exon_seq) + position_start_3prime
				else:
					cds_seq_multi_aligned += len(all_exons_seq[e_i])*"-"
					exon_i = len(all_exons_seq[e_i])

			alignement_multiple.write(">" + cds_id + "_" +gene_id + "\n")
			alignement_multiple.write(cds_seq_multi_aligned + "\n")					

	i = 1
	clusters = {}
	used = []
	while len(list_cds.keys()) > 0:
		cds_to_compare = list(list_cds.keys())[0]
		all_exon_infos = []
		exon_seq = list(all_cds_datas[cds_to_compare].keys())


		for e, primes in all_cds_datas[cds_to_compare].items():
			splice_5prime = all_cds_datas[cds_to_compare][e][3]		
			splice_3prime = all_cds_datas[cds_to_compare][e][4]
			all_exon_infos = [list_cds[cds_to_compare], splice_5prime, splice_3prime]

		all_exons_infos = all_cds_datas[cds_to_compare]
			
		exon_of_cds = list_cds[cds_to_compare]

		del list_cds[cds_to_compare]

		if not (cds_to_compare in used):
			used.append(cds_to_compare)

			cluster_id = "cluster" + str(i)
			i = i + 1
			clusters[cluster_id] = [cds_to_compare]

			for cds_id, exon_id in list_cds.items():
				flag = True
				all_exon_infos2 = all_cds_datas[cds_id]

				if list(all_cds_datas[cds_to_compare].keys()) == list(all_cds_datas[cds_id].keys()):
					for exon_id in all_cds_datas[cds_to_compare].keys():
						splice_5prime = all_cds_datas[cds_to_compare][exon_id][3]		
						splice_3prime = all_cds_datas[cds_to_compare][exon_id][4]						

						splice_5prime2 = all_cds_datas[cds_id][exon_id][3]		
						splice_3prime2 = all_cds_datas[cds_id][exon_id][4]						
						if splice_5prime != splice_5prime2 or splice_3prime or splice_3prime2:
							flag = False
							break

					if flag == True:
						clusters[cluster_id].append(cds_id)
						used.append(cds_id)

	for cluster_id, elts in clusters.items():
		cluster.write(">" + cluster_id + " :\n")
		cluster.write(', '.join(elts) + "\n")

def readTreeFromFile(path):
	fichier = open(path, "r")
	datas = fichier.readlines()[0]
	return str(datas)		


def main(gain, dup, lost, iteration, arg):
	k  = arg.k_nb_exons
	k1 = arg.k_ei_c
	k2 = arg.k_ei_c
	k3 = arg.k_t_c

	number_max_exon_transcript = 0
	number_max_exon_gene = 0
	prob_deletion = arg.k_indel
	prob_insertion = arg.k_indel
	
	tc6 = arg.random_selection
	tc5 = arg.alternative_five_prime
	tc3 = arg.alternative_three_prime
	tcsk = arg.exon_skipping
	tcme = arg.mutually_exclusive
	tcir = arg.intron_retention
	tctl = arg.transcript_loss
	k_intron = arg.k_intron

	guideTreeFile = arg.input_tree_file

	one_nt_dist = {}
	two_nt_dist = {}
	three_nt_dist = {}
	one_nt_dist_nt = []
	two_nt_dist_nt = {}
	three_nt_dist_nt = {}
	mean_exon_cds = 0
	sd_exon_cds = 0
	mean_exon_gene = 0
	sd_exon_gene = 0
	mean_exon_len = 0
	sd_exon_len = 0
	exon_lenght = 0
	exons_dict = {}
	exon_lenghts = []
	mean_len_introns = 0
	sd_len_intons = 0
	nt_inton_dist = []
	intron_dict = {}
	exon_loss_dup_gain = {}
	prob_loss = 0.4
	prob_gain = 0.5
	prob_dup = 0.1
	exon_to_gain = []
	mean_nb_indel = 1.5
	sd_len_indel = 1
	mean_len_indel = 5
	dict_indel_exons = {}
	dup_cost = 40
	#src = "output/"
	src2 = "_" + iteration + "_"
	tmp = guideTreeFile.split("/")[-1]
	dirName = tmp.split(".")[0]
	src = "Example/output/" + dirName + "/"
	try:
		os.mkdir("Example/output/" + dirName)		
		os.mkdir("Example/output/" + dirName + "/cds")
		os.mkdir("Example/output/" + dirName + "/genes")
		os.mkdir("Example/output/" + dirName + "/cds_gene")
		os.mkdir("Example/output/" + dirName + "/cluster")
		os.mkdir("Example/output/" + dirName + "/positions")
		os.mkdir("Example/output/" + dirName + "/pairwise_alignment")
		os.mkdir("Example/output/" + dirName + "/multiple_alignment")		

	except FileExistsError:
		pass
		print("Directory " , dirName ,  " already exists")

	guideTree = readTreeFromFile(guideTreeFile)

	tree = Tree(guideTree)

	nb_leaf = len(tree.get_leaf_names())

	codonMatrixFile = "Data/transitionCodonMatrix.txt"

	codonslist = ['TTT', 'TTC', 'TTA', 'TTG', 'TCT', 'TCC', 'TCA', 'TCG', 'TAT', 'TAC', 'TGT', 'TGC', 'TGG', 'CTT',
				  'CTC', 'CTA', 'CTG', 'CCT', 'CCC', 'CCA', 'CCG', 'CAT', 'CAC', 'CAA', 'CAG', 'CGT', 'CGC', 'CGA',
				  'CGG', 'ATT', 'ATC', 'ATA', 'ATG', 'ACT', 'ACC', 'ACA', 'ACG', 'AAT', 'AAC', 'AAA', 'AAG', 'AGT',
				  'AGC', 'AGA', 'AGG', 'GTT', 'GTC', 'GTA', 'GTG', 'GCT', 'GCC', 'GCA', 'GCG', 'GAT', 'GAC', 'GAA',
				  'GAG', 'GGT', 'GGC', 'GGA', 'GGG']   

	mean_len_introns, sd_len_intons, mean_exon_cds, sd_exon_cds, mean_exon_gene, sd_exon_gene, mean_exon_len, sd_exon_len, one_nt_dist,  two_nt_dist, three_nt_dist, nt_inton_dist = init(one_nt_dist, two_nt_dist, three_nt_dist)

	one_nt_dist_nt, two_nt_dist_nt, three_nt_dist_nt, nt_inton_dist_nt = generate_nt_for_random(one_nt_dist,
																								two_nt_dist,
																								three_nt_dist,
																								nt_inton_dist)
	#generation du nombre d'exon du géne racine et d'une CDS
	number_max_exon_transcript = generate_number_max_exon_transcript(mean_exon_cds, sd_exon_cds)
	number_max_exon_gene = generate_number_max_exon_gene(number_max_exon_transcript, k)


	#generation de la longuer des exons de chaque exon du gène racine, qui sont multiples de 3   
	for i in range(number_max_exon_gene):
		exon_len = generate_number_by_normal_law(mean_exon_len, sd_exon_len)
		exon_len = exon_len - exon_len % 3 + 3
		exon_lenghts.append(exon_len)


	end_exon = ""
	for i in range(len(exon_lenghts)):
		lenght = exon_lenghts[i]
		exon, end_exon = generate_exon(lenght, one_nt_dist_nt, two_nt_dist_nt, three_nt_dist_nt, end_exon)
		exons_dict["exon_" + str(i)] = exon
		#exon_loss_dup_gain est un dict qui pour chaque exon à chaque noeud de l'arbre contient tous les evns de son histoire d'évolution 
		exon_loss_dup_gain["exon_" + str(i)] = [0, [0, ""], 0]

	adjust_exons_extremities(exons_dict, end_exon)

	#generation d'un n*n exon qui seront utilisé lors des evns de gains d'exons
	for i in range(len(tree.get_leaf_names()) * len(tree.get_leaf_names())*10):		
		lenght = generate_number_by_normal_law(mean_exon_len, sd_exon_len)
		if lenght % 3 == 1:
			lenght += 2
		elif lenght % 3 == 2:
			lenght += 1
		exon_to_gain.append(generate_exon(lenght, one_nt_dist_nt, two_nt_dist_nt, three_nt_dist_nt, end_exon)[0])

	
	listKeyOrder = write_init(exons_dict, "exon")

	#generation de intron, les splices sites peuvent varier
	for i in range(nb_leaf * len(exon_lenghts) - 1):
		intron_lenght = generate_number_by_normal_law(mean_len_introns, sd_len_intons)
		intron = generate_intron(intron_lenght, nt_inton_dist_nt)
		intron_dict["intron_" + str(i)] = intron

	write_init(intron_dict, "intron")


	make_structure_evol(tree, exon_loss_dup_gain, prob_loss, prob_dup, prob_gain, exons_dict,
						exon_to_gain, listKeyOrder, k2)


	codon_matrix = initCodonMatrix(codonMatrixFile)

	codon_distribution = makeDistribution(codon_matrix, codonslist)

	make_substition_evol(tree, exons_dict, codon_distribution, listKeyOrder, k1)

	exon_to_gain = "".join(exon_to_gain)

	indel_exon_state = make_indel_evol(tree, exons_dict, codon_distribution, listKeyOrder, mean_nb_indel, mean_len_indel, sd_len_indel, exon_to_gain,
					dict_indel_exons,k1, prob_insertion, prob_deletion)


	exon_intron_of_genes, resulting_gene_exon, genes_seq= write_gene_in_file(tree, indel_exon_state, intron_dict, src, src2)
			


	make_protein_evol(tree, exon_intron_of_genes, resulting_gene_exon, tc5, tc3, tcsk, tcme, tc6, k3, intron_dict, src, src2, genes_seq)



if __name__ == "__main__":
	parser = build_arg_parser()
	arg = parser.parse_args()

	gain = arg.exon_i_change_gain
	dup = arg.exon_i_change_dup
	lost = arg.exon_i_change_loss

	cost = [gain, dup, lost]
	cmpt = 1
	
	n = arg.number_of_simulation

	for i in range(10):
		main(cost[0], cost[1], cost[2], "iteration_" + str(cmpt), arg) 
		cmpt += 1


