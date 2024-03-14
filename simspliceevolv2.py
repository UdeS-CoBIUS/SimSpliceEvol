""" SimSpliceEvol is a tool designed to simulate the evolution of sets of alternative transcripts along the branches of an input gene tree. In addition to traditional sequence evolution events, the simulation also incorporates events related to the evolution of gene exon-intron structures and alternative splicing. These events modify the sets of transcripts produced from genes. Data generated using SimSpliceEvol is valuable for testing spliced RNA sequence analysis methods, including spliced alignment of cDNA and genomic sequences, multiple cDNA alignment, identification of orthologous exons, splicing orthology inference, and transcript phylogeny inference. These tests are essential for methods that require knowledge of the real evolutionary relationships between the sequences.
> Usage:
======
    python3 simspliceevolv2.py [-h] -i INPUT_TREE_FILE [-it ITERATIONS]
                          [-dir_name DIRECTORY_NAME] [-eic_el EIC_EL]
                          [-eic_ed EIC_ED] [-eic_eg EIC_EG] [-c_i C_I]
                          [-c_d C_D] [-k_nb_exons K_NB_EXONS] [-k_eic K_EIC]
                          [-k_indel K_INDEL] [-k_tc K_TC]
                          [-tc_a5 ALTERNATIVE_FIVE_PRIME]
                          [-tc_a3 ALTERNATIVE_THREE_PRIME]
                          [-tc_es EXON_SKIPPING] [-tc_me MUTUALLY_EXCLUSIVE]
                          [-tc_ir INTRON_RETENTION] [-tc_tl TRANSCRIPT_LOSS]
> Reference:
======
    https://github.com/dondavy/SimSpliceEvol
"""
__authors__ = ("Wend Yam Donald Davy Ouedraogo")
__contact__ = ("wend.yam.donald.davy.usherbrooke.ca")
__copyright__ = "CoBIUS lab at Université de Sherbrooke, QC, CANADA"
__date__ = "2024-02-26"
__version__= "2.0.2"
__previousAuthor__ = ("Esaie Kuitche, PhD-2019")


import argparse
import random
import numpy as np
from ete3 import Tree, TreeStyle, NodeStyle, faces, TextFace
import itertools
import copy
import os
import math
import timeout_decorator
#import concurrent.futures
#import string
import pandas as pd
import random
import networkx as nx
#import matplotlib.pyplot as plt
#import matplotlib.colors as mcolors
#from networkx.drawing.nx_agraph import write_dot

time_transcripts_generation= 30

class colorsSimulation:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def build_arg_parser():
    '''The parser function'''
    parser = argparse.ArgumentParser(description="SimSpliceEvol program parameters")
    parser.add_argument('-i', '--input_tree_file', required=True,  help="input guide tree (default =Example/input/small.nw)", default="execution/inputs/small.nw")
    parser.add_argument('-dir_name', '--directory_name', default = '.', help="the output folder. By default, the current directory is used.")
    parser.add_argument('-it', '--iteration', default = '1', help="the output sub-folder and the name of the simulation.") 
    
    #relative frequencies of gain, duplication and loss of exons
    parser.add_argument('-eic_el', '--eic_el', default = 0.4, help="relative frequencies of exon-intron structure change by exon loss(default=0.1)")
    parser.add_argument('-eic_ed', '--eic_ed', default = 0.1,  help = "relative frequencies of exon-intron structure change by exon duplication(default=0.8)")
    parser.add_argument('-eic_eg', '--eic_eg', default = 0.5, help="relative frequencies of exon-intron structure change by exon gain(default=0.1)") 
    
    #relative frequencies of insertion and deletion events
    parser.add_argument('-c_i', '--c_i', default = 0.7, help="relative frequencies of insertion events(default=0.7)")
    parser.add_argument('-c_d', '--c_d', default = 0.3,  help = "relative frequencies of deletion events(default=0.3)")
    
    #multiplicative constance (depending on branches lengths)
    parser.add_argument('-k_nb_exons', '--k_nb_exons', default = 1.5, help="multiplicative constant for number of exons in gene (default =1.5)")
    parser.add_argument('-k_eic', '--k_eic',  help="multiplicative constant for exon-intron change (eic) rate  (default=25)",  default =0.5) 
    parser.add_argument('-k_indel', '--k_indel', help="multiplicative constant for codon indel rate (default = 5)", default = 0.5)
    parser.add_argument('-k_tc', '--k_tc', help="multiplicative constant for transcript change (default =5)", default = 3)
    
    # parameters for the internal nodes and leaves
    parser.add_argument('-tc_rs', '--random_selection', default =  0.0, help="relative frequence of random selection (default =1.0)")
    parser.add_argument('-tc_a5', '--alternative_five_prime', default =  0.5, help="relative frequence of alternative five prime in tc (default =0.1)")
    parser.add_argument('-tc_a3', '--alternative_three_prime', default = 0.5, help="relative frequence of alternative three prime in tc (default =0.1)")
    parser.add_argument('-tc_es', '--exon_skipping', default = 0.01, help="relative frequence of exon skipping in tc (default =0.5)")
    parser.add_argument('-tc_me', '--mutually_exclusive', default = 0.3, help="relative frequence of mutually exclusive in tc (default =0.15)")
    parser.add_argument('-tc_ir', '--intron_retention', default = 0.0, help="relative frequence of intron retention in tc (default =0.00)")
    parser.add_argument('-tc_tl', '--transcript_loss', default = 0.02, help="relative frequence of transcript loss in tc (default =0.4)")
        
    return parser

#######################################
##########  UTILITIES    ##############
#######################################

def random_generator_nucleotide(ONE_NT_DIST, TWO_NT_DIST, THREE_NT_DIST, NT_INTRON_DIST):
    '''RANDOM NUCLEOTIDE GENERATION'''
    all_one_nt_dist = list("".join([key * value for key, value in ONE_NT_DIST.items()]))
    random.shuffle(all_one_nt_dist)
    
    all_nt_intron_dist = list("".join([key * value for key, value in NT_INTRON_DIST.items()]))
    random.shuffle(all_nt_intron_dist)
    
    all_two_nt_dist = {}
    for k, v in TWO_NT_DIST.items():
        tmp_list = list("".join([key * value for key, value in v.items()]))
        random.shuffle(tmp_list)
        all_two_nt_dist[k] = tmp_list

    all_three_nt_dist = {}
    for k, v in THREE_NT_DIST.items():
        tmp_list = list("".join([key * value for key, value in v.items()]))
        random.shuffle(tmp_list)
        all_three_nt_dist[k] = tmp_list
    return all_one_nt_dist, all_two_nt_dist, all_three_nt_dist, all_nt_intron_dist

def random_number_generator_by_normal_law(mu, sigma):
    '''random_number_generator_by_normal_law'''
    val = int(np.random.normal(mu, sigma, 1)[0])
    if val > 0:
        return val
    else:
        return random_number_generator_by_normal_law(mu, sigma)
    
def random_generator_max_exon_transcript(MEAN_EXON_CDS, SD_EXON_CDS):
    '''RANDOM GENERATOR the number max of exons in transcripts'''
    return random_number_generator_by_normal_law(MEAN_EXON_CDS, SD_EXON_CDS)
    
def random_generator_max_exon_gene(NUMBER_MAX_EXON_TRANSCRIPT, K_NB_EXONS):
    '''RANDOM GENERATOR the number max of exons in genes'''
    if K_NB_EXONS >= 1:
        return int(round(K_NB_EXONS*NUMBER_MAX_EXON_TRANSCRIPT))
    else:
        return NUMBER_MAX_EXON_TRANSCRIPT

def exon_generator(length, all_one_nt_dist, all_two_nt_dist, all_three_nt_dist, STOP_CODONS, end_exon):
    """GENERATION EXONS SEQUENCES"""
    i = 0
    one = ""
    two = ""
    three = ""
    exon = ""
    if len(end_exon) == 1:
        while True:
            two = random.choice(all_two_nt_dist[end_exon])
            exon = two + random.choice(all_three_nt_dist[end_exon + two])
            if two+exon not in STOP_CODONS:
                break
    elif len(end_exon) == 2:
        while True:
            exon = random.choice(all_three_nt_dist[end_exon])
            if end_exon+exon not in STOP_CODONS:
                break
    while i<length:
        if i % 3 == 0:
            one = random.choice(all_one_nt_dist)
        elif i%3 == 1:
            two = one + random.choice(all_two_nt_dist[one])
        else:
            while True:
                three = two + random.choice(all_three_nt_dist[two])
                if three not in STOP_CODONS:
                    exon = exon + three
                    break
        i = i + 1
    if i % 3 == 0:
        return exon, ""
    elif i % 3 == 1:
        return exon+one, one
    else:
        return exon+two, two

def get_exons_extremities_adjusted(exons_dict, end_exon, STOP_CODONS):
    """ADJUSTEMENTS OF EXONS EXTREMITIES"""
    seq = exons_dict["exon_0"]

    if seq[:3] != "ATG":
        seq = "ATG" + seq
        exons_dict["exon_0"] = seq

    id_last = len(exons_dict.keys()) - 1
    seq = exons_dict["exon_" + str(id_last)]
    if len(end_exon) == 0:
        seq_tmp = seq[:-3] + random.choice(STOP_CODONS)
    else:
        seq_tmp = seq[:-len(end_exon)] + random.choice(STOP_CODONS)
    exons_dict["exon_" + str(id_last)] = seq_tmp

    if len(exons_dict.keys()) > 2:
        i = 1 
        while i < (len(exons_dict.keys())-1):
            if random.choice([True, False]):
                exons_dict["exon_" + str(i)] = "ATG" + exons_dict["exon_" + str(i)]
            i += 1
    return True

def write_init_structure(exons_dictionary, introns_dictionary):
    """Return the gene structure"""
    keys_exons = exons_dictionary.keys()
    keys_introns = introns_dictionary.keys()
    gene_structure = []

    for i in range(len(keys_exons)):
        name_exon = 'exon_{}'.format(i)
        if name_exon in keys_exons:
            gene_structure.append(name_exon)
            if i != len(keys_exons)-1:
                name_intron = 'intron_{}'.format(i)
                if name_intron in keys_introns:
                    gene_structure.append(name_intron)
    return gene_structure

def intron_generator(LENGTH, SPLICE_SITES, all_nt_intron_dist):
    i = 0
    splicesites = list(itertools.chain.from_iterable(SPLICE_SITES))
    splicesite = random.choice(splicesites)

    if len(splicesite)>0:
        splicesitetmp = splicesite.split("-")
        donnor = splicesitetmp[0]
        acceptor = splicesitetmp[1]
    else:
        donnor = random.choice(["A", "G", "C", "T"]) + "" + random.choice(["A", "G", "C", "T"])
        acceptor = random.choice(["A", "G", "C", "T"]) + "" + random.choice(["A", "G", "C", "T"])

    intron = donnor

    while i < (LENGTH - 4):
        intron += random.choice(all_nt_intron_dist)
        i += 1

    intron += acceptor   
    return intron

def make_distribution(CODON_MATRIX, CODONS_LIST):
    codons_dict = {}
    for i in range(len(CODONS_LIST)):
        l = []
        for j in range(len(CODONS_LIST)):
            l += [CODONS_LIST[j]] * int(CODON_MATRIX[i][j])
            codons_dict[CODONS_LIST[i]] = l
    codons_dict['TAG'] = ['TAG']
    codons_dict['TAA'] = ['TAA']
    codons_dict['TGA'] = ['TGA']
    return codons_dict

@timeout_decorator.timeout(time_transcripts_generation)
def transcripts_generation_ancestor(number, number_max_exons, exons_dict, tc_rs):
    exons_range = [int(str(_).split('_')[-1]) for _ in exons_dict.keys()]
    transcripts_already_picked = []
    transcripts_dict = {}
    if tc_rs == 0:
        exons_range.sort()
        transcripts_already_picked.append(exons_range)
    else:
        for transcript_number in range(number):
            exons_number_included = np.random.randint(1, number_max_exons)
            #choose exons
            not_found = True
            while (not_found):
                exons_got_choosen = random.sample(exons_range, exons_number_included)
                exons_got_choosen.sort()
                if exons_got_choosen not in transcripts_already_picked:
                    transcripts_already_picked.append(exons_got_choosen)
                    not_found = False
    
    for i, transcript_composition in enumerate(transcripts_already_picked):
        transcript_name = 'transcript_{}-{}'.format(0, i)
        transcripts_dict[transcript_name] = ['exon_{}'.format(numb) for numb in transcript_composition]

    return transcripts_dict

@timeout_decorator.timeout(time_transcripts_generation)
def transcripts_generation_nodes(number, number_max_exons, exons_list, current_list_of_transcript, gene_number, number_new_transcript):
    transcripts_already_picked = []
    for transcript_number in range(number):
        #choose exons
        not_found = True
        #print(number_max_exons, number_new_transcript)
        while (not_found):
            exons_number_included = np.random.randint(1, number_max_exons)
            exons_got_choosen = random.sample(exons_list, exons_number_included)
            exons_composition = [_ for _ in exons_list if _ in exons_got_choosen]
            if exons_composition not in transcripts_already_picked and exons_composition not in current_list_of_transcript:
                transcripts_already_picked.append(exons_composition)
                not_found = False
    transcripts_dict = {}
    for i, transcript_composition in enumerate(transcripts_already_picked):
        transcript_name = 'transcript_{}-{}'.format(gene_number, i+number_new_transcript)
        transcripts_dict[transcript_name] = transcript_composition

    return transcripts_dict

def choice_distinct(all_data, k):
    choices = []
    if k > len(all_data):
        k = len(all_data)
    while len(choices) < k:
        selection = random.choice(all_data)
        if selection not in choices:
            choices.append(selection)
    return choices

def placement_exons_in_evolution(event, index, new_exon_name, current_list, start_index, stop_index, exons_list):
    if new_exon_name not in current_list:
        new_list = current_list
        if event == 'dup':
            new_list.insert(index+1, new_exon_name)
            exons_list.append(new_exon_name)
        elif event == 'new':
            if start_index == stop_index-1:
                placement_exons_in_evolution('dup', start_index, new_exon_name, current_list, 'null', 'null', exons_list)
            else:
                choosed_place = random.randint(start_index+1, stop_index)
                new_list.insert(choosed_place, new_exon_name)
                exons_list.append(new_exon_name)
        return new_list, exons_list
    else:
        return current_list, exons_list

def update_transcripts_forest(current_forest, transcript, name_transcript, lca_events, event):
    #print(current_forest)
    is_new = True
    for i, tree in enumerate(current_forest):
        try:
            t = tree&transcript
            tmp_name = t.name
            t.name = t.name + '_AS_{}'.format(name_transcript)
            tmp_1 = Tree('{};'.format(name_transcript))
            tmp_2 = Tree('{};'.format(tmp_name))
            t.add_child(tmp_1)
            t.add_child(tmp_2)
            lca_events[t.name] = event
            is_new = False
        except:
            pass
            #raise ValueError('not found in the tree pos: {}'.format(i))
    if is_new:
        tmp_tree = Tree('{};'.format(transcript))
        current_forest.append(tmp_tree)
    
    return current_forest, lca_events
########################################
##########  SIMULATION    ##############
######################################## 
def simulation_program(tree, 
                       ancestral_gene_structure, exons_dict, introns_dict, transcripts_dict,
                       EIC_EL, EIC_ED, EIC_EG, C_I, C_D, K_EIC, K_INDEL,
                       TC_RS, TC_TL, TC_ES, TC_ME, TC_A3, TC_A5, TC_IR, K_TC,
                       all_one_nt_dist, all_two_nt_dist, all_three_nt_dist, STOP_CODONS,
                       MEAN_EXON_LEN, SD_EXON_LEN,
                       all_nt_intron_dist, MEAN_INTRON_LEN, SD_INTRON_LEN, SPLICE_SITES,
                       CODONS_MATRIX, CODONS_LIST
                       ):
    number_gain = 0
    number_new_transcript = len(transcripts_dict.keys())
    CODONS_DISTRIBUTION = make_distribution(CODONS_MATRIX, CODONS_LIST)
    CODONS_USED = []
    for codon in CODONS_DISTRIBUTION.keys():
        if (codon not in STOP_CODONS and codon !='ATG'):
            values_codons = CODONS_DISTRIBUTION[codon]
            for value in values_codons:
                CODONS_USED.append(value)
                
    # Transcripts evolution forest
    transcripts_names = list(transcripts_dict.keys())
    current_tr_trees = []
    lca_events = {}
    for transcript in transcripts_names:
        tmp_tree = Tree('{};'.format(transcript))
        current_tr_trees.append(tmp_tree)
    parent_transcripts_forest = current_tr_trees
    
    #Simulation begins ...
    for preorder_number, node in enumerate(tree.traverse('preorder')):
        if node.is_root():
            current_gene_structure = [_ for _ in ancestral_gene_structure if _.startswith('exon')]
            current_positions_of_all_exons_in_the_simulation = [_ for _ in ancestral_gene_structure if _.startswith('exon')]
            dict_tr = {}
            new_dict_transcripts_in_node = {}
            for transcript in transcripts_dict.keys():
                exons = transcripts_dict[transcript]
                dict_tr[transcript] = {}
                new_exons_up = []
                for exon in exons:
                    dict_tr[transcript][exon] = exons_dict[exon]
                    new_exons_up.append(exon)
                new_dict_transcripts_in_node[transcript] = new_exons_up

            # defined variables
            new_dict_transcripts_in_node_names = list(transcripts_dict.keys())
            new_dict_transcripts_in_node_sequences = dict_tr
            current_exons_dict = exons_dict
            #print('vvvvvvvvv')
            #print(new_dict_transcripts_in_node_names)
            #### ES
            
            #number of transcripts that have undergone exon skipping event
            NUMBER_TRANSCRIPTS_ES = int(math.ceil(TC_ES *  len(new_dict_transcripts_in_node_names)))
            #print(NUMBER_TRANSCRIPTS_ES, K_TC, TC_ES)
            #NUMBER_TRANSCRIPTS_ES = 2 #example
            
            # choose the transcripts affected
            transcripts_choosed = choice_distinct(new_dict_transcripts_in_node_names, NUMBER_TRANSCRIPTS_ES)
            
            #Apply the ES event on these transcripts
            for transcript in transcripts_choosed:
                number_exons = len(new_dict_transcripts_in_node[transcript])
                if number_exons == 1:
                    # transcript is deleted
                    pass
                else:
                    # which exon
                    ch_exon_pos = np.random.randint(0, number_exons-1)
                    exon_to_skip = new_dict_transcripts_in_node[transcript][ch_exon_pos]
                    # which transcript
                    name_transcript = 'transcript_{}-{}'.format(preorder_number, number_new_transcript)
                    number_new_transcript += 1
                    # m.a.j
                    new_dict_transcripts_in_node[name_transcript] = [exon for exon in new_dict_transcripts_in_node[transcript] if exon != exon_to_skip]
                    new_dict_transcripts_in_node_sequences[name_transcript] = {}
                    for exon in new_dict_transcripts_in_node[transcript]:
                        if exon != exon_to_skip:
                            new_dict_transcripts_in_node_sequences[name_transcript][exon] = new_dict_transcripts_in_node_sequences[transcript][exon]
                    new_dict_transcripts_in_node_names.append(name_transcript)
                    #update_phylogenies
                    parent_transcripts_forest, lca_events = update_transcripts_forest(parent_transcripts_forest, transcript, name_transcript, lca_events, 'es')
            
            #### ME
            
            #number of transcripts that have undergone a ME event
            NUMBER_TRANSCRIPTS_ME = int(math.ceil(TC_ME *  len(new_dict_transcripts_in_node_names)))
            #print(NUMBER_TRANSCRIPTS_ES, K_TC, TC_ES)
            #NUMBER_TRANSCRIPTS_ME = 2 #example
            
            # choose the transcripts affected
            transcripts_choosed = choice_distinct(new_dict_transcripts_in_node_names, NUMBER_TRANSCRIPTS_ME)
            
            # Apply the ME event on these transcripts
            for transcript in transcripts_choosed:
                number_exons = len(new_dict_transcripts_in_node[transcript])
                if number_exons < 4:
                    # cannot apply ME event on this transcript
                    pass
                else:
                    exons = new_dict_transcripts_in_node[transcript]
                    choosed_adj_sub = np.random.randint(1, len(exons)-2)
                    choosed_exons = [exons[choosed_adj_sub], exons[choosed_adj_sub+1]]
                    exons_compositions = []
                    isOk = True
                    for exon_me in choosed_exons:
                        tmp_list_transcript_desc = []
                        for exon in exons:
                            if exon != exon_me:
                                tmp_list_transcript_desc.append(exon)
                        for transcript_i in new_dict_transcripts_in_node.keys():
                            compo = new_dict_transcripts_in_node[transcript_i]
                            if compo == tmp_list_transcript_desc:
                                isOk = False
                                break
                        if not isOk:
                            break
                        else:
                            exons_compositions.append(tmp_list_transcript_desc)
                    if isOk:
                        for transcript_composition in exons_compositions:
                            # which transcript
                            name_transcript = 'transcript_{}-{}'.format(preorder_number, number_new_transcript)
                            number_new_transcript += 1
                            
                            # m.a.j
                            new_dict_transcripts_in_node[name_transcript] = transcript_composition
                            new_dict_transcripts_in_node_sequences[name_transcript] = {}
                            for exon in transcript_composition:
                                new_dict_transcripts_in_node_sequences[name_transcript][exon] = new_dict_transcripts_in_node_sequences[transcript][exon]
                            new_dict_transcripts_in_node_names.append(name_transcript)
                            #update_phylogenies
                            parent_transcripts_forest, lca_events = update_transcripts_forest(parent_transcripts_forest, transcript, name_transcript, lca_events, 'me')
            
            #### 5SS
            # number of transcripts affected
            NUMBER_TRANSCRIPTS_5SS = int(math.ceil(TC_A5 * len(new_dict_transcripts_in_node_names)))
            #NUMBER_TRANSCRIPTS_5SS = 2 #example
            #choose transcripts affected
            transcripts_choosed = choice_distinct(new_dict_transcripts_in_node_names, NUMBER_TRANSCRIPTS_5SS)
            # Apply 5SS event on these transcripts
            for transcript in transcripts_choosed:
                list_exons = new_dict_transcripts_in_node[transcript]
                number_exons = len(list_exons)
                if number_exons == 1:
                    # cannot apply on this transcript
                    pass
                else:
                    exon_five_prime = np.random.choice(list_exons[1:])
                    sequence_exon = new_dict_transcripts_in_node_sequences[transcript][exon_five_prime]
                    if len(sequence_exon) > 54:
                        range_to_choose = 9 # can be modifying
                        sequence_5_prime = sequence_exon[:range_to_choose+1]
                        if '*' in sequence_5_prime:
                            continue
                        else:
                            sequence_new_exon = ''.join(['-' for i in range(range_to_choose+1)]) + 'AG' +sequence_exon[range_to_choose+3:]
                            '''
                            sequence_new_exon = ''.join(['-' for i in range(range_to_choose+1)]) + 'AG' +sequence_exon[range_to_choose+3:]
                            pos_x = len(''.join(['-' for i in range(range_to_choose+1)])) + 1
                            pos_y = pos_x + 1
                            seq_tmp = list(sequence_exon)
                            seq_tmp[pos_x] = 'A'
                            seq_tmp[pos_y] = 'G'
                            '''
                            #print(len(sequence_new_exon), len(sequence_exon))
                            # which transcript
                            name_transcript = 'transcript_{}-{}'.format(preorder_number, number_new_transcript)
                            number_new_transcript += 1
                            # m.a.j
                            new_dict_transcripts_in_node[name_transcript] = new_dict_transcripts_in_node[transcript]
                            new_dict_transcripts_in_node_sequences[name_transcript] = {}
                            for exon in list_exons:
                                if exon == exon_five_prime:
                                    new_dict_transcripts_in_node_sequences[name_transcript][exon] = sequence_new_exon
                                    #current_exons_dict[exon] = ''.join(seq_tmp)
                                    current_exons_dict[exon] = sequence_new_exon
                                else:
                                    new_dict_transcripts_in_node_sequences[name_transcript][exon] = new_dict_transcripts_in_node_sequences[transcript][exon]
                            new_dict_transcripts_in_node_names.append(name_transcript)
                            #update_phylogenies
                            parent_transcripts_forest, lca_events = update_transcripts_forest(parent_transcripts_forest, transcript, name_transcript, lca_events, '5ss')
            
            
            #### 3SS
            # number of transcripts affected
            NUMBER_TRANSCRIPTS_3SS = int(math.ceil(TC_A3 * len(new_dict_transcripts_in_node_names)))
            #NUMBER_TRANSCRIPTS_3SS = 2 #example
            #choose transcripts affected
            transcripts_choosed = choice_distinct(new_dict_transcripts_in_node_names, NUMBER_TRANSCRIPTS_3SS)
            
                
            # Apply 3SS event on these transcripts
            for transcript in transcripts_choosed:
                list_exons = new_dict_transcripts_in_node[transcript]
                number_exons = len(list_exons)
                if number_exons == 1:
                    # cannot apply on this transcript
                    pass
                else:
                    exon_three_prime = np.random.choice(list_exons[:-1])
                    sequence_exon = new_dict_transcripts_in_node_sequences[transcript][exon_three_prime]
                    if len(sequence_exon) > 54:
                        range_to_choose = 9 # can be modifying
                        sequence_3_prime = sequence_exon[-range_to_choose:]
                        if '*' in sequence_3_prime:
                            continue
                        else:
                            sequence_new_exon = sequence_exon[:-range_to_choose-2] + 'GT' + ''.join(['-' for i in range(range_to_choose)])
                            '''
                            sequence_new_exon = sequence_exon[:-range_to_choose+2] + 'GT' + ''.join(['-' for i in range(range_to_choose)]) 
                            pos_x = len(sequence_exon[:-range_to_choose+2]) + 1
                            pos_y = pos_x + 1
                            seq_tmp = list(sequence_exon)
                            seq_tmp[pos_x] = 'G'
                            seq_tmp[pos_y] = 'T'
                            '''
                            #print(len(sequence_new_exon), len(sequence_exon))
                            # which transcript
                            name_transcript = 'transcript_{}-{}'.format(preorder_number, number_new_transcript)
                            number_new_transcript += 1
                            # m.a.j
                            new_dict_transcripts_in_node[name_transcript] = new_dict_transcripts_in_node[transcript]
                            new_dict_transcripts_in_node_sequences[name_transcript] = {}
                            for exon in list_exons:
                                if exon == exon_three_prime:
                                    new_dict_transcripts_in_node_sequences[name_transcript][exon] = sequence_new_exon
                                    #current_exons_dict[exon] = ''.join(seq_tmp)
                                    current_exons_dict[exon] = sequence_new_exon
                                else:
                                    new_dict_transcripts_in_node_sequences[name_transcript][exon] = new_dict_transcripts_in_node_sequences[transcript][exon]
                            new_dict_transcripts_in_node_names.append(name_transcript)
                            #update_phylogenies
                            parent_transcripts_forest, lca_events = update_transcripts_forest(parent_transcripts_forest, transcript, name_transcript, lca_events, '3ss')
            
            #### IR
            
            # number of transcripts affected
            NUMBER_TRANSCRIPTS_IR = int(math.ceil(TC_IR * len(new_dict_transcripts_in_node_names)))
            #NUMBER_TRANSCRIPTS_IR = 2 #example
            
            #choose transcripts_affected
            transcripts_choosed = choice_distinct(new_dict_transcripts_in_node_names, NUMBER_TRANSCRIPTS_IR)
            
            # Apply IR on these transcripts
            exons_news = []
            for transcript in transcripts_choosed:
                list_exons = new_dict_transcripts_in_node[transcript]
                number_exons = len(list_exons)
                
                if number_exons <= 2:
                    # cannot apply
                    pass
                else:
                    #intron retained
                    choosed_location = np.random.randint(1, number_exons-1)
                    pos_insertion_at_tr_level = choosed_location
                    exon_of_interest = list_exons[pos_insertion_at_tr_level]
                    length = int(random_number_generator_by_normal_law(MEAN_INTRON_LEN, SD_INTRON_LEN)/3.0)
                    new_exon_sequence = intron_generator(length, SPLICE_SITES, all_nt_intron_dist)
                    number_gain += 1
                    #pos_insertion = random.randint(1, len(current_gene_structure)-1)
                    #gene
                    pos_insertion = current_gene_structure.index(exon_of_interest)
                    start_exon = current_gene_structure[pos_insertion-1]
                    stop_exon = current_gene_structure[pos_insertion]
                    current_positions_of_all_exons_in_the_simulation, exons_new = placement_exons_in_evolution('new', 'null', "exon_#"+str(number_gain)+"-intron", current_positions_of_all_exons_in_the_simulation, current_positions_of_all_exons_in_the_simulation.index(start_exon), current_positions_of_all_exons_in_the_simulation.index(stop_exon), exons_new)
                    if "exon_#"+str(number_gain)+"-intron" in exons_new:
                        current_gene_structure.insert(pos_insertion, "exon_#"+str(number_gain)+"-intron")
                        current_exons_dict["exon_#"+str(number_gain)+"-intron"] = new_exon_sequence
                    
                    #transcript
                    # m.a.j
                    # which transcript
                    name_transcript = 'transcript_{}-{}'.format(preorder_number, number_new_transcript)
                    number_new_transcript += 1
                    update_composition= []
                    for exon in current_gene_structure:
                        if exon in new_dict_transcripts_in_node[transcript] or exon in exons_new:
                            update_composition.append(exon)
                    new_dict_transcripts_in_node[name_transcript] = update_composition
                    new_dict_transcripts_in_node_sequences[name_transcript] = {}
                    for exon in update_composition:
                        if exon in exons_new:
                            new_dict_transcripts_in_node_sequences[name_transcript][exon] = current_exons_dict[exon]
                        else:
                            new_dict_transcripts_in_node_sequences[name_transcript][exon] = new_dict_transcripts_in_node_sequences[transcript][exon]
                    new_dict_transcripts_in_node_names.append(name_transcript)
                    #update_phylogenies
                    parent_transcripts_forest, lca_events = update_transcripts_forest(parent_transcripts_forest, transcript, name_transcript, lca_events, 'ir')     
                    
            # introns generation
            new_node_structure = []
            current_introns_dict = {}
            for i, exon in enumerate(current_gene_structure):
                new_node_structure.append(exon)
                if i != len(current_gene_structure)-1:
                    intron_name = 'intron_{}'.format(i)
                    new_node_structure.append(intron_name)
                    length = random_number_generator_by_normal_law(MEAN_INTRON_LEN, SD_INTRON_LEN)
                    intron_sequence = intron_generator(length, SPLICE_SITES, all_nt_intron_dist)
                    current_introns_dict[intron_name] = intron_sequence   
            '''
            node.add_features(
                gene_name='gene_0',
                gene_structure=copy.deepcopy(ancestral_gene_structure),
                exons_dict=copy.deepcopy(exons_dict),
                introns_dict=copy.deepcopy(introns_dict),
                transcripts_dict=copy.deepcopy(transcripts_dict),
                transcripts_sequences_dict = copy.deepcopy(dict_tr),
                all_exons_in_the_simulation=[exon for exon in ancestral_gene_structure if exon.startswith('exon')]
                transcripts_forest=copy.deepcopy(current_tr_trees)
            )
            '''
            #print('vvvvvvvvv')
            #print(new_dict_transcripts_in_node_names)
            node.add_features(
                gene_name='gene_0',
                gene_structure=copy.deepcopy(new_node_structure),
                exons_dict=copy.deepcopy(current_exons_dict),
                introns_dict=copy.deepcopy(current_introns_dict),
                transcripts_dict=copy.deepcopy(new_dict_transcripts_in_node),
                transcripts_sequences_dict = copy.deepcopy(new_dict_transcripts_in_node_sequences),
                all_exons_in_the_simulation=copy.deepcopy(current_positions_of_all_exons_in_the_simulation)
                #transcripts_forest=copy.deepcopy(current_tr_trees)
            )
            #pass
        else:
            ###################################################################################
            ##############################  STRUCTURE EVOLUTION  ##############################
            ###################################################################################
            
            #gene structure evolution
            ##parent and global information
            current_positions_of_all_exons_in_the_simulation = copy.deepcopy(node.up.all_exons_in_the_simulation)
            parent_structure = copy.deepcopy(node.up.gene_structure)
            parent_transcripts = copy.deepcopy(node.up.transcripts_dict)
            parent_transcripts_sequences = copy.deepcopy(node.up.transcripts_sequences_dict)
            #parent_transcripts_forest = copy.deepcopy(node.up.transcripts_forest)
            ##current node information
            current_exons_dict = copy.deepcopy(node.up.exons_dict)
            CODON_SUBST_RATE = node.dist
            
            #rename node if node is an internal node otherwise, keep the node name
            if len(node.name) == 0:
                gene_name = 'gene_{}'.format(preorder_number)
            else:
                gene_name = node.name
                         
            ##loss-dup-gain (EXON LEVEL)
            ###loss event
            exons_parent = [_ for _ in parent_structure if _.startswith('exon')]
            nb_exons_parent = len(exons_parent)
            NB_LOSS = int(nb_exons_parent * EIC_EL * K_EIC * CODON_SUBST_RATE)
            
            #NB_LOSS = 2 #example
            #print(NB_LOSS)
            exons_to_delete = choice_distinct(exons_parent, NB_LOSS)
            current_gene_structure = [_ for _ in exons_parent if _ not in exons_to_delete]
            #print(current_gene_structure)
            for exon_del in exons_to_delete:
                current_exons_dict.pop(exon_del)
            
            ###duplication event
            NB_DUP = int(len(current_gene_structure) * EIC_ED * K_EIC * CODON_SUBST_RATE)
            #NB_DUP = 3 #example
            exons_to_duplicate = choice_distinct(current_gene_structure, NB_DUP)
            exons_dup = []
            for exon in exons_to_duplicate:
                name_exon_dup = exon+'-dup'
                current_positions_of_all_exons_in_the_simulation, exons_dup = placement_exons_in_evolution('dup', current_positions_of_all_exons_in_the_simulation.index(exon), name_exon_dup, current_positions_of_all_exons_in_the_simulation, 0, 0, exons_dup)

            updated_current_gene_structure = [_ for _ in current_positions_of_all_exons_in_the_simulation if _ in current_gene_structure or _ in exons_dup]
            current_gene_structure = updated_current_gene_structure   

            for exon in exons_dup:
                exon_original = exon[0:-4]
                sequence_exon_original = current_exons_dict[exon_original]
                current_exons_dict[exon] = sequence_exon_original.replace('*', '')

            ###gain event
            NB_GAIN = int(len(current_gene_structure) * EIC_EG * K_EIC * CODON_SUBST_RATE)
            #NB_GAIN = 3 #example
            exons_new = []
            for cmpt in range(NB_GAIN):
                if len(current_gene_structure) >= 3:
                    new_exon_len = random_number_generator_by_normal_law(MEAN_EXON_LEN, SD_EXON_LEN)
                    new_exon_len = new_exon_len - new_exon_len % 3 + 3
                    new_exon_sequence = exon_generator(new_exon_len, all_one_nt_dist, all_two_nt_dist, all_three_nt_dist, STOP_CODONS, "")[0]
                    # update current exons dict
                    number_gain += 1
                    pos_insertion = random.randint(1, len(current_gene_structure)-1)
                    start_exon = current_gene_structure[pos_insertion-1]
                    stop_exon = current_gene_structure[pos_insertion]
                    current_positions_of_all_exons_in_the_simulation, exons_new = placement_exons_in_evolution('new', 'null', "exon_#"+str(number_gain)+"-new", current_positions_of_all_exons_in_the_simulation, current_positions_of_all_exons_in_the_simulation.index(start_exon), current_positions_of_all_exons_in_the_simulation.index(stop_exon), exons_new)
                    if "exon_#"+str(number_gain)+"-new" in exons_new:
                        current_gene_structure.insert(pos_insertion, "exon_#"+str(number_gain)+"-new")
                        current_exons_dict["exon_#"+str(number_gain)+"-new"] = new_exon_sequence

            ################################################
            ##########  SEQUENCE EVOLUTION   ###############
            ################################################

            #CODON SUBSTITUTION
            NB_CODONS = int(len("".join([current_exons_dict[exon].replace('*','') for exon in current_gene_structure]))/3.0)
            nb_subst = int(NB_CODONS * CODON_SUBST_RATE * K_EIC)
            subt_exons_codons = [choice_distinct(current_gene_structure, 1)[0] for _ in range(nb_subst)]
            set_subt_exons_codons = set(subt_exons_codons)
            for exon_subt in list(set_subt_exons_codons):
                count_exon_subt = subt_exons_codons.count(exon_subt)
                exon_subt_sequence = current_exons_dict[exon_subt]
                length_exon_subt = len(exon_subt_sequence)
                already_picked_positions = [0, 1, 2, length_exon_subt-1, length_exon_subt-2, length_exon_subt-3]
                subt_positions = []
                for _ in range(count_exon_subt):
                    new_pos = choice_distinct(range(length_exon_subt), 1)[0]
                    if new_pos not in already_picked_positions:
                        if '*' not in list(exon_subt_sequence[new_pos:new_pos + 3]) and '-' not in list(exon_subt_sequence[new_pos:new_pos + 3]) :
                            subt_positions.append(new_pos)
                            already_picked_positions.extend([new_pos, new_pos+1, new_pos+2])
                for position in subt_positions:
                    codon_subtitute = random.choice(CODONS_DISTRIBUTION[exon_subt_sequence[position:position+3]])
                    exon_subt_sequence = exon_subt_sequence[:position] + codon_subtitute + exon_subt_sequence[position+3:]
                current_exons_dict[exon_subt] = exon_subt_sequence

            #INDEL EVOLUTION
            for exon in current_gene_structure:
                exon_sequence = current_exons_dict[exon]
                exon_length = len(exon_sequence)
                NUMBER_OF_CODONS = int(exon_length/3.0)
                # INSERTION
                if exon not in exons_parent:
                    
                    NB_CODONS_INSERTION = int((NUMBER_OF_CODONS-2) * K_INDEL * C_I * CODON_SUBST_RATE)
                    if NB_CODONS_INSERTION > 0:
                        range_start_codons = [_ for _ in range(3,exon_length-3,3)] #to avoid the insertion of a start/stop codon
                        codons_to_insert = choice_distinct(range_start_codons, NB_CODONS_INSERTION)
                        random_codons = []
                        for position in sorted(codons_to_insert):
                            random_codon = random.choices(CODONS_USED)
                            random_codons.append(random_codon[0])
                        inserted_sequence = []
                        for ind, position in enumerate(sorted(codons_to_insert)):
                            if ind == 0:
                                inserted_sequence.append(exon_sequence[0:position])
                                inserted_sequence.append(random_codons[ind])
                            elif ind == (len(sorted(codons_to_insert)) - 1):
                                inserted_sequence.append(exon_sequence[sorted(codons_to_insert)[ind-1]:position])
                                inserted_sequence.append(random_codons[ind])
                                inserted_sequence.append(exon_sequence[position:len(exon_sequence)-3])
                                inserted_sequence.append(exon_sequence[-3:])
                                
                            else:
                                inserted_sequence.append(exon_sequence[sorted(codons_to_insert)[ind-1]:position])
                                inserted_sequence.append(random_codons[ind])
                        current_exons_dict[exon] = ''.join(inserted_sequence)
                        
                # DELETION
                else:
                    NB_CODON_DELETIONS = int(NUMBER_OF_CODONS * K_INDEL * C_D * CODON_SUBST_RATE)
                    if NB_CODON_DELETIONS > 0:
                        range_start_codons = [_ for _ in range(0,exon_length,3)]
                        codons_to_delete = choice_distinct(range_start_codons, NB_CODON_DELETIONS)
                        nucleotides_to_delete = []
                        for codon_to_delete in codons_to_delete:
                            if '*' not in exon_sequence[codon_to_delete:codon_to_delete + 3]:
                                nucleotides_to_delete.append(codon_to_delete) # 1st nucleotide
                                nucleotides_to_delete.append(codon_to_delete + 1) # 2nd nucleotide
                                nucleotides_to_delete.append(codon_to_delete + 2) # 3rd nucleotide
                        new_exon_sequence = ''.join([exon_sequence[i] if i not in nucleotides_to_delete else '*' for i in range(exon_length)])
                        current_exons_dict[exon] = new_exon_sequence
                                  
            
            ##################################################
            ##########  TRANSCRIPT EVOLUTION   ###############
            ##################################################
            
            #print(parent_transcripts.keys())
            parent_transcripts_name = list(parent_transcripts.keys())
            #exons_in_node = [_ for _ in new_node_structure if _.startswith('exon')]
            transcripts_in_node = []
            transcripts_in_node_composition = []
            # CREATE CURRENT TRANSCRIPTS_DICT_SEQUENCES
            ## variables defined
            new_dict_transcripts_in_node_sequences = {}
            ## 
            #for transcript in new_dict_transcripts_in_node_names:
            #    transcript_composition = new_dict_transcripts_in_node[transcript]
            #    new_dict_transcripts_in_node_sequences[transcript] = {}
            #    for exon in transcript_composition:
            #        new_dict_transcripts_in_node_sequences[transcript][exon] = current_exons_dict[exon]     
            #print(current_gene_structure)
            for transcript in parent_transcripts_name:
                transcript_exons = parent_transcripts[transcript]
                isOk = True
                for exon in transcript_exons:
                    if exon not in current_gene_structure:
                        isOk = False
                        break
                if isOk:
                    transcript_name = 'transcript_{}-{}'.format(preorder_number, transcript.split('-')[-1])
                    transcripts_in_node.append(transcript_name)
                    transcripts_in_node_composition.append(transcript_exons)
                    new_dict_transcripts_in_node_sequences[transcript_name] = parent_transcripts_sequences[transcript]
                    # m.a.j trees
                    parent_transcripts_forest, lca_events = update_transcripts_forest(parent_transcripts_forest, transcript, transcript_name, lca_events, 'conservation')
            #print(transcripts_in_node)
            #break
            # Creation of transcripts

            NUMBER_OF_TRANSCRIPTS_GENERATED = random_number_generator_by_normal_law(3, 1)
        
            new_dict_transcripts_in_node = transcripts_generation_nodes(NUMBER_OF_TRANSCRIPTS_GENERATED, len(current_gene_structure), current_gene_structure, transcripts_in_node_composition, preorder_number, number_new_transcript)
            number_new_transcript += NUMBER_OF_TRANSCRIPTS_GENERATED
            for tr in new_dict_transcripts_in_node.keys():
                #update trees phylogenies
                parent_transcripts_forest, lca_events = update_transcripts_forest(parent_transcripts_forest, tr, tr, lca_events, 'creation')
                composition = new_dict_transcripts_in_node[tr]
                new_dict_transcripts_in_node_sequences[tr] = {}
                for exon in composition:
                    exon_sequence = current_exons_dict[exon]
                    new_dict_transcripts_in_node_sequences[tr][exon] = exon_sequence
            for i_n, tr in enumerate(transcripts_in_node):
                new_dict_transcripts_in_node[tr] = transcripts_in_node_composition[i_n]
            
                
            
            #variables defined
            #print(new_dict_transcripts_in_node.keys())
            new_dict_transcripts_in_node_names = list(new_dict_transcripts_in_node.keys())
            
            #print('*********************************************************')
            #print('*********************************************************')
            #print(list(new_dict_transcripts_in_node_sequences.keys()))
            #print(new_dict_transcripts_in_node_names)
            #print(list(new_dict_transcripts_in_node.keys()))
            
            
            #### LOSS OF TRANSCRIPTS
            NUMBER_TRANSCRIPTS_DELETED = int(len(new_dict_transcripts_in_node) * K_TC * TC_TL * CODON_SUBST_RATE)
            #print(NUMBER_TRANSCRIPTS_DELETED, K_TC, TC_TL)
            #NUMBER_TRANSCRIPTS_DELETED = 2 #example
            transcripts_to_delete = choice_distinct(new_dict_transcripts_in_node_names, NUMBER_TRANSCRIPTS_DELETED)
            update_transcripts_in_node = []
            for transcript in new_dict_transcripts_in_node_names:
                if transcript in transcripts_to_delete:
                    new_dict_transcripts_in_node.pop(transcript)
                    new_dict_transcripts_in_node_sequences.pop(transcript)
                if transcript not in transcripts_to_delete:
                    update_transcripts_in_node.append(transcript)
            new_dict_transcripts_in_node_names = update_transcripts_in_node
            #print('*********************************************************')
            #print('*********************************************************')
            #print(list(new_dict_transcripts_in_node_sequences.keys()))
            #print(new_dict_transcripts_in_node_names)
            #print(list(new_dict_transcripts_in_node.keys()))
            #break
            
            #### ES
            
            #number of transcripts that have undergone exon skipping event
            NUMBER_TRANSCRIPTS_ES = int(TC_ES * K_TC * CODON_SUBST_RATE * len(new_dict_transcripts_in_node_names))
            #print(NUMBER_TRANSCRIPTS_ES, K_TC, TC_ES)
            #NUMBER_TRANSCRIPTS_ES = 2 #example
            
            # choose the transcripts affected
            transcripts_choosed = choice_distinct(new_dict_transcripts_in_node_names, NUMBER_TRANSCRIPTS_ES)
            
            #Apply the ES event on these transcripts
            for transcript in transcripts_choosed:
                number_exons = len(new_dict_transcripts_in_node[transcript])
                if number_exons == 1:
                    # transcript is deleted
                    pass
                else:
                    # which exon
                    ch_exon_pos = np.random.randint(0, number_exons-1)
                    exon_to_skip = new_dict_transcripts_in_node[transcript][ch_exon_pos]
                    # which transcript
                    name_transcript = 'transcript_{}-{}'.format(preorder_number, number_new_transcript)
                    number_new_transcript += 1
                    # m.a.j
                    new_dict_transcripts_in_node[name_transcript] = [exon for exon in new_dict_transcripts_in_node[transcript] if exon != exon_to_skip]
                    new_dict_transcripts_in_node_sequences[name_transcript] = {}
                    for exon in new_dict_transcripts_in_node[transcript]:
                        if exon != exon_to_skip:
                            new_dict_transcripts_in_node_sequences[name_transcript][exon] = new_dict_transcripts_in_node_sequences[transcript][exon]
                    new_dict_transcripts_in_node_names.append(name_transcript)
                    #update_phylogenies
                    parent_transcripts_forest, lca_events = update_transcripts_forest(parent_transcripts_forest, transcript, name_transcript, lca_events, 'es')
            
                              
            #break
            
            #### ME
            
            #number of transcripts that have undergone a ME event
            NUMBER_TRANSCRIPTS_ME = int(TC_ME * K_TC * CODON_SUBST_RATE * len(new_dict_transcripts_in_node_names))
            #print(NUMBER_TRANSCRIPTS_ES, K_TC, TC_ES)
            #NUMBER_TRANSCRIPTS_ME = 2 #example
            
            # choose the transcripts affected
            transcripts_choosed = choice_distinct(new_dict_transcripts_in_node_names, NUMBER_TRANSCRIPTS_ME)
            
            # Apply the ME event on these transcripts
            for transcript in transcripts_choosed:
                number_exons = len(new_dict_transcripts_in_node[transcript])
                if number_exons < 4:
                    # cannot apply ME event on this transcript
                    pass
                else:
                    exons = new_dict_transcripts_in_node[transcript]
                    choosed_adj_sub = np.random.randint(1, len(exons)-2)
                    choosed_exons = [exons[choosed_adj_sub], exons[choosed_adj_sub+1]]
                    exons_compositions = []
                    isOk = True
                    for exon_me in choosed_exons:
                        tmp_list_transcript_desc = []
                        for exon in exons:
                            if exon != exon_me:
                                tmp_list_transcript_desc.append(exon)
                        for transcript_i in new_dict_transcripts_in_node.keys():
                            compo = new_dict_transcripts_in_node[transcript_i]
                            if compo == tmp_list_transcript_desc:
                                isOk = False
                                break
                        if not isOk:
                            break
                        else:
                            exons_compositions.append(tmp_list_transcript_desc)
                    if isOk:
                        for transcript_composition in exons_compositions:
                            # which transcript
                            name_transcript = 'transcript_{}-{}'.format(preorder_number, number_new_transcript)
                            number_new_transcript += 1
                            
                            # m.a.j
                            new_dict_transcripts_in_node[name_transcript] = transcript_composition
                            new_dict_transcripts_in_node_sequences[name_transcript] = {}
                            for exon in transcript_composition:
                                new_dict_transcripts_in_node_sequences[name_transcript][exon] = new_dict_transcripts_in_node_sequences[transcript][exon]
                            new_dict_transcripts_in_node_names.append(name_transcript)
                            #update_phylogenies
                            parent_transcripts_forest, lca_events = update_transcripts_forest(parent_transcripts_forest, transcript, name_transcript, lca_events, 'me')
            
            
                                 
            #### 5SS
            # number of transcripts affected
            NUMBER_TRANSCRIPTS_5SS = int(K_TC * CODON_SUBST_RATE * TC_A5 * len(new_dict_transcripts_in_node_names))
            #NUMBER_TRANSCRIPTS_5SS = 2 #example
            #choose transcripts affected
            transcripts_choosed = choice_distinct(new_dict_transcripts_in_node_names, NUMBER_TRANSCRIPTS_5SS)
            # Apply 5SS event on these transcripts
            for transcript in transcripts_choosed:
                list_exons = new_dict_transcripts_in_node[transcript]
                number_exons = len(list_exons)
                if number_exons == 1:
                    # cannot apply on this transcript
                    pass
                else:
                    exon_five_prime = np.random.choice(list_exons[1:])
                    sequence_exon = new_dict_transcripts_in_node_sequences[transcript][exon_five_prime]
                    if len(sequence_exon) > 54:
                        range_to_choose = 9 # can be modifying
                        sequence_5_prime = sequence_exon[:range_to_choose+1]
                        if '*' in sequence_5_prime:
                            continue
                        else:
                            sequence_new_exon = ''.join(['-' for i in range(range_to_choose+1)]) + 'AG' +sequence_exon[range_to_choose+3:]
                            #pos_x = len(''.join(['-' for i in range(range_to_choose+1)])) + 1
                            #pos_y = pos_x + 1
                            #seq_tmp = list(sequence_exon)
                            #seq_tmp[pos_x] = 'A'
                            #seq_tmp[pos_y] = 'G'
                            #print(len(sequence_new_exon), len(sequence_exon))
                            # which transcript
                            name_transcript = 'transcript_{}-{}'.format(preorder_number, number_new_transcript)
                            number_new_transcript += 1
                            # m.a.j
                            new_dict_transcripts_in_node[name_transcript] = new_dict_transcripts_in_node[transcript]
                            new_dict_transcripts_in_node_sequences[name_transcript] = {}
                            for exon in list_exons:
                                if exon == exon_five_prime:
                                    new_dict_transcripts_in_node_sequences[name_transcript][exon] = sequence_new_exon
                                    #current_exons_dict[exon] = ''.join(seq_tmp)
                                    current_exons_dict[exon] = sequence_new_exon
                                else:
                                    new_dict_transcripts_in_node_sequences[name_transcript][exon] = new_dict_transcripts_in_node_sequences[transcript][exon]
                            new_dict_transcripts_in_node_names.append(name_transcript)
                            #update_phylogenies
                            parent_transcripts_forest, lca_events = update_transcripts_forest(parent_transcripts_forest, transcript, name_transcript, lca_events, '5ss')
                                      
            #break           
                    
            #### 3SS
            # number of transcripts affected
            NUMBER_TRANSCRIPTS_3SS = int(K_TC * CODON_SUBST_RATE * TC_A3 * len(new_dict_transcripts_in_node_names))
            #NUMBER_TRANSCRIPTS_3SS = 2 #example
            #choose transcripts affected
            transcripts_choosed = choice_distinct(new_dict_transcripts_in_node_names, NUMBER_TRANSCRIPTS_3SS)
            
                
            # Apply 3SS event on these transcripts
            for transcript in transcripts_choosed:
                list_exons = new_dict_transcripts_in_node[transcript]
                number_exons = len(list_exons)
                if number_exons == 1:
                    # cannot apply on this transcript
                    pass
                else:
                    exon_three_prime = np.random.choice(list_exons[:-1])
                    sequence_exon = new_dict_transcripts_in_node_sequences[transcript][exon_three_prime]
                    if len(sequence_exon) > 54:
                        range_to_choose = 9 # can be modifying
                        sequence_3_prime = sequence_exon[-range_to_choose:]
                        if '*' in sequence_3_prime:
                            continue
                        else:
                            
                            sequence_new_exon = sequence_exon[:-range_to_choose-2] + 'GT' + ''.join(['-' for i in range(range_to_choose)]) 
                            #pos_x = len(sequence_exon[:-range_to_choose-2]) + 1
                            #pos_y = pos_x + 1
                            #seq_tmp = list(sequence_exon)
                            #seq_tmp[pos_x] = 'G'
                            #seq_tmp[pos_y] = 'T'
                            
                            # which transcript
                            name_transcript = 'transcript_{}-{}'.format(preorder_number, number_new_transcript)
                            number_new_transcript += 1
                            # m.a.j
                            new_dict_transcripts_in_node[name_transcript] = new_dict_transcripts_in_node[transcript]
                            new_dict_transcripts_in_node_sequences[name_transcript] = {}
                            for exon in list_exons:
                                if exon == exon_three_prime:
                                    new_dict_transcripts_in_node_sequences[name_transcript][exon] = sequence_new_exon
                                    #current_exons_dict[exon] = ''.join(seq_tmp)
                                    current_exons_dict[exon] = sequence_new_exon
                                else:
                                    new_dict_transcripts_in_node_sequences[name_transcript][exon] = new_dict_transcripts_in_node_sequences[transcript][exon]
                            new_dict_transcripts_in_node_names.append(name_transcript)
                            #update_phylogenies
                            parent_transcripts_forest, lca_events = update_transcripts_forest(parent_transcripts_forest, transcript, name_transcript, lca_events, '3ss')
                        
            
            #### IR
            
            # number of transcripts affected
            NUMBER_TRANSCRIPTS_IR = int(CODON_SUBST_RATE * K_TC * TC_IR * len(new_dict_transcripts_in_node_names))
            #NUMBER_TRANSCRIPTS_IR = 2 #example
            
            #choose transcripts_affected
            transcripts_choosed = choice_distinct(new_dict_transcripts_in_node_names, NUMBER_TRANSCRIPTS_IR)
            
            # Apply IR on these transcripts
            exons_news = []
            for transcript in transcripts_choosed:
                list_exons = new_dict_transcripts_in_node[transcript]
                number_exons = len(list_exons)
                
                if number_exons <= 2:
                    # cannot apply
                    pass
                else:
                    #intron retained
                    choosed_location = np.random.randint(1, number_exons-1)
                    pos_insertion_at_tr_level = choosed_location
                    exon_of_interest = list_exons[pos_insertion_at_tr_level]
                    length = int(random_number_generator_by_normal_law(MEAN_INTRON_LEN, SD_INTRON_LEN)/3.0)
                    new_exon_sequence = intron_generator(length, SPLICE_SITES, all_nt_intron_dist)
                    number_gain += 1
                    #pos_insertion = random.randint(1, len(current_gene_structure)-1)
                    #gene
                    pos_insertion = current_gene_structure.index(exon_of_interest)
                    start_exon = current_gene_structure[pos_insertion-1]
                    stop_exon = current_gene_structure[pos_insertion]
                    current_positions_of_all_exons_in_the_simulation, exons_new = placement_exons_in_evolution('new', 'null', "exon_#"+str(number_gain)+"-intron", current_positions_of_all_exons_in_the_simulation, current_positions_of_all_exons_in_the_simulation.index(start_exon), current_positions_of_all_exons_in_the_simulation.index(stop_exon), exons_new)
                    if "exon_#"+str(number_gain)+"-intron" in exons_new:
                        current_gene_structure.insert(pos_insertion, "exon_#"+str(number_gain)+"-intron")
                        current_exons_dict["exon_#"+str(number_gain)+"-intron"] = new_exon_sequence
                    
                    #transcript
                    # m.a.j
                    # which transcript
                    name_transcript = 'transcript_{}-{}'.format(preorder_number, number_new_transcript)
                    number_new_transcript += 1
                    update_composition= []
                    for exon in current_gene_structure:
                        if exon in new_dict_transcripts_in_node[transcript] or exon in exons_new:
                            update_composition.append(exon)
                    new_dict_transcripts_in_node[name_transcript] = update_composition
                    new_dict_transcripts_in_node_sequences[name_transcript] = {}
                    for exon in update_composition:
                        if exon in exons_new:
                            new_dict_transcripts_in_node_sequences[name_transcript][exon] = current_exons_dict[exon]
                        else:
                            new_dict_transcripts_in_node_sequences[name_transcript][exon] = new_dict_transcripts_in_node_sequences[transcript][exon]
                    new_dict_transcripts_in_node_names.append(name_transcript)
                    #update_phylogenies
                    parent_transcripts_forest, lca_events = update_transcripts_forest(parent_transcripts_forest, transcript, name_transcript, lca_events, 'ir')
                    
            # introns generation
            new_node_structure = []
            current_introns_dict = {}
            for i, exon in enumerate(current_gene_structure):
                new_node_structure.append(exon)
                if i != len(current_gene_structure)-1:
                    intron_name = 'intron_{}'.format(i)
                    new_node_structure.append(intron_name)
                    length = random_number_generator_by_normal_law(MEAN_INTRON_LEN, SD_INTRON_LEN)
                    intron_sequence = intron_generator(length, SPLICE_SITES, all_nt_intron_dist)
                    current_introns_dict[intron_name] = intron_sequence
            
            #m.a.j the forest in the parent
            #node.up.transcripts_forest = copy.deepcopy(parent_transcripts_forest)
   
            # m.a.j the gene structure
            node.add_features(
                gene_name=copy.deepcopy(gene_name),
                gene_structure=copy.deepcopy(new_node_structure),
                exons_dict=copy.deepcopy(current_exons_dict),
                introns_dict=copy.deepcopy(current_introns_dict),
                transcripts_dict=copy.deepcopy(new_dict_transcripts_in_node),
                transcripts_sequences_dict=copy.deepcopy(new_dict_transcripts_in_node_sequences),
                all_exons_in_the_simulation=copy.deepcopy(current_positions_of_all_exons_in_the_simulation)
                #transcripts_forest=copy.deepcopy(parent_transcripts_forest)
            )
    return tree, parent_transcripts_forest, lca_events

def draw_transcripts_phylogenies(ttree, all_transcripts, lca_events, tree, tree_name, source, iteration_name, df_transcripts_to_alignments):
    for n in ttree.traverse():
        n.add_features(all_transcripts_remove=all_transcripts,
                       lca_events_remove=lca_events,
                       df_transcripts_to_alignments_remove = df_transcripts_to_alignments)
    for n in ttree.traverse():
        n.dist = 1
        if n.is_leaf():
            nstyle = NodeStyle()
            nstyle['size'] = 15
            nstyle['shape'] = 'sphere'
            transcript = n.name
            #print(transcript)
            if transcript in all_transcripts:
                nstyle['fgcolor'] = 'gold'
            else:
                nstyle['fgcolor'] = 'lightgrey'
            n.set_style(nstyle)
        else:
            nstyle = NodeStyle()
            nstyle['size'] = 15
            nstyle['shape'] = 'sphere'
            node_name = n.name
            nstyle["hz_line_color"] = "black"
            event = lca_events[node_name]
            if event == 'ir': 
                nstyle['fgcolor'] = 'red'
            elif event == 'me':
                nstyle['fgcolor'] = 'orange' 
            elif event == '5ss':
                nstyle['fgcolor'] = 'violet'
            elif event == '3ss':
                nstyle['fgcolor'] = 'mediumblue' 
            elif event == 'es':
                nstyle['fgcolor'] = 'limegreen' 
            else:
                nstyle['fgcolor'] = 'black'
            n.set_style(nstyle)
    ts = TreeStyle()
    ts.scale = 120
    ts.branch_vertical_margin = 50
    ts.show_leaf_name = True
    ts.show_scale = False

    os_ttree_png = os.path.join(source, iteration_name, 'phylogenies', tree_name+'_w_msa_.png')
    os_ttree_svg = os.path.join(source, iteration_name, 'phylogenies', tree_name+'_w_msa_.svg')
    os_ttree_nwk = os.path.join(source, iteration_name, 'phylogenies', tree_name+'.nwk')
    ttree.render(os_ttree_png, w=183, units="mm", tree_style=ts)
    ttree.render(os_ttree_svg, w=300, units="mm", tree_style=ts)

    tss = TreeStyle()
    #tss.legend.add_face(TextFace("true alignment", fsize=5, ftype='courier'), column=0)
    tss.layout_fn= ColorCodedNode
    tss.scale = 120
    tss.branch_vertical_margin = 50
    tss.show_leaf_name = True
    tss.show_scale = False
    os_ttree_png = os.path.join(source, iteration_name, 'phylogenies', tree_name+'_msa_.png')
    os_ttree_svg = os.path.join(source, iteration_name, 'phylogenies', tree_name+'_msa_.svg')
    ttree.render(os_ttree_png, units="px", tree_style=tss)
    ttree.render(os_ttree_svg, units="mm", tree_style=tss)


    
    
    ttree.write(format=9, outfile=os_ttree_nwk)

def ColorCodedNode(n):
    all_transcripts = n.all_transcripts_remove
    df_transcripts_to_alignments = n.df_transcripts_to_alignments_remove
    lca_events = n.lca_events_remove
    n.dist = 1
    if n.is_leaf():
        nstyle = NodeStyle()
        nstyle['size'] = 15
        nstyle['shape'] = 'sphere'
        transcript = n.name
        #print(transcript)
        if transcript in all_transcripts:
            nstyle['fgcolor'] = 'gold'
        else:
            nstyle['fgcolor'] = 'lightgrey'
        

        if transcript in df_transcripts_to_alignments.keys():
            alignment_seq = df_transcripts_to_alignments[transcript]
            #faces.add_face(alignment_seq, 0, "aligned")
            seqFace = faces.SeqMotifFace(alignment_seq, gapcolor='black', fgcolor='red')
            n.add_face(seqFace, 0, "aligned")
        else:
            pass
        n.set_style(nstyle)
    else:
        nstyle = NodeStyle()
        nstyle['size'] = 15
        nstyle['shape'] = 'sphere'
        node_name = n.name
        nstyle["hz_line_color"] = "black"
        event = lca_events[node_name]
        if event == 'ir': 
            nstyle['fgcolor'] = 'red'
        elif event == 'me':
            nstyle['fgcolor'] = 'orange' 
        elif event == '5ss':
            nstyle['fgcolor'] = 'violet'
        elif event == '3ss':
            nstyle['fgcolor'] = 'mediumblue' 
        elif event == 'es':
            nstyle['fgcolor'] = 'limegreen' 
        else:
            nstyle['fgcolor'] = 'black'
        n.set_style(nstyle)
    '''
    os_ttree_png = os.path.join(source, iteration_name, 'phylogenies', tree_name+'.png')
    os_ttree_svg = os.path.join(source, iteration_name, 'phylogenies', tree_name+'.svg')
    os_ttree_nwk = os.path.join(source, iteration_name, 'phylogenies', tree_name+'.nwk')
    ttree.render(os_ttree_png, w=183, units="mm", tree_style=ts)
    ttree.render(os_ttree_svg, w=300, units="mm", tree_style=ts)
    ttree.write(format=9, outfile=os_ttree_nwk)
    '''     
    return True

def get_msa_and_clusters(tree, forest, lca_events, source, iteration_name):
    #print(tree)
    nodes_of_interest = [node for node in tree.traverse('preorder') if node.is_leaf()]
    node_of_interest = nodes_of_interest[-1]
    #columns
    all_exons_in_the_simulation = node_of_interest.all_exons_in_the_simulation
    #lines
    all_transcripts = []
    for node in nodes_of_interest:
        trs = list(node.transcripts_dict.keys())
        all_transcripts.extend(trs)
    #dataframe creation
    df = pd.DataFrame(columns=['id_transcript'] + all_exons_in_the_simulation, index=all_transcripts)
    
    #fill the df
    #print(df)
    all_transcripts_dict = []
    for exon in all_exons_in_the_simulation:
        for node in nodes_of_interest:
            transcripts_dict = node.transcripts_dict
            
            transcripts = list(transcripts_dict.keys())
            for transcript in transcripts:
                df.loc[transcript]['id_transcript'] = transcript
                exons_in_transcript = transcripts_dict[transcript]
                all_transcripts_dict.append({transcript:exons_in_transcript})
                if exon in exons_in_transcript:
                    exon_sequence = node.transcripts_sequences_dict[transcript][exon]
                    df.loc[transcript][exon] = exon_sequence
                else:
                    df.loc[transcript][exon] = 'nil'
     
    #print(df)
            
    lengths_exons = []
    exons_at_leaves = []
    for exon in all_exons_in_the_simulation:
        exons_in_column = [_ for _ in df[exon].values]
        lengths_exons_in_column = [len(_) for _ in exons_in_column if _ != 'nil']
        if len(lengths_exons_in_column) != 0:
            exons_at_leaves.append(exon)
        
            if all(lengths_exons_in_column):
                length = lengths_exons_in_column[0]
                lengths_exons.append(length)
                missed_seq = ''.join(['-' for _ in range(length)])
                df[exon] = df[exon].replace(['nil'], missed_seq)
                #print(length)
            else:
                raise('Errors')
    # retrieve the msa of transcripts
    df_transcripts_to_alignments = {}
    file = open('./{}/{}/{}/msa_transcripts.alg'.format(source, iteration_name, 'multiple_alignments'), 'w')
    for iter, row in df.iterrows():
        sequence = []
        for exon in exons_at_leaves:
            tmp_seq = row[exon]
            sequence.append(tmp_seq)
        sequence_alg = ''.join(sequence)
        df_transcripts_to_alignments[row.id_transcript] = sequence_alg
        file.write('>{}\n{}\n'.format(row.id_transcript, sequence_alg))
    file.close()
    
    groups_of_orthologs = []
    already_checked_transcripts = []
    for transcript_enum in range(len(all_transcripts_dict)):
        transcript_info = all_transcripts_dict[transcript_enum]
        transcript_id_keys = list(transcript_info.keys())
        #print(transcript_id_keys)
        transcript_id = transcript_id_keys[0]
        if transcript_id not in already_checked_transcripts:
            exons_id = transcript_info[transcript_id]
            same_structure_transcripts = set()
            same_structure_transcripts.add(transcript_id)
            already_checked_transcripts.append(transcript_id)
            for transcript_enum_ref in range(len(all_transcripts_dict)):
                transcript_info_ref = all_transcripts_dict[transcript_enum_ref]
                transcript_ref = list(transcript_info_ref.keys())[0]
                exons_ref = transcript_info_ref[transcript_ref]
                if transcript_ref != transcript_id and len(exons_ref) == len(exons_id):
                    is_the_same_structure = True
                    for exon_id in exons_id:
                        if exon_id not in exons_ref:
                            is_the_same_structure = False
                            break
                    if is_the_same_structure:
                        same_structure_transcripts.add(transcript_ref)
                        already_checked_transcripts.append(transcript_ref)
            groups_of_orthologs.append(same_structure_transcripts)
                                        
    # save phylogenies
    dict_forest = {}
    for enum_forest, ttree in enumerate(forest):
        leaves = [str(_.name) for _ in ttree.get_leaves()]
        if any([True if transcript_id in leaves else False for transcript_id in all_transcripts]):
            enum_forest += 1
            tree_name = 'phylo_{}'.format(enum_forest)
            for leaf in leaves:
                dict_forest[leaf] = tree_name
            if len(leaves) >= 2:
                draw_transcripts_phylogenies(ttree, all_transcripts, lca_events, tree, tree_name, source, iteration_name, df_transcripts_to_alignments)
        
    #save the clusters into file
    file_save_clusters = open('./{}/{}/{}/ortholog_groups.clusters'.format(source, iteration_name, 'clusters'), 'w')
    cluster_enum = 0
    is_there_not_cluster = False
    for group in groups_of_orthologs:
        #group = set()
        #for element in group_raw:
        #    if element in all_transcripts:
        #        group.add(element)
        if len(group) >= 2:
            cluster_enum += 1
            is_there_not_cluster = True
            transcript_zero = list(group)[0]
            phylo_name = dict_forest[transcript_zero]
            file_save_clusters.write('>cluster_{}:(see phylogenies\{})\n'.format(cluster_enum, phylo_name))
            for transcript in group:
                file_save_clusters.write('{}\n'.format(transcript))
            file_save_clusters.write('\n')
    if not is_there_not_cluster:
        file_save_clusters.write('Sorry! No ortholog groups for this simulation.')
    file_save_clusters.close()
    
    #print(df)
    return df

def write_data(df_transcripts, tree, source, iteration_name):
    #print(df_transcripts)
    file_pairwise = open("{}/{}/{}/pairwise_alignments.fasta".format(source, iteration_name, 'pairwise_alignments'), 'w')
    file_cdna = open("{}/{}/{}/cdna.fasta".format(source, iteration_name, 'cdna'), 'w')
    file_genes = open("{}/{}/{}/genes.fasta".format(source, iteration_name, 'genes'), 'w')
    file_transcripts = open("{}/{}/{}/transcripts.fasta".format(source, iteration_name, 'transcripts'), 'w')
    file_transcripts_to_gene = open("{}/{}/{}/mappings.txt".format(source, iteration_name, 'transcripts_to_gene'), 'w')
    file_exons_positions = open("{}/{}/{}/exons_positions.txt".format(source, iteration_name, 'exons_positions'), 'w')
    nodes = [node for node in tree.traverse('preorder') if node.is_leaf()]
    for node in nodes:
        gene_structure = node.gene_structure
        exons_dict = node.exons_dict
        gene = node.gene_name
        exons = [exon for exon in gene_structure if exon.startswith('exon')]
        transcripts = list(node.transcripts_dict.keys())
        introns_dict = node.introns_dict
        #columns = list(df_transcripts.columns)
        sequence_exons = ''
        sequence_gene = ''
        for struc in gene_structure:
            if struc.startswith('exon'):
                sequence_gene += exons_dict[struc]
            else:
                sequence_gene += introns_dict[struc]
        for exon in exons:
            sequence_exons += exons_dict[exon]
        file_cdna.write('>{}\n{}\n'.format(gene, sequence_exons.replace('-','').replace('*','')))
        file_genes.write('>{}\n{}\n'.format(gene, sequence_gene.replace('-','').replace('*','')))

        
        sequence_algs = {}
        for transcript in transcripts:
            exons_in_transcript = node.transcripts_dict[transcript]
            file_exons_positions.write('>{}:{}\n'.format(transcript, gene))
            start = 0
            start_g = 0
            seq_transcript = ''
            for struc in gene_structure:
                if struc.startswith('exon'):
                    exon = struc
                    #for exon in exons:
                
                    if exon not in exons_in_transcript:
                        seq = exons_dict[exon]
                        seq_missed = ''.join(['-' for _ in range(len(seq))])
                        seq_transcript += seq_missed
                        seq_exon_missed_without_gaps = seq_missed.replace('-','').replace('*','')
                        start_gene = start_g
                        end_gene = start_gene + len(seq_exon_missed_without_gaps)
                        start_g = end_gene
                    else:
                        seq = node.transcripts_sequences_dict[transcript][exon]
                        seq_transcript += seq
                        seq_exon_without_gaps = seq.replace('-','').replace('*','')
                        start_tr = start
                        end_tr = start_tr + len(seq_exon_without_gaps)
                        start = end_tr
                        # gene
                        start_gene = start_g
                        end_gene = start_gene + len(seq_exon_without_gaps)
                        start_g = end_gene
                        file_exons_positions.write('{}\t{}\t{}\t{}\n'.format(start_tr, end_tr, start_gene, end_gene))
                else:
                    intron_seq = introns_dict[struc]
                    start_gene = start_g
                    end_gene = start_gene + len(intron_seq)
                    start_g = end_gene
            file_exons_positions.write('\n')        
            sequence_algs[transcript] = seq_transcript
            file_pairwise.write('>{}\n{}\n>{}\n{}\n\n'.format(gene, sequence_exons, transcript, seq_transcript))
            file_transcripts.write('>{}\n{}\n'.format(transcript, seq_transcript.replace('*','').replace('-','')))
            file_transcripts_to_gene.write('>{}:{}\n'.format(transcript, gene))
    
    print(colorsSimulation.OKCYAN + '\n+++++\tExons Alignment Preview\t+++++\n' + colorsSimulation.ENDC)
    
    columns = [_ for _ in gene_structure if _.startswith('exon') and _ in df_transcripts.columns]
    print(columns)
    for iter, row in df_transcripts.iterrows():
        display_sequence = []
        id_tr = row.id_transcript
        for column in columns:
            data = row[column]
            if len(set(list(data))) == 1:
                display_sequence.append('---')
            else:
                display_sequence.append('🧬')
        cur_display_sequence = '\t'.join(display_sequence)
        display = f'{id_tr}\t\t{cur_display_sequence}'
        print(display)
            
    file_pairwise.close()
    file_cdna.close()
    file_transcripts.close()
    file_transcripts_to_gene.close()
    file_exons_positions.close()
    file_genes.close()
    return True

def create_output_folders(folder, iteration):
    folder_iteration = '{}/{}'.format(folder, iteration)
    folder_clusters = '{}/clusters'.format(folder_iteration)
    folder_pairwise_alignment = '{}/pairwise_alignments'.format(folder_iteration)
    folder_genes = '{}/genes'.format(folder_iteration)
    folder_cdna = '{}/cdna'.format(folder_iteration)
    folder_exons_positions = '{}/exons_positions'.format(folder_iteration)
    folder_transcripts_to_gene = '{}/transcripts_to_gene'.format(folder_iteration)
    folder_transcripts = '{}/transcripts'.format(folder_iteration)
    folder_relations = '{}/phylogenies'.format(folder_iteration)
    folder_msa = '{}/multiple_alignments'.format(folder_iteration)

    folders = [folder_clusters, folder_msa, folder_relations, folder_pairwise_alignment, folder_genes, folder_exons_positions, folder_transcripts_to_gene, folder_transcripts, folder_cdna]
    try:
        os.system('mkdir {}'.format(folder_iteration))
    except OSError:
            print("\t ❗ The output folder: {} already exists. The files inside will be modified".format(folder_iteration))
    for folder_path in folders:
        try:
            os.mkdir('{}'.format(folder_path))
            print("\t ✔️ The output folder: {} successfully created.".format(folder_path))
        except OSError:
            print("\t ❗ The output folder: {} already exists. The files inside will be modified".format(folder_path))
    return True

def simspliceevol(SRC, ITERATION_NAME, TREE_INPUT, K_NB_EXONS, K_INDEL, C_I, C_D, EIC_ED, EIC_EG, EIC_EL, K_EIC, K_TC, TC_RS, TC_A3, TC_A5, TC_ME, TC_ES, TC_IR, TC_TL):
    """The main function of the program.

    Args:
        inputs (_LIST_): [SRC, ITERATION_NAME, TREE_INPUT, ONE_NT_DIST, TWO_NT_DIST, THREE_NT_DIST, NT_INTRON_DIST, MEAN_EXON_CDS, MEAN_EXON_LEN, MEAN_INTRON_LEN, SD_EXON_CDS, SD_EXON_LEN, SD_INTRON_LEN, K_NB_EXONS, K_INDEL, C_I, C_D, EIC_ED, EIC_EG, EIC_EL, K_EIC, STOP_CODONS, SPLICE_SITES, CODONS_MATRIX, CODONS_LIST, K_TC, TC_RS, TC_A3, TC_A5, TC_ME, TC_ES, TC_IR, TC_TL]

    Returns:
        _ete3 Tree python object_: Tree python object from the ete3 library
    """
    
    #==========================================
    # CONSTANTES
    #==========================================
    CODONS_LIST = ['TTT', 'TTC', 'TTA', 'TTG', 'TCT', 'TCC', 'TCA', 'TCG', 'TAT', 'TAC', 'TGT', 'TGC', 'TGG', 'CTT',
				  'CTC', 'CTA', 'CTG', 'CCT', 'CCC', 'CCA', 'CCG', 'CAT', 'CAC', 'CAA', 'CAG', 'CGT', 'CGC', 'CGA',
				  'CGG', 'ATT', 'ATC', 'ATA', 'ATG', 'ACT', 'ACC', 'ACA', 'ACG', 'AAT', 'AAC', 'AAA', 'AAG', 'AGT',
				  'AGC', 'AGA', 'AGG', 'GTT', 'GTC', 'GTA', 'GTG', 'GCT', 'GCC', 'GCA', 'GCG', 'GAT', 'GAC', 'GAA',
				  'GAG', 'GGT', 'GGC', 'GGA', 'GGG']

    CODONS_MATRIX = [[650., 160.,  20.,  10.,  10.,   0.,   0.,   0.,  40.,   0.,  10.,
          0.,  10.,  20.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,  10.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,  10.,  10.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,  10.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.],
 [160., 670.,   0.,  10.,   0.,  10.,   0.,   0.,  10.,  40.,   0.,
         10.,  10.,   0.,  20.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,
         10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,  10.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,  10.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.],
 [ 20.,   0., 200., 190.,   0.,   0.,  10.,   0.,  10.,   0.,  10.,
          0.,   0., 120.,  10., 270.,  10.,   0.,   0.,   0.,   0.,   0.,
          0.,  10.,   0.,   0.,   0.,   0.,   0.,  20.,   0.,  50.,  10.,
          0.,   0.,  10.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,   0.,
          0.,  10.,   0.,  20.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.],
 [ 10.,  10., 190., 200.,   0.,   0.,  10.,  10.,   0.,   0.,   0.,
          0.,  10., 110.,  60., 170., 120.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,  10.,  30.,
          0.,   0.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         10.,  10.,   0.,  10.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.],
 [  10.,    0.,    0.,    0., -120.,  240.,  360.,  120.,   10.,
           0.,   20.,    0.,    0.,    0.,    0.,    0.,    0.,   30.,
           0.,    0.,    0.,   10.,    0.,   10.,    0.,   10.,    0.,
           0.,    0.,    0.,    0.,    0.,    0.,   50.,    0.,   20.,
           0.,   20.,    0.,    0.,    0.,  120.,    0.,    0.,    0.,
          10.,    0.,    0.,    0.,   40.,    0.,   20.,    0.,   10.,
           0.,    0.,    0.,   10.,    0.,    0.,    0.,    0.],
 [  0.,  10.,   0.,   0., 240., -80., 160., 350.,   0.,  10.,   0.,
         10.,   0.,   0.,   0.,   0.,   0.,   0.,  20.,   0.,   0.,   0.,
         10.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,  40.,   0.,  20.,   0.,  10.,   0.,   0.,   0., 110.,   0.,
          0.,   0.,   0.,   0.,   0.,  10.,  30.,  10.,  10.,   0.,  10.,
          0.,   0.,   0.,  10.,   0.,   0.,   0.],
 [  0.,   0.,  10.,  10., 360., 160., -80., 170.,   0.,   0.,  10.,
          0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,  30.,   0.,   0.,
          0.,  20.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,
         20.,   0.,  60.,  10.,  10.,   0.,  10.,   0.,  60.,  10.,  10.,
          0.,   0.,   0.,  10.,   0.,  20.,   0.,  40.,   0.,   0.,   0.,
         10.,   0.,   0.,   0.,  20.,   0.,   0.],
 [  0.,   0.,   0.,  10., 120., 350., 170., -50.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  20.,   0.,
          0.,   0.,  20.,   0.,   0.,   0.,  20.,   0.,   0.,   0.,  10.,
          0.,  20.,   0.,  50.,   0.,  10.,   0.,  10.,  10., 100.,   0.,
         10.,   0.,   0.,   0.,   0.,   0.,  20.,  10.,  50.,   0.,  10.,
          0.,  10.,   0.,   0.,   0.,  20.,   0.],
 [ 40.,  10.,  10.,   0.,  10.,   0.,   0.,   0., 580., 240.,  10.,
          0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  30.,
          0.,   0.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,
         10.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,  10.,   0.,   0.,
          0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.],
 [  0.,  40.,   0.,   0.,   0.,  10.,   0.,   0., 240., 630.,   0.,
         10.,  10.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         30.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.],
 [ 10.,   0.,  10.,   0.,  20.,   0.,  10.,   0.,  10.,   0., 250.,
        570.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  10.,
          0.,   0.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  20.,   0.,   0.,
          0.,  20.,   0.,  10.,   0.,  20.,   0.,  10.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.],
 [  0.,  10.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,  10., 570.,
        310.,   0.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,
          0.,   0.,  10.,   0.,  10.,   0.,  20.,   0.,  10.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.],
 [ 10.,  10.,   0.,  10.,   0.,   0.,   0.,   0.,  10.,  10.,   0.,
          0., 930.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.],
 [ 20.,   0., 120., 110.,   0.,   0.,   0.,   0.,   0.,   0.,  10.,
          0.,   0., 180., 160., 210.,  50.,  10.,   0.,   0.,   0.,  10.,
          0.,  10.,   0.,  10.,   0.,   0.,   0.,  30.,   0.,  10.,  10.,
         10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,  20.,   0.,  10.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.],
 [  0.,  20.,  10.,  60.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,
         10.,   0., 160., 350.,  90., 150.,   0.,  10.,   0.,   0.,   0.,
         10.,   0.,   0.,   0.,  10.,   0.,  10.,   0.,  40.,   0.,  10.,
          0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,  30.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.],
 [ 10.,   0., 270., 170.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,   0., 210.,  90.,  20.,  80.,   0.,   0.,  10.,   0.,   0.,
          0.,  10.,   0.,   0.,   0.,  10.,   0.,  10.,   0.,  40.,  20.,
          0.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,  10.,   0.,  20.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.],
 [  0.,  10.,  10., 120.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,  10.,  50., 150.,  80., 440.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,  10.,   0.,  10.,   0.,  10.,   0.,  20.,   0.,  20.,
          0.,   0.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,  10.,   0.,  20.,   0.,   0.,   0.,  20.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.],
 [  0.,   0.,   0.,   0.,  30.,   0.,  10.,   0.,   0.,   0.,   0.,
          0.,   0.,  10.,   0.,   0.,   0., 130., 230., 410., 100.,  10.,
          0.,   0.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         10.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,  10.,   0.,   0.,
          0.,  10.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,  10.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.],
 [  0.,   0.,   0.,   0.,   0.,  20.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,  10.,   0.,   0., 230., 150., 150., 330.,   0.,
         10.,   0.,  10.,   0.,  10.,   0.,  10.,   0.,   0.,   0.,   0.,
          0.,  10.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,  10.,   0.,
          0.,   0.,  10.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,  10.,
          0.,   0.,   0.,  10.,   0.,   0.,   0.],
 [  0.,   0.,   0.,   0.,   0.,   0.,  30.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,  10.,   0., 410., 150., 170., 150.,   0.,
          0.,  10.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,  10.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,  10.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,
         10.,   0.,   0.,   0.,  10.,   0.,   0.],
 [  0.,   0.,   0.,   0.,   0.,   0.,   0.,  20.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0., 100., 330., 150., 270.,   0.,
          0.,   0.,  10.,   0.,  10.,   0.,  10.,   0.,   0.,   0.,   0.,
          0.,  10.,   0.,  10.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,  10.,   0.,  10.,   0.,  20.,   0.,  10.,
          0.,  10.,   0.,   0.,   0.,  10.,   0.],
 [ 10.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,  30.,   0.,  10.,
          0.,   0.,  10.,   0.,   0.,   0.,  10.,   0.,   0.,   0., 400.,
        390.,  20.,  10.,  20.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,
         10.,   0.,   0.,   0.,  30.,   0.,  10.,   0.,  10.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.],
 [  0.,  10.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,  30.,   0.,
          0.,   0.,   0.,  10.,   0.,   0.,   0.,  10.,   0.,   0., 390.,
        360.,  10.,  30.,   0.,  20.,  10.,  20.,   0.,   0.,   0.,   0.,
          0.,  10.,   0.,   0.,   0.,  20.,   0.,  10.,   0.,  10.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,  10.,
          0.,  10.,   0.,  10.,   0.,   0.,   0.],
 [  0.,   0.,  10.,   0.,  10.,   0.,  20.,   0.,   0.,   0.,   0.,
          0.,   0.,  10.,   0.,  10.,   0.,   0.,   0.,  10.,   0.,  20.,
         10., 500., 180.,  10.,   0.,  50.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,  20.,   0.,  20.,   0.,  30.,   0.,  10.,   0.,  10.,
          0.,   0.,   0.,  10.,   0.,   0.,   0.,  10.,   0.,  10.,   0.,
         30.,   0.,   0.,   0.,  10.,   0.,   0.],
 [  0.,   0.,   0.,   0.,   0.,   0.,   0.,  20.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,  10.,   0.,  10.,   0.,  10.,  10.,
         30., 180., 470.,   0.,  20.,  10.,  50.,   0.,   0.,   0.,  10.,
          0.,  10.,   0.,  20.,   0.,  10.,   0.,  20.,   0.,  10.,   0.,
         10.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,  20.,   0.,  10.,
         10.,  30.,   0.,   0.,   0.,  10.,   0.],
 [   0.,    0.,    0.,    0.,   10.,    0.,    0.,    0.,   10.,
           0.,   10.,    0.,    0.,   10.,    0.,    0.,    0.,   10.,
           0.,    0.,    0.,   20.,    0.,   10.,    0., -170.,  270.,
         400.,  130.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
           0.,   10.,    0.,   20.,   10.,   10.,    0.,  170.,   70.,
           0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
           0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.],
 [  0.,   0.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,  10.,   0.,
         10.,   0.,   0.,  10.,   0.,  10.,   0.,  10.,   0.,  10.,   0.,
         20.,   0.,  20., 270., -40., 130., 360.,   0.,   0.,   0.,   0.,
          0.,  10.,   0.,   0.,   0.,  10.,   0.,  30.,   0.,  10.,  10.,
         40.,   0.,  10.,   0.,   0.,   0.,  20.,   0.,   0.,   0.,  10.,
          0.,  10.,   0.,  10.,   0.,   0.,   0.],
 [   0.,    0.,    0.,    0.,    0.,    0.,   10.,    0.,    0.,
           0.,    0.,    0.,    0.,    0.,    0.,   10.,    0.,    0.,
           0.,   10.,    0.,   10.,   10.,   50.,   10.,  400.,  130.,
        -460.,  310.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
           0.,    0.,    0.,   30.,    0.,   10.,    0.,  320.,  150.,
           0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
           0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.],
 [   0.,    0.,    0.,    0.,    0.,    0.,    0.,   20.,    0.,
           0.,    0.,    0.,   10.,    0.,   10.,    0.,   10.,    0.,
          10.,    0.,   10.,    0.,   20.,    0.,   50.,  130.,  360.,
         310., -310.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
          30.,    0.,    0.,    0.,   40.,    0.,   10.,   20.,  190.,
           0.,    0.,    0.,   10.,    0.,    0.,    0.,   30.,    0.,
           0.,    0.,   20.,    0.,   10.,    0.,   10.,    0.],
 [ 10.,   0.,  20.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,  30.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0., 480., 140., 160.,  10.,
         10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,  70.,   0.,  30.,  10.,  10.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.],
 [  0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,  40.,   0.,  20.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0., 140., 590.,  60.,  10.,
          0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,  10.,  70.,   0.,  30.,   0.,  10.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.],
 [ 10.,   0.,  50.,  10.,   0.,   0.,   0.,   0.,  10.,   0.,   0.,
          0.,   0.,  10.,   0.,  40.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0., 160.,  60., 430.,  20.,
          0.,   0.,  20.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,  10.,
          0.,  30.,   0., 110.,  10.,   0.,   0.,  10.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.],
 [ 10.,  10.,  10.,  30.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,
          0.,   0.,  10.,  10.,  20.,  20.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,  10.,   0.,   0.,   0.,   0.,  10.,  10.,  20., 750.,
          0.,   0.,  10.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         10.,   0.,   0.,  10.,  10.,   0.,   0.,  10.,  10.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.],
 [  0.,   0.,   0.,   0.,  50.,   0.,  20.,   0.,  10.,   0.,  10.,
          0.,   0.,  10.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,  10.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,
        140., 140., 320., 120.,  20.,   0.,  10.,   0.,  50.,  10.,   0.,
          0.,  30.,   0.,   0.,   0.,  20.,   0.,   0.,   0.,  10.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.],
 [  0.,   0.,   0.,   0.,   0.,  40.,   0.,  20.,   0.,   0.,   0.,
         10.,   0.,   0.,  10.,   0.,   0.,   0.,  10.,   0.,  10.,   0.,
         10.,   0.,  10.,   0.,  10.,   0.,   0.,   0.,  10.,   0.,   0.,
        140., 210.,  80., 280.,   0.,  10.,   0.,  10.,   0.,  50.,   0.,
          0.,   0.,  30.,   0.,   0.,   0.,  20.,   0.,  10.,   0.,  10.,
          0.,   0.,   0.,  10.,   0.,   0.,   0.],
 [  0.,   0.,  10.,   0.,  20.,   0.,  60.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,  10.,   0.,   0.,
          0.,  20.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  20.,  10.,
        320.,  80., 100., 190.,  10.,   0.,  20.,   0.,  20.,   0.,  10.,
         10.,   0.,   0.,  40.,   0.,   0.,   0.,  30.,   0.,   0.,   0.,
         10.,   0.,   0.,   0.,   0.,   0.,   0.],
 [  0.,   0.,   0.,  10.,   0.,  20.,  10.,  50.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,  10.,   0.,
          0.,   0.,  20.,   0.,   0.,   0.,  30.,   0.,   0.,   0.,  10.,
        120., 280., 190.,  90.,   0.,  10.,   0.,  10.,  10.,  30.,   0.,
         10.,   0.,   0.,   0.,  30.,   0.,   0.,   0.,  30.,   0.,   0.,
          0.,  10.,   0.,   0.,   0.,  10.,   0.],
 [  0.,   0.,   0.,   0.,  20.,   0.,  10.,   0.,  10.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,  30.,
          0.,  20.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         20.,   0.,  10.,   0., 500., 180.,  20.,  10.,  50.,  10.,  10.,
          0.,   0.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,  40.,   0.,
         10.,   0.,  10.,   0.,  10.,   0.,   0.],
 [  0.,   0.,   0.,   0.,   0.,  10.,   0.,  10.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,
         20.,   0.,  10.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,  10.,   0.,  10., 180., 610.,   0.,  10.,  10.,  40.,   0.,
         10.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,  30.,
          0.,   0.,   0.,  10.,   0.,   0.,   0.],
 [  0.,   0.,  10.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,  10.,
          0.,  30.,   0.,  20.,   0.,  30.,   0.,   0.,   0.,  10.,   0.,
         10.,   0.,  20.,   0.,  20.,   0., 540., 130.,  10.,   0.,  60.,
         20.,   0.,   0.,  10.,   0.,   0.,   0.,  10.,   0.,  10.,   0.,
         20.,   0.,   0.,   0.,  10.,   0.,   0.],
 [  0.,   0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,
         10.,   0.,  20.,  10.,  30.,   0.,  40.,   0.,   0.,   0.,   0.,
          0.,  10.,   0.,  10.,  10.,  10., 130., 570.,   0.,  10.,  10.,
         70.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,   0.,
          0.,  20.,   0.,   0.,   0.,  10.,   0.],
 [  0.,   0.,   0.,   0., 120.,   0.,  60.,  10.,  10.,   0.,  20.,
          0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,  10.,
          0.,  10.,   0.,  10.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,
         50.,   0.,  20.,  10.,  50.,  10.,  10.,   0., 160., 290.,  20.,
         10.,  10.,   0.,   0.,   0.,  20.,   0.,   0.,   0.,  20.,   0.,
         10.,   0.,  20.,   0.,  10.,  10.,   0.],
 [  0.,   0.,   0.,   0.,   0., 110.,  10., 100.,   0.,   0.,   0.,
         10.,   0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,
         10.,   0.,  10.,   0.,  10.,   0.,  10.,   0.,   0.,   0.,   0.,
         10.,  50.,   0.,  30.,  10.,  40.,   0.,  10., 290., 190.,   0.,
         10.,   0.,   0.,   0.,   0.,   0.,  20.,   0.,   0.,   0.,  20.,
          0.,  10.,   0.,  20.,   0.,  10.,   0.],
 [   0.,    0.,    0.,    0.,    0.,    0.,   10.,    0.,    0.,
           0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
           0.,    0.,    0.,    0.,    0.,   10.,    0.,  170.,   10.,
         320.,   20.,    0.,    0.,   10.,    0.,    0.,    0.,   10.,
           0.,   10.,    0.,   60.,   10.,   20.,    0., -130.,  440.,
           0.,    0.,   10.,    0.,    0.,    0.,    0.,    0.,    0.,
           0.,   10.,    0.,    0.,    0.,   10.,    0.,    0.],
 [   0.,    0.,    0.,   10.,    0.,    0.,    0.,   10.,    0.,
           0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
           0.,    0.,    0.,    0.,    0.,    0.,   10.,   70.,   40.,
         150.,  190.,    0.,    0.,    0.,   10.,    0.,    0.,   10.,
          10.,    0.,   10.,   20.,   70.,   10.,   10.,  440., -100.,
           0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
           0.,    0.,   10.,    0.,    0.,    0.,   20.,    0.],
 [ 10.,   0.,  10.,  10.,  10.,   0.,   0.,   0.,  10.,   0.,  20.,
          0.,   0.,  20.,   0.,  10.,   0.,  10.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.,  70.,  10.,  30.,   0.,
         30.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,   0.,
          0., 250., 110., 250.,  90.,  30.,   0.,  10.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.],
 [  0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,
         10.,   0.,   0.,  30.,   0.,  10.,   0.,  10.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,  70.,   0.,   0.,
          0.,  30.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
          0., 110., 410.,  60., 190.,   0.,  30.,   0.,  10.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.],
 [ 10.,   0.,  20.,  10.,   0.,   0.,  10.,   0.,   0.,   0.,  10.,
          0.,   0.,  10.,   0.,  20.,   0.,   0.,   0.,  10.,   0.,   0.,
          0.,  10.,   0.,   0.,   0.,   0.,   0.,  30.,   0., 110.,  10.,
          0.,   0.,  40.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,  10.,
          0., 250.,  60., 190., 140.,  10.,   0.,  30.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.],
 [  0.,   0.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         10.,   0.,   0.,   0.,   0.,  20.,   0.,   0.,   0.,  10.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,  10.,  10.,  30.,  10.,  10.,
          0.,   0.,   0.,  30.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,  90., 190., 140., 400.,   0.,  10.,   0.,  20.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.],
 [  0.,   0.,   0.,   0.,  40.,  10.,  20.,   0.,   0.,   0.,  20.,
          0.,   0.,  10.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,
         20.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,  20.,   0.,   0.,
          0.,  30.,   0.,  10.,   0., 360., 100., 240.,  50.,   0.,   0.,
         10.,   0.,  20.,   0.,  10.,   0.,   0.],
 [  0.,  10.,   0.,   0.,   0.,  30.,   0.,  20.,   0.,   0.,   0.,
         20.,   0.,   0.,  10.,   0.,   0.,   0.,  10.,   0.,  10.,   0.,
         10.,   0.,  10.,   0.,  20.,   0.,   0.,   0.,  10.,   0.,   0.,
          0.,  20.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,  20.,   0.,
          0.,   0.,  30.,   0.,  10., 100., 350.,  70., 190.,   0.,  10.,
          0.,  10.,   0.,  20.,   0.,   0.,   0.],
 [  0.,   0.,  10.,   0.,  20.,  10.,  40.,  10.,   0.,   0.,  10.,
          0.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,  10.,   0.,   0.,
          0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  10.,  10.,
          0.,   0.,  30.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,   0.,
          0.,  10.,   0.,  30.,   0., 240.,  70., 290., 120.,   0.,   0.,
         10.,   0.,  10.,   0.,  20.,  10.,   0.],
 [  0.,   0.,   0.,   0.,   0.,  10.,   0.,  50.,   0.,   0.,   0.,
         10.,   0.,   0.,   0.,   0.,  20.,   0.,   0.,   0.,  20.,   0.,
          0.,   0.,  20.,   0.,   0.,   0.,  30.,   0.,   0.,   0.,  10.,
          0.,  10.,   0.,  30.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,
          0.,   0.,  10.,   0.,  20.,  50., 190., 120., 340.,   0.,  10.,
          0.,  10.,   0.,  10.,   0.,  20.,   0.],
 [  0.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,  10.,
          0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         10.,   0.,   0.,   0.,  40.,   0.,  10.,   0.,  20.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0., 690., 140.,
         30.,  10.,  10.,   0.,   0.,   0.,   0.],
 [  0.,   0.,   0.,   0.,   0.,  10.,   0.,  10.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,  10.,   0.,
         10.,   0.,  10.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,  10.,   0.,   0.,   0.,  30.,   0.,   0.,   0.,  20.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,  10., 140., 660.,
         10.,  30.,   0.,  10.,   0.,   0.,   0.],
 [  0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,   0.,
          0.,  30.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,  10.,   0.,  10.,   0.,  20.,   0.,  10.,   0.,  10.,
          0.,   0.,   0.,   0.,   0.,  10.,   0.,  10.,   0.,  30.,  10.,
        680., 130.,   0.,   0.,  10.,   0.,   0.],
 [  0.,   0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,
         10.,   0.,  30.,   0.,  10.,   0.,  20.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,  10.,   0.,   0.,   0.,  20.,   0.,  10.,   0.,
         10.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,  10.,  10.,  30.,
        130., 660.,   0.,   0.,   0.,  10.,   0.],
 [  0.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,  20.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,  20.,   0.,  10.,   0.,  10.,   0.,
          0.,   0., 430., 140., 210., 140.,   0.],
 [  0.,   0.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,
         10.,   0.,   0.,   0.,  10.,   0.,  10.,   0.,   0.,   0.,   0.,
          0.,  10.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,  20.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,  20.,   0.,  10.,   0.,  10.,
          0.,   0., 140., 490.,  50., 190.,   0.],
 [  0.,   0.,   0.,   0.,   0.,   0.,  20.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,   0.,
          0.,  10.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,  10.,   0.,  10.,   0.,  10.,   0.,  10.,
          0.,   0.,   0.,   0.,   0.,  10.,   0.,  20.,   0.,   0.,   0.,
         10.,   0., 210.,  50., 360., 260.,   0.],
[  0.,   0.,   0.,   0.,   0.,   0.,   0.,  20.,   0.,   0.,   0.,
          0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,  10.,   0.,
          0.,   0.,  10.,   0.,   0.,   0.,  10.,   0.,   0.,   0.,   0.,
          0.,   0.,   0.,  10.,   0.,   0.,   0.,  10.,  10.,  10.,   0.,
         20.,   0.,   0.,   0.,   0.,   0.,   0.,  10.,  20.,   0.,   0.,
          0.,  10., 140., 190., 260., 260.,   0.],
[   0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
           0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
           0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
           0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
           0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
           0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
           0.,    0.,    0.,    0.,    0.,    0.,    0., 1000.]]
    
    ONE_NT_DIST = {'A': 24, 'C': 32, 'T': 19, 'G': 25}

    TWO_NT_DIST = {
		'A': {'A': 21, 'C': 28, 'T': 21, 'G': 30},
		'C': {'A': 27, 'C': 38, 'T': 21, 'G': 14},
		'T': {'A': 12, 'C': 35, 'T': 20, 'G': 33},
		'G': {'A': 29, 'C': 30, 'T': 19, 'G': 22}
	}

    THREE_NT_DIST = {
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
    
    NT_INTRON_DIST = {'A': 29, 'C': 20, 'T': 30, 'G': 21}
    MEAN_EXON_GENE = 0
    SD_EXON_GENE = 0
    MEAN_EXON_CDS = 9
    SD_EXON_CDS = 8
    MEAN_EXON_LEN = 170
    SD_EXON_LEN = 258
    MEAN_INTRON_LEN = 3760
    SD_INTRON_LEN = 20126
    
    #EIC_EL= 0.4
    #EIC_ED=0.1
    #EIC_EG=0.5
    #C_I = 0.4 # The relative frequencies of insertion event
    #C_D = 0.6 # The relative frequencies of deletion event
    
    STOP_CODONS = ['TAG', 'TAA', 'TGA']
    
    SPLICE_SITES = [["GT-AG"]*90,  ["GC-AG"]*7, [""]*3]
    
    ###########################
    
    tree = Tree(TREE_INPUT)
    
    print(tree)
    
    # Create output folders
    create_output_folders(SRC, ITERATION_NAME)
    
    all_one_nt_dist, all_two_nt_dist, all_three_nt_dist, all_nt_intron_dist = random_generator_nucleotide(ONE_NT_DIST, TWO_NT_DIST, THREE_NT_DIST, NT_INTRON_DIST)
    
    # Retrieve the maximum number of exons in transcripts
    NUMBER_MAX_EXON_TRANSCRIPT = random_generator_max_exon_transcript(MEAN_EXON_CDS, SD_EXON_CDS)
    
    # Retrieve the maximum number of exons in the gene
    NUMBER_MAX_EXON_GENE = random_generator_max_exon_gene(NUMBER_MAX_EXON_TRANSCRIPT, K_NB_EXONS)
    
    # Retrieve the lengths of exons sequences
    exons_lengths_not_adjusted = [random_number_generator_by_normal_law(MEAN_EXON_LEN, SD_EXON_LEN) for _ in range(NUMBER_MAX_EXON_GENE)]
    
    exons_lengths = [exon_length - exon_length%3 + 3 for exon_length in exons_lengths_not_adjusted]
    
    # Retrieve the sequences of exons
    end_exon = ""
    exons_dict = {}
    for i in range(len(exons_lengths)):
        length = exons_lengths[i]
        exon, end_exon = exon_generator(length, all_one_nt_dist, all_two_nt_dist, all_three_nt_dist, STOP_CODONS, end_exon)
        exons_dict["exon_" + str(i)] = exon
    get_exons_extremities_adjusted(exons_dict, end_exon, STOP_CODONS)
    
    # Retrieve the sequences of introns
    introns_dict = {}
    for i in range(len(exons_lengths)-1):
        length = random_number_generator_by_normal_law(MEAN_INTRON_LEN, SD_INTRON_LEN)
        intron = intron_generator(length, SPLICE_SITES, all_nt_intron_dist)
        introns_dict["intron_" + str(i)] = intron
    
    ancestral_gene_structure = write_init_structure(exons_dict, introns_dict)
    
    # Retrieve transcripts (exon-composition and their sequences)
    if TC_RS == 0.0:
        NUMBER_OF_TRANSCRIPTS_IN_LUCA = random_number_generator_by_normal_law(3, 1)
        transcripts_dict = transcripts_generation_ancestor(NUMBER_OF_TRANSCRIPTS_IN_LUCA, NUMBER_MAX_EXON_TRANSCRIPT, exons_dict, TC_RS)
    else:
        transcripts_dict = transcripts_generation_ancestor(0, 0, exons_dict, TC_RS)
        
    
    # simulation
    print(colorsSimulation.HEADER + '## ==== Simulation' + colorsSimulation.ENDC + colorsSimulation.OKGREEN + '\t<--->\t STARTED' + colorsSimulation.ENDC)
    tree_info, forest, lca_events = simulation_program(tree, 
                       ancestral_gene_structure, exons_dict, introns_dict, transcripts_dict,
                       EIC_EL, EIC_ED, EIC_EG, C_I, C_D, K_EIC, K_INDEL, 
                       TC_RS, TC_TL, TC_ES, TC_ME, TC_A3, TC_A5, TC_IR, K_TC,
                       all_one_nt_dist, all_two_nt_dist, all_three_nt_dist, STOP_CODONS,
                       MEAN_EXON_LEN, SD_EXON_LEN,
                       all_nt_intron_dist, MEAN_INTRON_LEN, SD_INTRON_LEN, SPLICE_SITES,
                       CODONS_MATRIX, CODONS_LIST
                       )

    # Retrieve multiple sequence alignment
    print(colorsSimulation.WARNING + '#### \t+++ RETRIEVING CLUSTERS and Phylogenies \t<--->\t STARTED'+ colorsSimulation.ENDC)
    df_transcripts = get_msa_and_clusters(tree_info, forest, lca_events, SRC, ITERATION_NAME)
    print(colorsSimulation.WARNING + '#### \t+++ RETRIEVING CLUSTERS and Phylogenies \t<--->\t SUCCESSFULLY  FINISHED'+ colorsSimulation.ENDC)
    
    print(colorsSimulation.OKCYAN +'###### \t+++ SAVING RESULTS \t<--->\t STARTED'+ colorsSimulation.WARNING)
    
    write_data(df_transcripts, tree_info, SRC, ITERATION_NAME)
    
    print(colorsSimulation.OKCYAN +'###### \t+++ SAVING RESULTS \t<--->\t SUCCESSFULLY FINISHED' + colorsSimulation.WARNING)
    
    print(colorsSimulation.HEADER + '## ==== Simulation' + colorsSimulation.ENDC + colorsSimulation.OKGREEN + '\t<--->\t SUCCESSFULLY FINISHED' + colorsSimulation.WARNING)
    
    return tree_info, df_transcripts

def successful_message():
    message = f"""
    🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟
    🌟                                                                                          🌟
    🌟   SimSpliceEvol was successfully executed! Check the results in the output directory     🌟
    🌟                                                                                          🌟
    🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟🌟
    """
    print(message)
    
def error_message():
    message = f"""
    ❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌
    ❌                                                                                                          ❌
    ❌   An Error occured! Check it out. Contact us for any questions by checking the github repository at      ❌
    ❌                             https://github.com/UdeS-CoBIUS/SimSpliceEvolEvol                             ❌
    ❌                                                                                                          ❌
    ❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌❌
    """
    print(message)

if __name__ == '__main__':
    
    #==========================================
    # PARAMETERS USER-DEFINED
    #==========================================
    
    parser = build_arg_parser()
    args = parser.parse_args()
    
    SRC = str(args.directory_name)
    ITERATION = str(args.iteration)
    
    EIC_EL= float(args.eic_el)
    EIC_ED = float(args.eic_ed)
    EIC_EG =float(args.eic_eg)
    
    C_I = float(args.c_i)
    C_D = float(args.c_d)
    
    TREE_INPUT = args.input_tree_file
    K_NB_EXONS = float(args.k_nb_exons)
    K_EIC = float(args.k_eic)
    K_INDEL = float(args.k_indel)
    K_TC = float(args.k_tc)
    
    TC_RS = float(args.random_selection)
    TC_A5 = float(args.alternative_five_prime)
    TC_A3 = float(args.alternative_three_prime)
    TC_ME = float(args.mutually_exclusive)
    TC_ES = float(args.exon_skipping)
    TC_IR = float(args.intron_retention)
    TC_TL = float(args.transcript_loss)
    
    #==========================================
    #ALGORITHM
    #==========================================
    
    #ITERATION = 1
    
    ITERATION_NAME = "_iteration_"+str(ITERATION)
    
    #inputs = [SRC, ITERATION_NAME, TREE_INPUT, K_NB_EXONS, K_INDEL, C_I, C_D, EIC_ED, EIC_EG, EIC_EL, K_EIC, K_TC, TC_RS, TC_A3, TC_A5, TC_ME, TC_ES, TC_IR, TC_TL]

    #simspliceevol(SRC, ITERATION_NAME, TREE_INPUT, K_NB_EXONS, K_INDEL, C_I, C_D, EIC_ED, EIC_EG, EIC_EL, K_EIC, K_TC, TC_RS, TC_A3, TC_A5, TC_ME, TC_ES, TC_IR, TC_TL)
    
    try:
        simspliceevol(SRC, ITERATION_NAME, TREE_INPUT, K_NB_EXONS, K_INDEL, C_I, C_D, EIC_ED, EIC_EG, EIC_EL, K_EIC, K_TC, TC_RS, TC_A3, TC_A5, TC_ME, TC_ES, TC_IR, TC_TL)
        successful_message()
    except:
        error_message()
        
