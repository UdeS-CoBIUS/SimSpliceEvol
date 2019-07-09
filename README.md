# SimSpliceEvol
SimSpliceEvol: Alternative splicing-aware simulation of biological sequence evolution

####Esaie Kuitche, Safa Jammali and Aïda Ouangraoua 
#####Contact Esaie.Kuitche.Kamela@USherbrooke.ca

##Requirements The program requires the following to be available

* python (3.7)
* argparse (https://pypi.python.org/pypi/argparse)
* ete3 (https://pypi.python.org/pypi/ete3/)
* argparse (https://pypi.python.org/pypi/argparse)

##Usage
```
usage: SimSpliceEvol.py [-h] [-ei_c_l EXON_I_CHANGE_LOSS]
                        [-ei_c_g EXON_I_CHANGE_GAIN]
                        [-ei_c_d EXON_I_CHANGE_DUP] [-tc_rs RANDOM_SELECTION]
                        [-tc_a5 ALTERNATIVE_FIVE_PRIME]
                        [-tc_a3 ALTERNATIVE_THREE_PRIME]
                        [-tc_ek EXON_SKIPPING] [-tc_me MUTUALLY_EXCLUSIVE]
                        [-tc_ir INTRON_RETENTION] [-tc_tl TRANSCRIPT_LOSS]
                        [-k_indel K_INDEL] [-k_ei_c K_EI_C] [-k_t_c K_T_C]
                        [-k_intron K_INTRON] [-k_nb_exons K_NB_EXONS]
                        [-n NUMBER_OF_SIMULATION] [-i INPUT_TREE_FILE]

SimSpliceEvol program parameters

optional arguments:
  -h, --help            show this help message and exit
  -ei_c_l EXON_I_CHANGE_LOSS, --exon_i_change_loss EXON_I_CHANGE_LOSS
                        relative frequence of exon loss
  -ei_c_g EXON_I_CHANGE_GAIN, --exon_i_change_gain EXON_I_CHANGE_GAIN
                        relative frequence of exon gain
  -ei_c_d EXON_I_CHANGE_DUP, --exon_i_change_dup EXON_I_CHANGE_DUP
                        relative frequence of exon duplication
  -tc_rs RANDOM_SELECTION, --random_selection RANDOM_SELECTION
                        relative frequence of random selection
  -tc_a5 ALTERNATIVE_FIVE_PRIME, --alternative_five_prime ALTERNATIVE_FIVE_PRIME
                        relative frequence of alternative five prime
  -tc_a3 ALTERNATIVE_THREE_PRIME, --alternative_three_prime ALTERNATIVE_THREE_PRIME
                        relative frequence of alternative three prime
  -tc_ek EXON_SKIPPING, --exon_skipping EXON_SKIPPING
                        relative frequence of exon skipping
  -tc_me MUTUALLY_EXCLUSIVE, --mutually_exclusive MUTUALLY_EXCLUSIVE
                        relative frequence of mutually exclusive
  -tc_ir INTRON_RETENTION, --intron_retention INTRON_RETENTION
                        relative frequence of intron retention
  -tc_tl TRANSCRIPT_LOSS, --transcript_loss TRANSCRIPT_LOSS
                        relative frequence of transcript loss
  -k_indel K_INDEL, --k_indel K_INDEL
                        user-defined constant
  -k_ei_c K_EI_C, --k_ei_c K_EI_C
                        user-defined constant
  -k_t_c K_T_C, --k_t_c K_T_C
                        user-defined constant
  -k_intron K_INTRON, --k_intron K_INTRON
                        user-defined constant
  -k_nb_exons K_NB_EXONS, --k_nb_exons K_NB_EXONS
                        user-defined constant
  -n NUMBER_OF_SIMULATION, --number_of_simulation NUMBER_OF_SIMULATION
                        Number of simulation wanted
  -i INPUT_TREE_FILE, --input_tree_file INPUT_TREE_FILE
                        input guide Tree

```

##Input files
###example of the input guide tree
Example
```
(ggal:0.129452,(mdom:0.116455,((mmus:0.114942,hsap:0.0962098)1:0.0001,btaus:0.11365)1:0.0143738)1:0.0100642);

```
##Ouptut
* output/cds fasta file of CDS
* output/cds_gene mapping file of CDS with gene
* output/cluster orthology group of CDS
* output/genes fasta file of genes
* output/multiple_alignment multiple alignment of all CDS
* output/pairwise_alignment alignment of gene against its CDS
* output/positions positions of each exon in its CDS and its the gene

##Running SimSpliceEvol on an example
```
python3.7 src/SimSpliceEvol.py
```
