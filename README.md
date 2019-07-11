# SimSpliceEvol
SimSpliceEvol: Alternative splicing-aware simulation of biological sequence evolution

### Esaie Kuitche, Safa Jammali and Aïda Ouangraoua 

### Contact Esaie.Kuitche.Kamela@USherbrooke.ca

### Requirements The program requires the following to be available

* python (3.7)
* ete3 (https://pypi.python.org/pypi/ete3/)
* argparse (https://pypi.python.org/pypi/argparse)

### Usage
```
usage: SimSpliceEvol.py [-h] [-i INPUT_TREE_FILE] [-n NUMBER_OF_SIMULATION]
                        [-k_indel K_INDEL] [-k_eic K_EIC] [-k_tc K_TC]
                        [-k_intron K_INTRON] [-k_nb_exons K_NB_EXONS]
                        [-eic_l EXON_I_CHANGE_LOSS]
                        [-eic_g EXON_I_CHANGE_GAIN] [-eic_d EXON_I_CHANGE_DUP]
                        [-tc_rs RANDOM_SELECTION]
                        [-tc_a5 ALTERNATIVE_FIVE_PRIME]
                        [-tc_a3 ALTERNATIVE_THREE_PRIME]
                        [-tc_ek EXON_SKIPPING] [-tc_me MUTUALLY_EXCLUSIVE]
                        [-tc_ir INTRON_RETENTION] [-tc_tl TRANSCRIPT_LOSS]

SimSpliceEvol program parameters

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_TREE_FILE, --input_tree_file INPUT_TREE_FILE
                        input guide tree (default = Example/input/large.nw)
  -n NUMBER_OF_SIMULATION, --number_of_simulation NUMBER_OF_SIMULATION
                        number of simulations
  -k_indel K_INDEL, --k_indel K_INDEL
                        multiplicative constant for codon indel rate (default=0.5)
  -k_eic K_EIC, --k_eic K_EIC
                        multiplicative constant for exon-intron change (eic) rate (default =5)
  -k_tc K_TC, --k_tc K_TC
                        multiplicative constant for transcript change (tc) rate (default =5)
  -k_intron K_INTRON, --k_intron K_INTRON
                        multiplicative constant for substitution rate in intron (default =1.5)
  -k_nb_exons K_NB_EXONS, --k_nb_exons K_NB_EXONS
                        multiplicative constant for number of exons in gene (default =1.5)
  -eic_l EXON_I_CHANGE_LOSS, --exon_i_change_loss EXON_I_CHANGE_LOSS
                        relative frequence of exon loss in eic (default = 0.4)
  -eic_g EXON_I_CHANGE_GAIN, --exon_i_change_gain EXON_I_CHANGE_GAIN
                        relative frequence of exon gain in eic (default = 0.5)
  -eic_d EXON_I_CHANGE_DUP, --exon_i_change_dup EXON_I_CHANGE_DUP
                        relative frequence of exon duplication in eic (default = 0.1)
  -tc_rs RANDOM_SELECTION, --random_selection RANDOM_SELECTION
                        relative frequence of random selection in tc (default
                        =0.5)
  -tc_a5 ALTERNATIVE_FIVE_PRIME, --alternative_five_prime ALTERNATIVE_FIVE_PRIME
                        relative frequence of alternative five prime in tc
                        (default =0.1)
  -tc_a3 ALTERNATIVE_THREE_PRIME, --alternative_three_prime ALTERNATIVE_THREE_PRIME
                        relative frequence of alternative three prime in tc
                        (default =0.1)
  -tc_ek EXON_SKIPPING, --exon_skipping EXON_SKIPPING
                        relative frequence of exon skipping in tc (default
                        =0.1)
  -tc_me MUTUALLY_EXCLUSIVE, --mutually_exclusive MUTUALLY_EXCLUSIVE
                        relative frequence of mutually exclusive in tc
                        (default =0.1)
  -tc_ir INTRON_RETENTION, --intron_retention INTRON_RETENTION
                        relative frequence of intron retention in tc (default
                        =0.05)
  -tc_tl TRANSCRIPT_LOSS, --transcript_loss TRANSCRIPT_LOSS
                        relative frequence of transcript loss in tc (default
                        =0.4)
```

## Input files

### example of input guide tree (examples in Examples/input/<filename>.nw)

```
(ggal:0.129452,(mdom:0.116455,((mmus:0.114942,hsap:0.0962098)1:0.0001,btaus:0.11365)1:0.0143738)1:0.0100642);

```
## Ouptut directory and files
* output/<filename>/cds: fasta file of CDS
* output/<filename>/cds_gene: mapping file for CDS and genes
* output/<filename>/cluster: splicing orthology groups of CDS
* output/<filename>/genes: fasta file of gene sequences
* output/<filename>/multiple_alignment: multiple alignment of all CDS
* output/<filename>/pairwise_alignment: pairwise splice alignment of all genes against all CDS
* output/<filename>/positions: exon composition of each CDS with location in CDS and gene

## Running SimSpliceEvol on an example
```
python3.7 src/SimSpliceEvol.py
```
