#!/bin/bash 

### This program takes a list of sequences and labels them based on a comparison to a wild type sequence 
./label.sh ../data/library.lst  GGACGCTTCATATAATCCTAATGATATGGTTTGGGAGTTTCTACCAAGAGCCTTAAACTCTTGATTATGAAGTGAA > ../data/library.fasta


### For example 
# wt-  GGACGCTTCATATAATCCTAATGATATGGTTTGGGAGTTTCTACCAAGAGCCTTAAACTCTTGATTATGAAGTGAA
# seq- CGACGCTTCATATAATCCTAATGATATGGTTTGGGAGTTTCTACCAAGAGCCTTAAACTCTTGATTATGAAGTGAA

# The sequence seq would be labeled as

# >G1C
#CGACGCTTCATATAATCCTAATGATATGGTTTGGGAGTTTCTACCAAGAGCCTTAAACTCTTGATTATGAAGTGAA

