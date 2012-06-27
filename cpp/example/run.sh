#!/bin/bash

../build/Debug/demos/makeAln Phix_S1_L001_R1_001.fastq Phix_S1_L001_R2_001.fastq library.fa > therm1.dat

### This creates a dataset containing
#CATC CGGATAGAAAAGAAACAACAACAACAAC 6 -1 ATTAATTCTTTAATATAAACTATCCGTTCG
#multiplex_id,seq_id,pos,quality,align_frag
