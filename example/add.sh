#!/bin/bash

../cpp/build/Debug/demos/makeAln -c AAAGAAACAACAACAACAAC -s ./data/Phix_S1_L001_R1_001.fastq -e ./data/Phix_S1_L001_R2_001.fastq -l ./data/library.fasta -b ./data/barcodes.fasta -n 9 -o ./data/add4lib.dat -d 0

