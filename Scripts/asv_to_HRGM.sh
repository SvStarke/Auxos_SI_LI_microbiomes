#!/bin/bash

makeblastdb -in Data/16S_rRNA_library_filtered_edit.fasta -dbtype nucl
blastn -query Data/asv_seqs.fna -db Data/16S_rRNA_library_filtered_edit.fasta -perc_identity 95 -qcov_hsp_perc 95 -outfmt 6 > Output/asvs_to_HRGM.m8

