# Comparative Analysis of Amino Acid Auxotrophies and Peptidase Profiles in Non-Dysbiotic and Dysbiotic Small Intestinal Microbiomes

## Data preprocessing

### What are the relative abundances of HRGM genomes?

1. DADA2 workflow for ASV inference and quantification

   ```sh
   Rscript Scripts/DADA2_SIBO.R
   Rscript Scripts/DADA2_Reimagine.R
   ```

2. Mapping of 16S sequences (i.e., ASVs) to 16S rRNA genes from HRGM genomes.

   ```sh
   Scripts/./asv_to_HRGM.sh
   Scripts/./asv_to_HRGM_SIBO.sh
   ```

   

### Which protein encoding genes in HRGM genomes putatively have extracellular proteolytic/peptidolytic activity?

#### Prediction of genes with peptidolytic activity

1. Get the MEROPS database files `merops_scan.lib` (fasta seqs) and `domain.sql` (SQL database with all MEROPS Codes and the info if these are peptidases or inhibitors.

   https://www.ebi.ac.uk/merops/

   > The file "merops_scan.lib" is a subset of pepunit.lib containing the sequences used for the MEROPS batch Blast. It contains a non-redundant library of protein sequences in FastA format of the peptidase units for all the family type examples and peptidase/inhibitor holotypes. The library can be searched by use of FastA without further modification, but must be converted and indexed for BLAST searches. (Description from the MEROPS website).

1. Filter `domain.sql` to get a data table with only peptidases and omitting inhibitors

   ```sh
   Rscript Scripts/n_prepare_MeropsDomainTable.R
   ```

2. Filter MEROPS Database "Scan Sequences" (=non-redundant library of protein sequences in FastA format of the peptidase units for all the <u>family type examples</u> and peptidase/inhibitor holotypes) and retain all sequences that are listed in the table `Data/domain2.csv.gz`, which is an output of (2).

   ```sh
   Rscript Scripts/n_filterMeropsScanSeqs.R
   ```

3. Search for HRGM genes with putative peptidase activity. 

   ```sh
   makeblastdb -in merops_scan_filt.lib -dbtype prot
   blastp -query combined.faa -db merops_scan_filt.lib -evalue 0.01 -num_threads 12 -outfmt 6 > merops_HRGM.m8
   gzip merops_HRGM.m8
   ```

   where the file `merops_scan_filt.lib` is an output of (3), the file `combined.faa` is a fasta file with all protein sequences of all HRGM genomes.

#### Prediction of putative peptidases with extracellular catalytic activity

All protein-coding genes were processed with SignalP-6.0 using the following script/call, to predict signal peptides that control protein secretion and translocation.

```sh
#!/bin/bash

signalp6 --fastafile $sample.faa --organism other --output_dir $sample --format txt --mode fast --model_dir signalp-6-package/models/
```

In this script, the variable `$sample ` is the respective HRGM genome ID e.g.: "HRGM_Genome_0003".



Finally, per genome, the intersect between sequences with a signal peptide for secretion/translocation and the sequences with putative peptidolytic activity were determined. Only those sequences (putatively extracellular & putativly peptidolytic activity) were considered in all following analysis (below).


## Final data analysis and visualisation

#### Initialize all required data in R

```R
source("n_initData.R")
```

#### Analyse REIMAGE cohort data

```R
source("n_reimagine.R")
```

#### Analyse SIBO cohort data

```R
source("n_sibo.R")
```