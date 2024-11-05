# K-MER Analysis #




K-mer_proportion is a PERL software suite designed to analyze K-mer representativity in coding sequence datasets, supporting both genomic and genic data analysis. This package was developed to investigate the APOBEC3 mutational footprint in human virus genomes, as described in this study (https://doi.org/10.1371/journal.ppat.1008718).

In brief, we estimate K-mer representation in synthetic genomes by calculating the ratio of observed versus expected K-mer counts, following the approach described by Werman et al. (https://doi.org/10.1093/ve/vev015). A K-mer is defined as a set of sequences sharing a common motif. For instance, the NTC K-mer includes the sequences ATC, CTC, GTC, and TTC. The observed count of K-mers is directly calculated from the synthetic genome, while the expected count is estimated by averaging K-mer counts across one thousand random shuffles of the sequences. This null model, based on one thousand iterations, provides a baseline representation of random K-mer occurrence in sequences with the same nucleotide composition.

The resulting observed-to-expected K-mer ratio indicates the representation of each K-mer in the analyzed sequence. A ratio < 0 suggests underrepresentation of the K-mer, while a ratio around 0 indicates no significant deviation from expected representation. 

The data repository contains FASTA files of coding sequences (CDS) for human viruses, along with their corresponding GenBank annotation files.

## GENOMIC K-MER analysis ##
###  Synthetic coding genome
cds_ncds_split.pl

For each genome, a "synthetic coding genome" is generated by concatenating coding sequences.

###  Input sequence split to reduce counter analysis 
cds_ncds_split.02.pl

The synthetic coding genome FASTA file is split to facilitate parallel processing.

###  K-mer ratio determination
Trikmer_analyser_0.5_global1000.pl

This script analyzes the preformatted files to determine K-mer ratios.


## GENIC K-MER analysis ##


###  Input sequence split to reduce counter analysis 
cds_ncds_split.03.pl

The synthetic coding genome FASTA file is split to facilitate parallel processing.

###  K-mer ratio determination
Trikmer_analyser_0.5_genic1000.pl

This script analyzes the preformatted files to determine K-mer ratios specific to genic regions.

