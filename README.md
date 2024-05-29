# Least-mutated-sequence

-------------
1. Background

MutationSiteFinder v.1.0 and MutationSiteCounter v.1.0 are two corresponding programs of my manuscript entitled "A Method for Calculating the Least Mutated Sequence in DNA Alignment Based on Point Mutation Sites" (on bioRxiv, [https://doi.org/10.1101/2023.11.14.567125](https://www.biorxiv.org/content/10.1101/2023.11.14.567125v1)).

They are Perl programs which are designed to calculate the least mutated sequence in a DNA Alignment file.

MutationSiteFinder v.1.0 takes a DNA Alignment file in FASTA format and output a matrix file.

MutationSiteCounter v.1.0 takes a matrix file and output a sorted result.

-------------
2. Author

Libing Shen 

email: shenlibing@ihup.org.cn

-------------
3. Usage

There are two steps to calculate the least mutated sequence in a DNA Alignment file.

First, run perl MutationSiteFinder DNA_alignment_file.
Second, run perl MutationSiteCounter matrix_file.

Example: 

(Under prompt sign: > or $ or #)

perl MutationSiteFinder sample_data.fa

perl MutationSiteCounter mutation_site_matrix.txt


-------------
4. License

Libing Shen @ International Human Phenome Institutes (Shanghai)

-------------


5. Related preprint publications

The preprint manuscript: https://zenodo.org/records/10645490.

The preprint manuscript: https://zenodo.org/records/10703955.


