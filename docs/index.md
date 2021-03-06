# PESCA

PESCA (proportion of population-specific and shared causal variants) is a
software package for estimating proportions of population-specific and shared
causal variants from GWAS summary statistics.

## Software requirement

### For compiling and running PESCA

* g++ (version 4.7.2)
* [boost C++ library (1.55)](https://www.boost.org/)
* BLAS library

To compile the PESCA program, please `cd` into `src/`, and then type `make`.

## Data

### Genome partition

We provide definitions of regions that are approximately independent in
EAS and EUR [here](https://github.com/huwenboshi/pesca/tree/master/partition).

### Reference panel

Reference panels should be in [PLINK format](https://www.cog-genomics.org/plink/2.0/input#bed).

We provide 1000 Genomes reference panel for EAS and EUR [here](https://ucla.box.com/s/o3k9nhosmtjgjufvjjfcshbjynkko8j5).
All SNPs in this reference panel have minor allele frequency greater than 5%
in both populations, and have been LD pruned with \\(R^2\\) threshold of 0.95.

## Contact

* **Huwenbo Shi**: hshi_AT_hsph.harvard.edu
* **Kathy Burch**: kathy.s.burch_AT_gmail.com
