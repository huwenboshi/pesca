# Estimating genome-wide proportion of population-specific and shared causal variants

This page describes the pipeline to estimate the genome-wide prior
probabilities of a SNP to be causal in one population (i.e.
population-specific) or in both populations (i.e. shared).

### Recommended pipeline

#### Step 1: estimate genome-wide heritability of the trait

This step can be performend using [LDSC](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation).

#### Step 2: extract summary statistics and LD of each region for each chromosome

In this step, one

* Intersects the GWAS summary statistics with the reference panel. We provide
reference panels in PLINK format for East Asians and Europeans
[here](https://ucla.box.com/s/o3k9nhosmtjgjufvjjfcshbjynkko8j5). This reference
panel contains LD pruned SNPs (\\(R^2 > 0.95\\)) with MAF greater than 5% in
both populations.

* Splits GWAS summary statistics data into regions for each chromosome. Format
of the GWAS summary statistics data used by PESCA can be found
[here](https://huwenboshi.github.io/pesca/input_format/#gwas-summary-statistics-data).

* Computes the corresponding LD matrix of each region. We provide region
definitions that are approximately LD-independent across East Asians
and Europeans [here](https://github.com/huwenboshi/pesca/tree/master/partition).
Format of LD matrix can be found [here](https://huwenboshi.github.io/pesca/input_format/#ld-matrix).

We recommend using the following directory organization.

```nohighlight
<trait name>/
    chr1/
        chr1_region1_zscore_pop1.txt
        chr1_region1_zscore_pop2.txt
        chr1_region1_ld_pop1.txt
        chr1_region1_ld_pop2.txt
        ...
        chr1_region100_zscore_pop1.txt
        chr1_region100_zscore_pop2.txt
        chr1_region100_ld_pop1.txt
        chr1_region100_ld_pop2.txt
        ...
        chr1_zscore_list_pop1.txt
        chr1_zscore_list_pop2.txt
        chr1_ld_list_pop1.txt
        chr1_ld_list_pop2.txt
    chr2/
        ...
    ...
    chr22/
        ...
```

* `chr1_region1_zscore_pop1.txt` and `chr1_region2_zscore_pop2.txt` are
the GWAS summary statistics data files in chromosome 1 region 1 for
population 1 and 2 respectively
(see [here](https://huwenboshi.github.io/pesca/input_format/#gwas-summary-statistics-data)
for format).

* `chr1_region1_ld_pop1.txt` and `chr1_region1_ld_pop2.txt`
are the corresponding LD matrices (see [here](https://huwenboshi.github.io/pesca/input_format/#ld-matrix)
for format).

* `chr1_zscore_list_pop1.txt` is a list of paths to all the `chr1_region*_zscore_pop1.txt`
files, i.e. GWAS summary statistics data in all regions on chromosome 1 for
population 1. `chr1_zscore_list_pop2.txt` is a list of paths to all the
`chr1_region*_zscore_pop2.txt` files.

* `chr1_ld_list_pop1.txt` is a list of paths to all the `chr1_region*_ld_pop1.txt`
files, i.e. LD matrices for all regions on chromosome 1 for
population 1. `chr1_ld_list_pop2.txt` is a list of paths to all the
`chr1_region*_ld_pop2.txt` files.

#### Step 3: run PESCA for each chromosome

Command for running PESCA

PESCA estimates genome-wide proportions of population-specific and shared
causal variants via the following command.

```nohighlight
<directory to pesca>/pesca \
    --mode fit \
    --zscore1 <list of GWAS summary statistics files for population 1> \
    --zscore2 <list of GWAS summary statistics files for population 2> \
    --ld1 <list of LD files for population 1> \
    --ld2 <list of LD files for population 2> \
    --nburn <number of MCMC burn-ins, 5000 by default> \
    --nsample <number of MCMC samples, 5000 by default> \
    --lambda <shrinkage parameter> \
    --sigmasq1 $sigmasq1 \
    --sigmasq2  $sigmasq2 \
    --totnsnp <total number of SNPs across all chromosomes> \
    --max_iter <maximum number of EM iterations> \
    --out <output file name>
```

Here are the meaning of the flags:

* `--model` tells PESCA whether to estimate genome-wide prior or per-SNP
posterior probabilities. There are two options for this flag `fit` and `post`.
For estimating per-SNP posterior see
[here](https://huwenboshi.github.io/pesca/posterior/).

* `--zscore1` specifies a text file containing a list of paths to GWAS
summary statistics for population 1, one path per line. Typically, there
should be one such text file for each chromosome. And each line in the text
file corresponds to GWAS summary statistics of one region on that chromosome.
See [here](https://huwenboshi.github.io/pesca/input_format/#gwas-summary-statistics-data)
for format of GWAS summary statistics data for PESCA.

* `--zscore2` specifies a text file containing a list of paths to GWAS
summary statistics for population 2, one path per line. Typically, there
should be one such text file for each chromosome. And each line in the text
file corresponds to GWAS summary statistics of one region on that chromosome.
See [here](https://huwenboshi.github.io/pesca/input_format/#gwas-summary-statistics-data)
for format of GWAS summary statistics data for PESCA.

* `--ld1` specifies a text file containing a list of paths to LD matrices
for population 1, one path per line. Typically, there should be one such
text file for each chromosome. And each line in the text file corresponds
to LD matrix of one region on that chromosome. Additionally, each line
in the text file should correspond to the same region listed by the
`--zscore1` flag.

* `--ld2` specifies a text file containing a list of paths to LD matrices
for population 2, one path per line. Typically, there should be one such
text file for each chromosome. And each line in the text file corresponds
to LD matrix of one region on that chromosome. Additionally, each line
in the text file should correspond to the same region listed by the
`--zscore2` flag.

* `--nburn` specifies the number of burn-ins for the MCMC. The default is 5000.

* `--nsample` specifies the number of samples for the MCMC. The default is 5000.

* `--sigmasq1` specifies genome-wide SNP-heritability (e.g. estimated by
LDSC) for population 1 multiplied by sample size of the GWAS in population 1.

* `--sigmasq2` specifies genome-wide SNP-heritability (e.g. estimated by
LDSC) for population 2 multiplied by sample size of the GWAS in population 2.

* `--totnsnp` specifies total number of SNPs across all chromosomes.

* `--max_iter` specifies maximum number of EM iterations. The default is 200.

* `--out` specifies the output file name.
