# Estimating genome-wide proportion of population-specific and shared causal variants

This page describes the pipeline to estimate the genome-wide prior
probabilities of a SNP to be causal in one population (i.e.
population-specific) or in both populations (i.e. shared).

### Recommended pipeline

#### Step 1 -- estimate heritability of the trait

This step can be performend using [LDSC](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation).

#### Step 2 -- extract summary statistics and LD of each region

#### Step 3 -- run PESCA

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
