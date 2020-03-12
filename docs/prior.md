# Estimating genome-wide proportion of population-specific and shared causal variants

This page describes the pipeline to estimate the genome-wide prior
probabilities of a SNP to be causal in one population (i.e.
population-specific) or in both populations (i.e. shared).

### Recommended pipeline

#### Step 0: estimate genome-wide heritability of the trait

This step can be performend using [LDSC](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation).

<div style="background-color:rgba(230, 230, 250, 1.0);">
( <b>Note</b>: This step should be performed using all SNPs without any
MAF filtering or LD pruning, except for those required by LDSC. )
</div>

#### Step 1: extract summary statistics and LD of each region for each chromosome

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

#### Step 2: run PESCA for each chromosome

PESCA estimates genome-wide proportions of population-specific and shared
causal variants via the following command. This step could be run for each
chromosome in parallel.

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

<div style="background-color:rgba(230, 230, 250, 1.0);">
( <b>Note</b>: This is the number of SNPs after intersecting with the
reference panels, which have been MAF filtered and LD pruned. )
</div>

* `--sigmasq1` specifies genome-wide SNP-heritability (e.g. estimated by
LDSC) for population 1 multiplied by sample size of the GWAS in population 1.

* `--sigmasq2` specifies genome-wide SNP-heritability (e.g. estimated by
LDSC) for population 2 multiplied by sample size of the GWAS in population 2.

* `--totnsnp` specifies total number of SNPs across all chromosomes.

* `--max_iter` specifies maximum number of EM iterations. The default is 200.

* `--out` specifies the output file name.

The following is an example output, `<output file name>.log`.
```text
iter    nsnp    q00     q01     q10     q11      f00     f01       f10       f11
0       9960    9389.8  190.067 190.067 190.067  0       -3.9      -3.9      3.9
1       9960    9433.58 177.374 184.488 164.561  0       -3.97377  -3.93445  3.85947
2       9960    9466.98 164.664 187.16  141.2    0       -4.05166  -3.9236   3.76987
3       9960    9478.17 164.541 173.529 143.756  0       -4.05359  -4.0004   3.86536
4       9960    9502.39 164.518 157.391 135.704  0       -4.05628  -4.10056  3.90802
5       9960    9524.52 154.285 154.576 126.616  0       -4.12283  -4.12094  3.9233
...     ...     ...     ...     ...     ...      ...     ...       ...       ...
```
* `iter` is the iteration number of the EM algorithm.

* `nsnp` is the total number of SNPs.

* `q00` is the number of SNPs that are causal in neither population.

* `q01` is the number of SNPs that are causal in population 1.

* `q10` is the number of SNPs that are causal in population 2.

* `q11` is the number of SNPs that are causal in both populations.

* `f00`, `f01`, and `f10`, and `f11` are parameters of the multivariate Bernoulli
(MVB) model.

We use results from the last EM iteration as the final estimates.

#### Step 3: aggregate results from all chromosomes

After step 2, there should be one `<output file name>.log` for each
chromosome.


