# Estimating posterior probabilities of each SNP being population-specific and shared

This page describes the pipeline to estimate the posterior probabilities of a
SNP being population-specific and shared.

### Step 0: estimate genome-wide prior probabilities

This step should be performend following steps described
[here](https://huwenboshi.github.io/pesca/prior/).

### Step 1: estimate per-SNP posterior probabilities

PESCA estimates per-SNP posterior probabilities of each SNP beling
population-specific and shared via the following command. This step could be
run for each chromosome in parallel.

```text
<directory to pesca>/pesca \
    --mode post \
    --f00 <MVB f00 parameter from step 0> \
    --f01 <MVB f00 parameter from step 0> \
    --f10 <MVB f00 parameter from step 0> \
    --f11 <MVB f00 parameter from step 0> \
    --zscore1 <list of GWAS summary statistics files for population 1> \
    --zscore2 <list of GWAS summary statistics files for population 2> \
    --ld1 <list of LD files for population 1> \
    --ld2 <list of LD files for population 2> \
    --nburn <number of MCMC burn-ins, 5000 by default> \
    --nsample <number of MCMC samples, 5000 by default> \
    --lambda <shrinkage parameter> \
    --sigmasq1 <sample size times heritability in population 1> \
    --sigmasq2 <sample size times heritability in population 1> \
    --totnsnp <total number of SNPs across all chromosomes> \
    --max_iter <number of independent MCMC chains> \
    --out <output file name>_chr<chromosome #>
```

Here are the meaning of the flags:

* `--mode` tells PESCA whether to estimate genome-wide prior or per-SNP
posterior probabilities. There are two options for this flag `fit` and `post`.
`post` tells PESCA to estimate posteriors. For estimating priors see 
[here](https://huwenboshi.github.io/pesca/priors/).

* `f00`, `f01`, `f10`, `f11` are multivariate Bernoulli (MVB) model parameters.
This estimates should be obtained from step 0.

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

<div style="background-color:rgba(230, 230, 250, 1.0);">
( <b>Note</b>: This is the number of SNPs after intersecting with the
reference panels, which have been MAF filtered and LD pruned. )
</div>

* `--max_iter` specifies the number of independent MCMC chains to estimate
posteriors. PESCA uses the average posteriors as the final posteriors. We
recommend setting this value to 20.

* `--out` specifies the output file name.

The following is an example output, `<output file name>_chr<chromosome #>.txt`.

```text
SNP          BP          A0_GWAS1   A1_GWAS1    Z_GWAS1     A0_GWAS2  A1_GWAS2   Z_GWAS2     P_GWAS1_ONLY    P_GWAS2_ONLY    P_BOTH
rs12567000   101242593   C          T           -0.669797   C         T          -1.20845    0.00013         0               0.0736
rs12404364   101255697   T          C           1.40208     T         C          1.45765     0               0.00016         0.10239
rs12043207   101302813   A          G           -0.720678   A         G          0.356864    0               0               0.05533
rs7541397    101322866   G          A           -2.03966    G         A          -3.05884    0.01002         0.00356         0.4613
rs3737581    101341083   G          A           -0.879881   G         A          0.411767    0.0005          0               0.08021
...          ...         ...        ...         ...         ...       ...        ...         ...             ...             ...
```

* `SNP` is the ID of the SNP.

* `BP` is the base pair position of the SNP.

* `A0_GWAS1` is the non-effect allele of GWAS in population 1.

* `A1_GWAS1` is the effect allele of GWAS in population 1.

* `Z_GWAS1` is the Z-score of GWAS in population 1.

* `A0_GWAS2` is the non-effect allele of GWAS in population 2.

* `A1_GWAS2` is the effect allele of GWAS in population 2.

* `Z_GWAS2` is the Z-score of GWAS in population 2.

* `P_GWAS1_ONLY` is the posterior probability that the SNP is causal only in
population 1.

* `P_GWAS2_ONLY` is the posterior probability that the SNP is causal only in
population 2.

* `P_BOTH` is the posterior probability that the SNP is causal in both
populations.
