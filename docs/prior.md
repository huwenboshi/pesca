# Estimating genome-wide proportion of population-specific and shared causal variants

This page describes the pipeline to estimate the genome-wide prior
probabilities of a SNP to be causal in one population (i.e.
population-specific) or in both populations (i.e. shared).

### Typical command

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
    --totnsnp <total number of SNPs in the list of GWAS summary statistics files> \
    --max_iter <maximum number of EM iterations> \
    --out <output file name>
```

Here are the meaning of the flags:

* `--model` tells PESCA whether to estimate genome-wide prior or per-SNP
posterior probabilities. There are two options for this flag `fit` and `post`.
For estimating per-SNP posterior see
[here](https://huwenboshi.github.io/pesca/posterior/).
