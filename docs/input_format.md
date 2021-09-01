# Input Format

This page describes file formats of GWAS summary statistics data and LD matrix
used by PESCA.

## GWAS summary statistics data

GWAS summary statistics data should have the following columns, ordered as
listed below.

* `SNP`: rs ID of the SNP (e.g. rs62442)
* `BP`: Base pair position of the SNP. 
* `A1`: Effect allele (i.e. trait increasing allele) of the SNP. The sign of
Z-score is with respect to this allele.
* `A2`: The other allele of the SNP.
* `Z`: Z-score of the SNP.
* `N`: Sample size.

The following is an example of GWAS summary statistics file.
```text
SNP         BP         A2      A1      Z           N 
rs131529    16870162   A       G       0.985754    50000
rs131538    16871137   G       A       1.382624    50000
rs131545    16871951   G       A       1.304881    50000
rs140375    16873251   C       A       1.494041    50000
rs131566    16879530   G       A       0.981185    50000
...         ...        ...     ...     ...         ...
```

There should be a GWAS summary statistics file for each non-overlapping
region.

## GWAS summary statistics data list file

GWAS summary statistics data list file contains a list of paths to GWAS
summary statistics data split by region. Typically, there should be one such
list file for each chromosome and each population. A lsit of GWAS summary
statistics files can be specified to PESCA through the `--zscore1` and
`--zscore2` (one flag for each population) flags.

The following is an example of GWAS summary statistics data list file.

```text
/<path>/<to>/<GWAS>/<summary>/<data>/region1_zscore_pop1.txt
/<path>/<to>/<GWAS>/<summary>/<data>/region2_zscore_pop1.txt
                            ...
/<path>/<to>/<GWAS>/<summary>/<data>/region100_zscore_pop1.txt
                            ...
```

## LD matrix

The LD matrix file is a text file that stores the LD matrix of a region. Each
LD matrix file should correspond to a GWAS summary statistics file of a
region. There should be M lines and M columns in the LD matrix file, where M
is the number of SNPs in the region.

The following is an example of a LD matrix file.
```text
1 0.964359 -0.0862651 -0.0846509 -0.0339536 ...
0.964359 1 -0.0655564 -0.0446302 -0.0588997 ...
-0.0862651 -0.0655564 1 0.859188 0.0562084  ...
-0.0846509 -0.0446302 0.859188 1 0.048621   ...
-0.0339536 -0.0588997 0.0562084 0.048621 1  ...
...        ...        ...       ...      ...
```

## LD matrix list file

LD matrix list file contains a list of paths to LD matrices split by region.
Typically, there should be one such list file for each chromosome and each
population. LD matrix list files can be specified to PESCA through the
`--ld1` and `--ld2` (one flag for each population) flags.

The following is an example of LD matrixlist file.

```text
/<path>/<to>/<GWAS>/<summary>/<data>/region1_ld_pop1.txt
/<path>/<to>/<GWAS>/<summary>/<data>/region2_ld_pop1.txt
                            ...
/<path>/<to>/<GWAS>/<summary>/<data>/region100_ld_pop1.txt
                            ...
```

LD matrix list files can be specified to PESCA through the `--ld1` and
`--ld2` flags (one flag for each population).

<div style="background-color:rgba(230, 230, 250, 1.0);">
( <b>Note</b>: Each line in the list of LD matrix file should correpond to a
line in the list of GWAS summary statistics file. )
</div>
