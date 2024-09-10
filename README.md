# Gentianella-calcis
Presented here are scripts used in data processing and analysis for "Population genomic data from threatened New Zealand *Gentianella calcis* (Gentianaceae) challenges conservation prioritisation."

## *Stacks* parameter trials without replicates
As per the suggestion of [Paris et al. 2017](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12775) we performed intial parameter trials in Stacks v2.64 using a representative subset of 12 samples.

## *Stacks* parameter trials with replicates
In our parameter trials, the number of loci continued to increase up to *M* 6 in our trials, necessitating an alternative approach to optimise parameters. We used the method of [Mastretta-Yanes et al. 2014](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12291) and performed additional trials using available replicate pairs of the previously used subset.

## Estimation of SNP error
We estimated SNP error (defined as the number of differences in SNPs between sample replicate pairs divided by the total number of SNPs) using modified code from Mastretta-Yanes et al. to compare alongside number of loci produced for each parameter combination in *Stacks* to see which produced the greatest number of *-R* 80 loci while also reducing SNP error.

## Identify paralogs
As paralogs may represent an additional source of error in RAD datasets, we used *HDplot* to SNPs with the expected proportion of heterozygous individuals (H) greater than 5 and a “read-ratio deviation” (D) less than -4 or greater than 4 as putative paralogs. Any RAD locus containing one of these SNPs was exlcuded from further analysis.

## Population genetic analyses in R
We used various R packages including *adegenet* v2.1.10, *ggplot2*, *StAMMP* v1.6.3, *ade4* v.1.7-20, *dartR* v2.9.7, and *vcfR* v1.13.0 to import, manipulate, analyise and visualise VCF files produced from *Stacks*.

## Putative outlier detection and visualisation
We used *pcadapt* v4.4.3 to identify putatively neutral and outlier SNPs in a version of our dataset. These were separated into their own genlight objects using *dartR* and separate analysis of each undertaken to compare their patterns.
