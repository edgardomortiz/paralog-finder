# paralog-finder
Detects and blacklists paralog RAD loci analyzed in Stacks, based on the McKinney 2017 method (doi:10.1111/1755-0998.12613)

## Description
These scripts allow the identification of paralog RAD loci based on the method of McKinney 2017 (doi:10.1111/1755-0998.12613). However, we introduced some modifications and also made some additions. The main modification we made is the way in which the percentage of heterozygote individuals is calculated in datasets with varying degrees of missing data per locus. McKinney et al. (2017) divide the number of heterozygote individuals in the locus over the total number of individuals in the dataset, their results were not skewed because they either simulated data (0% missing data) or use very stringent filters in empirical data (at most 10% missing data per locus). Our calculations consider instead,the number of individuals present at each locus, since typical datasets usually have many loci represented by only a few individuals.

As an addition to their method, we also provide a script that creates a blacklist of paralog loci that can be used in `Stacks` to run the `populations` module and exclude the paralogs from the calculations and the output matrices for downstream analyses in other software.

Briefly, the scripts must be run in three steps:
- **Step 1:** Calculate loci statistics from the VCF file produced by Stacks with the script `HDplot_process_vcf.py`
- **Step 2:** Plot the statistics and decide limits for maximum heterozigosity (H), maximum and minimum read ratio deviation (D), and minimum of samples per locus with the R script `HDplot_graphs.R`
- **Step 3:** Provide the parameters decided by the inspection of the graphs to the script `blacklist_paralogs.py` to create the blacklist of paralog loci and a whitelist of singleton loci.

## Requirements:
- Python 2.7 with the packages: `scipy`, `pandas`, `numpy`, and `statsmodels`, these are easy to install if you use [`miniconda`](https://conda.io/miniconda.html).
- R with the CRAN packages: `ggplot2`, `ggExtra`, and `argparse`.
