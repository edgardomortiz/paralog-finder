#!/usr/bin/env python2
# -*- coding: utf-8 -*-


'''
Script based on McKinney et al. 2017's HDplot_python.py at:
http://datadryad.org/bitstream/handle/10255/dryad.123460/HDplot_python.py?sequence=1

Our modification calculates percentage of heterozygotes per locus
dividing the number of of heterozygote individuals over the number
of individuals present at the locus (unlike McKinney's which always
divides over the total number of individuals in the dataset, which
didn't skew their results because they had at most 10% missing data
per locus).

We also integrate the separate python script called vcf_to_depth.py at:
http://datadryad.org/bitstream/handle/10255/dryad.123461/vcf_to_depth.py?sequence=1
into the main script.
'''


__author__      = "Edgardo M. Ortiz"
__credits__     = "Juan D. Palacio-Mejía"
__version__     = "1.0"
__email__       = "e.ortiz.v@gmail.com"
__date__        = "2017-10-21"


import scipy.stats
import pandas as pd
import numpy as np
import statsmodels
import os
import gzip
import argparse


def vcf_to_allele_depth(vcf_file, out_file):
    if vcf_file.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open
    with opener(vcf_file) as INFILE:
        with open(out_file, 'w') as OUTFILE:
            header_lines = 0
            for line in INFILE:
                if line.startswith('##'):
                    header_lines += 1
                elif line.startswith('#CHROM'):
                    print 'skipped {} header lines'.format(header_lines)
                    header = line.strip().strip('#').split('\t')
                    inds = header[9:]
                    print 'found {} individuals'.format(len(inds))
                else:
                    tabs = line.split('\t')
                    contig = tabs[0]
                    pos = tabs[1]
                    locus_ID = tabs[2]
                    genotypes = tabs[9:]
                    depth_a_of_ind = dict()
                    depth_b_of_ind = dict()
                    for gen_idx, gen in enumerate(genotypes):
                        if gen.split(':')[0] in ['1/0', '0/1']: # if het
                            depth_a_of_ind[inds[gen_idx]] = int(gen.split(':')[2].split(',')[0])
                            depth_b_of_ind[inds[gen_idx]] = int(gen.split(':')[2].split(',')[1])
                    sum_a = sum(depth_a_of_ind.values())
                    sum_b = sum(depth_b_of_ind.values())
                    num_hets = len(depth_b_of_ind.values())

                    # Our modification: get actual number of samples in the locus
                    num_samples = int(tabs[7].split(";")[0].replace("NS=",""))
                    
                    if sum_a+sum_b > 0:
                        OUTFILE.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(contig, pos, locus_ID, sum_a, sum_b, np.float(sum_a)/(sum_a+sum_b), num_hets, num_samples))


def main():
	parser = argparse.ArgumentParser(description="Processes a VCF input file produced by Stacks for plotting heterozigosity and read depth deviations using the method of McKinney et al. 2017 (doi: 10.1111/1755-0998.12613)")
	parser.add_argument("-i", "--input", action="store", dest="filename", required=True,
		help="Name of VCF input file, must have read depth per allele in each individual (Stacks format)")
	args = parser.parse_args()

	filename = args.filename
	depths_file = args.filename.split(".")[0]+".depths"

	# LOAD VCF FILE AND EXTRACT SEQUENCE COUNTS FROM HETEROZYGOUS INDIVIDUALS
	vcf_to_allele_depth(vcf_file=filename, out_file = depths_file)

	# LOAD COUNT DATA INTO DATAFRAME
	depths = pd.read_csv(depths_file, sep = '\t', header = None)
	depths.columns = ['contig', 'pos', 'locus_ID', 'depth_a' , 'depth_b', 'ratio', 'num_hets', 'num_samples']
	depths.head()

	#SUM READ COUNTS PER LOCUS
	depths['total_depth'] = depths['depth_a'] + depths['depth_b']
	depths['depth_per_het'] = depths['total_depth']/[np.float(xx) for xx in depths['num_hets']]
	depths.head()

	# CALCULATE HETEROZYGOSITY
	depths['hetPerc']=depths['num_hets']/depths['num_samples']
	depths.head()

	# CALCULATE EXPECTED STANDARD DEVIATION BASED ON BINOMIAL DISTRIBUTION
	depths['std'] = scipy.stats.binom(n = depths['total_depth'], p = .5).std()
	depths.head()

	# CALCULATE Z-SCORE BASED ON STANDARD DEVIATION
	depths['z'] = -(depths['total_depth']/2. - depths['depth_a'])/ depths['std']
	depths.head()

	# WRITE OUTPUT FILE CONTAINING DEPTH AND BIAS INFORMATION
	depths.to_csv(depths_file+"Bias", sep="\t")

	# REMOVE INTERMEDIATE FILE
	os.remove(depths_file)


if __name__ == "__main__":
    main()
