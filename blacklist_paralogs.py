#!/usr/bin/env python2
# -*- coding: utf-8 -*-


'''
After running HDplot (McKinney et al. 2017, Molecular Ecology Resources) to
identify paralogs by determining the "read depth ratio (D)" and "percentage
of heterozygotes" we need to parse the output and extract the list of loci
that must be excluded from the Stack analysis (i.e., create a blacklist)
'''


__author__      = "Edgardo M. Ortiz"
__credits__     = "Juan D. Palacio-Mejía"
__version__     = "1.0"
__email__       = "e.ortiz.v@gmail.com"
__date__        = "2017-10-15"


import sys
import argparse


def main():
	parser = argparse.ArgumentParser(description="Creates blacklist of paralog loci (and whitelist of singletons) based on the method McKinney et al. 2017 (doi: 10.1111/1755-0998.12613)")
	parser.add_argument("-i", "--input", action="store", dest="filename", required=True,
		help="Name of the input file, the table with extension .depthsBias produced by HDplot.py")
	parser.add_argument("--maxH", action="store", dest="max_hetPerc", type=float, default=0.6,
		help="Maximum proportion of heterozygotes in a locus, default=0.6 taken from McKinney et al. 2016")
	parser.add_argument("--minN", action="store", dest="min_num_samples", type=int, default=1,
		help="Minimum number of samples in locus, default=1")
	parser.add_argument("--minD", action="store", dest="min_z", type=int, default=-7,
		help="Lower limit of read ratio deviation (D), default=-7 taken from McKinney et al. 2016")
	parser.add_argument("--maxD", action="store", dest="max_z", type=int, default=7,
		help="Upper limit of read ratio deviation (D), default=7 taken from McKinney et al. 2016")
	args = parser.parse_args()

	filename = args.filename
	max_hetPerc = args.max_hetPerc
	min_z = args.min_z
	max_z = args.max_z
	min_num_samples = args.min_num_samples

	# Blacklist filename will be the same as input file, replacing extension by "_paralogs.blacklist"
	blacklist = filename.split(".")[0]+"_paralogs.blacklist"

	# Whitelist filename will be the same as input file, replacing extension by "_singletons.whitelist"
	whitelist = filename.split(".")[0]+"_singletons.whitelist"

	paralogs = []
	singletons = []

	print "Retaining loci with at least "+str(min_num_samples)+" sample and with proportion of heterozygotes ≤ "+str(max_hetPerc)+" and D between "+str(min_z)+" and "+str(max_z)

	# Process file
	with open(filename) as dbias:
		for line in dbias:

			# Skip first line
			if not line.startswith("\t"):

				# Split line in tabs
				tabs = line.strip("\n").split("\t")

				# Check minimum number of samples in locus
				if tabs[8] >= min_num_samples:

					# Check maximum percentage of heterozygotes
					if float(tabs[11]) <= max_hetPerc:

						# Check D is within range:
						if float(tabs[13]) >= min_z and float(tabs[13]) <= max_z:

							# If true, append locus name to singletons list
							singletons.append(tabs[3].split("_")[0])

						# If D not within range append to list of paralogs
						else:
							paralogs.append(tabs[3].split("_")[0])

					# If percentage of heterozygotes above threshold append to paralogs list
					else:
						paralogs.append(tabs[3].split("_")[0])

				# If minimum number of sample in locus is not reached
				else:
					paralogs.append(tabs[3].split("_")[0])


	# Eliminate duplicate loci names from both lists, some SNPs may send loci to
	# both lists, therefore we have to be sure that names don't repeat between them
	paralogs = set(paralogs)
	print str(len(paralogs))+" loci written to blacklist of paralogs"
	singletons = list(set(singletons) - set(paralogs)) # If locus in list of paralogs eliminate it from list of singletons
	print str(len(singletons))+" loci written to whitelist of singletons"

	# Write both lists
	out_blacklist = open(blacklist, "w")
	out_whitelist = open(whitelist, "w")

	out_blacklist.write("\n".join(paralogs))
	out_blacklist.close()

	out_whitelist.write("\n".join(singletons))
	out_whitelist.close()


if __name__ == "__main__":
    main()


