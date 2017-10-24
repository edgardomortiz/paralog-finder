#!/usr/bin/env Rscript

# Script based on McKinney et al. 2017's HDplot.R at:
# http://datadryad.org/bitstream/handle/10255/dryad.123458/HDplot.R?sequence=1
# Re-adapted to parse arguments and to create additional graphs for
# Number of samples in locus (N) vs Heterozigosity (H) or Read ratio deviation (D)


# Author: Edgardo M. Ortiz
# Date: 2017-10-21


library(ggplot2)
library(ggExtra)
suppressMessages(library(argparse))

# Create parser object
parser <- ArgumentParser()
parser$add_argument("-i", "--input", type="character", required=TRUE,
                    help="Name of the file with extension .depthsBias")
parser$add_argument("--minD", type="integer", default=0,
                    help="Minimum D to display in graph, default = min from data")
parser$add_argument("--maxD", type="integer", default=0,
                    help="Maximum D to display in graph, default = max from data")
parser$add_argument("--transp", type="double", default=0.01, 
                    help="Transparency for points in graph, default = 0.01")
args <- parser$parse_args() 

# HDplot
HDplotData<-read.delim(args$input,header=TRUE)
cat(nrow(HDplotData), "SNPs loaded from", args$input, "\n")

# D vs Heterozigosity graphs
if ( args$minD != 0 | args$maxD != 0 ) {
  minD <- floor(args$minD)
  maxD <- ceiling(args$maxD)
  if (maxD-minD >= 20) {
    ticks_D <- round((maxD-minD)/20, digits=0)
  } else {ticks_D <- round((maxD-minD)/20, digits=2)
  }
  HDplotData <- subset(HDplotData, z>=args$minD & z<=args$maxD)
} else {
  minD <- floor(min(HDplotData$z, na.rm=TRUE))
  maxD <- ceiling(max(HDplotData$z, na.rm=TRUE))
  if (maxD-minD >= 20) {
    ticks_D <- round((maxD-minD)/20, digits=0)
  } else {ticks_D <- round((maxD-minD)/20, digits=2)
  }
}

minH <- 0
maxH <- round(max(HDplotData$hetPerc, na.rm=TRUE), 1)
if (maxH >= 0.5) {
  ticks_H <- round((maxH-minH)/20, digits=2)
} else {
  ticks_H <- round((maxH-minH)/10, digits=4)
}

minN <- min(HDplotData$num_samples, na.rm=TRUE)
maxN <- max(HDplotData$num_samples, na.rm=TRUE)
if (maxN-minN == 1) {
  seq_N <- c(minN,maxN)
} else {
  seq_N <- seq(minN,maxN,round((maxN-minN)/20, digits=0))
}

if (minD != maxD) {
  # Read ratio deviation (D) vs. heterozigosity
  DvH <- ggplot(HDplotData, aes(hetPerc,z)) + geom_point(alpha=args$transp, size=0.8) + theme_minimal() +
    scale_x_continuous(name="Proportion of heterozygotes (H)", breaks=seq(minH,maxH,ticks_H)) +
    scale_y_continuous(name="Read ratio deviation (D)", breaks=seq(minD,maxD,ticks_D)) 
  DvH <- ggMarginal(DvH, type="boxplot", size=20, fill='#FFFFFF', outlier.size=0.8)
  ggsave("graph1_D_vs_Het.png", plot=DvH, width = 8, height = 8)
  message("graph1_D_vs_Het.png saved")
} else {
  cat("There is no variation in read radio deviation, graph 1 won't be plotted. D=", maxD)
}

if (minN != maxN) {
  # Number of samples in locus vs. heterozigosity
  NvH <- ggplot(HDplotData, aes(hetPerc,num_samples)) + geom_point(alpha=args$transp, size=0.8) + theme_minimal() +
    scale_x_continuous(name="Proportion of heterozygotes (H)", breaks=seq(minH,maxH,ticks_H)) +
    scale_y_continuous(name="Number of samples in locus", breaks=seq_N)
  NvH <- ggMarginal(NvH, type="boxplot", size=20, fill='#FFFFFF', outlier.size=0.8)
  ggsave("graph2_NumSamplesInLocus_vs_Het.png", plot=NvH, width = 8, height = 8)
  message("graph2_NumSamplesInLocus_vs_Het.png saved")

  # Read ratio deviation (D) vs. number of samples in locus
  DvN <- ggplot(HDplotData, aes(num_samples,z)) + geom_point(alpha=args$transp, size=0.8) + theme_minimal() +
    scale_x_continuous(name="Number of samples in locus", breaks=seq_N) +
    scale_y_continuous(name="Read ratio deviation (D)", breaks=seq(minD,maxD,ticks_D))
  DvN <- ggMarginal(DvN, type="boxplot", size=20, fill='#FFFFFF', outlier.size=0.8)
  ggsave("graph3_D_vs_NumSamplesInLocus.png", plot=DvN, width = 8, height = 8)
  message("graph3_D_vs_NumSamplesInLocus.png saved")
} else {
  cat("There is no variation in locus depth, graphs 3 and 2 won't be plotted.")
}

# Allelic ratio vs. heterozigosity
RvH <- ggplot(HDplotData, aes(hetPerc,ratio)) + geom_point(alpha=args$transp, size=0.8) + theme_minimal() +
  scale_x_continuous(name="Proportion of heterozygotes (H)", breaks=seq(minH,maxH,ticks_H)) +
  scale_y_continuous(name="Allelic ratio (a:b)", breaks=c(0.25, 0.5, 0.75), labels=c("1:3", "1:1", "3:1"))
RvH <- ggMarginal(RvH, type="boxplot", size=20, fill='#FFFFFF', outlier.size=0.8)
ggsave("graph4_AlleleRatio_vs_Het.png", plot=RvH, width = 8, height = 8)
message("graph4_AlleleRatio_vs_Het.png saved")

