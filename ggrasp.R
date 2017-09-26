#!/usr/bin/env Rscript

###############################################################################
#                                                                             #
#       Copyright (c) 2017 J. Craig Venter Institute.                         #
#       All rights reserved.                                                  #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.    #
#                                                                             #
###############################################################################
###############################################################################

library(getopt); #reading in parameters
library(ggrasp);
library(methods);

#copied from sample_medioids.R
params=c(
  "input", "i", 1, "character",
  "type", "t", 2, "character",
  "output", "o", 1, "character",
  "method", "m", 2, "character",
  "offset", "d", 2, "numeric", #optional conversion of similarity to distance matrix using offset provided
  "verbose", "v", 2, "character",
  "layout", "l", 2, "character",
  "names", "n", 2, "character",
  "clusters", "c", 2, "character",
  "ranks", "r", 2, "character",
  "threshold", "h", 2, "character",
	"keep", "k", 2, "numeric",
	"zscore", "z", 2, "character",
	"minlambda", "x", 2, "character",
	"metainfo", "a", 2, "character",
	"start","s", 2, "numeric",
	"end","e", 2, "numeric",
	"package","p", 2, "character",
	"writetable", "0", "2", "numeric",
	"writepseudo", "1", "2", "numeric",
	"writetree", "2", "2", "numeric",
	"writetrimtree", "3", "2", "numeric",
	"writeitol", "4", "2", "numeric",
	"plothist", "5", "2", "numeric",
	"plotgmm", "6", "2", "numeric",
	"plottree", "7", "2", "numeric",
	"plotclust", "8", "2", "numeric",
	"plottrim", "9", "2", "numeric",
	"plottype", "y", "2", "character"
	
)


opt=getopt(spec=matrix(params, ncol=4, byrow=TRUE))

script_name=unlist(strsplit(commandArgs(FALSE)[4],"=")[1])[2]

usage = paste (
  "\nUsage:\n\n", script_name, "\n\n",
  "\n# Input Variables\n",
  "\t-i, --input <input file name for tab delimited distance matrix with row and column headers, newick file, or aligned multiple fasta file>\n",
  "\t-t, --type <type of the tree format to makes. Default is 'complete'>",
  "\t-d, --offset <optional offset to convert a similarity matrix to a distance matrix>\n",
  "\t-m, --method  <optional method to use for hclust such as complete, single or average or nj to perform neighbor-joining>\n",
  "\t-v, --verbose <optional verbal threshold: 1 = no output, 2 = text only, 3 = graphics in pdf. Default = 3>\n",
  "\t-n, --names <optional name file to label the exterior of the tree. Multiple files can be included with a comma between>\n",
  "\t-r, --ranks <optional input file with ranks of the genes for medoid selection. Ranks should be from 1 (best) to n (worst) in two columns with the genome ID in the first>\n",
  "\n# Clustering Variables\n",
  "\t-c, --clusters <optional number of clusters to generate>\n",
  "\t-h, --threshold <optional threshold to use to cluster the genomes>\n",
  "\t-k, --keep <optional flag to keep all reference clusters as sub-and show subclusters>\n",
  "\t-z, --zscore <optional value to combine all gaussian mixtures within this z-score. Default = 1. Use 0 to turn off>\n",
  "\t-a, --metadata <optional input file containing metadata for the genomes. It should be tab-deleniated>\n",
  "\t-x, --minlambda <optional value used to discard a Gaussian distribution if it has less than this percentage of the total distribution. Default is 0.005>\n",
  "\t-s, --start <optional value with the minimum number of gmm to examine. Default is 2>\n",
  "\t-e, --end <optional value with the maximum number of gmm to examine. Default is 10>\n",
  "\t-p, --packages <optional value with the package to use for the Gaussian mixture model estimate. Default is bgmm. \n",	
  "\n# Output Variables\n",
  "\t-o, --output <output file name for html table, pdf plots, cluster and txt files>\n",

  "\t-l, --layout <optional phylogeny layout such as rectangular, circular, radial and slanted. Default = circular>\n", 
 "\t-0, --writetable <optional flag used to save the cluster table>\n",	
  "\t-1, --writepseudo <optional flag used to save the psuedo fasta lsit>\n",	
  "\t-2, --writetree <optional flag used to save the cluster table>\n",	
  "\t-3, --writetree <optional flag used to save the cluster table>\n",	
  "\t-4, --writetree <optional flag used to save the cluster table>\n",	 
  "\t-5, --plothist <optional flag used to save the histogram with the plot>\n",	
  "\t-6, --plotgmm <optional flag used to save the histogram with the Gaussian distribution overlayed>\n",	
  "\t-7, --plottree <optional flag used to save the full tree plot>\n",	
  "\t-8, --plotclust <optional flag used to save the full tree plot with the clusters highlighted>\n",	
  "\t-9, --plottrim <optional flag used to save the trimmed tree with only the medoids>\n",	
  "\t-y, --plottype <optional flag used to provide the plot format>\n",	

  "This script will select sample genome medoids by selecting a treshhold using Gaussian Mixture Models.\n",
  "\n",
  "\n");

is.binary = is.binary.tree;
if(!length(opt$input))
{
  
  cat(usage)
  cat(length(opt))
  q(status=-1)
}
if (!length(opt$output))
{
	cat("No output file name inputed. Using default \"output\"\n\n");
	OutputFileName = "Output";
}else
{
	OutputFileName <- opt$output;
}

if (!is.null(opt$metainfo))
{
	if (file.exist(opt$metainfo))
	{
		meta.info <- read.table(opt$metainfo, sep="\t", header=T, stringsAsFactors=F, row.names = 1);
	}else
	{ 
		cat(paste("\nCannot find metainfo file ", opt$metainfo, ". Skipping this step\n\n", sep=""));
	}
}


medioid_run_str = "";


dist_con <- opt$offset;
in.file = opt$input;
in.format = "";
if (!(is.null(opt$type)))
{
	in.format = opt$type;
	if (!(tolower(opt$type) %in% c("matrix", "tree", "fasta")))
	{
		cat(paste("Format ", in.format, " is not recognized. Allowing program to select...\n", sep=""));
		in.format = "";
	}
}


verbosity = 3;
verb_lev_str = c("No text or PDF", "Text Only", "Both Text and PDF Graphics")

z.limit.input = 1.0;
min.lambda = 0.005;
if (!is.null(opt$z.remove))
{
	z.limit.input = as.numeric(opt$zscore);
}
if (!is.null(opt$min.lambda))
{
	min.lambda = as.numeric(opt$minlambda);
}

if (!is.null(opt$verbosity)) {
  if (is.integer(as.integer(opt$verbosity)))
  {
    verbosity = as.integer(opt$verbosity);
    if (verbosity >0 & verbosity < 4)
    {
      cat("\nSetting output verbosity to", verb_lev_str[verbosity], "\n\n")
    }else
    {
      varbosity = 3;
      cat("\nIncorrected Verbosity setting (1,2,3 accepted.) Output verbosity remains at", verb_lev_str[verbosity], "\n\n")
    }
  }else
  {
    cat("\nIncorrected Verbosity setting (1,2,3 accepted.) Output verbosity remains at", verb_lev_str[verbosity], "\n\n")
    
  }
}



offset.val = 0
if (!is.null(opt$offset)) {
  if (verbosity > 1) {
    cat("\nConverting similarity matrix to distance matrix using threshold: ", opt$offset, "\n\n"); }
 offset.val = opt$offset
}



accept_hclust <- c("complete", "average", "single", "nj")
hclust_m = "complete"
if (!is.null(opt$method))
{
	if (opt$method %in% accept_hclust)
    {
      hclust_m <- opt$method;
      if (verbosity > 1)
      {
        cat("\nUsing", hclust_m,"to make phylogeny\n\n");
      }
    }else
    {
      cat("\nSuggested unacceptible phylogenetic method...\n\nQuiting..\n");
      q()
    }
  }
RankNameFile = "";
if (!is.null(opt$ranks))
{
  RankNameFile = opt$ranks;
}
gg.1 <- ggrasp.load(file=in.file, file.format=in.format, offset=as.numeric(offset.val), tree.method=hclust_m, rank.file=RankNameFile)
if(is.null(opt$input_newick_file))
{
	write(gg.1@phy, file=paste(OutputFileName, "tree", "nwk", sep="."))
}

layout = "circular";
if (!is.null(opt$layout)) {
  layout = opt$layout;
  if (verbosity > 1) {
    cat("\nSetting phylogeny plot layout to: ", layout, "\n\n"); }
}

if (!is.null(opt$names))
{
  name_file <- opt$names;
  name_file = InputNameFile;
  if (verbosity > 1)
  {
    if (grepl(",", name_file))
    {
      cat("\nAdding files ", name_file, " to label exterior of tree\n\n")
    }else
    {
      cat("\nAdding file ", name_file, " to label exterior of tree\n\n")
      
    }
  }
}

accept.package = c("bgmm","mixtools")
run.type.in = "bgmm"
if(!is.null(opt$package))
{
	if (!(opt$package %in% accept.package))
	{
		cat(paste(opt$package, " is not found. Default bgmm is being used...\n\n", sep=""));
		
	}else
	{
		run.type.in = opt$package;
	}
}
g.start = 2;
g.end = 10;
if (!is.null(opt$start))
{
	g.start = opt$start;
}

if (!is.null(opt$end))
{
	g.end = opt$end;
}
z.score = 1;
if (!is.null(opt$zscore))
{
	z.score = opt$zscore;
}
if(is.null(opt$cutoff_num) && is.null(opt$threshold))
{
		gg.2 <- ggrasp.cluster(gg.1, gmm.start=g.start, gmm.max=g.end, z.limit=z.score, run.type=run.type.in)

  if (verbosity > 1)
  {
    cat(paste("\n", as.character(length(gg.2@gmm$mu))," Gaussian Mixture Models identified...\n\n", sep=""));
 
    cat(paste("\nMedoid Threshold set at ", as.character(gg.2@h),"\n\nClustering the Genomes with threshold\n\n", sep="")); 
  }

  
}else
{
  if (is.null(opt$threshold))
  {	
	if (verbosity > 1)
	{
		cat(paste("Getting the top ", opt$clusters, " clades by walking down the tree..\n\n"));
	}
	clusters = as.numeric(opt$clusters)
	gg.2 <- ggrasp.cluster(gg.1, num.clusters = clusters)
    if (verbosity > 1)
	{
		cat(paste("\nMedoid Threshold set at ", as.character(gg.2@h),"\n\nClustering the Genomes with threshold\n\n", sep="")); 
	}
  }else
  {
	thrsh = as.numeric(opt$threshold);

	gg.2 <- ggrasp.cluster(gg.1, threshold = thrsh)
	if (verbosity > 1)
	{
		cat(paste(as.character(max(gg.2@clusters)), " clusters made using threshold of " , thrsh,"\n\n", sep="")); 
	}
  }
}

plottype = "pdf"
if (!is.null(opt$plottype))
{
	opt$plottype = tolower(opt$plottype);
	if (!(opt$plottype %in% c("pdf", "png")))
	{
		cat(paste(opt$plottype, " is not recognized. Default value is being used"));
	}
	else
	{
		plottype = opt$plottype;
	}
	cat(plottype)
}


if(!is.null(opt$writetable))
{
	
	file.name <- paste(OutputFileName, "table", "txt", sep=".");
	write(capture.output(print(gg.2, "table", opt$writetable)), file.name);
}

if(!is.null(opt$writetree))
{
	
	file.name <- paste(OutputFileName, "tree", "txt", sep=".");
	write(capture.output(print(gg.1, "tree", opt$writetree)), file.name);
}

if(!is.null(opt$writetrimtree))
{
	
	file.name <- paste(OutputFileName, "tree", "trim", "txt", sep=".");
	write(capture.output(print(gg.2, "tree", opt$writetree)), file.name);
}

if(!is.null(opt$writeitol))
{
	
	file.name <- paste(OutputFileName, "itol", "txt", sep=".");
	ggrasp.write(gg.2, "tree", file=file.name);
}

if(!is.null(opt$writepseudo))
{
	
	file.name <- paste(OutputFileName, "pfasta", "txt", sep=".");
	write(capture.output(print(gg.2, "list", opt$writetable)), file.name);
}

file.name <- paste(OutputFileName, "medoids", "txt", sep=".");
write(gg.2@medoids, file.name);

if(!is.null(opt$plothist))
{
	
	file.name <- paste(OutputFileName, "hist", plottype, sep=".");
	ggsave(plot(gg.2, "hist"), file = file.name, device=plottype);
}

if(!is.null(opt$plottree))
{
	
	file.name <- paste(OutputFileName, "tree", plottype, sep=".");
	ggsave(plot(gg.1, "tree"), file = file.name, device=plottype);
}

if(!is.null(opt$plotclust))
{
	
	file.name <- paste(OutputFileName, "tree", plottype, sep=".");
	ggsave(plot(gg.2, "tree"), file = file.name, device=plottype);
}

if(!is.null(opt$plotgmm))
{
	
	file.name <- paste(OutputFileName, "gmm", plottype, sep=".");
	ggsave(plot(gg.2, "gmm"), file = file.name, device=plottype);
}



if(!is.null(opt$plottrim))
{
	
	file.name <- paste(OutputFileName, "trimmed", plottype, sep=".");
	ggsave(plot(gg.2, "trimmed"), file = file.name, device=plottype);
}
if (verbosity > 1)
{
  cat(paste("Finished\n\n\n")); 
 }