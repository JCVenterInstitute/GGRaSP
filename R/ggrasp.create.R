

#' @import ggplot2
#' @import mixtools
#' @import ape
#' @import bgmm
#' @import colorspace
#' @import methods
#' @importFrom grDevices col2rgb rainbow
#' @importFrom graphics strwidth
#' @importFrom stats as.dist dnorm hclust pchisq
#' @importFrom utils read.delim read.table write.table
#' @title ggrasp.load
#' @description ggrasp.load() initializes a class GGRaSP object from a file containing either a tree, a distance matrix or a multi-fasta alignment. The returned object can subsequently be clustered using ggrasp.cluster().
#' 
#' 
#' @param file File containing the tree, matrix or sequence alignment used to initialize the ggrasp object. Required.
#' @param file.format The format the file is in, with tree, fasta and matrix accepted. If not given the program makes a guess.
#' @param rank.file File containing the ranks of genomes in a tab-delineated file with the genome in column 1 and the rank in column 2. The rank is a non-negative number.
#' @param offset Numeric representing a perfect match. Default is 0.
#' @param tree.method The method used to make the tree from a distance matrix. "Complete" (Default), "Average", "Single", and "nj" (Neighbor Joining) are currently available. 
#' 
#' @return Returns a class GGRaSP variable
#' @examples #The following data is from Chavda et al 2016 which phylotyped Enterobacter genomes
#' # Our example uses the data underpinning the tree shown in Figure 2
#' # Also included is a ranking file to prioritize closed Enterobactor genomes
#' 
#' library(ggrasp);
#' tree.file <- system.file("extdata", "Enter.kSNP.tree", package="ggrasp")
#' rank.file.in <- system.file("extdata", "Enter.kSNP.ranks", package="ggrasp")
#' Enter.tree <- ggrasp.load(tree.file, file.format = "tree", rank.file = rank.file.in);
#'
#' # Other options include loading by fasta file:
#' fasta.file <- system.file("extdata", "Enter.kSNP2.fasta", package="ggrasp")
#' rank.file.in <- system.file("extdata", "Enter.kSNP.ranks", package="ggrasp")
#' Enter.tree <- ggrasp.load(fasta.file, file.format = "fasta", rank.file =rank.file.in)
#'
#' # and by distance matrix. Since this distance matrix is actually percent identity,
#' # we will us an offset of 100
#' mat.file <- system.file("extdata", "Enter.ANI.mat", package="ggrasp")
#' rank.file.in <- system.file("extdata", "Enter.kSNP.ranks", package="ggrasp")
#' Enter.in <- ggrasp.load(mat.file, file.format = "matrix", rank.file =rank.file.in, offset = 100)
#'
#' # Use summary() to examine the data loaded
#' summary(Enter.in)
#'
#' #Use plot() to see the tree
#' \donttest{plot(Enter.in)}
#' @export



ggrasp.load = function(file, file.format, rank.file, offset, tree.method = "complete")
 {
	
requireNamespace("ggplot2");
requireNamespace("mixtools");
requireNamespace("ape");
requireNamespace("bgmm");
requireNamespace("colorspace");
requireNamespace("methods");
	phy = NULL;
	if (!exists("is.binary", mode="function"))
	{
		is.binary <- is.binary.tree;
	}
	if (missing(file))
	{
		cat("Please enter in a file");
		return();
	}
	if (!file.exists(file))
	{
		cat(paste("Cannot find file ", file ,"\n Please use an existing file\n", sep=""));
		return();
	}
	if (!missing(file.format))
	{
		file.format = tolower(file.format);
		if ( file.format != "tree" & file.format != "fasta" & file.format != "matrix")
		{
			cat(paste("Cannot recognize format ", file.format, ". Ignoring...", sep=""));
			file.format = NULL;
		
		}
	}
	if (missing(file.format))
	{
		file.format = NULL;
	}	
	
	if(is.null(file.format))
	{
		try.1 <- scan(file, what = "character");
		good = 0;
		if (substring(try.1[1], 0, 1) == "(")
		{
			phy <- read.tree(text = paste(try.1, collapse=""));
			if (exists("phy"))
			{
				phy = .root_midpoint(phy);
				if (is.binary(phy)==FALSE)
				{
					phy = multi2di(phy);
				}
				distance_matrix = cophenetic.phylo(phy);
				good = 1;
				file.format="done";
			}
		}
		if (substring(try.1[1], 0, 1) == ">")
		{
			fasta.in <- read.FASTA(file)
			if (exists("fasta.in"))
			{
				distance_matrix <- as.matrix(dist.dna(fasta.in, model="raw", ));
				good = 1;
				file.format="done";
			}
		}
		if (good == 0)
		{
			file.format = "matrix";
		}
	}
	if (file.format == "matrix")
	{
		distance_dataframe <- read.delim(file, sep="\t", header=TRUE, check.names=FALSE, comment.char="", quote="", row.names=1, strip.white=TRUE);
		if (exists("distance_dataframe"))
		{
			distance_matrix <- as.matrix(distance_dataframe);
			if (!missing(offset))
			{
				distance_matrix = offset - distance_matrix; 
			}	
		}
		else
		{
			cat("File not in matrix format. Please check and fix...\n\n");
			return();
		}
	}
	if (file.format == "fasta")
	{
		fasta.in <- read.FASTA(file)
		if (exists("fasta.in"))
		{
			distance_matrix <- as.matrix(dist.dna(fasta.in, model="raw"));
		}
		else
		{
			cat("File not recognized as fasta. Please check and fix...\n\n");
			return();
		}
	}		
	
	if (file.format == "tree")
	{
		phy = read.tree(file);
		if (exists("phy"))
		{
			phy = .root_midpoint(phy);
			if (is.binary(phy)==FALSE)
			{
				phy = multi2di(phy);
			}
			distance_matrix = cophenetic.phylo(phy);
		}
		else
		{
			cat("File not recognized as a newick tree. Please check and fix...\n\n");
			return();
		}
	}
	dist_mat <- as.dist(distance_matrix);
	accept_hclust <- c("complete", "average", "single", "nj")
	tree.method = tolower(tree.method)
	cat(paste(tree.method, "\n"));
	if (is.null(phy))
	{
		if (tree.method %in% accept_hclust)
		{
			hclust_m <-tree.method;
			if (tree.method == "nj")
			{
				phy <- nj(dist_mat);
				phy = .root_midpoint(phy);
				if (is.binary(phy)==FALSE)
				{
					phy = multi2di(phy);
				}
				
			}else
			{
				cat(summary(dist_mat));
				hc_phy <- hclust(dist_mat, method=tree.method)
				phy <- as.phylo(hc_phy);
				phy = .root_midpoint(phy);
				if (is.binary(phy)==FALSE)
				{
					phy = multi2di(phy);
				}
				
			}   
		}else
		{
			cat("\nSuggested unacceptible phylogenetic method...\n\nQuiting..\n");
			return()
		}
	}
	if (!missing(rank.file))
	{
		if (file.exists(rank.file))
		{
			rank.list <- read.table(rank.file, sep="\t");
			rank.lst = rank.list$V2;
			names(rank.lst) = rank.list$V1;
			if (sum(!names(rank.lst) %in% rownames(distance_matrix))>0) 
			{
				cat("Names in provided rank file do not match... Discarding all those that do not match...\n\n");
				cat("The non-matching names are:\n");
				cat(paste(names(rank.lst)[!names(rank.lst) %in% rownames(distance_matrix)], collapse="\n"))
			    rank.lst <- rank.lst[names(rank.lst) %in% rownames(distance_matrix)]
			}
			oth.names <- rownames(distance_matrix)[!rownames(distance_matrix) %in% names(rank.lst)];
			if (length(oth.names) > 0)
			{
				oth.ids <- rep(max(rank.lst)+1, length(oth.names)); 
				names(oth.ids) <- oth.names;
				rank.lst <- c(rank.lst, oth.ids)
			}
   
		}
	}
	if (!exists("rank.lst"))
	{
		rank.lst = vector(mode="numeric", length = 0);
	}
	phy.df = "";
	if (exists("phy"))
	{
		phy.df = write.tree(phy);
	}
	rownames(distance_matrix) = gsub(" ",  "_", rownames(distance_matrix));
	colnames(distance_matrix) = gsub(" ",  "_", colnames(distance_matrix));
	
	ggrasp.1 <- new("ggrasp", dist.mat = distance_matrix, rank = rank.lst, phy = phy.df, cluster=vector("numeric"), h = 0, medoids = vector("numeric"), gmm = data.frame(), gmm.orig = data.frame())
	return(ggrasp.1);
}




#' An S4 class representing the GGRaSP data and output
#'
#' @slot dist.mst The distance matrix showing the distances between different genomes
#' @slot phy The phylogenetic tree in newick format
#' @slot rank The ranks of the respective genomes with lower getting higher priority in being called as a medoid
#' @slot cluster A vector giving the numeric cluster ID for each genome
#' @slot h The threshold variable used to make the clusters
#' @slot medoids A vector giving the medoid for each cluster
#' @slot gmm A data.frame containing all the gaussian distributions used to find the threshold when available
#' @slot gmm.orig A data.frame containing all the gaussian distributions prior to cleaning. Used to recalculate the threshold when needed 

 
setClass("ggrasp", representation(dist.mat="matrix", phy = "character", rank="vector", cluster="vector", h = "numeric", medoids="vector", gmm="data.frame", gmm.orig="data.frame"));

setMethod("show", "ggrasp", function(object)summary(object))



#
.root_midpoint = function(phy1)
{
  #unroots the tree if rooted
  con.dist = cophenetic.phylo(phy1)
 
 if (is.rooted(phy1))
  {
    phy1 <- unroot(phy1)
  }
  #finds the two tips with the longest point between
  con.dist = cophenetic.phylo(phy1)
  tmp_n = rownames(con.dist)[rowSums(con.dist==max(con.dist))>0];
  node1 <- tmp_n[1];
	node2 <- tmp_n[con.dist[tmp_n[1], tmp_n] == max(con.dist)][1];
  #gets the node that is the most recent common ancestor
  mrca.tmp = getMRCA(phy1, c(node1, node2))
  #n1 is the integer that is the location of node 1 in the ape phylogeny tips
  n1 <- node1;
  if (node1 %in% phy1$tip.label)
  {
    n1 = (1:length(phy1$tip.label))[phy1$tip.label %in% node1];
  }
  #n2 is the same for node 2
  n2 <- node2;
  if (node2 %in% phy1$tip.label)
  {
    n2 = (1:length(phy1$tip.label))[phy1$tip.label %in% node2];
  }
  #something is not working... so quitting
  if (!(n1 %in% phy1$edge[,2]) | !(n2 %in% phy1$edge[,2]) | is.null(mrca.tmp))
  {
	cat(paste("Not working here.. ", node1, " ", node2, "\n\n", sep=""));
    q();
  }
  
  #Gets the distance between tip1 and the ancestor
  dst1 = phy1$edge.length[phy1$edge[,2]==n1];
  curr_nd = phy1$edge[phy1$edge[,2]==n1,1];
  while(curr_nd != mrca.tmp)
  {
    
    dst1 = dst1+phy1$edge.length[phy1$edge[,2]==curr_nd];
    curr_nd = phy1$edge[phy1$edge[,2]==curr_nd,1];
  }
  
  dst2 = phy1$edge.length[phy1$edge[,2]==n2];
  curr_nd = phy1$edge[phy1$edge[,2]==n2,1];
  while(curr_nd != mrca.tmp)
  {
    dst2 = dst2+phy1$edge.length[phy1$edge[,2]==curr_nd];
    curr_nd = phy1$edge[phy1$edge[,2]==curr_nd,1];
  }
  diff = abs(dst2 - dst1)/2;
  tot = abs(dst2 + dst1)/2
  #print(c(tot, diff))
  if (dst2 > dst1)
  {
    dst_f = phy1$edge.length[phy1$edge[,2]==n2];
    curr_nd = n2;
    while(dst_f < tot)
    {
      curr_nd = phy1$edge[phy1$edge[,2]==curr_nd,1];
      dst_f = dst_f+phy1$edge.length[phy1$edge[,2]==curr_nd];
    }
    dst_f = tot-(dst_f - phy1$edge.length[phy1$edge[,2]==curr_nd]);
  }
  else
  {
    
    dst_f = phy1$edge.length[phy1$edge[,2]==n1];
    curr_nd = n1;
    while(dst_f < tot)
    {
      curr_nd = phy1$edge[phy1$edge[,2]==curr_nd,1];
      dst_f = dst_f+phy1$edge.length[phy1$edge[,2]==curr_nd];
    }
    dst_f = tot-(dst_f - phy1$edge.length[phy1$edge[,2]==curr_nd]);  
  }
  tmp.tree <- read.tree(text="(a:0.1,b:0.1);");
  tmp2.tree = bind.tree(x = phy1, y = tmp.tree, where=curr_nd, position=dst_f);
  tmp2.tree = root(tmp2.tree, "a")
  tmp2.tree = drop.tip(tmp2.tree, c("a","b"));
  return(tmp2.tree);
}



