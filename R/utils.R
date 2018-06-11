
#' @title ggrasp.write
#' @description writes formatted information from a class GGRaSP object to a file. Multiple output options are available. 
#' 
#' @param x ggrasp-class object to be written
#' @param type Format of the data printed, either "tree" (New Hampshire extended style), "table" where the medoids or representative are shown in a table format, "list" where the information is shown in a pseudo-fasta format, or "itol" which prints out a file that can be loaded into the itol phylogeny viewer (http://itol.embl.de) which will color the clades of the different clusters
#' @param rank.level Maximum level of the rank to show. Ignored pre-clustering. After clustering, 0 will show only the medoids, -1 will show all values independent of rank, and any value >= 1 will show all the genomes less than or equal to that rank (including medoids). Default is 0 (only the medoids)
#' @param file String pointing to file where the data will be saved. If no file is given, the result will be printed out on the screen.
#'
#' @examples
#' #Getting the ggrasp object
#' Enter.tree <- ggrasp.load(system.file("extdata", "Enter.kSNP.tree", package="ggrasp"), 
#' file.format = "tree", rank.file =system.file("extdata", "Enter.kSNP.ranks", package="ggrasp"));
#' \donttest{Enter.tree.cluster <- ggrasp.cluster(Enter.tree)}
#'
#' #Default examples: using the initizalized ggrasp object will 
#' #write the newick tree string to "tree.nwk"
#' ggrasp.write(Enter.tree, type="tree", file=file.path(tempdir(), "tree.nwk"));
#' 
#' # Using the clustered ggrasp object will write a text file with the clusters saved as an ITOL clade
#' # In conjecture with the phylogeny, this is readable by 
#' # ITOL web phylogeny visualizer (http://itol.embl.de/) 
#' \donttest{ggrasp.write(Enter.tree.cluster, type="itol", file=file.path(tempdir(), "tree.itol.clade.txt"));}
#
#' @export

ggrasp.write <- function(x, file, type, rank.level)
{
	if (missing(x))
	{
		return();
	}
	write.out = "";
	
	if (length(x@medoids) == 0)
	{
		write.out <- (x@phy);
	}
	else
	{
		if (missing(type))
		{
			write.out <- (write.table(data.frame(medoid = x@medoids, cluster = sapply(x@medoids, function(y) {x@cluster[names(x@cluster)==y]})), sep="\t", quote=F, row.names=F));
	
		}
		else
		{
		type = tolower(type);
		if (type == "itol")
		{
			if (length(x@medoids) > 0)
			{
				itol = paste("TREE_COLORS", "SEPARATOR TAB", "DATA", sep="\n");
				col1 = as.data.frame(t(col2rgb(rainbow_hcl(max(x@cluster), start= 0.5)))); 
				cols = paste("rgba(", paste(col1$red, col1$green, col1$blue, "0.5", sep=","), ")", sep="");
				m.1 <- mrca(read.tree(text=x@phy))
				for (i in 1:max(x@cluster))
				{
					n.1 <- names(x@cluster)[x@cluster ==i]
					if (length(n.1) > 1)
					{
						m.1[rownames(m.1) %in% n.1, colnames(m.1) %in% n.1]
						m2 = min(m.1[rownames(m.1) %in% n.1, colnames(m.1) %in% n.1][1, 2:length(n.1)])
						out = paste(n.1[1], "|",  n.1[m.1[rownames(m.1) %in% n.1, colnames(m.1) %in% n.1][1,]==m2][1], sep="");
					}
					else
					{
						out = n.1[1];
					}
					out1 = paste(out, "range", cols[i], paste("cluster", i), sep="\t");
					out2 = paste(out, "clade", cols[i], "normal", "10", sep="\t");
			
					itol = paste(itol, out1, out2, sep="\n");
				}
			
				write.out = itol;
			}
		}
		if (type=="tree")
		{
			if (length(x@rank) > 0)
			{
				if (missing(rank.level))
				{
					write.out <- (.print.tree(x));
				}
				else
				{
					if (rank.level==0)
					{
						write.out <-(.print.tree(x));
					}
					else
					{
						if (rank.level == -1)
						{
							write.out <-(x@phy);
				
						}
						else
						{
							write.out <-(.print.tree.all(x, rank.level));
						}
					}
				}
			}
			else
			{
				if (missing(rank.level))
				{
					write.out <-(.print.tree(x));
				}
				else
				{
					if (rank.level==0)
					{
						write.out <-(.print.tree(x));
				
					}
					else
					{
						if (rank.level == -1)
						{
							write.out <-(x@phy);
				
						}
						else
						{
							if (rank.level >0)
							{
								write.out <-(.print.tree(x));

							}
						}
					}
				}
			}
		}
		if (type=="table")
		{
			if (length(x@rank) == 0)
			{
				write.out <-(.print.table(x, 0));
				
			}
			else
			{
				if (missing(rank.level))
				{
					write.out <-(.print.table(x, 0));
				}
				else
				{
					if(!is.numeric(rank.level) | rank.level==0)
					{
						write.out <-(.print.table(x,0));
				
					}
					else
					{
						if (rank.level == -1)
						{
							write.out <-(write.table(data.frame(genome = names(x@cluster), cluster = x@cluster)));
				
						}
						else
						{
							write.out <-(.print.table(x, rank.level));
						}
					}
				}
			}
		}
		if (type=="list")
		{
			if (length(x@rank)==0)
			{
				write.out <-(.write.centroid.fasta(x, 0));
				
			}
			else
			{
				if (missing(rank.level))
				{
					write.out <-(.write.centroid.fasta(x, 0));
				}
				else
				{
					if(!is.numeric(rank.level) | rank.level==0)
					{
						write.out <-(.write.centroid.fasta(x, 0));
					}
					else
					{
						if (rank.level == -1)
						{
							write.out <-(write.table(data.frame(genome = names(x@cluster), cluster = x@cluster)));
				
						}
						else
						{
							write.out <-(.write.centroid.fasta(x, rank.level));
						}
					}
				}
			}
		}
		}
	}
	if (missing(file))
	{
		cat(write.out)
		return();
	}
	else
	{
		write(write.out, file = file);
	}
}

#' @title print.ggrasp
#' @description prints formatted information from a class GGRaSP object. Multiple output options are available. 
#' 
#' @param x ggrasp-class object to be printed
#' @param type Format of the data printed, either "tree" (new hampshire extended style), "table" where the medoids or representative are shown in a table format, or "list" where the information is shown in a pseudo-fasta format
#' @param rank.level Maximum level of the rank to show. Ignored pre-clustering. After clustering, 0 will show only the medoids, -1 will show all values independent of rank, and any value >= 1 will show all the genomes less than or equal to that rank (including medoids). Default is 0 (only the medoids)
#' @param ... ignored
#'
#' @examples
#' #Getting the ggrasp object
#' Enter.tree <- ggrasp.load(system.file("extdata", "Enter.kSNP.tree", package="ggrasp"), 
#' file.format = "tree", rank.file =system.file("extdata", "Enter.kSNP.ranks", package="ggrasp"));
#' \donttest{Enter.tree.cluster <- ggrasp.cluster(Enter.tree)}
#'
#' #Default examples: using the initizalized ggrasp object will print the newick tree string 
#' print(Enter.tree);
#' 
#' # Using the clustered ggrasp object will print the medoids and their respective clusters
#' \donttest{print(Enter.tree.cluster)}
#
#' #Below are examples of using different output formats and rank levels
#' \donttest{print(Enter.tree.cluster, "tree")}
# Making a table with all the medoids and top-ranked genomes that are non-medoids
#' \donttest{print(Enter.tree.cluster, "table", 1)}
# Making a table with all genomes independent of rank
#' \donttest{print(Enter.tree.cluster, "table", 0)}

#' @method print ggrasp
#' @export
print.ggrasp <- function(x, type, rank.level, ...){
	if (missing(x))
	{
		return();
	}
	if (length(x@medoids) == 0)
	{
		print(x@phy);
	}
	else
	{
		if (missing(type))
		{
			print(write.table(data.frame(medoid = x@medoids, cluster = sapply(x@medoids, function(y) {x@cluster[names(x@cluster)==y]})), sep="\t", quote=F, row.names=F));
	
		}
		else
		{
		type = tolower(type);
		if (type=="tree")
		{
			if (length(x@rank) > 0)
			{
				if (missing(rank.level))
				{
					print(.print.tree(x));
				}
				else
				{
					if (rank.level==0)
					{
						print(.print.tree(x));
					}
					else
					{
						if (rank.level == -1)
						{
							print(x@phy);
				
						}
						else
						{
							print(.print.tree.all(x, rank.level));
						}
					}
				}
			}
			else
			{
				if (missing(rank.level))
				{
					print(.print.tree(x));
				}
				else
				{
					if (rank.level==0)
					{
						print(.print.tree(x));
				
					}
					else
					{
						if (rank.level == -1)
						{
							print(x@phy);
				
						}
						else
						{
							if (rank.level >0)
							{
								print(.print.tree(x));

							}
						}
					}
				}
			}
		}
		if (type=="table")
		{
			if (length(x@rank) == 0)
			{
				print(.print.table(x, 0));
				
			}
			else
			{
				if (missing(rank.level))
				{
					print(.print.table(x, 0));
				}
				else
				{
					if(!is.numeric(rank.level) | rank.level==0)
					{
						print(.print.table(x,0));
				
					}
					else
					{
						if (rank.level == -1)
						{
							print(write.table(data.frame(genome = names(x@cluster), cluster = x@cluster)));
				
						}
						else
						{
							print(.print.table(x, rank.level));
						}
					}
				}
			}
		}
		if (type=="list")
		{
			if (length(x@rank)==0)
			{
				print(.write.centroid.fasta(x, 0));
				
			}
			else
			{
				if (missing(rank.level))
				{
					print(.write.centroid.fasta(x, 0));
				}
				else
				{
					if(!is.numeric(rank.level) | rank.level==0)
					{
						print(.write.centroid.fasta(x, 0));
					}
					else
					{
						if (rank.level == -1)
						{
							print(write.table(data.frame(genome = names(x@cluster), cluster = x@cluster)));
				
						}
						else
						{
							print(.write.centroid.fasta(x, rank.level));
						}
					}
				}
			}
		}
		}
	}
}




#' @title plot.ggrasp
#' @description plots a class GGRaSP variable either the full tree, a reduced tree, or the gmm.
#' 
#' @param x ggrasp-class object to be plotted
#' @param type Type of plot generated, either "tree" (the full tree with the clusters shown as grouped leaves), "reduced" (tree with only the medoids shown), "hist" (histogram of the distribution of the distances) or "gmm" (histogram of the distances overlayed with the Gaussian distributions)
#' @param layout The layout style of the tree, either "circular" (default), "radial", "slanted", "linear" or "rectangular" ("linear" or "rectangular" are the same).
#' @param ... ignored
#'
#' @return A ggplot object containing the plot. It can be printed to standard output or saved using ggsave.
#' @method plot ggrasp
#'
#' @export



plot.ggrasp <- function(x, type, layout = "circular", ...)
{
	if (length(x@medoids) == 0)
	{
		if (missing(type))
		{
			return(.plot_hist(as.vector(x@dist.mat), (max(x@dist.mat)-min(x@dist.mat))/50));
		}
		type = tolower(type);
		if (type == "tree")
		{
			phy.print = read.tree(text=x@phy);
			return(.plot_tree(phy.print, method=layout));
		}
		if (type=="histogram")
		{
			return(.plot_hist(as.vector(x@dist.mat), (max(x@dist.mat)-min(x@dist.mat))/50));
		}
	}
	if (missing(type))
	{
		phy.print = read.tree(text=x@phy);
		return(.plot_tree(phy.print, x@cluster, x@medoids, layout, x@h));
	}
	type = tolower(type);
	if (type=="tree")
	{
		phy.print = read.tree(text=x@phy);
		return(.plot_tree(phy.print, x@cluster, x@medoids, layout, x@h));
	}
	if (type=="gmm")
	{
		stp_val = (max(x@dist.mat)-min(x@dist.mat))/50;
		for (i in 1:length(x@gmm$sigma))
		{
			if (stp_val > x@gmm$sigma[i])
			{
				stp_val =  x@gmm$sigma[i];
			}
		}
		return(.plot_gmm(x@gmm, as.vector(x@dist.mat), stp_val, x@h));
	}
	if (type=="histogram")
	{
		return(.plot_hist(as.vector(x@dist.mat), (max(x@dist.mat)-min(x@dist.mat))/50));
	}
	if (type=="trimmed")
	{
		phy.print = read.tree(text= x@phy);
		phy.new = drop.tip(phy.print, phy.print$tip.label[!phy.print$tip.label %in% x@medoids]);
		phy.new[phy.new$tip];
		return(.plot_tree(phy.new, layout));
	}
}

#' @title summary.ggrasp
#' @description prints a summary of the class GGRaSP variable. Output includes medoids and cutoff value after the clustering
#' @param object ggrasp-class object
#' @param ... ignored
#' @method summary ggrasp

#' 
#' @export

 
summary.ggrasp <- function(object, ...)
{
	cat(paste("Contains ", nrow(object@dist.mat), " genomes", sep=""));	
	if (length(object@rank) > 0)
	{
		cat(" and is ranked\n");
	}
	else
	{
		cat("\n");
	}
	if (length(object@cluster) > 0)
	{
		cat(paste("\nContains ", max(object@cluster), " clusters using threshold at ", object@h, "\n", sep=""));
	}
}

#Interior functons: not callable except by GGRaSP high level calls

#Creates a ggplot2 phylogeny, by calling .plot_tree. Redundant
.draw_tree = function(phy1, clusters, medoids, layout="rectangular")
{
  graph_tree <- .plot_tree(phy1, clusters, medoids, method = layout)
  return(graph_tree)
}




#Goes through the clusters and returns all non-representitve genomes with a rank at or below the keep.rank variable
#The other variables are generated from the GGRaSP S4 class variables
#Maybe change to return a list?
.get.sub.cnts = function(medoids, clusters, rank.lst, keep.rank)
{
	ret.cnts <- data.frame(id = c(), cnts = c());
	for (i in 1:length(medoids))
	{
		if (is.null(keep.rank))
		{
			ret.cnts = rbind(ret.cnts, data.frame(id = medoids[i], cnts = sum(clusters==i)));
		}
		else
		{
			for (j in names(clusters)[clusters == i & !names(clusters) %in% medoids])
			{
				print(j)
			}
		}
	}

}

#Makes a string with the pseudo fasta list of the clusters. Medoids are shown as the header.
#x is a GGRaSP class variable. If keep.rank is given, all non-representative genomes at or below
#this rank are printed with a double >> as a flag. The non-representative ranked genomes are printed
#first, followed by the rest of the genomes.
.write.centroid.fasta = function(x, keep.rank)
{
  out.str = "";
 
  for (i in 1:length(x@medoids))
  {
	tmp ="";
    if (length(x@rank) > 0)
    {
      tmp = paste("\t", x@rank[names(x@rank) == x@medoids[i]], sep="")
    }
	tmp2 = sum(x@cluster ==i);
    out.str = paste(out.str, ">", x@medoids[i], tmp, "\t", tmp2, "\n", sep="");
	out.str2 = "";
	out.str3 = "";
    for (j in names(x@cluster)[x@cluster == i & !names(x@cluster) %in% x@medoids])
    {
      tmp = "";
		add = "";

      if (length(x@rank)>0)
      {
        tmp = paste("\t", x@rank[names(x@rank) == j], sep="")
		if (!is.null(keep.rank)) 
		{
			if (x@rank[names(x@rank) == j] <= keep.rank)
			{		
				add = ">>";
      			out.str3 = paste(out.str3, add, j, tmp, "\t0\n", sep="");		
			}else
			{
				out.str2 = paste(out.str2, add, j, tmp, "\t0\n", sep="");
			}
		}else
		{
			out.str2 = paste(out.str2, add, j, tmp, "\t0\n", sep="");
		}      
      }
	  else
		{
			out.str2 = paste(out.str2, add, j, tmp, "\t0\n", sep="");
		}
	}
	out.str = paste(out.str, out.str2, out.str3, sep="");
  }
  return(out.str);
}

#Creates a ggplot2 style object containing the tree. Tree is midpoint rooted.
#Can take a while- recommend using ITOL on the printed phylogeny.
.draw_tree = function(phy1)
{
	#
  if (!is.rooted(phy1))
  {
    con.dist <- cophenetic.phylo(phy1);
    tmp <- rownames(con.dist)[rowSums(con.dist==max(con.dist))>0];
    phy1 = .root_midpoint(phy1)
  }
  if (!is.binary(phy1))
  {
	phy1 <- multi2di(phy1);
  }
  tmp_lst = lapply(1:length(phy1$tip.label), function(x)
  {
    data.frame(x1 = 0, y1 = 0, x2 = 0, y2 = 0, type = "node", make = "both", is.tip = T, label=phy1$tip.label[x], id = x);
  });
  tmp_r = 1:length(phy1$tip.label);
  nd = node.depth(phy1, method=2);
  nd.mrca = mrca(phy1);
  for (i in 2:max(nd))
  {
    nds.tmp = (1:length(nd))[nd == i]
    hits = rep(0, length(tmp_r)); 
    for (j in nds.tmp)
    {
      tmp_nd = phy1$edge[phy1$edge[,1]==j,2];
      hits[tmp_r ==tmp_nd[1]]= 1;
      hits[tmp_r ==tmp_nd[2]]= 1;
    }
    tmp_lst2 = lapply(nds.tmp, function(x)
    {
      tmp_nd = phy1$edge[phy1$edge[,1]==x,2];
      if (sum(rowSums(nd.mrca==tmp_nd[1])>0) > sum(rowSums(nd.mrca==tmp_nd[2])>0))
      {
        swap.nd = tmp_nd[1]; tmp_nd[1] = tmp_nd[2]; tmp_nd[2] = swap.nd;
      }
      tmp1 = tmp_lst[[(1:length(tmp_r))[tmp_r ==tmp_nd[1]]]];
      dst1 = phy1$edge.length[phy1$edge[,2]== tmp_nd[1]];
      x1_n = tmp1$x1[tmp1$id == tmp_nd[1] & tmp1$type=="node"];
      y1_n = tmp1$y1[tmp1$id == tmp_nd[1] & tmp1$type=="node"];
      tmp1$x1 = tmp1$x1 + dst1;
      tmp1$x2 = tmp1$x2 + dst1;
      tmp2 = tmp_lst[[(1:length(tmp_r))[tmp_r == tmp_nd[2]]]];
      dst2 = phy1$edge.length[phy1$edge[,2] == tmp_nd[2]];
      x2_n = tmp2$x1[tmp2$id == tmp_nd[2] & tmp2$type=="node"];
      y2_all = max(c(tmp1$y1, tmp1$y2))+1;
      y2_n = tmp2$y1[tmp2$id == tmp_nd[2] & tmp2$type=="node"] +y2_all;
      tmp2$x1 = tmp2$x1 + dst2;
      tmp2$x2 = tmp2$x2 + dst2;
      tmp2$y1 = tmp2$y1 + y2_all;
      tmp2$y2 = tmp2$y2 + y2_all;
      
      tmp1 = rbind(tmp1,
                   tmp2,
                   data.frame(x1 = 0, y1 = y1_n, x2 = dst1, y2 = y1_n, type = "edge", make = "rect", is.tip = F, label=NA, id = x),
                   data.frame(x1 = 0, y1 = y2_n, x2 = dst2, y2 = y2_n, type = "edge", make = "rect", is.tip = F, label=NA, id = x),
                   data.frame(x1 = 0, y1 = (y1_n+y2_n)/2, x2 = dst1, y2 = y1_n, type = "edge", make = "slant", is.tip = F, label=NA, id = x),
                   data.frame(x1 = 0, y1 = (y1_n+y2_n)/2, x2 = dst2, y2 = y2_n, type = "edge", make = "slant", is.tip = F, label=NA, id = x),
                   data.frame(x1 = 0, y1 = y1_n, x2 = 0, y2 = y2_n, type = "edge", make = "rect", is.tip = F, label=NA, id = x),
                   data.frame(x1 = dst1/2, y1 = min(tmp1$y1), x2 = max(tmp1$x2), y2 = max(tmp1$y2), type="box", make = "both",is.tip = F, label=tmp1$label[tmp1$type=="node" & tmp1$id == tmp_nd[1]], id = tmp_nd[1]),
                   data.frame(x1 = dst2/2, y1 = min(tmp2$y1), x2 = max(tmp2$x2), y2 = max(tmp2$y2), type="box", make = "both",is.tip = F, label=tmp2$label[tmp2$type=="node" & tmp2$id == tmp_nd[2]], id = tmp_nd[2]),
                   data.frame(x1 = 0, y1 = (y1_n+y2_n)/2, x2 = 0, y2 = (y1_n+y2_n)/2, type="node", make = "both", is.tip=F, label=as.character(x), id=x));
    });
    tmp_r = c(tmp_r[hits==0], nds.tmp);
    tmp_lst=c(tmp_lst[hits==0], tmp_lst2);
    # print(c(length(tmp_r), length(tmp_lst), sum(hits ==1)))
    
  }
  if (length(tmp_lst)==1)
  {
    return(tmp_lst[[1]])
  }
  return(tmp_lst)
}

.make_arc_text = function(string.text, x_plot_max, x_plot_max_add, y_plot_max)
{
  str.tmp = strsplit(string.text, ",")[[1]];
  ret.df = NULL;
  for (i in 1:length(str.tmp))
  {
    x1 <- strsplit(str.tmp[i], "/")[[1]];
    x2 = x1[length(x1)];
    df.x = -.75 - 2*(length(str.tmp)-i);
    y.kern.len = 10;
    arc.text.df <- data.frame(chr = strsplit(x2, "")[[1]], x.1 = rep(x_plot_max+df.x*x_plot_max_add, nchar(x2)), y.1 = c(0, sapply(1:(nchar(x2)-1), function(x) {strwidth(paste(strsplit(x2, "")[[1]][1:x], collapse=""), "figure") * y.kern.len})));
    y.mx <- max(arc.text.df$y.1)/2
    arc.text.df$y.1 = arc.text.df$y.1 - y.mx
    arc.text.df$y.1[arc.text.df$y.1 <0] = arc.text.df$y.1[arc.text.df$y.1<0] + y_plot_max
    ang3 = -1*(arc.text.df$y.1)/y_plot_max * 360;
    ang3[ang3 < 0] = 360 + ang3[ang3 < 0];
    ang3[ang3 > 90 & ang3 < 270] = ang3[ang3 > 90 & ang3 < 270] -180;
    ang3[ang3 < 0] = 360 + ang3[ang3 < 0];
    arc.text.df$ang = ang3;
    if (is.null(ret.df))
    {
      ret.df = arc.text.df;
    }
    else
    {
      ret.df = rbind(ret.df, arc.text.df)
    }
  }
  return(ret.df);
}


.plot_gmm = function(gmm, gmm.x, stp = 0.1, thresh= NULL){
   count = NULL
  n1 <- as.vector(gmm.x);
  
  df_area <- data.frame(
    dist.id = as.vector(sapply(1:length(gmm$mu), function(x) {sapply(seq(min(n1), max(n1), by=stp), function(y){x})})),
    loc = as.vector(sapply(1:length(gmm$mu), function(x) {sapply(seq(min(n1), max(n1), by=stp), function(y){y})})),
    den = as.vector(sapply(1:length(gmm$mu), function(x) {sapply(seq(min(n1), max(n1), by=stp), function(y){dnorm(x = y, mean=gmm$mu[x], sd=gmm$sigma[x])*gmm$lambda[x]})})));
  
  
  grph_gmm <- ggplot2::ggplot()+ggplot2::geom_histogram(aes_string(x = "n1", y= "..count.."), binwidth= stp)+with(df_area, geom_area(data=df_area, aes(x = loc, y = den*length(n1)*stp, fill=as.factor(dist.id), alpha = 0.5), position="identity"));
  grph_gmm <- grph_gmm + labs(x = "Percent Difference", y = "Number of Genome Pairs", fill="Gaussian Mixture Model")+ guides(alpha="none");
  if (!is.null(thresh))
  {
    grph_gmm <- grph_gmm + geom_segment(aes(x = thresh, y = 0, xend = thresh, yend = (max(df_area$den)*stp*length(n1))), linetype=2);
  }
  return(grph_gmm);
}



.plot_hist = function(gmm.x, stp = 0.1)
{
  count = NULL
  n1 <- as.vector(gmm.x);

  
  
  grph_gmm <- ggplot()+geom_histogram(aes_string(x = "n1", y= "..count.."), binwidth= stp);
  grph_gmm <- grph_gmm + labs(x = "Percent Difference", y = "Number of Genome Pairs")

  return(grph_gmm);
}


.get_name_line = function(phy_df, name_file=NULL, x_max = 0, num=NULL, id=NULL)
{
  phy.names <- substring(phy_df$label[phy_df$is.tip==TRUE],0);
  y_plot_max = max(c(phy_df$y1, phy_df$y2));
  if (!is.null(name_file))
  {
    name.list <- read.table(name_file, sep="\t", stringsAsFactors = F);
    name.df <- NULL; 
    strt = 0; strt.val = name.list$V2[name.list$V1 == phy.names[1]];
    rnbw <- rainbow(length(unique(name.list$V2)));
    num.hit = 1;
    curr = 0;
    names(rnbw) = substring(unique(name.list$V2),0);
    for (i in 2:nrow(name.list))
    {
      curr = curr + 1;
      curr.val <- name.list$V2[name.list$V1 == phy.names[i]]
      if (curr.val != strt.val)
      {
        if (is.null(name.df))
        {
          name.df = data.frame(y1 = max(0,strt - 0.33), x1 = x_max, x2 = x_max, y2 = curr - 0.66, id =  strt.val, color=rnbw[names(rnbw) == strt.val], num=num, name = id, num.hit = num.hit);
        }else
        {
          name.df = rbind(name.df, data.frame(y1 = max(0,strt - 0.33), x1 = x_max, x2 = x_max, y2 = curr - 0.66, id = strt.val, color=rnbw[names(rnbw) == strt.val], num=num, name = id, num.hit = num.hit))
        }
        strt.val <- curr.val;
        strt = curr;
        num.hit = 0;
      }
      num.hit = num.hit + 1;
    }
    if (is.null(name.df))
    {
      
      name.df = data.frame(y1 = strt - 0.33, x1 = x_max, x2 = x_max, y2 = curr, id = strt.val, color=rnbw[names(rnbw) == strt.val], num=num, name = id, num.hit = num.hit);
    }else
    {
      name.df = rbind(name.df, data.frame(y1 = strt - 0.33, x1 = x_max, x2 = x_max, y2 = curr, id = strt.val, color=rnbw[names(rnbw) == strt.val], num=num, name = id, num.hit = num.hit))
    }
    ang2[ang2 < 0] = 360 + ang2[ang2 < 0];
    ang2 = 90-(name.df$y1 + name.df$y2)/2/y_plot_max * 360;
    hjst = rep(-0.1, nrow(name.df));
    ang2[ang2 < 0] = 360 + ang2[ang2 < 0];
    hjst[ang2 > 90 & ang2 < 270] = 1.1
    ang2[ang2 > 90 & ang2 < 270] = ang2[ang2 > 90 & ang2 < 270] -180;
    ang2[ang2 < 0] = 360 + ang2[ang2 < 0];
    name.df$ang2 = ang2;
    name.df$hjst = hjst;
    return(list(name.df, rnbw)) 
  }
}

.plot_tree = function(phy, clusters=NULL, medoids=NULL, method = "circular", thresh = NULL, name_file=NULL, ring_labels=NULL)
{
  if (!is.rooted(phy))
  {
    phy = .root_midpoint(phy)
  }
  is.ultra = is.ultrametric(phy);
  phy_df <- .draw_tree(phy);
  make_type = "rect";
  dfr = 10;
  if (method == "radial" | method == "slanted")
  {
    make_type = "slant";
  }
  y_plot_max = sum(phy_df$is.tip)#+sum(phy_df$is.tip)/40;
  x_plot_max_add = max(node.depth.edgelength(phy))/dfr;
  x_plot_max = max(node.depth.edgelength(phy)) + max(node.depth.edgelength(phy))/dfr;
   p <- ggplot2::ggplot() + with(phy_df, ggplot2::geom_segment(data=phy_df[phy_df$type == "edge" & phy_df$make==make_type,], aes(x = x1, y= y1, xend=x2, yend=y2)));
  if (method == "circular" | method == "radial")
  {
    p <- p + coord_polar(theta="y")
  }
  p <- p + theme(axis.text.x=element_blank(), axis.text.y=element_blank())+labs(x="", y="") + scale_y_continuous(limits = c(0, y_plot_max));
  if (is.null(clusters) & is.null(medoids))
  {
    ang = 90-(phy_df$y1[phy_df$type=="node" & phy_df$is.tip] + phy_df$y2[phy_df$type=="node" & phy_df$is.tip])/2/y_plot_max * 360;
    hjst = rep(-0.1, sum(phy_df$type == "node" & phy_df$is.tip));
    ang[ang < 0] = 360 + ang[ang < 0];
    hjst[ang > 90 & ang < 270] = 1.1
    ang[ang > 90 & ang < 270] = ang[ang > 90 & ang < 270] -180;
    ang[ang < 0] = 360 + ang[ang < 0];
  
	tip.lab.df <- data.frame(phy_df[phy_df$type == "node" & phy_df$is.tip,], ang = ang, hjst = hjst);
	if (method == "circular" | method == "radial")
	{
		sz = 300/nrow(tip.lab.df);
		p <- p + with(tip.lab.df,  ggplot2::geom_text(data = tip.lab.df, aes(x = x1, y = y1, label = label, angle = ang, hjust = hjst), size=sz))
	}else
	{
		sz = 200/nrow(tip.lab.df);
		 p <- p +with(tip.lab.df, ggplot2::geom_text(data = tip.lab.df, aes(x = x1, y = y1, label = label, hjust = -0.1),  size=sz)+ guides(alpha=F, fill=F, size=F))
	}
	p <- p + xlim(c(0, max(phy_df$x1)+max(phy_df$x1)/20));
  }
  df.col = NULL;
  if (!is.null(clusters) & !is.null(medoids))
  {
    rn <- rainbow(1+sum(sapply(unique(clusters), function(x) {sum(clusters==x)})>1));
	if (is.ultra & !is.null(thresh))
	{
		cat("Drawing ultrametric line at threshold...\n\n")
		p <- p + ggplot2::geom_segment(aes(x = max(phy_df$x2)-thresh, y = 0, xend = max(phy_df$x2)-thresh, yend = max(phy_df$y2)+1), linetype=2)
	}
    cnt = 1;
    clust_df = NULL;
    for (i in 1:max(clusters))
    {
    
      if (sum(clusters==i) > 1)
      {
        tmp <- getMRCA(phy, names(clusters)[clusters==i]);
        cnt = cnt +1;
        ang = 90-(phy_df$y1[phy_df$type=="box" & phy_df$id %in% tmp] + phy_df$y2[phy_df$type=="box" & phy_df$id %in% tmp])/2/y_plot_max * 360;
        hjst = -0.1;
        if (ang < 0)
        {
          ang = ang + 360;
          if (hjst  == 1.1)
          {
            hjst = -0.1;
          }
        }
        if (ang > 90 & ang < 270)
        {
          ang = ang -180
          hjst = 1.1
          if (ang < 0)
          {
            ang = ang + 360
          }
        }
        if (is.null(clust_df))
        {
          clust_df <- data.frame(phy_df[phy_df$type=="box" & phy_df$id %in% tmp,], label_med = as.character(medoids[i]), col = tmp, ang = ang, hjst = hjst, cnt = sum(clusters==i));
          #df.col <- data.frame(id = tmp, col = rn[cnt])
        }else
        {
          clust_df <- rbind(clust_df, data.frame(phy_df[phy_df$type=="box" & phy_df$id %in% tmp,], label_med = as.character(medoids[i]), col = tmp, ang = ang, hjst = hjst, cnt = sum(clusters==i)));
          #df.col <- rbind(df.col, data.frame(id = tmp, col = rn[cnt]))
        }
      }
    }

	if (is.null(name_file))
    {
      p <- p + with(clust_df,  ggplot2::geom_rect(data=clust_df, aes(xmin = x1, ymin = y1, xmax=x2, ymax = y2, alpha = 0.5), fill="grey50"))
    }
    if (!is.ultra & !is.null(thresh))
	{
		p <- p +  with(clust_df, ggplot2::geom_segment(data=clust_df, aes(x = x1 - thresh, y = y1, xend = x1-thresh, yend = y2), linetype=2))
	}

	if (!is.null(name_file))
    {
      name.l1 <- strsplit(name_file, ",");
      name.df = NULL;
      for (j in 1:length(name.l1[[1]]))
      {
        if (is.null(name.df))
        {
          tmp.lst <- .get_name_line(phy_df, name.l1[[1]][j], x_plot_max, paste("Row ", as.character(j), sep="", collapse=""), paste("Row ", as.character(j),":", name.l1[[1]][j], collapse="", sep=""));
          
          name.df <-tmp.lst[[1]];
          names(tmp.lst[[2]]) <- paste("Row ", as.character(j), ": ", names(tmp.lst[[2]]), sep="");
          df.col2 <- tmp.lst[[2]];
          x_plot_max = x_plot_max + 2 * x_plot_max_add;
        }
        else{
          tmp.lst <- .get_name_line(phy_df, name.l1[[1]][j], x_plot_max, paste("Row ", as.character(j), sep="", collapse=""), paste("Row ", as.character(j),":", name.l1[[1]][j], collapse="", sep=""))
          
          name.df <-rbind(name.df, tmp.lst[[1]]);
          names(tmp.lst[[2]]) <- paste("Row ", as.character(j), ": ", names(tmp.lst[[2]]), sep="");
          df.col2 <- c(df.col2, tmp.lst[[2]])
          x_plot_max = x_plot_max + 2*x_plot_max_add;
          
        }
      }
      leg.list = list();
      if (!is.null(name.df))
      {
        p_old = p;
        #boxes for 
        
        p <- p + with(name.df, ggplot2::geom_rect(data=name.df[name.df$num.hit >1,], aes(xmin = x1, ymin = y1, xmax= x2+x_plot_max_add *0.9, ymax= y2, fill =  num), alpha = 0.5));
        #labels for all
        
        p <- p + with(name.df, ggplot2::geom_text(data=name.df, aes(x = x1, y = (y1+y2)/2, label=paste(toupper(substring(id, 0, 1)), substring(id,2,3), sep=""), angle=ang2, hjust=hjst2), size=2.5, fontface="bold"));
        #ring labels

		if (is.null(ring_labels))
        {
          arc.text.df = .make_arc_text(name_file, x_plot_max, x_plot_max_add, y_plot_max);
        }else
        {
          labels.l1 <- strsplit(ring_labels, ",");
          if (length(labels.l1[[1]]) == length(name.l1[[1]]))
          {
            arc.text.df = .make_arc_text(ring_labels, x_plot_max, x_plot_max_add, y_plot_max);
          }else
          {
            cat("Ring labels have a different length than number of ring. Defaulting to the file names...");
            arc.text.df = .make_arc_text(name_file, x_plot_max, x_plot_max_add, y_plot_max);
            
            
          }

        }
        p <- p + with(arc.text.df, ggplot2::geom_text(data=arc.text.df, aes(x = x.1, y = y.1, label=substring(chr,0), xend=x.1, yend=y.1, angle = ang), size=4, hjust = 0, vjust = 0, fontface="bold"));
        p <- p + ggplot2::scale_fill_manual(label=unique(name.df$num), values = rainbow(length(unique(name.df$num))))
        p <- p + with(clust_df, ggplot2::geom_rect(data=clust_df, aes(xmin = x1, ymin = y1, xmax=x_plot_max, ymax = y2, alpha = 0.5), fill="grey50"))

      }
    }
    p <- p +with(clust_df, ggplot2::geom_segment(data=clust_df, aes(x = x_plot_max, y = (y1), xend=x_plot_max, yend = y2)))
    
    if (method == "circular" | method == "radial")
    {
      p <- p + with(clust_df, ggplot2::geom_text(data=clust_df, aes(x = x_plot_max, y = (y1+y2)/2, label=label_med, hjust =hjst, angle = ang)));
    }else
    {
      p <- p + with(clust_df, ggplot2::geom_text(data=clust_df, aes(x = x_plot_max, y = (y1+y2)/2, label=label_med, hjust =-0.1)));
    }
    p <- p+ ggplot2::guides(alpha=F, fill=F)+labs(fill="Outside Ring Categories", color="Outside Ring Derivation Files")

  }

  p <- p + ggplot2::theme(axis.line = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank() , panel.border=element_blank(), panel.background = element_blank(), axis.ticks.y = element_blank(), axis.ticks.x = element_blank())

  return(p);
}



.print.table <- function(ggrasp.data, rank.level)
{
	out = data.frame();
	if (rank.level == 0)
	{
		for (med.id in ggrasp.data@medoids)
		{
			clust.id = ggrasp.data@cluster[names(ggrasp.data@cluster) == med.id];
			med.num = sum(ggrasp.data@cluster==clust.id);
			med.list = paste(names(ggrasp.data@cluster)[ggrasp.data@cluster == clust.id], collapse=",");
			out = paste(out, paste(med.id, clust.id,  med.num, med.list, sep="\t"), sep="\n");
		}
	}
	if (rank.level > 0)
	{
		near.neigh = .find_nearest_neighbor(ggrasp.data, keep.rank = rank.level);
		for (med.id in unique(near.neigh[,2]))
		{
			clust.id = ggrasp.data@cluster[names(ggrasp.data@cluster) == med.id];
			
			med.num = sum(near.neigh[,2]==med.id);
			med.list = paste(near.neigh[near.neigh[,2]==med.id,1], collapse=",");
			out = paste(out, paste(med.id, clust.id,  med.num, med.list, sep="\t"), sep="\n");
		}
	}
	return(out);
}

.print.tree <- function(ggrasp.data)
{
	phy2 <- read.tree(text = ggrasp.data@phy);
	phy2 = drop.tip(phy2, phy2$tip.label[!phy2$tip.label %in% ggrasp.data@medoids]);
	for (tip.id in phy2$tip.label)
	{
		med.num = sum(ggrasp.data@cluster == ggrasp.data@cluster[names(ggrasp.data@cluster)==tip.id]);
		med.id = "";
		if(length(ggrasp.data@medoids) > 0 & med.num > 0)
		{
			med.id = ggrasp.data@medoids[med.num];
		}
		add.tip.info = "";
		phy2$tip.label[phy2$tip.label == tip.id] = paste(tip.id, "[&&NHX!count=", med.num, ";medoid=", med.id,";]", sep="");
	}
	phy2.str <- write.tree(phy2);
	phy2.str <- gsub( "!", ":", phy2.str);
	return(phy2.str);
}

.print.tree.all <- function(ggrasp.data, rank.level)
{
	phy2 <- read.tree(text = ggrasp.data@phy);
	phy2 = drop.tip(phy2, phy2$tip.label[phy2$tip.label %in% names(ggrasp.data@rank)[ggrasp.data@rank > rank.level] & !(phy2$tip.label %in% ggrasp.data@medoids)]);
	near.neigh = .find_nearest_neighbor(ggrasp.data, keep.rank = rank.level);
	for (tip.id in phy2$tip.label)
	{
		med.num = sum(near.neigh[,2]==tip.id);
		add.tip.info = "";
		phy2$tip.label[phy2$tip.label == tip.id] = paste(tip.id, "[&&NHX!count=", med.num,"]", sep="");
	}
	phy2.str <- write.tree(phy2);
	phy2.str <- gsub( "!", ":", phy2.str);
	return(phy2.str);
}

.find_nearest_neighbor = function(x, keep.rank)
{
	good.list <- names(x@rank)[x@rank <= keep.rank | names(x@rank) %in% x@medoids];
	match.list <- names(x@rank)[x@rank > keep.rank & !(names(x@rank) %in% x@medoids)];
	ret.neigh <- cbind(match.list, rep(NA, length(match.list)));
	for (i in match.list)
	{
		ret.neigh[ret.neigh[,1]==i,2] = 
		colnames(x@dist.mat)[colnames(x@dist.mat) %in% good.list][x@dist.mat[rownames(x@dist.mat) ==i, 
		colnames(x@dist.mat) %in% good.list] == min(x@dist.mat[rownames(x@dist.mat) == i, colnames(x@dist.mat) %in% good.list])][1];
	}
	return(ret.neigh);
}
	