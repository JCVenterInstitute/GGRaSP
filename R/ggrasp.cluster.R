#' @title ggrasp.cluster
#' @description ggrasp.cluster() clusters the genomes in a GGRaSP class variable and assigns the most representative genome in each cluster after accounting for rank as a medoid.
#' 
#' 
#' @param ggrasp.data  Required. If neither a threshold or a num.cluster is given, a mixed model of Gaussian distributions is used to estimate a threshold to use the cluster.
#' @param threshold The threshold used to cluster together all genomes within this distance.
#' @param num.clusters Create this number of clusters independent of the cluster. 
#' @param gmm.start Number of Gaussian distributions to start the examination. Must be at least 2 and not greater than the gmm.max.  
#' @param gmm.max Maximum number of Gaussian distributions to examine. Has to be at least 2. 10 is the default
#' @param z.limit All Gaussian distributions with means within this number of standard deviations will be reduced to only the larger distribution. Defaults to 1. Set to 0 to keep all non-overlapping distributions.
#' @param min.lambda All Gaussian distributions with lambda value (proportion of the total distribution) below this value are removed before calculating the threshold. Default is 0.005. Set to 0 to keep all. 
#' @param run.type String giving the package to use to get the mixture model. Currently "bgmm" (default) and mixtools" are implemented.
#' @param left.dist Number giving the number Gaussian distribution model immediately to the left of the threshold used. 1 is the default. Only value between 1 and k-1 where k is the total number of number of Gaussian distributions. 
#' 
#' @return Returns a class GGRaSP variable with the clusters and medoids assigned. In cases where the Gaussian Mixture Model was used to estimate the cutoff threshold, the descriptive values of the different distributions is also stored
#' @examples #The following data is from Chavda et al 2016 which phylotyped Enterobacter genomes
#' # Our example uses the data underpinning the tree shown in Figure 2
#' # Also included is a ranking file to prioritize closed Enterobactor genomes
#'
#' #Loading the tree 
#' library(ggrasp);
#' tree.file <- system.file("extdata", "Enter.kSNP.tree", package="ggrasp")
#' rank.file.in <- system.file("extdata", "Enter.kSNP.ranks", package="ggrasp")
#' Enter.tree <- ggrasp.load(tree.file, file.format = "tree", rank.file = rank.file.in);
#'
#' #Clustering the tree using a threshold estimated by Gaussian Mixture Models (GMMs)
#' \donttest{Enter.tree.cluster <- ggrasp.cluster(Enter.tree)}
#'
#'
#' #Use print to get a list of the medoids selected
#' \donttest{print(Enter.tree.cluster)}
#'
#' #Re-clustering the tree using a threshold estimated by the GMMs but without the distribution
#' #cleaning (i.e. removing the overlapping and low count distributions)
#' \donttestEnter.tree.reclust <- ggrasp.recluster(Enter.tree.cluster, z.limit=0, min.lambda = 0)}
#'
#' #Use plot to examine the tree with the clusters highlighted and the medoid genome names on the edge
#' \donttest{plot(Enter.tree.cluster)}
#'
#' #Additional printing and plotting options are availible with plot() and print(). 
#' #For more information refer to ?plot.ggrasp and ?print.ggrasp
#' @export



ggrasp.cluster = function(ggrasp.data, threshold, num.clusters, z.limit = 1, gmm.start = 2, gmm.max = 10, min.lambda=0.005, run.type="bgmm", left.dist = 1)
{
	if(missing("ggrasp.data") || class(ggrasp.data) != "ggrasp")
	{
		cat("Required class GGRaSP data to be performed. Please run ggrasp.load() first\n");
		return();
	}
	if (!missing("threshold"))
	{
		if (is.numeric(threshold) == FALSE)
		{
			cat("Incorrect threshold variable... Not using\n\n");
			threshold = NA;
		}
	}
	if (!missing("num.clusters"))
	{
		if (is.numeric(num.clusters) == FALSE)
		{
			cat("Incorrect threshold variable... Not using\n\n");
			num.clusters = NA;
		}
	}
	if ((missing("threshold") ) & (missing("num.clusters") ))
	{
		if (!run.type %in% c("bgmm", "mixtools"))
		{
			cat("Unable to find run mode. Using bgmm instead...\n\n");
			run.type = "bgmm";
		}
		if (gmm.start < 2 || gmm.start > gmm.max)
		{
			cat("Number of starting GMM out of range. Using default of 2\n\n");
			gmm.start = 2;
		}
		gmm.first <- .make_gmm(as.dist(ggrasp.data@dist.mat), gmm.start, gmm.max, run.type);
		if ("mModel" %in% class(gmm.first))
		{
			gmm.first.2 <- data.frame(mu = gmm.first$mu[,1], sigma = gmm.first$cvar[,,1], lambda = gmm.first$pi);
		}
		else
		{
			gmm.first.2 <- data.frame(mu = gmm.first$mu, sigma = gmm.first$sigma, lambda = gmm.first$lambda)

		}
		
		#ks.val <- ks.test(as.vector(as.dist(ggrasp.data@dist.mat)), test.in)
		#cat(paste("KS Value (Goodness of fit) :", ks.val, "\n",sep=""));
		gmm.orig = data.frame(mu = gmm.first.2$mu, sigma = gmm.first.2$sigma, lambda = gmm.first.2$lambda);
		
		gmm.best <- .remove.gmm(gmm.first.2, z.limit, min.lambda);
		if (left.dist < 1 || left.dist > length(gmm.best$mu) -1)
		{
			cat("Gaussian Mixture Model Order Threshold provided (", left.dist, ") is outside the range. Defaulting to 1...\n\n");
			left.dist = 1;
		}
		thresh.hold <- .find_trough(gmm.best, -7, num =left.dist);
		cat(paste("Cutoff:", thresh.hold, "\n", sep=""));
		gmm.out <- data.frame(mu = gmm.best$mu, sigma = gmm.best$sigma, lambda = gmm.best$lambda);
		if (is.null(thresh.hold))
		{
			cat("Cannot find threshold... Not making clusters..\n\n");
			return(ggrasp.data);
		}
		clusters <- .get_tree_clusters(read.tree(text= ggrasp.data@phy), thresh.hold)
	
	}else
	{
		if (!missing("threshold") & !missing("num.clusters"))
		{
			cat("Both threshold and number of clades selected. Threshold will be used...\n\n")
			clusters <- .get_tree_clusters(read.tree(text = ggrasp.data@phy), threshold)
			thresh.hold = threshold; 
		}

		if (missing("threshold") & !missing("num.clusters"))
		{
			clusters = .get_tree_clusters_with_num(read.tree(text=ggrasp.data@phy), num.clusters);
			thresh.hold = max(sapply(1:max(clusters), function(x) {max(ggrasp.data@dist.mat[rownames(ggrasp.data@dist.mat) %in% names(clusters)[clusters==x], rownames(ggrasp.data@dist.mat) %in% names(clusters)[clusters==x]])})); 
		}
		if (!missing("threshold") & missing("num.clusters"))
		{
			clusters <- .get_tree_clusters(read.tree(text = ggrasp.data@phy), threshold)
			thresh.hold = threshold; 
		}
	}
	ggrasp.data@h <- thresh.hold;
	ggrasp.data@cluster <- clusters;
	ggrasp.data@medoids <- sapply((1:max(clusters)), .clust.medoid, ggrasp.data@dist.mat, clusters, ggrasp.data@rank)
	if (exists("gmm.out"))
	{
		ggrasp.data@gmm <- gmm.out;
	}
	if (exists("gmm.orig"))
	{
		ggrasp.data@gmm.orig <- gmm.orig;
	}	
	return(ggrasp.data)
}


#' @title ggrasp.addRanks 
#' @description adds a rank file to a GGRaSP object. If clusters have been defined, the medoids will be re-defined
#' 
#' @param x the GGRaSP object for which the ranks will be added.
#' @param rank.file string pointing to a file containing the ranks 
#'
#' @return A GGRaSP object where the ranks have been entirely redefined with the ranks in rank.file
#' 
#' @export

ggrasp.addRanks = function(x, rank.file)
{
	if (file.exists(rank.file))
	{
		rank.list <- read.table(rank.file, sep="\t");
		rank.lst = rank.list$V2;
		names(rank.lst) = rank.list$V1;
		if (sum(!names(rank.lst) %in% rownames(x@dist.mat))>0) 
		{
			cat("Names in provided rank file do not match... Discarding all those that do not match...\n\n");
			rank.lst <- rank.lst[names(rank.lst) %in% rownames(x@dist.mat)]
		}
		oth.names <- rownames(x@dist.mat)[!rownames(x@dist.mat) %in% names(rank.lst)];
		if (length(oth.names) > 0)
		{
			oth.ids <- rep(max(rank.lst)+1, length(oth.names)); 
			names(oth.ids) <- oth.names;
			rank.lst <- c(rank.lst, oth.ids)
		}
		x@rank = rank.lst;
		if (length(x@clusters)> 0)
		{
			x@medoids <- sapply((1:max(x@clusters)), .clust.medoid, x@dist.mat, x@clusters, x@rank)
	
		}
	}
	else
	{
		cat(paste("Cannot find file ", rank.file, ". Ignoring...\n", sep="")); 
	}
}

#' @title ggrasp.recluster 
#' @description recalculates a threshold and the resulting cluster using the previously defined Gaussian Mixture Model and provided threshold-determining factors. Requires the ggrasp.cluster to already have run 
#' 
#' @param x the GGRaSP object for which the ranks will be added.
#' @param z.limit All Gaussian distributions with means within this number of standard deviations will be reduced to only the larger distribution. Defaults to 1. Set to 0 to keep all non-overlapping distributions.
#' @param min.lambda All Gaussian distributions with lambda value (proportion of the total distribution) below this value are removed before calculating the threshold. Default is 0.005. Set to 0 to keep all. 
#' @param left.dist Number giving the number Gaussian distribution model immediately to the left of the threshold used. 1 is the default. Only value between 1 and k-1 where k is the total number of number of Gaussian distributions. 
#'
#' @return A GGRaSP object with the recalculated thresholds and the medoids using a previously generated GMM 
#' @examples #The following data is from Chavda et al 2016 which phylotyped Enterobacter genomes
#' # Our example uses the data underpinning the tree shown in Figure 2
#'
#' #Loading the tree 
#' library(ggrasp);
#' tree.file <- system.file("extdata", "Enter.kSNP.tree", package="ggrasp")
#' Enter.tree <- ggrasp.load(tree.file, file.format = "tree");
#'
#' #Clustering the tree using a threshold estimated by Gaussian Mixture Models (GMMs)
#' \donttest{Enter.tree.cluster <- ggrasp.cluster(Enter.tree)}
#'
#'
#' #Use print to get a list of the medoids selected
#' \donttest{print(Enter.tree.cluster)}
#'
#' #Re-clustering the tree using a threshold estimated by the GMMs but without the distribution
#' #cleaning (i.e. removing the overlapping and low count distributions)
#' \donttest{Enter.tree.reclust <- ggrasp.recluster(Enter.tree.cluster, z.limit=0, min.lambda = 0)}
#'
#
#' @export

ggrasp.recluster = function(x, z.limit=1, min.lambda=0.005, left.dist = 1)
{
	if (missing("x") || class(x) != "ggrasp")
	{
		cat("");
		return();
	}
	if (length(x@gmm.orig) == 0)
	{
		cat("No Gaussian Mixture Model Found for GGRaSP value. Running ggrasp.cluster using the provided variables\n\n");
		x = ggrasp.cluster(x, z.limit = z.limit, min.lambda = min.lambda, left.dist = left.dist);
	}
	else
	{
		gmm.best <- .remove.gmm(x@gmm.orig, z.limit, min.lambda);
		if (left.dist < 1 || left.dist > length(gmm.best$mu) -1)
		{
			cat("Gaussian Mixture Model Order Threshold provided (", left.dist, ") is outside the range. Defaulting to 1...\n\n");
			left.dist = 1;
		}
		thresh.hold <- .find_trough(gmm.best, -7, num =left.dist);
		cat(paste("Cutoff:", thresh.hold, "\n", sep=""));
		gmm.out <- data.frame(mu = gmm.best$mu, sigma = gmm.best$sigma, lambda = gmm.best$lambda);
		if (is.null(thresh.hold))
		{
			cat("Cannot find threshold... Not making clusters..\n\n");
			return(x);
		}
		clusters <- .get_tree_clusters(read.tree(text= x@phy), thresh.hold);
		x@h <- thresh.hold;
		x@cluster <- clusters;
		x@medoids <- sapply((1:max(clusters)), .clust.medoid, x@dist.mat, clusters, x@rank)
		if (exists("gmm.out"))
		{
			x@gmm <- gmm.out;
		}
	
	}
	return(x)
}	

#Interior function used to divide a phylogentic tree (in ape format) into
#subtrees using by cuting at the value given. Requires a rooted tree
.get_tree_clusters = function(nj_tree, thrsh)
{
  if (class(nj_tree) != "phylo")
  {
    if(class(nj_tree) == "hclust")
    {
      nj_tree = as.phylo(nj_tree);
    }
    else
    {
      cat("\n\nThis function only works with ape defined phylogenetic tree structures...")
      return(NULL);
    }
  }
  if (!is.binary.tree(nj_tree))
  {
    cat("\n\nThis function only works with binary trees...\n\nMaking a binary tree by randomly resolving tree...\n\n");
	nj_tree = multi2di(nj_tree);
   # return(NULL);
  }
  #tmp.clusters is the storage for the clusters. Each tip is assigned to a cluster
  tmp.clusters = 1:length(nj_tree$tip.label);
  #NA is used as the default unassigned state and each cluster name is the tip label
  tmp.clusters[] = NA;
  names(tmp.clusters) = nj_tree$tip.label;
  
  nj.mrca = mrca(nj_tree)
  #con.dist is the total distance along the edges between the tips (here the genes)
  con.dist <- cophenetic.phylo(nj_tree);
  #need a counting variable for the clusters
  clust.num = 1;
  #cycle through all the interior nodes of the phylogeny 
  for (i in min(nj_tree$edge[,1]):max(nj_tree$edge[,1]))
  {
    #getting all the children node...
    tmp.grp <- nj_tree$edge[nj_tree$edge[,1]==i,2];
    if (length(tmp.grp) == 3) 
    {
      #First case, node has 3 children (and therefore no ancestor) and so all the branches are children
      #id.grpN are the tips associated with each of the branches
      id.grp1 = rownames(nj.mrca)[rowSums(nj.mrca == tmp.grp[1])>0];
      id.grp2 = rownames(nj.mrca)[rowSums(nj.mrca == tmp.grp[2])>0];
      id.grp3 = rownames(nj.mrca)[rowSums(nj.mrca == tmp.grp[3])>0];
     }else
    {
      #First case, node has 2 children (and therefore 1 ancestor) and so the branches are children + ancestor
      id.grp1 = rownames(nj.mrca)[rowSums(nj.mrca == tmp.grp[1])>0];
      id.grp2 = rownames(nj.mrca)[rowSums(nj.mrca == tmp.grp[2])>0];
      id.grp3 = rownames(nj.mrca)[!rownames(nj.mrca) %in% c(id.grp1, id.grp2)];                                                   
    }
    #Making sure that there are tips 
    if (length(id.grp1) > 0 & length(id.grp2) > 0 & length(id.grp3) > 0)
    {
      #currently using maximum distance- maybe add in other methods to look at the distances
      #dist is the distances between each of the clusters
      dist.3 <- c(max(con.dist[rownames(con.dist) %in% id.grp1, colnames(con.dist) %in% id.grp2]),
                  max(con.dist[rownames(con.dist) %in% id.grp1, colnames(con.dist) %in% id.grp3]),
                  max(con.dist[rownames(con.dist) %in% id.grp2, colnames(con.dist) %in% id.grp3]));    
      #Only want where one branch vs branch comparison is less than the threshold- this means that this could be a cap
      if(sum(dist.3 < thrsh) ==1)
        {
          #First choose which branch pair we are looking at
          comb.num.1 = id.grp1; comb.num.2 = id.grp2;
          if (dist.3[2]< thrsh)
          {
            comb.num.2 = id.grp3;
          }
          if (dist.3[3]< thrsh)
          {
            comb.num.1 = id.grp3;
          }
          #Check to see if the tips are already part of a cluster
          tmp1 <- tmp.clusters[!is.na(tmp.clusters) & names(tmp.clusters) %in% comb.num.1];
          tmp2 <- tmp.clusters[!is.na(tmp.clusters) & names(tmp.clusters) %in% comb.num.2];
          #Option 1: they are not. Assign to the lowest blank number, incriment the counter
          if (length(tmp1) + length(tmp2) == 0)
          {
            tmp.clusters[names(tmp.clusters) %in% comb.num.1] = clust.num;
            tmp.clusters[names(tmp.clusters) %in% comb.num.2] = clust.num;
            clust.num = clust.num + 1;
          }
          else
          {
            #case 2: tip group 1 is already named; use it to name tip group 2
            if (length(unique(tmp1))==1 & length(tmp2) ==0)
            {
              tmp.clusters[names(tmp.clusters) %in% comb.num.2] = unique(tmp1)[1];
            }else
            {
              #case 2: tip group 2 is already named; use it to name tip group 1
              if (length(unique(tmp2))==1 & length(tmp1) ==0)
              {
                tmp.clusters[names(tmp.clusters) %in% comb.num.1] = unique(tmp2)[1];
              }else
              {
                #case 4: both are named and are collapse cluster to the lowest number
                tmp.clusters[names(tmp.clusters) %in% comb.num.1] = min(unique(tmp2)[1], unique(tmp1)[1]);
                tmp.clusters[names(tmp.clusters) %in% comb.num.2] = min(unique(tmp2)[1], unique(tmp1)[1]);
                
              }
            }
            
          }
       
      }    
      }
  }
  #Fill in all the singleton clusters
  tmp.clusters[is.na(tmp.clusters)]= clust.num -1 + (1:sum(is.na(tmp.clusters)))
  
  #Make sure that there are no empty cluster numbers
  tmp.ord <- order(unique(tmp.clusters))
  tmp.clusters2 = tmp.clusters;
  for (i in 1:length(tmp.ord))
  {
    tmp.clusters2[tmp.clusters==tmp.ord[i]] = i;
  }
  return(tmp.clusters2)
}


.clust.medoid = function(i, pdistmat, pclusters, rank.sc) {
  #Medoids are the cluster member that has the minimal distance 
  #Any rank trimming is done prior to this function
  ind = rownames(pdistmat) %in% names(pclusters)[(pclusters == i)];
  if (is.null(rank.sc) | length(rank.sc)==0)
  {
    rank.sc <- rep(1, nrow(pdistmat));
    names(rank.sc) <- rownames(pdistmat); 
  }
  rank.sc = rank.sc[rownames(pdistmat)];
  if (sum(rank.sc[ind] == min(rank.sc[ind]))==1)
  {
    return(names(rank.sc)[ind][rank.sc[ind] == min(rank.sc[ind])]);
  }
  if (sum(ind) <= 1) {
    return (rownames(pdistmat)[ind]) ## medoid of a single object is the object
  } else {
    return(names(which.min(rowSums( pdistmat[ind & rank.sc == min(rank.sc[ind]), ind] ))))
  }
}

#Returns the N clusters with the nodes farther
.get_tree_clusters_with_num = function(nj.tree, num_clusters)
{
	#current cluster count. Starts at one, goes to num_clusters provided
  curr.clust = 1;
  #gets the node that's the most recent common ancestor for each genome pair
  nj.mrca = mrca(nj.tree); 
  #Stores the cluster membership of each genome
  tmp.clusters = rep(1, length(nj.tree$tip.label));
  names(tmp.clusters) = nj.tree$tip.label;
  
  #Distance matrix for the tree
  d.mat = cophenetic.phylo(nj.tree);
  nj.node.depth <- rep(0, nj.tree$Nnode);
  
  tmp <- length(nj.tree$tip.label);
  #Takes advanges of the ape system for label 0:(# of Nodes) are the tip and the subsequent are the internal nodes
  for (i in (nj.tree$Nnode + tmp):(tmp+1))
  {
	x2 <- nj.tree$edge[nj.tree$edge[,1]==i,2];
	x3 <- nj.tree$edge.length[nj.tree$edge[,1]==i]
	if (x2[1] > tmp)
	{
		x3[1] = x3[1] + nj.node.depth[x2[1]-tmp];
	}
	if (x2[2] > tmp)
	{
		x3[2] = x3[2] + nj.node.depth[x2[2]-tmp];
	}
	nj.node.depth[i - tmp] = max(x3);
  }
  pos = 1;
  nj.node.height <- length(nj.tree$tip.label) + order(nj.node.depth, decreasing=T)
  curr.clust = curr.clust + 1;
  while(curr.clust <= num_clusters)
  {
    tmp.grp <- nj.tree$edge[nj.tree$edge[,1]==nj.node.height[pos],2];
    if (length(tmp.grp) > 1)
    {
      for (i in 2:length(tmp.grp))
      {
        id.grp = rownames(nj.mrca)[rowSums(nj.mrca == tmp.grp[i])>0];
        tmp.clusters[names(tmp.clusters) %in% id.grp] = curr.clust;
        curr.clust = curr.clust + 1;
      }
    }
    pos = pos + 1;
  }
  return(tmp.clusters);
}

#Makes the gaussian mixture models
.make_gmm <- function(mat, start.iter, max.iter = 10, run.type = "bgmm")
{
  #Make the matrix into a vector so it can be plotted/used to make a 
  n1 <- as.vector(mat);
  #i is the number of guassian distributions
  i = start.iter;
  cat(paste("Run with # ", i, " Gaussian Distribution\n", sep=""));
  if (run.type == "mixtools")
  {
	old <- normalmixEM(n1, k = i)
	}
	else
	{
		old <- unsupervised(n1, k = i)
	}
  #the previous gaussian mixture model
  orig <- old;
  #keep on looking at the number of gaussian mixture models until either:
  #1. the newest gaussian mixture model is not a better fit than the previous using log-likelihood test
  #2. the log likelihood is positive (this is to prevent over-fitting in practice)
  if (start.iter < max.iter)
  {
  for (i in (start.iter+1):max.iter)
  {
	cat(paste("Run with # ", i, " Gaussian Distribution\n", sep=""));
    
	if (run.type == "mixtools")
	{
		# mu = seq(min(n1), max(n1), (max(n1)-min(n1))/(i-1))
		new <- normalmixEM(n1, k = i, verb = F, maxrestarts=5);
		if(class(new) != "mixEM" )
		{
			cat("Normal Mixture Failed")
			break;
		}
		j1 <- new$loglik - old$loglik;
	}
	else
	{
		new <- unsupervised(n1, k = i)

		j1 <- new$likelihood - old$likelihood;
	}
	if (!is.null(j1))
	{
		j <- 1- pchisq(2 * (j1), 1);
		if (j > 0.05)
		{
			break;
		}
		old <- new;
	}
  }
  }
  #Need to add a grph of this here
  cat("Done\n");
  return(old)
  
}

#Removes any gmms that is either too close to other or insufficient numbers
#Does not try to combine them; just discards the lower one
.remove.gmm = function(gmm, z.limit = 1.0, min.lambda = 0.005)
{
  i = 1;
  gmm <- gmm[gmm$lambda >= min.lambda,];
  gmm.ord <- order(gmm$lambda, decreasing=T);
  n1 <- unlist(sapply(1:(length(gmm$mu)-1), function(x) { sapply((x+1):length(gmm$mu), function(y) {abs(gmm$mu[x]- gmm$mu[y])/max(gmm$sigma[x], gmm$sigma[y])})}));
  n1.1 <- unlist(sapply(1:(length(gmm$mu)-1), function(x) { sapply((x+1):length(gmm$mu), function(y) {x})}))
  n1.2 <- unlist(sapply(1:(length(gmm$mu)-1), function(x) { sapply((x+1):length(gmm$mu), function(y) {y})}))
  while(min(n1) < z.limit & length(gmm$mu) > 2)
  {
    n2 <- n1.1[n1==min(n1)];
    n3 <- n1.2[n1==min(n1)];
    if (gmm$lambda[n2] > gmm$lambda[n3])
    {
      gmm = gmm[-n3,];
    }
    else
    {
      gmm = gmm[-n2,];
    }
    n1 <- unlist(sapply(1:(length(gmm$mu)-1), function(x) { sapply((x+1):length(gmm$mu), function(y) {abs(gmm$mu[x]- gmm$mu[y])/max(gmm$sigma[x], gmm$sigma[y])})}));
    n1.1 <- unlist(sapply(1:(length(gmm$mu)-1), function(x) { sapply((x+1):length(gmm$mu), function(y) {x})}))
    n1.2 <- unlist(sapply(1:(length(gmm$mu)-1), function(x) { sapply((x+1):length(gmm$mu), function(y) {y})}))

  }
  return(gmm);  
}

#Uses modified newtonian method to get the minimal point between 2 gaussian distributions
.find_trough <- function(gmm, sigdig=-7, num=1)
{
  all.min = NULL;
  for (j in 2:length(gmm$mu))
  {
	n1 <-floor(log10(gmm$mu[order(gmm$mu)][j] - gmm$mu[order(gmm$mu)][j-1]))-1
	st = gmm$mu[order(gmm$mu)][j-1];
	en = gmm$mu[order(gmm$mu)][j];
	if (is.null(n1)) 
	{
		cat("Clustering Failed...\n\n");
		cat(paste(gmm$mu));
		cat("\n");
		return(NULL);
	}
	#use ne
	for (i in seq(n1, sigdig, by = -1))
	{
		d1 <- dnorm(seq(st,en,by=(10**i)), mean=gmm$mu[order(gmm$mu)][j-1], sd=gmm$sigma[order(gmm$mu)][j-1])*gmm$lambda[order(gmm$mu)][j-1]+dnorm(seq(st, en, by=(10**i)), mean=gmm$mu[order(gmm$mu)][j], sd=gmm$sigma[order(gmm$mu)][j])*gmm$lambda[order(gmm$mu)][j];
		nm1 <- order(d1)[1];
		old_st = st;
		st = max(old_st, old_st + (10**i) * (nm1-1))
		ret = old_st + (10**i) * nm1
		en = min(old_st + (10**i) * (nm1+1), en)
	}
	if (is.null(all.min))
	{
		all.min = ret;
	}else
	{
		all.min = c(all.min, ret);
	}
	}
  return(sort(all.min)[num]);
}
