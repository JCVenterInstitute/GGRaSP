# GGRaSP

GGRaSP (Gaussian Genome Representative Selector with Prioritization) is an R-package which can generate and return a reprentative set of genomes from a large group of genomes with a defined relationship. The reprentative set is select either using user-defined cutoff or cluster number values or else *de novo* calculates the clusters based on modeling the genome relationships with a Gaussian Mixture Model. The default value returned is the list of reprentative genomes, but the package also allows for multiple outputs including text files, plots, and trees. To allow for high-throughput analysis, we have included an Rscript file that can run GGRaSP from the command line (though it does require GGRaSP to be installed to the default R location).

## Installing GGRaSP

GGRaSP can be installed using the install_git command in devtools. However this will not install the command line script, which needs to be installed seperately as shown below. The example data is also pulled for demonstration purposes.
```
R
>library(devtools)
>library(git2r)
>install_git("https://github.com/JCVenterInstitute/GGRaSP.git", branch="master", credentials = cred_user_pass(username, password))
>quit()
```

if this is not working, it is also the options to download the ggrasp_1.0.tar.gz file from github and install ggrasp from this
```
git clone https://github.com/JCVenterInstitute/GGRaSP.git
R
>library(devtools)
>install("./GGRaSP/")
```

getting the individual files
```
git fetch https://github.com/JCVenterInstitute/GGRaSP/
From https://github.com/JCVenterInstitute/GGRaSP
 * branch            HEAD       -> FETCH_HEAD

git checkout -m FETCH_HEAD -- ggrasp.R

git checkout -m FETCH_HEAD -- examples/
```

## Using GGRaSP

GGRaSP can be used in two ways. The first uses the R-console to use the GGRaSP functions to load the genomes, cluster the genomes, and report the reprentative genomes. The second uses the command line Rscript. Below we will show the default version of running both with the example data from Chavdra et al 2016.

### Using GGRaSP in the R console

GGRaSP is centered around two primary functions to load and analyze the data, with three primary output variables that allow for analysis and reporting of the clustering. For more detailed descriptions of all the GGRaSP R functions, please examine the vigenettes.

To start with, simply load the library and genome-relationship file, here the Chavda ANI similarity matrix provided in the examples file. An offset of 100 is used to transform the simularity matrix to a distance matrix and uses a complete hclust (equivilant to an UPMGA) to make the phylogeny.
```
>library(ggrasp)
># The file is Enter.ANI.mat 
>enter.in.ggrasp <- ggrasp.load(system.file("extdata", "Enter.ANI.mat", package="ggrasp"), file.format="matrix", tree.method="complete", offset=100)
```
Now that the the simularity matrix is loaded in, cluster it using the default values, checkout the cutoff and number of distributions using the summary variable, and visualize the Gaussian Distributions:
```
>enter.cluster <- ggrasp.cluster(enter.in.ggrasp);
>enter.cluster
>plot(enter.cluster, "gmm")
```

### Running GGRaSP on the command line

GGRaSP also includes an Rscript program that can be run on the command line. Since GGRaSP does allow for multiple parameter changes, the script can take in multiple parameters for (i) input, (ii) clustering, and (iii) outputing. It also can run a simplified default version. For example, the same analysis as above (load in a similarty matrix, make a UPMGA tree, default GMM cluster, print the medoids and save the gmm plot) can be run below:
```
ggrasp.R -i ./examples/Enter.ANI.mat -d 100 -o ./enter.test.out --plotgmm
```

The output files will be: enter.test.out.gmm.pdf (gmm plot) and: enter.out.medoids.txt (medoid list).

A list of all the available flags for the scripts are as follows:

```
# Input Variables
		-i, --input <input file name for tab delimited distance matrix with row and column headers, newick file, or aligned multiple fasta file>
		-t, --type <type of the tree format to makes. Default is 'complete'>    
		-d, --offset <optional offset to convert a similarity matrix to a distance matrix>
		-m, --method  <optional method to use for hclust such as complete, single or average or nj to perform neighbor-joining>
		-v, --verbose <optional verbal threshold: 1 = no output, 2 = text only, 3 = graphics in pdf. Default = 3>
		-n, --names <optional name file to label the exterior of the tree. Multiple files can be included with a comma between>
		-r, --ranks <optional input file with ranks of the genes for medoid selection. Ranks should be from 1 (best) to n (worst) in two columns with the genome ID in the first>

# Clustering Variables
		-c, --clusters <optional number of clusters to generate>
		-h, --threshold <optional threshold to use to cluster the genomes>
		-k, --keep <optional flag to keep all reference clusters as sub-and show subclusters>
		-z, --zscore <optional value to combine all gaussian mixtures within this z-score. Default = 1. Use 0 to turn off>
		-x, --minlambda <optional value used to discard a Gaussian distribution if it has less than this percentage of the total distribution. Default is 0.005>
		-s, --start <optional value with the minimum number of gmm to examine. Default is 2>
		-e, --end <optional value with the maximum number of gmm to examine. Default is 10>
		-p, --packages <optional value with the package to use for the Gaussian mixture model estimate. Default is bgmm. Mixtools also can be used>

# Output Variables (medoid list is always printed)
		-o, --output <output file name for html table, pdf plots, cluster and txt files>
		-y, --plottype <optional flag used to provide the plot format>
		-l, --layout <optional phylogeny layout such as rectangular, circular, radial and slanted. Default = circular>
		-a, --metadata <optional input file containing metadata for the genomes to plot on the ring. It should be tab-deleniated>
      
		-0, --writetable <optional flag used to save the cluster table>
		-1, --writepseudo <optional flag used to save the psuedo fasta list>
		-2, --writetree <optional flag used to save the newick style full tree >
		-3, --writetrimtree <optional flag used to save the cluster table>
		-4, --writeitol <optional flag used to save an itol style file to view the cluster information>
		
		-5, --plothist <optional flag used to save the histogram with the plot>
		-6, --plotgmm <optional flag used to save the histogram with the Gaussian distribution overlayed>
		-7, --plottree <optional flag used to save the full tree plot>
		-8, --plotclust <optional flag used to save the full tree plot with the clusters highlighted>
		-9, --plottrim <optional flag used to save the trimmed tree with only the medoids>
		This script will select sample genome medoids by selecting a treshhold using Gaussian Mixture Models.
```

##Citing GGRaSP

When using GGRaSP in your analysis, please cite GGRaSP as follows:
```
	GGRaSP: a R-package for selecting representative genomes using Gaussian mixture models. Thomas H Clarke, Lauren M Brinkac, Granger Sutton, and Derrick E Fouts. Bioinformatics, bty300, https://doi.org/10.1093/bioinformatics/bty300


```