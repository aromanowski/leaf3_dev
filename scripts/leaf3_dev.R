############################################################################
# leaf3_dev.R                                                              #
# --------------                                                           #
# This script analyses L3 gene expression of A. thaliana plants grown      #
# under WL (LD 12 h : 12 h) or 2 EoD FR treatments (since 6 or 18 D.A.S.)  #
# at 3 developmental timepoints (13, 16 and 20 D.A.S.).                    #
# It was used to analyse the data described in Romanowski A. et al., 2021  #
#                                                                          #
############################################################################


### Notes
# Many lines of codes are commented. Uncomment them as needed for your usage

###############################################
#             Requires                        #
###############################################
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install()}

if(!require(ASpli)) BiocManager::install("ASpli")
if(!require(GenomicFeatures)) BiocManager::install("GenomicFeatures")
if(!require(AnnotationDbi)) BiocManager::install("AnnotationDbi")
if(!require(BiocParallel)) BiocManager::install("BiocParallel")

if(!require(devtools)) {
  install.packages("devtools")
  devtools::install_github("kassambara/ggpubr")}
if(!require(RColorBrewer)) install.packages("RcolorBrewer")
if(!require(gplots)) install.packages("gplots")
if(!require(dplyr)) install.packages("dplyr")
if(!require(tidyverse)) install.packages("tidyverse")

if(!require(org.At.tair.db)) BiocManager::install("org.At.tair.db")
if(!require(GO.db)) BiocManager::install("GO.db")
if(!require(GOstats)) install.packages("GOstats")

if(!require(pathview)) BiocManager::install("pathview")
if(!require(KEGGREST)) BiocManager::install("KEGGREST")
if(!require(clusterProfiler)) BiocManager::install("clusterProfiler")

###############################################
#             Includes                        #
###############################################
library(ASpli)
library(BiocParallel)
library(GenomicFeatures)
library(AnnotationDbi)
library(gplots)
library(ggpubr)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
# Other libraries are loaded in the required section

###############################################
#             Begin Analysis                  #
###############################################

# set base directory  ### modify accordingly ###
basedir <- "c:/leaf3_dev/comp"
setwd(basedir)

# load auxiliary functions
source("../scripts/common_functions.R")

# Make the file genome sqlite from the transcriptome, using AtRTD2 transcriptome genes.gtf
# This GTF was edited so that it had the correct chromosome naming
# genome<-makeTxDbFromGFF("/data/BioData/references/Arabidopsis_thaliana/AtRTD2/AtRTD2_edited.gtf", format = "gtf") 
# saveDb(genome, file="genomeAtRTD2.sqlite")
# features<-binGenome(genome)
# save(features, file="featuresAtRTD2.sqlite")

# Load the genome
# genome <- loadDb(file="genomeAtRTD2.sqlite")
# load(file="featuresAtRTD2.RData")

###########################################################
# Summary of AtRTD2 features:                             #
#                                                         #
# * Number of extracted Genes = 34212                     #
# * Number of extracted Exon Bins = 238662                #
# * Number of extracted intron bins = 178027              #
# * Number of extracted trascripts = 82190                #
# * Number of extracted junctions = 151944                #
# * Number of AS bins (not include external) = 41863      #
# * Number of AS bins (include external) = 41941          #
# * Classified as:                                        #
#   ES bins = 1686	(4%)                                  #
#   IR bins = 13033	(31%)                                 #
#   Alt5'ss bins = 4244	(10%)                             #
#   Alt3'ss bins = 7683	(18%)                             #
#   Multiple AS bins = 15217	(36%)                       #
#   classified as:                                        #
#                 ES bins = 1627	(11%)                   #
#                 IR bins = 5060	(33%)                   #
#                 Alt5'ss bins = 2941	(19%)               #
#           			Alt3'ss bins = 5001	(33%)               #
###########################################################

#############################################################
# Create a target file with description of the experiment   #
# Example:                                                  #
# sample	        bam	                        condition     #
# Pr_D13_WL_01	  01_Pr_D13_WL_102_val.bam	  ctrl          #
# Pr_D13_WL_02	  02_Pr_D13_WL_103_val.bam	  ctrl          #
# Pr_D13_FR06_01	03_Pr_D13_FR06_104_val.bam	Pr_D13_FR06   #
# Pr_D13_FR06_02	04_Pr_D13_FR06_106_val.bam	Pr_D13_FR06   #
# B_D16_WL_01	    05_B_D16_WL_207_val.bam     B_D16_WL      #
# B_D16_WL_02	    06_B_D16_WL_208_val.bam     B_D16_WL      #
# B_D16_FR06_01	  09_B_D16_FR06_211_val.bam	  B_D16_FR06    #
# B_D16_FR06_02	  10_B_D16_FR06_212_val.bam	  B_D16_FR06    #
# B_D20_WL_01	    15_B_D20_WL_119_val.bam     B_D20_WL      #
# B_D20_WL_02	    16_B_D20_WL_120_val.bam	    B_D20_WL      #
# B_D20_FR06_01	  19_B_D20_FR06_123_val.bam	  B_D20_FR06    #
# B_D20_FR06_02 	20_B_D20_FR06_124_val.bam	  B_D20_FR06    #
# B_D20_FR18_01	  25_B_D20_FR18_231_val.bam	  B_D20_FR18    #
# B_D20_FR18_02	  26_B_D20_FR18_232_val.bam	  B_D20_FR18    #
#############################################################

# Load the leaf3_dev targets and bam files
targets <- read.table("targets.txt", stringsAsFactors = FALSE, header = TRUE, row.names = 1)
# bam <- loadBAM(targets,cores=6) # load BAMS
# save(bam, file="bam.RData") # Save BAM data as an RData file
# load(file="bam.RData") # Load BAM data from .RData file

# Get count data and save it as an .Rdata file
# l is read length from the experiment. Here it is set to 150L
# max intron size is conserved in different plant families (Yan et al., 2013 - doi: 10.1007/s11427-013-4540-y)
# counts <- readCounts(features, bam, cores=7, readLength = 150, targets=targets, maxISize = 5000)
# save(counts, file="counts_leaf3_dev.Rdata") # Save raw counts data

##########################################
#         Leaf Dev (Raw reads)           #
##########################################

# Change to project directory
setwd(basedir)

# Load target files
targets <- read.table("targets.txt", stringsAsFactors = FALSE, header = TRUE, row.names = 1)

# Load Counts data
load(file ="counts_leaf3_dev.Rdata")

# Get raw counts and define groups for EdgeR
cg <- countsg(counts)
lev <- c("ctrl",
        "Pr_D13_FR06",
        "B_D16_WL",
        "B_D16_FR06",
        "B_D20_WL",
        "B_D20_FR06",
        "B_D20_FR18")
group <- factor(targets$condition, levels = lev)

# Create the design matrix
design <- model.matrix(~0+group)
colnames(design) <- levels(group)

# Filter and keep genes with rd > 0.05
# Genes that have very low counts across all the libraries should be removed prior to downstream analysis. This is justified on both biological
# and statistical grounds. From biological point of view, a gene must be expressed at some minimal level before it is likely to be translated 
# into a protein or to be considered biologically important. From a statistical point of view, genes with consistently low counts are very unlikely 
# be assessed as significantly DE because low counts do not provide enough statistical evidence for a reliable judgment to be made. Such genes can 
# therefore be removed from the analysis without any loss of information.
df<-filterByRd(cg, targets = targets, min = 0.05, type="any")

# Create DGEList element with raw counts
# This object is easy to use as it can be manipulated like an ordinary list in R, and it can also be subsetted like a matrix. The main components of
# a DGEList object are a matrix of read counts, sample information in the data.frame format and optional gene annotation. We enter the counts into a
# DGEList object using the function DGEList in edgeR:
y <- DGEList(counts= df[,setdiff(colnames(df), c("symbol", "locus_overlap", "gene_coordinates","start","end","length","effective_length"))] ,
             group=group)
colnames(y)=rownames(targets) # Set Column names
y$samples$lib.size = colSums(y$counts) # Recalculate library sizes

# Filter to remove low counts
# Genes that have very low counts across all the libraries should be removed prior to downstream 
# analysis. This is justified on both biological and statistical grounds. From biological point of
# view, a gene must be expressed at some minimal level before it is likely to be translated into a
# protein or to be considered biologically important. From a statistical point of view, genes with
# consistently low counts are very unlikely be assessed as significantly DE because low counts do
# not provide enough statistical evidence for a reliable judgement to be made. Such genes can
# therefore be removed from the analysis without any loss of information.
# 
# As a rule of thumb, we require that a gene have a count of at least 10 or so in at least some 
# libraries before it is considered to be expressed in the study. In practice, the filtering is
# actually based on count-per-million (CPM) values so as to avoid favoring genes that are expressed
# in larger libraries over those expressed in smaller libraries. In edgeR, the filtering can be 
# accomplished using a function that takes into account the library sizes and the experimental 
# design:
keep <- filterByExpr(y, design)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]


# Plot the library sizes
# The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names
par(mar=c(9,4.1,4.1,2.1))
barplot(y$samples$lib.size, names=colnames(y), las=2)
# Add a title to the plot
title("Barplot of library sizes")

# Save the library sizes plot
png(file="library_sizes.png",    # create PNG for the MDS        
    width = 8*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)        # smaller font size))
par(mar=c(9,4.1,4.1,2.1))
barplot(y$samples$lib.size, names=colnames(y), las=2)
# Add a title to the plot
title("Barplot of library sizes")
dev.off()


# Count data is not normally distributed, so if we want to examine the distributions of the raw counts we need to log the counts. 
# We'll use box plots to check the distribution of the read counts on the log2 scale. 
# We can use the cpm function to get log2 counts per million, which are corrected for the different library sizes. 
# The cpm function also adds a small offset to avoid taking log of zero.
logcounts <- cpm(y,log=TRUE)

# Create a histogram of the log transformed data to see if it follows a normal distribution
hist(logcounts)

png(file="histogram_logCPM.png",    # create PNG for the MDS        
    width = 8*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)        # smaller font size))
par(mar=c(9,4.1,4.1,2.1))
hist(logcounts)
dev.off()

# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

png(filename = (paste0("logCPM_Boxplot",".png")),    # create PNG for the MDS        
    width = 10*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)
par(mar=c(9,4.1,4.1,2.1))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")
dev.off()

# Correlation analysis of replicate samples
dir.create("correlation")
df_logcounts <- data.frame(logcounts)
for (i in (0:((length(colnames(df_logcounts))/2)-1)))
{
  print(2*i+1)
  print(colnames(df_logcounts)[2*i+1])
  print(2*i+1+1)
  print(colnames(df_logcounts)[2*i+1+1])
  print(" - ")
  gg <- ggscatter(data.frame(df_logcounts[,(2*i+1):(2*i+2)]),
                  x = colnames(df_logcounts)[2*i+1],
                  y = colnames(df_logcounts)[2*i+2],
                  color = "red",
                  add = "reg.line",
                  add.params = list(color = "black", fill = "lightgray"),
                  conf.int = TRUE,
                  cor.coef = TRUE,
                  cor.coeff.args = list(method = "pearson", label.sep = "\n"),
                  xlab = colnames(df_logcounts)[2*i+1],
                  ylab = colnames(df_logcounts)[2*i+2]
  )
  ggarrange(gg) %>% 
    ggexport(filename = paste0("correlation/", i, " - log2 ", colnames(df_logcounts)[2*i+1], " vs ", colnames(df_logcounts)[2*i+2], ".png"),
             pointsize = 12,
             height = 4.5 * 600,
             width = 4.5 * 600,
             res = 600)
}
dev.off()

# Create a Multi Dimensional Scaling Plot
# The RNA samples can be clustered in two dimensions using multi-dimensional scaling (MDS) plots. This is both an analysis step and a quality
# control step to explore the overall differences between the expression profiles of the different samples.
# In the MDS plot, the distance between each pair of samples can be interpreted as the leading log-fold change between the samples for the genes
# that best distinguish that pair of samples. By default, leading fold-change is defined as the root-mean-square of the largest 500 log2-fold changes
# between that pair of samples.

# We specify the option to let us plot only one plot
par(mfrow=c(1,1))
# Check number and order of the samples
levels(factor(targets$condition))
## Choose different colours for each sample type
col.cell <- c("violet", "purple", "blue", "cyan", "green", "dark gray","orange", "pink", "magenta", "red", "dark red", "brown", "gray", "black")[factor(targets$condition)]
data.frame(targets$condition,col.cell)

# MDS with sample type colouring
plotMDS(y,col=col.cell)

png(filename = (paste0("MDS_plot",".png")),    # create PNG for the MDS        
    width = 10*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)        # smaller font size)
plotMDS(y,col=col.cell)
title("MDS by Sample")
dev.off()

# Hierarchical clustering with heatmaps
# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
head(var_genes)
# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)
# Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

## Get some nice colours for the variable genes heatmap
mypalette <- brewer.pal(7,"Spectral")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for sample variable
col.cell <- c("violet", "purple", "blue", "cyan", "green", "dark gray","orange", "pink", "magenta", "red", "dark red", "brown", "gray", "black")[factor(targets$condition)]

# Plot the heatmap
library(gplots)
par(mar=c(9,4.1,4.1,2.1))
heatmap.2(highly_variable_lcpm, 
          col=rev(morecols(50)),
          trace="column", 
          main="Variability of genes across samples",
          ColSideColors=col.cell,scale="row",
          cexRow=1,
          cexCol=1,
          margins=c(7,6),
          srtCol=45,
          labRow=FALSE)

# Save the heatmap
png(file="Top500_var_genes.heatmap.png",    # create PNG for the MDS        
    width = 10*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)        # smaller font size))
par(mar=c(9,4.1,4.1,2.1))
heatmap.2(highly_variable_lcpm, 
          col=rev(morecols(50)),
          trace="column", 
          main="Variability of genes across samples",
          ColSideColors=col.cell,scale="row",
          cexRow=1,
          cexCol=1,
          margins=c(7,6),
          srtCol=45,
          labRow=FALSE)
dev.off()

########################################################################
#         Leaf Dev (Fitted values and normalized expression values     #
########################################################################

# Create design model to define comparison groups
# Linear modeling and differential expression analysis in edgeR requires a design matrix to be specified. The design matrix records which
# treatment conditions were applied to each samples, and it also defines how the experimental effects are parametrized in the linear models.
#
# Similar example from Limma (see user guide), with a timecourse:
# lev <- c("wt.0hr","wt.6hr","wt.24hr","mu.0hr","mu.6hr","mu.24hr")
# f <- factor(targets$Target, levels=lev)
# design <- model.matrix(~0+f)
# colnames(design) <- lev
# fit <- lmFit(eset, design)

lev <- c("ctrl",
        "Pr_D13_FR06",
        "B_D16_WL",
        "B_D16_FR06",
        "B_D20_WL",
        "B_D20_FR06",
        "B_D20_FR18")
group <- factor(targets$condition, levels = lev)
design <- model.matrix(~0+group, data = y$samples)
colnames(design) <- unique(targets$condition)
save(design, file = "design.RData")

# Calculate normalization factors by library
# Normalization by trimmed mean of M values (TMM) (Robinson and Oshlack 2010) is performed by using the calcNormFactors function, which returns
# the DGEList argument with only the  norm.factors changed. It calculates a set of normalization factors, one for each sample, to eliminate 
# composition biases between libraries. The product of these factors and the library sizes defines the effective library size, which replaces 
# the original library size in all downstream analyses.
y<-calcNormFactors(y)

# Check before and after TMM normalization effect
# The expression profiles of individual samples can be explored more closely with mean-difference (MD) plots. An MD plot visualizes the library
# size-adjusted log-fold change between two libraries (the difference) against the average log-expression across those libraries (the mean).
par(mfrow=c(2,2)) # plot 2 by 2 graphs in the same layout
plotMD(logcounts,column = 7) # before normalization
abline(h=0,col="blue")
plotMD(logcounts,column = 8) # before normalization
abline(h=0,col="blue")
plotMD(y,column = 7) # after TMM normalization
abline(h=0,col="red")
plotMD(y,column = 8) # after TMM normalization
abline(h=0,col="red")

# Normalization graphs (one per sample side by side)
pdf("blade_samples_before_and_after_TMM_normalization.pdf",
    paper = "a4r")
for (i in seq(from = 0, to = length(rownames(targets)), by=8)) {
  par(mfrow=c(2,4), mar=c(9,4.1,4.1,2.1))
  
  print(i)
  for (j in seq(c(0:7))) {
    if ((i+j) <= length(rownames(targets))) {  
      print(paste0("j =", j, " and i=", i+j))
      plotMD(logcounts,column = i+j) # before normalization
      abline(h=0,col="blue")
      plotMD(y,column = i+j) # after TMM normalization
      abline(h=0,col="red")
    }
  }
  if (i > length(rownames(targets))) { i = length(rownames(targets))}
  
}
dev.off()

# Estimate dispersions:
# EdgeR estimates an empirical Bayes moderated dispersion for each individual gene. It also estimates a common dispersion, which is a
# global dispersion estimate averaged over all genes, and a trended dispersion where the dispersion of a gene is predicted from its abundance
y<-estimateDisp(y, design = design, robust = TRUE)
save(y, file = "y.RData")

# This returns a DGEList object with additional components (common.dispersion,  trended.dispersion and tagwise.dispersion) added to hold the 
# estimated dispersions. Here robust=TRUE has been used to protect the empirical Bayes estimates against the possibility of outlier genes with 
# exceptionally large or small individual dispersions (Phipson et al. 2016).

# Plot the dispersions:
# The vertical axis of the plotBCV plot shows square-root dispersion, also known as biological coefficient of variation (BCV) 
# (McCarthy, Chen, and Smyth 2012). For RNA-seq studies, the NB dispersions tend to be higher for genes with very low counts. 
# The dispersion trend tends to decrease smoothly with abundance and to asymptotic to a constant value for genes with larger counts.
par(mfrow=c(1,1)) # Plot only one graph
plotBCV(y) # Plot the dispersion of the data

# Save the dispersion plot
png(file="blade_data_dispersions.png",    # create PNG for the MDS        
    width = 9*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)        # smaller font size))
plotBCV(y) # Plot the dispersion of the data
dev.off()

# Calculate fitted.values counts (pseudo counts) and perform ANOVA like function (Quasi Likelihood F-test)

# The NB model can be extended with quasi-likelihood (QL) methods to account for gene-specific variability from both biological and technical
# sources (Lund et al. 2012; Lun, Chen, and Smyth 2016). Under the QL framework, the NB dispersion trend is used to describe the overall biological 
# variability across all genes, and gene-specific variability above and below the overall level is picked up by the QL dispersion. In the QL approach,
# the individual (tagwise) NB dispersions are not used.
# The estimation of QL dispersions is performed using the glmQLFit function.
# Setting robust=TRUE in glmQLFit is usually recommended (Phipson et al. 2016). This allows gene-specific prior df estimates, with lower values for
# outlier genes and higher values for the main body of genes. This reduces the chance of getting false positives from genes with extremely high or low
# raw dispersions, while at the same time increasing statistical power to detect differential expression for the main body of genes.
fit <- glmQLFit(y, design = design, robust = TRUE)
save(fit, file = "fit.RData")

# This returns a DGEGLM object with the estimated values of the GLM coefficients for each gene. It also contains a number of empirical Bayes (EB)
# statistics including the QL dispersion trend, the squeezed QL dispersion estimates and the prior degrees of freedom (df). The QL dispersions can be
# visualized by plotQLDisp function.
plotQLDisp(fit)


# Save the QL dispersion plot
png(file="QL_dispersions.png",    # create PNG for the MDS        
    width = 9*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)        # smaller font size))
plotQLDisp(fit) # Plot the QL dispersion of the data
dev.off()

# We use QL F-tests instead of the more usual likelihood ratio tests (LRT) as they give stricter error rate control by accounting for the uncertainty
# in dispersion estimation
qlf <- glmQLFTest(fit)

fdr<-p.adjust(qlf$table$PValue, method="BH") # Calculate FDR
tt<-cbind(qlf$table, fdr) 
final<-list(qlf, tt)
names(final)=c("qlf", "full")

# Genes that pass read density filter and likelihood test
dim(final$qlf)
genes_qlf <- dim(final$qlf)[1]

# Calculate normalized logCPM data
logcounts2 <- cpm(y, normalized.lib.sizes = TRUE, log=TRUE)

# Create a histogram to see if the transformed data follows a normal distribution
hist(logcounts2)

# save the plot
png(file="histogram_logCPM_normalised.png",    # create PNG for the MDS        
    width = 9*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)        # smaller font size))
hist(logcounts2)
dev.off()


# Calculate normalized counts (cpm)
norm_counts <- cpm(y, normalized.lib.sizes = TRUE) 

# Plot some example genes
par(mfrow=c(2,2), mar=c(9,4.1,4.1,2.1))
barplot(norm_counts["AT4G16780", ], las = 2, ylab = "Normalised counts (CPM)") # Example value
title("ATHB2 - AT4G16780")
barplot(norm_counts["AT2G46970", ], las = 2, ylab = "Normalised counts (CPM)") # Example value
title("PIL1 - AT2G46970")
barplot(norm_counts["AT1G02340", ], las = 2, ylab = "Normalised counts (CPM)") # Example value
title("HFR1 - AT1G02340")
barplot(norm_counts["AT1G70560", ], las = 2, ylab = "Normalised counts (CPM)") # Example value
title("TAA1 - AT1G70560")

# Plot example shade response genes
png(file="gene_controls_shade_responders.png",    # create PNG for the MDS        
    width = 9*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)        # smaller font size))
par(mfrow=c(2,2), mar=c(9,4.1,4.1,2.1))
barplot(norm_counts["AT4G16780", ], las = 2, ylab = "Normalised counts (CPM)") # Example value
title("ATHB2 - AT4G16780")
barplot(norm_counts["AT2G46970", ], las = 2, ylab = "Normalised counts (CPM)") # Example value
title("PIL1 - AT2G46970")
barplot(norm_counts["AT1G02340", ], las = 2, ylab = "Normalised counts (CPM)") # Example value
title("HFR1 - AT1G02340")
barplot(norm_counts["AT1G70560", ], las = 2, ylab = "Normalised counts (CPM)") # Example value
title("TAA1 - AT1G70560")
dev.off()

# Check distributions of samples using boxplots
par(mfrow = c(1,1))
boxplot(logcounts2, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts2),col="blue")
title("Boxplots of logCPMs (normalised)")

png(filename = (paste0("logCPM_normalized_Boxplot",".png")),    # create PNG for the MDS        
    width = 10*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)
par(mar=c(9,4.1,4.1,2.1))
boxplot(logcounts2, xlab="", ylab="Log2 counts per million",las=2)
# Add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts2),col="blue")
title("Boxplots of logCPMs (normalised)")
dev.off()

# Correlation analysis of replicate samples
library(ggpubr)
dir.create("correlation")
df_logcounts2 <- data.frame(logcounts2)
for (i in (0:((length(colnames(df_logcounts2))/2)-1)))
{
  print(2*i+1)
  print(colnames(df_logcounts2)[2*i+1])
  print(2*i+1+1)
  print(colnames(df_logcounts2)[2*i+1+1])
  print(" - ")
  gg <- ggscatter(data.frame(df_logcounts2[,(2*i+1):(2*i+2)]),
                  x = colnames(df_logcounts2)[2*i+1],
                  y = colnames(df_logcounts2)[2*i+2],
                  color = "red",
                  add = "reg.line",
                  add.params = list(color = "black", fill = "lightgray"),
                  conf.int = TRUE,
                  cor.coef = TRUE,
                  cor.coeff.args = list(method = "pearson", label.sep = "\n"),
                  xlab = colnames(df_logcounts2)[2*i+1],
                  ylab = colnames(df_logcounts2)[2*i+2]
  )
  ggarrange(gg) %>% 
    ggexport(filename = paste0("correlation/", i, " - normalised log2 ", colnames(df_logcounts2)[2*i+1], " vs ", colnames(df_logcounts2)[2*i+2], ".png"),
             pointsize = 12,
             height = 4.5 * 600,
             width = 4.5 * 600,
             res = 600)
}
dev.off()


# Create a Multi Dimensional Scaling Plot of the normalized data
## Let's choose colours for the different samples
col.cell <- c("violet", "purple", "blue", "cyan", "green", "dark gray","orange", "pink", "magenta", "red", "dark red", "brown", "gray", "black")[factor(targets$condition)]
data.frame(targets$condition,col.cell)
plotMDS(y, col=col.cell) # Plot MDS

# Save the MDS plot to a png file
png(filename = (paste0("MDS_plot_normalized_data",".png")),    # create PNG for the MDS        
    width = 10*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)        # smaller font size)
# We specify the option to let us plot only one plot
par(mfrow=c(1,1))
plotMDS(y,col=col.cell)
title("MDS by Sample (normalized data)")
dev.off()

# Hierarchical clustering with heatmaps
# We estimate the variance for each row in the logcounts matrix (normalized values)
var_genes <- apply(logcounts2, 1, var)
head(var_genes)
# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)
# Subset logcounts2 matrix
highly_variable_lcpm <- logcounts2[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

## Get some nicer colours
library(RColorBrewer)
mypalette <- brewer.pal(7,"Spectral")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for sample variable
col.cell <- c("violet", "purple", "blue", "cyan", "green", "dark gray","orange", "pink", "magenta", "red", "dark red", "brown", "gray", "black")[factor(targets$condition)]

# Plot the heatmap
library(gplots)
par(mar=c(9,4.1,4.1,2.1))
heatmap.2(highly_variable_lcpm, 
          col=rev(morecols(50)),
          trace="column", 
          main="Variability of genes across samples (normalized data)",
          ColSideColors=col.cell,scale="row",
          cexRow=1,
          cexCol=1,
          margins=c(7,6),
          srtCol=45,
          labRow=FALSE)

# Save the heatmap
png(file="Top500_var_genes_normalized_data.heatmap.png",    # create PNG for the MDS        
    width = 10*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)        # smaller font size))
par(mar=c(9,4.1,4.1,2.1))
heatmap.2(highly_variable_lcpm, 
          col=rev(morecols(50)),
          trace="column", 
          main="Variability of genes across samples (normalized data)",
          ColSideColors=col.cell,scale="row",
          cexRow=1,
          cexCol=1,
          margins=c(7,6),
          srtCol=45,
          labRow=FALSE)
dev.off()

# Save all count data
write.table(qlf$fitted.values, file="blade.genes.fitted.values.tab", sep="\t", row.names = TRUE, col.names = NA, na = "NA")
write.table(cg, file="blade.genes.raw.values.tab", sep="\t", row.names = TRUE, col.names = NA, na = "NA")
write.table(df, file="blade.genes.filtered.values.tab", sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
write.table(logcounts, file="blade.genes.logCPM.unnormalised.values.tab", sep="\t", row.names = TRUE, col.names = NA, na = "NA")
write.table(logcounts2, file="blade.genes.logCPM.normalised.values.tab", sep="\t", row.names = TRUE, col.names = NA, na = "NA")
write.table(norm_counts, file="blade.genes.CPM.normalised.values.tab", sep="\t", row.names = TRUE, col.names = NA, na = "NA")

# Add annotation data to normalised counts and logCPM data
md <- read.csv(file = "araport11_metadata_17-11-2018_and_curator_summary_31-12-2018.txt", na.strings=c("", "NA"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
logcounts2.annot <- getAnnnot(logcounts2, md)
norm_counts.annot <- getAnnnot(norm_counts, md)
write.table(logcounts2.annot, file="annot.blade.genes.logCPM.normalised.values.tab", sep="\t", row.names = NA, col.names = TRUE, na = "NA")
write.table(norm_counts.annot, file="annot.blade.genes.CPM.normalised.values.tab", sep="\t", row.names = NA, col.names = TRUE, na = "NA")

# Rearrange normalised counts (CPM) to make the data more readable and user friendly. Then save the rearranged CPM data.
df1 <- norm_counts
colnames(df1) <- rownames(targets)

df2 <- cbind(df1[,1], #Pr_D13_WL
             df1[,2],
             df1[,5], #B_D16_WL
             df1[,6],
             df1[,9], #B_D20_WL
             df1[,10],
             
             df1[,3], #Pr_D13_FR06
             df1[,4],
             df1[,7], #B_D16_FR06
             df1[,8],
             df1[,11], #B_D20_FR06
             df1[,12],
             
             df1[,13], #B_D20_FR18
             df1[,14])

rownames(df2) <- rownames(df1)
colnames(df2) <- c(colnames(df1)[1], #Pr_D13_WL
                   colnames(df1)[2],
                   colnames(df1)[5], #B_D16_WL
                   colnames(df1)[6],
                   colnames(df1)[9], #B_D20_WL
                   colnames(df1)[10],
                   
                   colnames(df1)[3], #Pr_D13_FR06
                   colnames(df1)[4],
                   colnames(df1)[7], #B_D16_FR06
                   colnames(df1)[8],
                   colnames(df1)[11], #B_D20_FR06
                   colnames(df1)[12],
                   
                   colnames(df1)[13], #B_D20_FR18
                   colnames(df1)[14])
write.table(df2, file="blade.genes.CPM.normalised.values.rearranged.tab", sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")

# Rearrange log normalised counts (logCPM) to make the data more readable and user friendly. Then save the rearranged logCPM data.
logcounts3 <- cbind(logcounts2[,1], #Pr_D13_WL
                    logcounts2[,2],
                    logcounts2[,5], #B_D16_WL
                    logcounts2[,6],
                    logcounts2[,9], #B_D20_WL
                    logcounts2[,10],
                    
                    logcounts2[,3], #Pr_D13_FR06
                    logcounts2[,4],
                    logcounts2[,7], #B_D16_FR06
                    logcounts2[,8],
                    logcounts2[,11], #B_D20_FR06
                    logcounts2[,12],
                    
                    logcounts2[,13], #B_D20_FR18
                    logcounts2[,14])

rownames(logcounts3) <- rownames(logcounts2)
colnames(logcounts3) <- c(colnames(logcounts2)[1], #Pr_D13_WL
                          colnames(logcounts2)[2],
                          colnames(logcounts2)[5], #B_D16_WL
                          colnames(logcounts2)[6],
                          colnames(logcounts2)[9], #B_D20_WL
                          colnames(logcounts2)[10],
                          
                          colnames(logcounts2)[3], #Pr_D13_FR06
                          colnames(logcounts2)[4],
                          colnames(logcounts2)[7], #B_D16_FR06
                          colnames(logcounts2)[8],
                          colnames(logcounts2)[11], #B_D20_FR06
                          colnames(logcounts2)[12],
                          
                          colnames(logcounts2)[13], #B_D20_FR18
                          colnames(logcounts2)[14])
write.table(logcounts3, file="blade.genes.logCPM.normalised.values.rearranged.tab", sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")

# Add annotations to rearranged data
df2.annot <- getAnnnot(df2, md)
logcounts3.annot <- getAnnnot(logcounts3, md)
write.table(df2.annot, file="annot.blade.genes.CPM.normalised.values.rearranged.tab", sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
write.table(logcounts3.annot, file="annot.blade.genes.logCPM.normalised.values.rearranged.tab", sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")

#################################
#       Sample comparisons      #
#################################

# The QL framework provides more accurate type I error rate control than an exact test, as it accounts for the uncertainty of the dispersion estimates.
# In contrast, the exact test assumes that the estimated dispersion is the true value, which can result in some inaccuracy. (The "exact" refers to the 
# fact that the p-value is calculated exactly rather than relying on approximations; however, this only results in exact type I error control when the 
# estimated and true dispersions are the same.) For this reason, I prefer using the QL methods whenever I apply edgeR.
# 
# The QL methods (and GLM methods) are also more flexible with respect to the experimental design. For example, if you got a second batch of samples, 
# all you would need to do in a GLM framework would be to change the design matrix, while the exact test methods can't handle an extra blocking factor
# for the batch.
#
# In summary, while both of the methods will work for your data set, the QL F-test is probably the better choice. There are some situations where the 
# QL F-test doesn't work well - for example, if you don't have replicates, you'd have to supply a fixed dispersion, which defeats the whole point of 
# modelling estimation uncertainty. Another situation is where the dispersions are very large and the counts are very small, whereby some of the 
# approximations in the QL framework seem to fail. In such cases, I usually switch to the LRT rather than using the exact test, for the reasons of 
# experimental flexibility that I mentioned above.
# Source: https://support.bioconductor.org/p/84291/

setwd(basedir)
dir.create("DE")
setwd(paste0(basedir,"/DE", collapse = NULL))
md <- read.csv(file = "../araport11_metadata_17-11-2018_and_curator_summary_31-12-2018.txt", na.strings=c("", "NA"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
md <- md[,!names(md) %in% c("Curator_summary_.31.12.2018.","Description_.17.11.2018.")]
#load("/data/projects/leaf_dev/comp/blade/fit.RData")

# Set cut-off values
p.value <- 0.05
fdr.gen <- 0.1

# Test for DE genes of primordia WL (ctrl) vs FR06
et <- glmQLFTest(fit, contrast=makeContrasts("Pr_D13_FR06-ctrl",levels=design))
topTags(et)
fdr.gen <- p.adjust(et$table$PValue, method="BH")
et_merge<-cbind(et$table[,setdiff(colnames(et$table),c("logCPM"))], fdr.gen)
# Get A. thaliana best hit and its description for each gene
# First column of metadata file must contain Gene IDs
et_merge <- getAnnnot(et_merge, md = md)
et_merge <- et_merge[order(et_merge$logFC, et_merge$PValue),] # Oder by logFC and p value
top_DE<-which(et_merge$PValue < p.value & et_merge$fdr.gen < fdr.gen)
length(top_DE)
detags <- rownames(y)[top_DE]
png(file="smear_00_Pr_D13_FR06_vs_WL.png",    # create PNG for the MDS        
    width = 5*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)        # smaller font size))
plotSmear(et, de.tags=detags)
title("Primordia Day 13 FR06 vs WL")
dev.off()
write.table(et_merge[top_DE,], file = "00-primordiaDE.tab", sep = "\t", row.names = FALSE)
write.table(et_merge, file = "00-primordiaDE_full.tab", sep = "\t", row.names = FALSE)

# Test for DE genes of blade d16 WL (ctrl) vs FR06
rm(et, et_merge, fdr.gen, top_DE)
et <- glmQLFTest(fit, contrast=makeContrasts("B_D16_FR06-B_D16_WL",levels=design))
topTags(et)
fdr.gen <- p.adjust(et$table$PValue, method="BH")
et_merge<-cbind(et$table[,setdiff(colnames(et$table),c("logCPM"))], fdr.gen)
# Get A. thaliana best hit and its description for each gene
# First column of metadata file must contain Gene IDs
et_merge <- getAnnnot(et_merge, md = md)
et_merge <- et_merge[order(et_merge$logFC, et_merge$PValue),] # Oder by logFC and p value
top_DE<-which(et_merge$PValue < p.value & et_merge$fdr.gen < fdr.gen)
length(top_DE)
detags <- rownames(y)[top_DE]
png(file="smear_01_B_D16_FR06_vs_WL.png",    # create PNG for the MDS        
    width = 5*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)        # smaller font size))
plotSmear(et, de.tags=detags)
title("Blade Day 16 FR06 vs WL")
dev.off()
write.table(et_merge[top_DE,], file = "01-blade_d16_FR06_vs_WLDE.tab", sep = "\t", row.names = FALSE)
write.table(et_merge, file = "01-blade_d16_FR06_vs_WLDE_full.tab", sep = "\t", row.names = FALSE)

# Test for DE genes of blade d20 WL (ctrl) vs FR06
rm(et, et_merge, fdr.gen, top_DE)
et <- glmQLFTest(fit, contrast=makeContrasts("B_D20_FR06-B_D20_WL",levels=design))
topTags(et)
fdr.gen <- p.adjust(et$table$PValue, method="BH")
et_merge<-cbind(et$table[,setdiff(colnames(et$table),c("logCPM"))], fdr.gen)
# Get A. thaliana best hit and its description for each gene
# First column of metadata file must contain Gene IDs
et_merge <- getAnnnot(et_merge, md = md)
et_merge <- et_merge[order(et_merge$logFC, et_merge$PValue),] # Oder by logFC and p value
top_DE<-which(et_merge$PValue < p.value & et_merge$fdr.gen < fdr.gen)
length(top_DE)
detags <- rownames(y)[top_DE]
png(file="smear_02_B_D20_FR06_vs_WL.png",    # create PNG for the MDS        
    width = 5*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)        # smaller font size))
plotSmear(et, de.tags=detags)
title("Blade Day 20 FR06 vs WL")
dev.off()
write.table(et_merge[top_DE,], file = "02-blade_d20_FR06_vs_WLDE.tab", sep = "\t", row.names = FALSE)
write.table(et_merge, file = "02-blade_d20_FR06_vs_WLDE_full.tab", sep = "\t", row.names = FALSE)

# Test for DE genes of blade d20 WL (ctrl) vs FR18
rm(et, et_merge, fdr.gen, top_DE)
et <- glmQLFTest(fit, contrast=makeContrasts("B_D20_FR18-B_D20_WL",levels=design))
topTags(et)
fdr.gen <- p.adjust(et$table$PValue, method="BH")
et_merge<-cbind(et$table[,setdiff(colnames(et$table),c("logCPM"))], fdr.gen)
# Get A. thaliana best hit and its description for each gene
# First column of metadata file must contain Gene IDs
et_merge <- getAnnnot(et_merge, md = md)
et_merge <- et_merge[order(et_merge$logFC, et_merge$PValue),] # Oder by logFC and p value
top_DE<-which(et_merge$PValue < p.value & et_merge$fdr.gen < fdr.gen)
length(top_DE)
detags <- rownames(y)[top_DE]
png(file="smear_03_B_D20_FR18_vs_WL.png",    # create PNG for the MDS        
    width = 5*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)        # smaller font size))
plotSmear(et, de.tags=detags)
title("Blade Day 20 FR18 vs WL")
dev.off()
write.table(et_merge[top_DE,], file = "03-blade_d20_FR18_vs_WLDE.tab", sep = "\t", row.names = FALSE)
write.table(et_merge, file = "03-blade_d20_FR18_vs_WLDE_full.tab", sep = "\t", row.names = FALSE)

# Test for DE genes of blade d20 FR06 (ctrl) vs FR18
rm(et, et_merge, fdr.gen, top_DE)
et <- glmQLFTest(fit, contrast=makeContrasts("B_D20_FR18-B_D20_FR06",levels=design))
topTags(et)
fdr.gen <- p.adjust(et$table$PValue, method="BH")
et_merge<-cbind(et$table[,setdiff(colnames(et$table),c("logCPM"))], fdr.gen)
# Get A. thaliana best hit and its description for each gene
# First column of metadata file must contain Gene IDs
et_merge <- getAnnnot(et_merge, md = md)
et_merge <- et_merge[order(et_merge$logFC, et_merge$PValue),] # Oder by logFC and p value
top_DE<-which(et_merge$PValue < p.value & et_merge$fdr.gen < fdr.gen)
length(top_DE)
detags <- rownames(y)[top_DE]
png(file="smear_04_B_D20_FR18_vs_FR06.png",    # create PNG for the MDS        
    width = 5*600,        # 10 x 600 pixels
    height = 5*600,
    res = 600,            # 300 pixels per inch
    pointsize = 8)        # smaller font size))
plotSmear(et, de.tags=detags)
title("Blade Day 20 FR18 vs FR06")
dev.off()
write.table(et_merge[top_DE,], file = "04-blade_d20_FR18_vs_FR06DE.tab", sep = "\t", row.names = FALSE)
write.table(et_merge, file = "04-blade_d20_FR18_vs_FR06DE_full.tab", sep = "\t", row.names = FALSE)

#################################
#           GO-Terms            #
#################################

## Requires
if(!require(org.At.tair.db)) BiocManager::install("org.At.tair.db")
if(!require(GO.db)) BiocManager::install("GO.db")
if(!require(GOstats)) BiocManager::install("GOstats")

## Includes
library(org.At.tair.db)
library(GO.db)
library(GOstats) # Requires RSQLite v2.1.0 to work correctly in R3.6.3

## Start GO analysis

# Add project name
project <- paste0(basedir,"/DE/")
sub_project = "GO"

setwd(basedir)
dir.create(project)
dir.create(paste0(project, "/", sub_project))

set.seed(12345)
universe <- read.csv2(file = paste0(basedir,"/blade.genes.CPM.normalised.values.tab"),
                      sep = "\t",
                      header = TRUE,
                      dec = ".")
colnames(universe)[1] <- "Gene_ID"
rownames(universe) <- universe$Gene_ID
universe <- rownames(universe)

# Set cut-off values
logFC <- 0.58
fdr.gen <- 0.05

# Get GO terms for 00-primordiaDE
project <- paste0(basedir,"/DE/")
sub_project = "GO"
all_genes <- read_tsv(file = paste0(basedir,"/DE/00-primordiaDE.tab"))
index_ALL <- which(abs(all_genes$logFC) > logFC & all_genes$fdr.gen < fdr.gen) #3046 - 2410
index_UP <- which(all_genes$logFC > logFC & all_genes$fdr.gen < fdr.gen) #1190 - 923
index_DOWN <- which(all_genes$logFC < -logFC & all_genes$fdr.gen < fdr.gen) #1856 - 1487
dir.create(paste0(project, "/", sub_project, "/00-primordiaDE"))
universe <- universe
PValue = 1 # cut-off set to 1 to retrieve the full list of GO terms
fdr = 1 # cut-off set to 1 to retrieve the full list of GO terms
annotation = "org.At.tair"
sub_project = "GO/00-primordiaDE"
all = TRUE
subset <- all_genes$Gene_ID[index_UP]
description = "primordia_d13_FR06_vs_WL_UP"
writeGSEA(subset, universe, PValue = PValue, fdr, annotation, ontology = "BP", description, project, sub_project, all)
writeGSEA(subset, universe, PValue = PValue, fdr, annotation, ontology = "MF", description, project, sub_project, all)
writeGSEA(subset, universe, PValue = PValue, fdr, annotation, ontology = "CC", description, project, sub_project, all)
subset <- all_genes$Gene_ID[index_DOWN]
description = "primordia_d13_FR06_vs_WL_DOWN"
writeGSEA(subset, universe, PValue = PValue, fdr, annotation, ontology = "BP", description, project, sub_project, all)
writeGSEA(subset, universe, PValue = PValue, fdr, annotation, ontology = "MF", description, project, sub_project, all)
writeGSEA(subset, universe, PValue = PValue, fdr, annotation, ontology = "CC", description, project, sub_project, all)

# Get GO terms for 01-blade_d16_FR06_vs_WLDE
setwd(basedir)
project <- paste0(basedir,"/DE/")
sub_project = "GO"
all_genes <- read_tsv(file = paste0(basedir,"/DE/01-blade_d16_FR06_vs_WLDE.tab"))
index_ALL <- which(abs(all_genes$logFC) > logFC & all_genes$fdr.gen < fdr.gen) #2069 - 1570
index_UP <- which(all_genes$logFC > logFC & all_genes$fdr.gen < fdr.gen) #824 - 608
index_DOWN <- which(all_genes$logFC < -logFC & all_genes$fdr.gen < fdr.gen) #1245 - 962
dir.create(paste0(project, "/", sub_project, "/01-blade_d16_FR06_vs_WLDE"))
universe <- universe
PValue = 1
fdr = 1
annotation = "org.At.tair"
project <- paste0(basedir,"/DE/")
all = TRUE
sub_project = "GO/01-blade_d16_FR06_vs_WLDE"
subset <- all_genes$Gene_ID[index_UP]
description = "blade_d16_FR06_vs_WLDE_UP"
writeGSEA(subset, universe, PValue = 1, fdr, annotation, ontology = "BP", description, project, sub_project, all)
writeGSEA(subset, universe, PValue = 1, fdr, annotation, ontology = "MF", description, project, sub_project, all)
writeGSEA(subset, universe, PValue = 1, fdr, annotation, ontology = "CC", description, project, sub_project, all)
subset <- all_genes$Gene_ID[index_DOWN]
description = "blade_d16_FR06_vs_WLDE_DOWN"
writeGSEA(subset, universe, PValue = 1, fdr, annotation, ontology = "BP", description, project, sub_project, all)
writeGSEA(subset, universe, PValue = 1, fdr, annotation, ontology = "MF", description, project, sub_project, all)
writeGSEA(subset, universe, PValue = 1, fdr, annotation, ontology = "CC", description, project, sub_project, all)

# Get GO terms for 02-blade_d20_FR06_vs_WLDE
setwd(basedir)
project <- paste0(basedir,"/DE/")
sub_project = "GO"
all_genes <- read_tsv(file = paste0(basedir,"/DE/02-blade_d20_FR06_vs_WLDE.tab"))
index_ALL <- which(abs(all_genes$logFC) > logFC & all_genes$fdr.gen < fdr.gen) #2529 - 1962
index_UP <- which(all_genes$logFC > logFC & all_genes$fdr.gen < fdr.gen) #1142 - 865
index_DOWN <- which(all_genes$logFC < -logFC & all_genes$fdr.gen < fdr.gen) #1387 - 1097
dir.create(paste0(project, "/", sub_project, "/02-blade_d20_FR06_vs_WLDE"))
universe <- universe
PValue = 1
fdr = 1
annotation = "org.At.tair"
project <- paste0(basedir,"/DE/")
all = TRUE
sub_project = "GO/02-blade_d20_FR06_vs_WLDE"
subset <- all_genes$Gene_ID[index_UP]
description = "blade_d20_FR06_vs_WLDE_UP"
writeGSEA(subset, universe, PValue = 1, fdr, annotation, ontology = "BP", description, project, sub_project, all)
writeGSEA(subset, universe, PValue = 1, fdr, annotation, ontology = "MF", description, project, sub_project, all)
writeGSEA(subset, universe, PValue = 1, fdr, annotation, ontology = "CC", description, project, sub_project, all)
subset <- all_genes$Gene_ID[index_DOWN]
description = "blade_d20_FR06_vs_WLDE_DOWN"
writeGSEA(subset, universe, PValue = 1, fdr, annotation, ontology = "BP", description, project, sub_project, all)
writeGSEA(subset, universe, PValue = 1, fdr, annotation, ontology = "MF", description, project, sub_project, all)
writeGSEA(subset, universe, PValue = 1, fdr, annotation, ontology = "CC", description, project, sub_project, all)


# Get GO terms for 03-blade_d20_FR18_vs_WLDE
setwd(basedir)
project <- paste0(basedir,"/DE/")
sub_project = "GO"
all_genes <- read_tsv(file = paste0(basedir,"/DE/03-blade_d20_FR18_vs_WLDE.tab"))
index_ALL <- which(abs(all_genes$logFC) > logFC & all_genes$fdr.gen < fdr.gen) #3011 - 2440
index_UP <- which(all_genes$logFC > logFC & all_genes$fdr.gen < fdr.gen) #1478 - 1212
index_DOWN <- which(all_genes$logFC < -logFC & all_genes$fdr.gen < fdr.gen) #1533 - 1228
dir.create(paste0(project, "/", sub_project, "/03-blade_d20_FR18_vs_WLDE"))
universe <- universe
PValue = 1
fdr = 1
annotation = "org.At.tair"
project <- paste0(basedir,"/DE/")
all = TRUE
sub_project = "GO/03-blade_d20_FR18_vs_WLDE"
subset <- all_genes$Gene_ID[index_UP]
description = "blade_d20_FR18_vs_WLDE_UP"
writeGSEA(subset, universe, PValue = 1, fdr, annotation, ontology = "BP", description, project, sub_project, all)
writeGSEA(subset, universe, PValue = 1, fdr, annotation, ontology = "MF", description, project, sub_project, all)
writeGSEA(subset, universe, PValue = 1, fdr, annotation, ontology = "CC", description, project, sub_project, all)
subset <- all_genes$Gene_ID[index_DOWN]
description = "blade_d20_FR18_vs_WLDE_DOWN"
writeGSEA(subset, universe, PValue = 1, fdr, annotation, ontology = "BP", description, project, sub_project, all)
writeGSEA(subset, universe, PValue = 1, fdr, annotation, ontology = "MF", description, project, sub_project, all)
writeGSEA(subset, universe, PValue = 1, fdr, annotation, ontology = "CC", description, project, sub_project, all)


# Get GO terms for 04-blade_d20_FR18_vs_FR06DE
setwd(basedir)
project <- paste0(basedir,"/DE/")
sub_project = "GO"
all_genes <- read_tsv(file = paste0(basedir,"/DE/04-blade_d20_FR18_vs_FR06DE.tab"))
index_ALL <- which(abs(all_genes$logFC) > logFC & all_genes$fdr.gen < fdr.gen) #26 - 7
index_UP <- which(all_genes$logFC > logFC & all_genes$fdr.gen < fdr.gen) #22 - 6
index_DOWN <- which(all_genes$logFC < -logFC & all_genes$fdr.gen < fdr.gen) #4 - 1
dir.create(paste0(project, "/", sub_project, "/04-blade_d20_FR18_vs_FR06DE"))
universe <- universe
PValue = 1
fdr = 1
annotation = "org.At.tair"
project <- paste0(basedir,"/DE/")
all = TRUE
sub_project = "GO/04-blade_d20_FR18_vs_FR06DE"
subset <- all_genes$Gene_ID[index_UP]
description = "blade_d20_FR18_vs_FR06DE_UP"
writeGSEA(subset, universe, PValue = 1, fdr, annotation, ontology = "BP", description, project, sub_project, all)
writeGSEA(subset, universe, PValue = 1, fdr, annotation, ontology = "MF", description, project, sub_project, all)
writeGSEA(subset, universe, PValue = 1, fdr, annotation, ontology = "CC", description, project, sub_project, all)
subset <- all_genes$Gene_ID[index_DOWN]
description = "blade_d20_FR18_vs_FR06DE_DOWN"
writeGSEA(subset, universe, PValue = 1, fdr, annotation, ontology = "BP", description, project, sub_project, all)
writeGSEA(subset, universe, PValue = 1, fdr, annotation, ontology = "MF", description, project, sub_project, all)
writeGSEA(subset, universe, PValue = 1, fdr, annotation, ontology = "CC", description, project, sub_project, all)

#######################################################
#                KEGG Pathway analysis                #
#######################################################

## Requires
if(!require(pathview)) BiocManager::install("pathview")
if(!require(KEGGREST)) BiocManager::install("KEGGREST")
if(!require(clusterProfiler)) BiocManager::install("clusterProfiler")

## Includes
library(pathview) # load library for KEGG pathway visualization
library(KEGGREST) # get names of pathways to visualize
library(clusterProfiler) # Load library for KEGG enrichment analysis

pathways <- keggList("pathway", "ath")

# Pull all genes for each pathway
pathway.codes <- sub("path:", "", names(pathways)) 
genes.by.pathway <- sapply(pathway.codes,
                           function(pwid){
                             pw <- keggGet(pwid)
                             if (is.null(pw[[1]]$GENE)) return(NA)
                             pw2 <- pw[[1]]$GENE[c(TRUE,FALSE)] # may need to modify this to c(FALSE, TRUE) for other organisms
                             pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
                             return(pw2)
                           }
)


# KEGG enrichment analysis will be performed by ClusterProfiler
# adapted from https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html#clusterprofiler

#Set wd
setwd(basedir)

# Load DE data
gene.data.d13 <- read.delim("DE/00-primordiaDE_full.tab", row.names = 1)
gene.data.d16 <- read.delim("DE/01-blade_d16_FR06_vs_WLDE_full.tab", row.names = 1)
gene.data.d20 <- read.delim("DE/02-blade_d20_FR06_vs_WLDE_full.tab", row.names = 1)
gene.data.d20.lt <- read.delim("DE/03-blade_d20_FR18_vs_WLDE_full.tab", row.names = 1)

# Set cut-off values
logFC <- 0.58
fdr.gen <- 0.05
p.value <- 0.05

# Create dir
setwd(basedir)
dir.create(paste0(basedir,"/DE/KEGG"))
universe <- rownames(gene.data.d13)

# Get enrichment for primordia d13 dataset
gene.data <- data.frame(gene.data.d13$logFC, row.names = rownames(gene.data.d13))

# Down
genes <- gene.data.d13$PValue[which((gene.data.d13$logFC < -logFC) & gene.data.d13$fdr.gen < fdr.gen)] #1856 - 1487
names(genes) <- row.names(gene.data.d13)[which((gene.data.d13$logFC < -logFC) & gene.data.d13$fdr.gen < fdr.gen)]
genes_KEGGresult_down<-pathwayEnrichment(genes, genes.by.pathway, pathways, universe)
drawKEGG(genes_KEGGresult_down[which(genes_KEGGresult_down$p.value < p.value),], "d13.down", paste0(basedir,"/DE/KEGG/d13_down"), gene.data, pathways)
write.table(genes_KEGGresult_down[which(genes_KEGGresult_down$p.value < p.value),], file = paste0(basedir,"/DE/KEGG/d13_down","/KEGG_pathways.tab"), sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
genes_in_pathways <- getGenesInPathway(genes, genes.by.pathway, pathways = as.character(genes_KEGGresult_down$pathway.code[which(genes_KEGGresult_down$p.value < p.value)]))
write.table(genes_in_pathways, file=paste0(basedir,"/DE/KEGG/d13_down","/KEGG_pathways_genes.tab"))
KEGG.table <- data.frame()
KEGG.table <- cbind(day = "d13",
                    condition = "down",
                    genes_KEGGresult_down[which(genes_KEGGresult_down$p.value < p.value),])
rm(genes_KEGGresult_down)

kegg.cp <- enrichKEGG(gene = names(genes), organism = 'ath', universe = universe) # alternative method
kegg.cp.results<-kegg.cp@result
write.table(kegg.cp.results[which(kegg.cp.results$pvalue < p.value),], file = paste0(basedir,"/DE/KEGG/d13_down","/KEGG_clusterProfiler_pathways.tab"), sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
rm(kegg.cp, kegg.cp.results)

# Up
genes<-gene.data.d13$PValue[which((gene.data.d13$logFC > logFC) & gene.data.d13$fdr.gen < fdr.gen)] #1190 - 923
names(genes) <- row.names(gene.data.d13)[which((gene.data.d13$logFC > logFC) & gene.data.d13$fdr.gen < fdr.gen)]
genes_KEGGresult_up<-pathwayEnrichment(genes, genes.by.pathway, pathways, universe)
drawKEGG(genes_KEGGresult_up[which(genes_KEGGresult_up$p.value < p.value),], "d13.up", paste0(basedir,"/DE/KEGG/d13_up"), gene.data, pathways)
write.table(genes_KEGGresult_up[which(genes_KEGGresult_up$p.value < p.value),], file = paste0(basedir,"/DE/KEGG/d13_up","/KEGG_pathways.tab"), sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
genes_in_pathways <- getGenesInPathway(genes, genes.by.pathway, pathways = as.character(genes_KEGGresult_up$pathway.code[which(genes_KEGGresult_up$p.value < p.value)]))
write.table(genes_in_pathways, file=paste0(basedir,"/DE/KEGG/d13_up","/KEGG_pathways_genes.tab"))
drawKEGGpathway("04075", "d13", paste0(basedir,"/DE/KEGG/plant_signaling"), gene.data, pathways)
KEGG.table <- rbind(KEGG.table,
                    cbind(day = "d13",
                          condition = "up",
                          genes_KEGGresult_up[which(genes_KEGGresult_up$p.value < p.value),])
)
rm(genes_KEGGresult_up)
kegg.cp <- enrichKEGG(gene = names(genes), organism = 'ath', universe = universe) # alternative method
kegg.cp.results<-kegg.cp@result
write.table(kegg.cp.results[which(kegg.cp.results$pvalue < p.value),], file = paste0(basedir,"/DE/KEGG/d13_up","/KEGG_clusterProfiler_pathways.tab"), sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
rm(kegg.cp, kegg.cp.results)

# Get enrichment for d16 dataset
gene.data <- data.frame(gene.data.d16$logFC, row.names = rownames(gene.data.d16))
# Down
genes <- gene.data.d16$PValue[which((gene.data.d16$logFC < -logFC) & gene.data.d16$fdr.gen < fdr.gen)] #1245 - 962
names(genes) <- row.names(gene.data.d16)[which((gene.data.d16$logFC < -logFC) & gene.data.d16$fdr.gen < fdr.gen)]
genes_KEGGresult_down<-pathwayEnrichment(genes, genes.by.pathway, pathways, universe)
drawKEGG(genes_KEGGresult_down[which(genes_KEGGresult_down$p.value < p.value),], "d16.down", paste0(basedir,"/DE/KEGG/d16_down"), gene.data, pathways)
write.table(genes_KEGGresult_down[which(genes_KEGGresult_down$p.value < p.value),], file = paste0(basedir,"/DE/KEGG/d16_down","/KEGG_pathways.tab"), sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
genes_in_pathways <- getGenesInPathway(genes, genes.by.pathway, pathways = as.character(genes_KEGGresult_down$pathway.code[which(genes_KEGGresult_down$p.value < p.value)]))
write.table(genes_in_pathways, file=paste0(basedir,"/DE/KEGG/d16_down","/KEGG_pathways_genes.tab"))
KEGG.table <- rbind(KEGG.table,
                    cbind(day = "d16",
                          condition = "down",
                          genes_KEGGresult_down[which(genes_KEGGresult_down$p.value < p.value),])
)
rm(genes_KEGGresult_down)
kegg.cp <- enrichKEGG(gene = names(genes), organism = 'ath', universe = universe) # alternative method
kegg.cp.results<-kegg.cp@result
write.table(kegg.cp.results[which(kegg.cp.results$pvalue < p.value),], file = paste0(basedir,"/DE/KEGG/d16_down","/KEGG_clusterProfiler_pathways.tab"), sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
rm(kegg.cp, kegg.cp.results)

# Up
genes<-gene.data.d16$PValue[which((gene.data.d16$logFC > logFC) & gene.data.d16$fdr.gen < fdr.gen)] #824 - 608
names(genes) <- row.names(gene.data.d16)[which((gene.data.d16$logFC > logFC) & gene.data.d16$fdr.gen < fdr.gen)]
genes_KEGGresult_up<-pathwayEnrichment(genes, genes.by.pathway, pathways, universe)
drawKEGG(genes_KEGGresult_up[which(genes_KEGGresult_up$p.value < p.value),], "d16.up", paste0(basedir,"/DE/KEGG/d16_up"), gene.data, pathways)
write.table(genes_KEGGresult_up[which(genes_KEGGresult_up$p.value < p.value),], file = paste0(basedir,"/DE/KEGG/d16_up","/KEGG_pathways.tab"), sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
genes_in_pathways <- getGenesInPathway(genes, genes.by.pathway, pathways = as.character(genes_KEGGresult_up$pathway.code[which(genes_KEGGresult_up$p.value < p.value)]))
write.table(genes_in_pathways, file=paste0(basedir,"/DE/KEGG/d16_up","/KEGG_pathways_genes.tab"))
drawKEGGpathway("04075", "d16", paste0(basedir,"/DE/KEGG/plant_signaling"), gene.data, pathways)
KEGG.table <- rbind(KEGG.table,
                    cbind(day = "d16",
                          condition = "up",
                          genes_KEGGresult_up[which(genes_KEGGresult_up$p.value < p.value),])
)
rm(genes_KEGGresult_up)
kegg.cp <- enrichKEGG(gene = names(genes), organism = 'ath', universe = universe) # alternative method
kegg.cp.results<-kegg.cp@result
write.table(kegg.cp.results[which(kegg.cp.results$pvalue < p.value),],file = paste0(basedir,"/DE/KEGG/d16_up","/KEGG_clusterProfiler_pathways.tab"), sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
rm(kegg.cp, kegg.cp.results)

# Get enrichment for d20 dataset
gene.data <- data.frame(gene.data.d20$logFC, row.names = rownames(gene.data.d20))
# Down
genes <- gene.data.d20$PValue[which((gene.data.d20$logFC < -logFC) & gene.data.d20$fdr.gen < fdr.gen)] #1387 - 1097
names(genes) <- row.names(gene.data.d20)[which((gene.data.d20$logFC < -logFC) & gene.data.d20$fdr.gen < fdr.gen)]
genes_KEGGresult_down<-pathwayEnrichment(genes, genes.by.pathway, pathways, universe)
drawKEGG(genes_KEGGresult_down[which(genes_KEGGresult_down$p.value < p.value),], "d20.down", paste0(basedir,"/DE/KEGG/d20_down"), gene.data, pathways)
write.table(genes_KEGGresult_down[which(genes_KEGGresult_down$p.value < p.value),], file = paste0(basedir,"/DE/KEGG/d20_down","/KEGG_pathways.tab"), sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
genes_in_pathways <- getGenesInPathway(genes, genes.by.pathway, pathways = as.character(genes_KEGGresult_down$pathway.code[which(genes_KEGGresult_down$p.value < p.value)]))
write.table(genes_in_pathways, file=paste0(basedir,"/DE/KEGG/d20_down","/KEGG_pathways_genes.tab"))
KEGG.table <- rbind(KEGG.table,
                    cbind(day = "d20",
                          condition = "down",
                          genes_KEGGresult_down[which(genes_KEGGresult_down$p.value < p.value),])
)
rm(genes_KEGGresult_down)
kegg.cp <- enrichKEGG(gene = names(genes), organism = 'ath', universe = universe) # alternative method
kegg.cp.results<-kegg.cp@result
write.table(kegg.cp.results[which(kegg.cp.results$pvalue < p.value),], file = paste0(basedir,"/DE/KEGG/d20_down","/KEGG_clusterProfiler_pathways.tab"), sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
rm(kegg.cp, kegg.cp.results)

# Up
genes<-gene.data.d20$PValue[which((gene.data.d20$logFC > logFC) & gene.data.d20$fdr.gen < fdr.gen)] #1142 - 865
names(genes) <- row.names(gene.data.d20)[which((gene.data.d20$logFC > logFC) & gene.data.d20$fdr.gen < fdr.gen)]
genes_KEGGresult_up<-pathwayEnrichment(genes, genes.by.pathway, pathways, universe)
drawKEGG(genes_KEGGresult_up[which(genes_KEGGresult_up$p.value < p.value),], "d20.up", paste0(basedir,"/DE/KEGG/d20_up"), gene.data, pathways)
write.table(genes_KEGGresult_up[which(genes_KEGGresult_up$p.value < p.value),], file = paste0(basedir,"/DE/KEGG/d20_up","/KEGG_pathways.tab"), sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
genes_in_pathways <- getGenesInPathway(genes, genes.by.pathway, pathways = as.character(genes_KEGGresult_up$pathway.code[which(genes_KEGGresult_up$p.value < p.value)]))
write.table(genes_in_pathways, file=paste0(basedir,"/DE/KEGG/d20_up","/KEGG_pathways_genes.tab"))
drawKEGGpathway("04075", "d20", paste0(basedir,"/DE/KEGG/plant_signaling"), gene.data, pathways)
KEGG.table <- rbind(KEGG.table,
                    cbind(day = "d20",
                          condition = "up",
                          genes_KEGGresult_up[which(genes_KEGGresult_up$p.value < p.value),])
)
rm(genes_KEGGresult_up)
kegg.cp <- enrichKEGG(gene = names(genes), organism = 'ath', universe = universe) # alternative method
kegg.cp.results<-kegg.cp@result
write.table(kegg.cp.results[which(kegg.cp.results$pvalue < p.value),], file = paste0(basedir,"/DE/KEGG/d20_up","/KEGG_clusterProfiler_pathways.tab"), sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
rm(kegg.cp, kegg.cp.results)

# Get enrichment for d20.lt dataset
gene.data <- data.frame(gene.data.d20.lt$logFC, row.names = rownames(gene.data.d20.lt))
# Down
genes <- gene.data.d20.lt$PValue[which((gene.data.d20.lt$logFC < -logFC) & gene.data.d20.lt$fdr.gen < fdr.gen)] #1533 - 1228
names(genes) <- row.names(gene.data.d20.lt)[which((gene.data.d20.lt$logFC < -logFC) & gene.data.d20.lt$fdr.gen < fdr.gen)]
genes_KEGGresult_down<-pathwayEnrichment(genes, genes.by.pathway, pathways, universe)
drawKEGG(genes_KEGGresult_down[which(genes_KEGGresult_down$p.value < p.value),], "d20.lt.down", paste0(basedir,"/DE/KEGG/d20.lt_down"), gene.data, pathways)
write.table(genes_KEGGresult_down[which(genes_KEGGresult_down$p.value < p.value),], file = paste0(basedir,"/DE/KEGG/d20.lt_down","/KEGG_pathways.tab"), sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
genes_in_pathways <- getGenesInPathway(genes, genes.by.pathway, pathways = as.character(genes_KEGGresult_down$pathway.code[which(genes_KEGGresult_down$p.value < p.value)]))
write.table(genes_in_pathways, file=paste0(basedir,"/DE/KEGG/d20.lt_down","/KEGG_pathways_genes.tab"))
KEGG.table <- rbind(KEGG.table,
                    cbind(day = "d20.lt",
                          condition = "down",
                          genes_KEGGresult_down[which(genes_KEGGresult_down$p.value < p.value),])
)
rm(genes_KEGGresult_down)
kegg.cp <- enrichKEGG(gene = names(genes), organism = 'ath', universe = universe) # alternative method
kegg.cp.results<-kegg.cp@result
write.table(kegg.cp.results[which(kegg.cp.results$pvalue < p.value),], file = paste0(basedir,"/DE/KEGG/d20.lt_down","/KEGG_clusterProfiler_pathways.tab"), sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
rm(kegg.cp, kegg.cp.results)

# Up
genes<-gene.data.d20.lt$PValue[which((gene.data.d20.lt$logFC > logFC) & gene.data.d20.lt$fdr.gen < fdr.gen)] #1478 - 1212
names(genes) <- row.names(gene.data.d20.lt)[which((gene.data.d20.lt$logFC > logFC) & gene.data.d20.lt$fdr.gen < fdr.gen)]
genes_KEGGresult_up<-pathwayEnrichment(genes, genes.by.pathway, pathways, universe)
drawKEGG(genes_KEGGresult_up[which(genes_KEGGresult_up$p.value < p.value),], "d20.lt.up", paste0(basedir,"/DE/KEGG/d20.lt_up"), gene.data, pathways)
write.table(genes_KEGGresult_up[which(genes_KEGGresult_up$p.value < p.value),], file = paste0(basedir,"/DE/KEGG/d20.lt_up","/KEGG_pathways.tab"), sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
genes_in_pathways <- getGenesInPathway(genes, genes.by.pathway, pathways = as.character(genes_KEGGresult_up$pathway.code[which(genes_KEGGresult_up$p.value < p.value)]))
write.table(genes_in_pathways, file=paste0(basedir,"/DE/KEGG/d20.lt_up","/KEGG_pathways_genes.tab"))
drawKEGGpathway("04075", "d20.lt", paste0(basedir,"/DE/KEGG/plant_signaling"), gene.data, pathways)
KEGG.table <- rbind(KEGG.table,
                    cbind(day = "d20.lt",
                          condition = "up",
                          genes_KEGGresult_up[which(genes_KEGGresult_up$p.value < p.value),])
)
rm(genes_KEGGresult_up)
kegg.cp <- enrichKEGG(gene = names(genes), organism = 'ath', universe = universe) # alternative method
kegg.cp.results<-kegg.cp@result
write.table(kegg.cp.results[which(kegg.cp.results$pvalue < p.value),], file = paste0(basedir,"/DE/KEGG/d20.lt_up","/KEGG_clusterProfiler_pathways.tab"), sep="\t", row.names = TRUE, col.names = TRUE, na = "NA")
rm(kegg.cp, kegg.cp.results)

setwd(paste0(basedir,"/DE/KEGG/"))

# Compile down KEGG pathways
KEGG.d13 <- read.table("d13_down/KEGG_pathways.tab")
KEGG.d16 <- read.table("d16_down/KEGG_pathways.tab")
KEGG.d20 <- read.table("d20_down/KEGG_pathways.tab")
KEGG.d20.lt <- read.table("d20.lt_down/KEGG_pathways.tab")
KEGG.d13$timepoint <- "d13"
KEGG.d16$timepoint <- "d16"
KEGG.d20$timepoint <- "d2006"
KEGG.d20.lt$timepoint <- "d2018"

KEGG.full.timepoints <- rbind(KEGG.d13,
                              KEGG.d16,
                              KEGG.d20,
                              KEGG.d20.lt)

write_tsv(KEGG.full.timepoints, file = "KEGG_pathways_down.tab")

# Compile up KEGG pathways
KEGG.d13 <- read.table("d13_up/KEGG_pathways.tab")
KEGG.d16 <- read.table("d16_up/KEGG_pathways.tab")
KEGG.d20 <- read.table("d20_up/KEGG_pathways.tab")
KEGG.d20.lt <- read.table("d20.lt_up/KEGG_pathways.tab")
KEGG.d13$timepoint <- "d13"
KEGG.d16$timepoint <- "d16"
KEGG.d20$timepoint <- "d2006"
KEGG.d20.lt$timepoint <- "d2018"

KEGG.full.timepoints <- rbind(KEGG.d13,
                              KEGG.d16,
                              KEGG.d20,
                              KEGG.d20.lt)

write_tsv(KEGG.full.timepoints, file = "KEGG_pathways_up.tab")

# Compile down KEGG pathways
KEGG.d13 <- read.table("d13_down/KEGG_clusterProfiler_pathways.tab")
KEGG.d16 <- read.table("d16_down/KEGG_clusterProfiler_pathways.tab")
KEGG.d20 <- read.table("d20_down/KEGG_clusterProfiler_pathways.tab")
KEGG.d20.lt <- read.table("d20.lt_down/KEGG_clusterProfiler_pathways.tab")
KEGG.d13$timepoint <- "d13"
KEGG.d16$timepoint <- "d16"
KEGG.d20$timepoint <- "d2006"
KEGG.d20.lt$timepoint <- "d2018"

KEGG.full.timepoints <- rbind(KEGG.d13,
                              KEGG.d16,
                              KEGG.d20,
                              KEGG.d20.lt)

write_tsv(KEGG.full.timepoints, path = "KEGG_pathways_clusterProfiler_down.tab")

# Compile up KEGG pathways
KEGG.d13 <- read.table("d13_up/KEGG_clusterProfiler_pathways.tab")
KEGG.d16 <- read.table("d16_up/KEGG_clusterProfiler_pathways.tab")
KEGG.d20 <- read.table("d20_up/KEGG_clusterProfiler_pathways.tab")
KEGG.d20.lt <- read.table("d20.lt_up/KEGG_clusterProfiler_pathways.tab")
KEGG.d13$timepoint <- "d13"
KEGG.d16$timepoint <- "d16"
KEGG.d20$timepoint <- "d2006"
KEGG.d20.lt$timepoint <- "d2018"

KEGG.full.timepoints <- rbind(KEGG.d13,
                              KEGG.d16,
                              KEGG.d20,
                              KEGG.d20.lt)

write_tsv(KEGG.full.timepoints, path = "KEGG_pathways_clusterProfiler_up.tab")


# Get plant signaling pathway of combined logFC of every timepoint
setwd(paste0(basedir,"/DE"))
gene.data.d13 <- read.delim(file="00-primordiaDE_full.tab", sep="\t", stringsAsFactors = FALSE, header = TRUE, dec= ".")
gene.data.d16 <- read.delim(file="01-blade_d16_FR06_vs_WLDE_full.tab", sep="\t", stringsAsFactors = FALSE, header = TRUE, dec= ".")
gene.data.d20 <- read.delim(file="02-blade_d20_FR06_vs_WLDE_full.tab", sep="\t", stringsAsFactors = FALSE, header = TRUE, dec= ".")
gene.data.d20.lt <- read.delim(file="03-blade_d20_FR18_vs_WLDE_full.tab", sep="\t", stringsAsFactors = FALSE, header = TRUE, dec= ".")
degs.d13 <- gene.data.d13[which(abs(gene.data.d13$logFC) > logFC & gene.data.d13$PValue < p.value & gene.data.d13$fdr.gen < fdr.gen),] # 3046 - 2410
degs.d16 <- gene.data.d16[which(abs(gene.data.d16$logFC) > logFC & gene.data.d16$PValue < p.value & gene.data.d16$fdr.gen < fdr.gen),] # 2069 - 1570
degs.d20 <- gene.data.d20[which(abs(gene.data.d20$logFC) > logFC & gene.data.d20$PValue < p.value & gene.data.d20$fdr.gen < fdr.gen),] # 2529 - 1962
degs.d20.lt <- gene.data.d20.lt[which(abs(gene.data.d20.lt$logFC) > logFC & gene.data.d20.lt$PValue < p.value & gene.data.d20.lt$fdr.gen < fdr.gen),] # 3011 - 2440
all.degs <- unique(c(degs.d13$Gene_ID, degs.d16$Gene_ID, degs.d20$Gene_ID, degs.d20.lt$Gene_ID)) #6357

gene.data.tp.combined <- data.frame(Gene_ID = all.degs,
                                    logFC.d13 = gene.data.d13$logFC[match(all.degs,gene.data.d13$Gene_ID)],
                                    logFC.d16 = gene.data.d16$logFC[match(all.degs,gene.data.d16$Gene_ID)],
                                    logFC.d20 = gene.data.d20$logFC[match(all.degs,gene.data.d20$Gene_ID)],
                                    logFC.d20.lt = gene.data.d20.lt$logFC[match(all.degs,gene.data.d20.lt$Gene_ID)])

rownames(gene.data.tp.combined) <- gene.data.tp.combined$Gene_ID
gene.data.tp.combined <- gene.data.tp.combined[,-1]

# draw plant hormone signaling for combined datapoints
drawKEGGpathway("04075", "combined_tp", paste0(basedir,"/DE/KEGG/plant_signaling"), gene.data.tp.combined, pathways)
setwd(paste0(basedir,"/DE/KEGG/"))
## Transcription
drawKEGGpathway("03020", "combined_tp", paste0(basedir,"/DE/KEGG/transcription"), gene.data.tp.combined, pathways)  # RNA polymerase
drawKEGGpathway("03022", "combined_tp", paste0(basedir,"/DE/KEGG/transcription"), gene.data.tp.combined, pathways)  # Basal transcription factors
drawKEGGpathway("03040", "combined_tp", paste0(basedir,"/DE/KEGG/transcription"), gene.data.tp.combined, pathways)  # Spliceosome
setwd(paste0(basedir,"/DE/KEGG/"))
## Translation
drawKEGGpathway("03010", "combined_tp", paste0(basedir,"/DE/KEGG/translation"), gene.data.tp.combined, pathways)  #Ribosome
drawKEGGpathway("00970", "combined_tp", paste0(basedir,"/DE/KEGG/translation"), gene.data.tp.combined, pathways)  #Aminoacyl-tRNA biosynthesis
drawKEGGpathway("03013", "combined_tp", paste0(basedir,"/DE/KEGG/translation"), gene.data.tp.combined, pathways)  #RNA transport
drawKEGGpathway("03015", "combined_tp", paste0(basedir,"/DE/KEGG/translation"), gene.data.tp.combined, pathways)  #mRNA surveillance pathway
drawKEGGpathway("03008", "combined_tp", paste0(basedir,"/DE/KEGG/translation"), gene.data.tp.combined, pathways)  #Ribosome biogenesis in eukaryotes
setwd(paste0(basedir,"/DE/KEGG/"))
## Replication and repair
drawKEGGpathway("03030", "combined_tp", paste0(basedir,"/DE/KEGG/replication_and_repair"), gene.data.tp.combined, pathways)  #DNA replication
drawKEGGpathway("03410", "combined_tp", paste0(basedir,"/DE/KEGG/replication_and_repair"), gene.data.tp.combined, pathways)  #Base excision repair
drawKEGGpathway("03420", "combined_tp", paste0(basedir,"/DE/KEGG/replication_and_repair"), gene.data.tp.combined, pathways)  #Nucleotide excision repair
drawKEGGpathway("03430", "combined_tp", paste0(basedir,"/DE/KEGG/replication_and_repair"), gene.data.tp.combined, pathways)  #Mismatch repair
drawKEGGpathway("03440", "combined_tp", paste0(basedir,"/DE/KEGG/replication_and_repair"), gene.data.tp.combined, pathways)  #Homologous recombination
drawKEGGpathway("03450", "combined_tp", paste0(basedir,"/DE/KEGG/replication_and_repair"), gene.data.tp.combined, pathways)  #Non-homologous end-joining
setwd(paste0(basedir,"/DE/KEGG/"))
## Nucleotide metabolism
drawKEGGpathway("00230", "combined_tp", paste0(basedir,"/DE/KEGG/nucleotide_metabolism"), gene.data.tp.combined, pathways)  #Purine metabolism
drawKEGGpathway("00240", "combined_tp", paste0(basedir,"/DE/KEGG/nucleotide_metabolism"), gene.data.tp.combined, pathways)  #Pyrimidine metabolism
setwd(paste0(basedir,"/DE/KEGG/"))
## Energy metabolism
drawKEGGpathway("00190", "combined_tp", paste0(basedir,"/DE/KEGG/energy_metabolism"), gene.data.tp.combined, pathways)  #  Oxidative phosphorylation
drawKEGGpathway("00195", "combined_tp", paste0(basedir,"/DE/KEGG/energy_metabolism"), gene.data.tp.combined, pathways)  #  Photosynthesis
drawKEGGpathway("00196", "combined_tp", paste0(basedir,"/DE/KEGG/energy_metabolism"), gene.data.tp.combined, pathways)  #  Photosynthesis - antenna proteins
drawKEGGpathway("00710", "combined_tp", paste0(basedir,"/DE/KEGG/energy_metabolism"), gene.data.tp.combined, pathways)  #  Carbon fixation in photosynthetic organisms
drawKEGGpathway("00910", "combined_tp", paste0(basedir,"/DE/KEGG/energy_metabolism"), gene.data.tp.combined, pathways)  #  Nitrogen metabolism
drawKEGGpathway("00920", "combined_tp", paste0(basedir,"/DE/KEGG/energy_metabolism"), gene.data.tp.combined, pathways)  #  Sulfur metabolism
setwd(paste0(basedir,"/DE/KEGG/"))
## Signal transduction
drawKEGGpathway("04016", "combined_tp", paste0(basedir,"/DE/KEGG/signaling"), gene.data.tp.combined, pathways)  #  MAPK signaling pathway - plant
drawKEGGpathway("04070", "combined_tp", paste0(basedir,"/DE/KEGG/signaling"), gene.data.tp.combined, pathways)  #  Phosphatidylinositol signaling system
drawKEGGpathway("04075", "combined_tp", paste0(basedir,"/DE/KEGG/signaling"), gene.data.tp.combined, pathways)  #  Plant hormone signal transduction
setwd(paste0(basedir,"/DE/KEGG/"))
## Transport and catabolism
drawKEGGpathway("04144", "combined_tp", paste0(basedir,"/DE/KEGG/transport_and_catabolism"), gene.data.tp.combined, pathways)  #  Endocytosis
drawKEGGpathway("04145", "combined_tp", paste0(basedir,"/DE/KEGG/transport_and_catabolism"), gene.data.tp.combined, pathways)  #  Phagosome
drawKEGGpathway("04146", "combined_tp", paste0(basedir,"/DE/KEGG/transport_and_catabolism"), gene.data.tp.combined, pathways)  #  Peroxisome
drawKEGGpathway("04136", "combined_tp", paste0(basedir,"/DE/KEGG/transport_and_catabolism"), gene.data.tp.combined, pathways)  #  Autophagy - other
setwd(paste0(basedir,"/DE/KEGG/"))


#########################################################################################
#                                                                                       #
#                                       Plots Section                                   #
#                                                                                       #
#########################################################################################

## Requires
if(!require(gplots)) install.packages("gplots")
if(!require(dplyr)) install.packages("dplyr")
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(devtools)) {
  install.packages("devtools")
  devtools::install_github("kassambara/ggpubr")}

## Includes
library(gplots)
library(dplyr)
library(tidyverse)
library(ggpubr)

#######################################
# Tidyverse plots for the manuscript  #
# Cell cycle                          #
#######################################

#Set wd
setwd(basedir)
dir.create("gene_plots")
dir.create("gene_plots/cell_cycle")

# Load metadata
md <- read_tsv(file = "araport11_metadata_17-11-2018_and_curator_summary_31-12-2018.txt")

# Create descriptive factors
treatment <- factor(c("white light", "early EoD FR", "late EoD FR"), levels = c("white light", "early EoD FR", "late EoD FR"))
day = factor(c("d13", "d16", "d20"), levels = c("d13", "d16", "d20"))

# Load DE data
gene.de.data.d13 <- read_tsv("DE/00-primordiaDE_full.tab")
gene.de.data.d16 <- read_tsv("DE/01-blade_d16_FR06_vs_WLDE_full.tab")
gene.de.data.d20 <- read_tsv("DE/02-blade_d20_FR06_vs_WLDE_full.tab")
gene.de.data.d20.lt <- read_tsv("DE/03-blade_d20_FR18_vs_WLDE_full.tab")

# Add treatment and day information to the data
gene.de.data.d13 <- gene.de.data.d13[,c(1,8:11)] %>% mutate(treatment = treatment[2], day = day[1])
gene.de.data.d16 <- gene.de.data.d16[,c(1,8:11)] %>% mutate(treatment = treatment[2], day = day[2])
gene.de.data.d20 <- gene.de.data.d20[,c(1,8:11)] %>% mutate(treatment = treatment[2], day = day[3])
gene.de.data.d20.lt <- gene.de.data.d20.lt[,c(1,8:11)] %>% mutate(treatment = treatment[3], day = day[3])

# Generate combined data
tp.de.data <- rbind(gene.de.data.d13, gene.de.data.d16, gene.de.data.d20, gene.de.data.d20.lt)

# Load logCPM expression data
gene.logCPM.data <- read_tsv("blade.genes.logCPM.normalised.values.tab")
colnames(gene.logCPM.data)[1] <- "Gene_ID"

# Reconstruct day by day data and add treatment and day information to the data
gene.logCPM.WL.d13.data <- rbind(gene.logCPM.data[,c(1,2)] %>% dplyr::rename(logCPM = Pr_D13_WL_01) %>% mutate(treatment = treatment[1], day = day[1])
                                 ,gene.logCPM.data[,c(1,3)] %>% dplyr::rename(logCPM = Pr_D13_WL_02) %>% mutate(treatment = treatment[1], day = day[1])
)

gene.logCPM.WL.d16.data <- rbind(gene.logCPM.data[,c(1,6)] %>% dplyr::rename(logCPM = B_D16_WL_01) %>% mutate(treatment = treatment[1], day = day[2])
                                 ,gene.logCPM.data[,c(1,7)] %>% dplyr::rename(logCPM = B_D16_WL_02) %>% mutate(treatment = treatment[1], day = day[2])
)
gene.logCPM.WL.d20.data <- rbind(gene.logCPM.data[,c(1,10)] %>% dplyr::rename(logCPM = B_D20_WL_01) %>% mutate(treatment = treatment[1], day = day[3])
                                 ,gene.logCPM.data[,c(1,11)] %>% dplyr::rename(logCPM = B_D20_WL_02) %>% mutate(treatment = treatment[1], day = day[3])
)
gene.logCPM.et.d13.data <- rbind(gene.logCPM.data[,c(1,4)] %>% dplyr::rename(logCPM = Pr_D13_FR06_01) %>% mutate(treatment = treatment[2], day = day[1])
                                 ,gene.logCPM.data[,c(1,5)] %>% dplyr::rename(logCPM = Pr_D13_FR06_02) %>% mutate(treatment = treatment[2], day = day[1])
)
gene.logCPM.et.d16.data <- rbind(gene.logCPM.data[,c(1,8)] %>% dplyr::rename(logCPM = B_D16_FR06_01) %>% mutate(treatment = treatment[2], day = day[2])
                                 ,gene.logCPM.data[,c(1,9)] %>% dplyr::rename(logCPM = B_D16_FR06_02) %>% mutate(treatment = treatment[2], day = day[2])
)
gene.logCPM.et.d20.data <- rbind(gene.logCPM.data[,c(1,12)] %>% dplyr::rename(logCPM = B_D20_FR06_01) %>% mutate(treatment = treatment[2], day = day[3])
                                 ,gene.logCPM.data[,c(1,13)] %>% dplyr::rename(logCPM = B_D20_FR06_02) %>% mutate(treatment = treatment[2], day = day[3])
)
gene.logCPM.lt.d20.data <- rbind(gene.logCPM.data[,c(1,14)] %>% dplyr::rename(logCPM = B_D20_FR18_01) %>% mutate(treatment = treatment[3], day = day[3])
                                 ,gene.logCPM.data[,c(1,15)] %>% dplyr::rename(logCPM = B_D20_FR18_02) %>% mutate(treatment = treatment[3], day = day[3])
)

# Generate combined ts data
gene.logCPM.data <- rbind(gene.logCPM.WL.d13.data,
                          gene.logCPM.WL.d16.data,
                          gene.logCPM.WL.d20.data,
                          gene.logCPM.et.d13.data,
                          gene.logCPM.et.d16.data,
                          gene.logCPM.et.d20.data,
                          gene.logCPM.lt.d20.data)

# Load Core Cell Cycle DE data
core.cc.de.data <- read_tsv("DE/selected_genes/cell_cycle_genes_for_manuscript.txt")

### Plot 1: Line plot LogFC in EoD FR1 and EoD FR2 of core cell cycle genes
# Select genes to plot logFC data
gene_to_plot <- tp.de.data %>%
  filter(Gene_ID %in% core.cc.de.data$Gene_ID)

gene_to_plot <- gene_to_plot %>% mutate(group = sapply(gene_to_plot$Gene_ID, function(x) toString(core.cc.de.data$Subtype[which(core.cc.de.data$Gene_ID == x)])))
gene_to_plot <- gene_to_plot %>% mutate(symbol = sapply(gene_to_plot$Gene_ID, function(x) toString(core.cc.de.data$Symbol[which(core.cc.de.data$Gene_ID == x)])))

# for (i in unique(gene_to_plot$Gene_ID)) {
#   gene_to_plot$logFC[gene_to_plot$Gene_ID == i] <- scale(gene_to_plot$logFC[gene_to_plot$Gene_ID == i], center = TRUE, scale = TRUE)
# }

# summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval
tgc <- summarySE(gene_to_plot, measurevar="logFC", groupvars=c("group", "treatment", "day"))

'%!in%' <- function(x,y)!('%in%'(x,y))

tgc1 <- tgc %>%
  filter(group %!in% c("F-Box", "SIM/SMR", "SKIP"))
tgc1 <- tgc1 %>%
  filter(treatment == "early EoD FR")

tgc2 <- tgc %>%
  filter(group %in% c("F-Box", "SIM/SMR", "SKIP"))
tgc2 <- tgc2 %>%
  filter(treatment == "early EoD FR")



geneplot1 <- ggplot(data = tgc1, aes(x = day, y = logFC, color = group, group = group)) +
  geom_line(size=1.5) +
  geom_point(size=2) +
  # scale_color_manual(values=c("#A00000", "#f08080")) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 45, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(colour = "grey20", size = 12),
        text = element_text(size = 14),
        strip.text = element_text(size = 14))+
  xlab("D.A.S.") + ylab(expression('log'[2]*'FC')) +
  scale_x_discrete(breaks = unique(gene_to_plot$day), labels = levels(day)) 

geneplot2 <- ggplot(data = tgc2, aes(x = day, y = logFC, color = group, group = group)) +
  geom_line(size=1.5) +
  geom_point(size=2) +
  # scale_color_manual(values=c("#A00000", "#f08080")) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 45, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(colour = "grey20", size = 12),
        text = element_text(size = 14),
        strip.text = element_text(size = 14))+
  xlab("D.A.S.") + ylab(expression('log'[2]*'FC')) +
  scale_x_discrete(breaks = unique(gene_to_plot$day), labels = levels(day)) 

ggsave("gene_plots/EoD_FR_cell_cycle_groups.tiff",
       plot = ggarrange(geneplot1, geneplot2, ncol = 1, nrow = 2),
       device = "tiff",
       scale = 1,
       width = 26,
       height = 24,
       units = "cm",
       dpi = 600,
       limitsize = TRUE)

# Load Core Cell Cycle DE data
core.cc.de.data <- read_tsv("DE/selected_genes/cell_cycle_genes_for_manuscript.txt")

# Plot core cell cycle genes
gene_to_plot <- gene.logCPM.data %>%
  filter(Gene_ID %in% core.cc.de.data$Gene_ID)

gene_to_plot <- gene_to_plot %>% mutate(group = sapply(gene_to_plot$Gene_ID, function(x) toString(core.cc.de.data$Subtype[which(core.cc.de.data$Gene_ID == x)])))
gene_to_plot <- gene_to_plot %>% mutate(symbol = sapply(gene_to_plot$Gene_ID, function(x) toString(core.cc.de.data$Symbol[which(core.cc.de.data$Gene_ID == x)])))

for (i in unique(gene_to_plot$Gene_ID)) {
  gene_to_plot$logCPM[gene_to_plot$Gene_ID == i] <- scale(gene_to_plot$logCPM[gene_to_plot$Gene_ID == i], center = TRUE, scale = TRUE)
}

### Plot 3: LogCPM in WL, EoD FR1 and EoD FR2 of core cell cycle genes (3 facets)
# summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval
tgc <- summarySE(gene_to_plot, measurevar="logCPM", groupvars=c("Gene_ID","symbol","treatment", "day"))

multigeneplot3 <- ggplot(data = tgc, aes(x = day, y = logCPM, color = Gene_ID, group = Gene_ID)) +
  geom_errorbar(aes(ymin=logCPM-se, ymax=logCPM+se), width=.1) +
  geom_line() +
  geom_point() +
  # ylim(0,(max(tgc$logCPM)+max(tgc$se)+1)) +
  scale_color_hue(labels = c(unique(lapply(tgc$Gene_ID, function(x) toString(md$Symbol[which(md$Gene_ID == x)]))))) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(colour = "grey20", size = 12),
        text = element_text(size = 14))+
  xlab("D.A.S.") + ylab("Scaled normalised expression (logCPM)") +
  scale_x_discrete(breaks = unique(gene_to_plot$day), labels = unique(gene_to_plot$day)) +
  ggtitle("Cell Cycle Genes",
          subtitle = "Leaf blade development - expression values") +
  scale_fill_discrete(name = "Symbol", labels = tgc$symbol) +
  facet_wrap(~treatment)


multigeneplot3

### Plot 4: LogCPM in WL, EoD FR1 and EoD FR2 of core cell cycle genes (one plot per gene)
rm(tgc)
tgc <- summarySE(gene_to_plot, measurevar="logCPM", groupvars=c("Gene_ID","symbol","treatment", "day"))

multigeneplot4 <- ggplot(data = tgc, aes(x = day, y = logCPM, color = treatment, group = treatment)) +
  geom_errorbar(aes(ymin=logCPM-se, ymax=logCPM+se), width=.1) +
  geom_line(size=1) +
  geom_point(size=2) +
  scale_color_manual(values=c("#a9a9a9", "#A00000", "#f08080")) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(colour = "grey20", size = 12),
        text = element_text(size = 14))+
  xlab("D.A.S.") + ylab("Scaled normalised expression (logCPM)") +
  scale_x_discrete(breaks = unique(gene_to_plot$day), labels = unique(gene_to_plot$day)) +
  ggtitle("Cell Cycle Genes",
          subtitle = "Leaf blade development - expression values") +
  scale_fill_discrete(name = "Symbol", labels = tgc$symbol) +
  facet_wrap(~symbol, ncol = 8)

multigeneplot4

### Plot 5: Line plot of average LogCPM of WL, EoD FR1 and EoD FR2 of core cell cycle genes
rm(tgc)
tgc <- summarySE(gene_to_plot, measurevar="logCPM", groupvars=c("treatment", "day"))

geneplot5 <- ggplot(data = tgc, aes(x = day, y = logCPM, color = treatment, group = treatment)) +
  geom_errorbar(aes(ymin=logCPM-se, ymax=logCPM+se), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=2) +
  scale_color_manual(values=c("#a9a9a9", "#A00000", "#f08080")) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(colour = "grey20", size = 12),
        text = element_text(size = 14))+
  xlab("D.A.S.") + ylab("Scaled normalised expression (logCPM)") +
  scale_x_discrete(breaks = unique(gene_to_plot$day), labels = unique(gene_to_plot$day)) +
  ggtitle("Core Cell Cycle Genes",
          subtitle = "Leaf blade development - expression values")

geneplot5

### Plot 5: Line plot of average LogCPM per group of WL, EoD FR1 and EoD FR2 of core cell cycle genes (one plot per gene group)
rm(tgc)
tgc <- summarySE(gene_to_plot, measurevar="logCPM", groupvars=c("treatment", "day", "group"))

multigeneplot6 <- ggplot(data = tgc, aes(x = day, y = logCPM, color = treatment, group = treatment)) +
  geom_errorbar(aes(ymin=logCPM-se, ymax=logCPM+se), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=2) +
  scale_color_manual(values=c("#a9a9a9", "#A00000", "#f08080")) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(colour = "grey20", size = 12),
        text = element_text(size = 14),
        strip.text = element_text(size = 14))+
  xlab("D.A.S.") + ylab("Scaled normalised expression (logCPM)") +
  scale_x_discrete(breaks = unique(gene_to_plot$day), labels = unique(gene_to_plot$day)) +
  ggtitle("Core Cell Cycle Genes",
          subtitle = "Leaf blade development - expression values")+
  facet_wrap(~group, ncol = 5)

multigeneplot6

##########################################################################
# Re-use previous section to plot TFs and hormone Light signalling genes #
##########################################################################
# PIL2, HFR1, HEC1, PAR1 and PAR2 / BEE1, BEE3, BIM1-3, BES1 and BEH2 / BBX6, 17, 21, 23, 27, 28, and 29 / CIB1, SPT, ABI5
# List of TF signalling genes

# Main figure TFs
gene.list.de.data <- data.frame(Gene_ID = c(md$Gene_ID[which(md$Symbol == "PIL1")],
                                            md$Gene_ID[which(md$Symbol == "PIL2")],
                                            md$Gene_ID[which(md$Symbol == "HFR1")],
                                            md$Gene_ID[which(md$Symbol == "PAR1")],
                                            md$Gene_ID[which(md$Symbol == "PAR2")],
                                            md$Gene_ID[which(md$Symbol == "HB-2")],
                                            md$Gene_ID[which(md$Symbol == "HAT2")],
                                            md$Gene_ID[which(md$Symbol == "CIB1")],
                                            md$Gene_ID[which(md$Symbol == "SPT")]
)
)



# BBX6 - COL5, BBX17, BBX21, BBX23, BBX27, BBX28, BBX29, NF-YA2, NF-YA3, NF-YA5, NF-YA7, NF-YA8, NF-YA9,NF-YB6, NF-YB8, NF-YC12
# For supp.TFs
gene.list2.de.data <- data.frame(Gene_ID = c(md$Gene_ID[which(md$Symbol == "BBX6 - COL5")],
                                             md$Gene_ID[which(md$Symbol == "BBX17")],
                                             md$Gene_ID[which(md$Symbol == "BBX21")],
                                             md$Gene_ID[which(md$Symbol == "BBX23")],
                                             md$Gene_ID[which(md$Symbol == "BBX27")],
                                             md$Gene_ID[which(md$Symbol == "BBX28")],
                                             md$Gene_ID[which(md$Symbol == "BBX29")],
                                             md$Gene_ID[which(md$Symbol == "NF-YA2")],
                                             md$Gene_ID[which(md$Symbol == "NF-YA3")],
                                             md$Gene_ID[which(md$Symbol == "NF-YA5")],
                                             md$Gene_ID[which(md$Symbol == "NF-YA7")],
                                             md$Gene_ID[which(md$Symbol == "NF-YA8")],
                                             md$Gene_ID[which(md$Symbol == "NF-YA9")],
                                             md$Gene_ID[which(md$Symbol == "NF-YB6")],
                                             md$Gene_ID[which(md$Symbol == "NF-YB8")],
                                             md$Gene_ID[which(md$Symbol == "NF-YC12")]
)
)

# SAUR23, SAUR24, ARF18, IAA19, IAA29, GH3.3, WES1  
# SAUR9, SAUR19, SAUR20, SAUR22, SAUR29, SAUR63, SAUR68, ARF10, MP, IAA2, IAA34, IAA7, IAA5, BRU6, DFL2, DFL1 

# Main figure hormones
gene.list3.de.data <- data.frame(Gene_ID = c(md$Gene_ID[which(md$Symbol == "SAUR23")],
                                             md$Gene_ID[which(md$Symbol == "SAUR24")],
                                             md$Gene_ID[which(md$Symbol == "ARF18")],
                                             md$Gene_ID[which(md$Symbol == "IAA19")],
                                             md$Gene_ID[which(md$Symbol == "GH3.3")],
                                             md$Gene_ID[which(md$Symbol == "WES1")],
                                             md$Gene_ID[which(md$Symbol == "SAUR9")],
                                             md$Gene_ID[which(md$Symbol == "SAUR19")],
                                             md$Gene_ID[which(md$Symbol == "SAUR20")],
                                             md$Gene_ID[which(md$Symbol == "SAUR22")],
                                             md$Gene_ID[which(md$Symbol == "SAUR29")],
                                             md$Gene_ID[which(md$Symbol == "SAUR63")],
                                             md$Gene_ID[which(md$Symbol == "SAUR68")],
                                             md$Gene_ID[which(md$Symbol == "ARF10")],
                                             md$Gene_ID[which(md$Symbol == "MP")],
                                             md$Gene_ID[which(md$Symbol == "IAA2")],
                                             md$Gene_ID[which(md$Symbol == "IAA34")],
                                             md$Gene_ID[which(md$Symbol == "IAA7")],
                                             md$Gene_ID[which(md$Symbol == "IAA5")],
                                             md$Gene_ID[which(md$Symbol == "BRU6")],
                                             md$Gene_ID[which(md$Symbol == "DFL2")],
                                             md$Gene_ID[which(md$Symbol == "DFL1")]
)
)

# AT4G10320, AT5G49030, AT4G17300, AT3G02760, AT1G72550, AT1G66530, AT1G25350, AT2G31170, AT1G29870, AT3G48110, AT4G26870, AT4G33760, AT3G55400,
# AT5G16715, AT5G22800, AT3G11710, AT1G09620, AT4G04350, AT1G17960, AT2G04842, AT2G25840, AT1G66520, AT1G28350, AT3G62120, AT5G52520, AT5G26710
# AT5G64050, AT4G32915, AT1G11870 
# Supp. fig tRNAcyl
gene.list4.de.data <- data.frame(Gene_ID = c("AT4G10320", "AT5G49030", "AT4G17300", "AT3G02760", "AT1G72550", "AT1G66530",
                                             "AT1G25350", "AT2G31170", "AT1G29870", "AT3G48110", "AT4G26870", "AT4G33760",
                                             "AT3G55400", "AT5G16715", "AT5G22800", "AT3G11710", "AT1G09620", "AT4G04350",
                                             "AT1G17960", "AT2G04842", "AT2G25840", "AT1G66520", "AT1G28350", "AT3G62120",
                                             "AT5G52520", "AT5G26710", "AT5G64050", "AT4G32915", "AT1G11870"
)
)

# Auxin signalling
# ARF1, MP, ARF10, ARF11, ARF18, IAA1, IAA2, IAA5, IAA7, AXR3, IAA19, IAA29, IAA34
gene.list5.de.data <- data.frame(Gene_ID = c(md$Gene_ID[which(md$Symbol == "ARF1")],
                                             md$Gene_ID[which(md$Symbol == "MP")],
                                             md$Gene_ID[which(md$Symbol == "ARF10")],
                                             md$Gene_ID[which(md$Symbol == "ARF11")],
                                             md$Gene_ID[which(md$Symbol == "ARF18")],
                                             md$Gene_ID[which(md$Symbol == "IAA1")],
                                             md$Gene_ID[which(md$Symbol == "IAA2")],
                                             md$Gene_ID[which(md$Symbol == "IAA5")],
                                             md$Gene_ID[which(md$Symbol == "IAA7")],
                                             md$Gene_ID[which(md$Symbol == "AXR3")],
                                             md$Gene_ID[which(md$Symbol == "IAA19")],
                                             md$Gene_ID[which(md$Symbol == "IAA29")],
                                             md$Gene_ID[which(md$Symbol == "IAA34")]
)
)



### Plot main fig TF genes
gene_to_plot <- gene.logCPM.data %>%
  filter(Gene_ID %in% gene.list.de.data$Gene_ID)
gene_to_plot <- gene_to_plot %>% mutate(symbol = sapply(gene_to_plot$Gene_ID, function(x) toString(md$Symbol[which(md$Gene_ID == x)])))
for (i in unique(gene_to_plot$Gene_ID)) {
  gene_to_plot$logCPM[gene_to_plot$Gene_ID == i] <- scale(gene_to_plot$logCPM[gene_to_plot$Gene_ID == i], center = TRUE, scale = TRUE)
}

### Multigene Plot: LogCPM in WL, EoD FR1 and EoD FR2 of core cell cycle genes (one plot per gene)
rm(tgc)
tgc <- summarySE(gene_to_plot, measurevar="logCPM", groupvars=c("Gene_ID","symbol","treatment", "day"))
multigeneplot <- ggplot(data = tgc, aes(x = day, y = logCPM, color = treatment, group = treatment)) +
  geom_errorbar(aes(ymin=logCPM-se, ymax=logCPM+se), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=2) +
  scale_color_manual(values=c("#a9a9a9", "#A00000", "#f08080")) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 45, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(colour = "grey20", size = 12),
        text = element_text(size = 14),
        strip.text = element_text(size = 14))+
  xlab("D.A.S.") + ylab("Scaled normalised expression (logCPM)") +
  scale_x_discrete(breaks = unique(gene_to_plot$day), labels = unique(gene_to_plot$day)) +
  ggtitle("EoD FR responsive TFs",
          subtitle = "Leaf blade development - expression values") +
  scale_fill_discrete(name = "Symbol", labels = tgc$symbol) +
  facet_wrap(factor(symbol, levels  = c("PIL1", "PIL2", "HFR1",
                                        "PAR1", "PAR2", "HB-2",
                                        "HAT2", "CIB1", "SPT")
  ) ~ ., ncol = 3)

multigeneplot
ggsave("gene_plots/EoD_FR_responsive_TFs.tiff",
       plot = last_plot(),
       device = "tiff",
       scale = 1,
       width = 24,
       height = 14,
       units = "cm",
       dpi = 600,
       limitsize = TRUE)

### Plot sup fig TF genes
gene_to_plot <- gene.logCPM.data %>%
  filter(Gene_ID %in% gene.list2.de.data$Gene_ID)
gene_to_plot <- gene_to_plot %>% mutate(symbol = sapply(gene_to_plot$Gene_ID, function(x) toString(md$Symbol[which(md$Gene_ID == x)])))
for (i in unique(gene_to_plot$Gene_ID)) {
  gene_to_plot$logCPM[gene_to_plot$Gene_ID == i] <- scale(gene_to_plot$logCPM[gene_to_plot$Gene_ID == i], center = TRUE, scale = TRUE)
}

### Multigene Plot: LogCPM in WL, EoD FR1 and EoD FR2 of core cell cycle genes (one plot per gene)
rm(tgc)
tgc <- summarySE(gene_to_plot, measurevar="logCPM", groupvars=c("Gene_ID","symbol","treatment", "day"))
multigeneplot <- ggplot(data = tgc, aes(x = day, y = logCPM, color = treatment, group = treatment)) +
  geom_errorbar(aes(ymin=logCPM-se, ymax=logCPM+se), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=2) +
  scale_color_manual(values=c("#a9a9a9", "#A00000", "#f08080")) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 45, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(colour = "grey20", size = 12),
        text = element_text(size = 14))+
  xlab("D.A.S.") + ylab("Scaled normalised expression (logCPM)") +
  scale_x_discrete(breaks = unique(gene_to_plot$day), labels = unique(gene_to_plot$day)) +
  ggtitle("Shade responsive TFs",
          subtitle = "Leaf blade development - expression values") +
  scale_fill_discrete(name = "Symbol", labels = tgc$symbol) +
  facet_wrap(factor(symbol, levels  = c("BBX6 - COL5", "BBX17", "BBX21", "BBX23", "BBX27", "BBX28",
                                        "BBX29", "NF-YA2", "NF-YA3", "NF-YA5", "NF-YA7", "NF-YA8",
                                        "NF-YA9", "NF-YB6", "NF-YB8", "NF-YC12")
  ) ~ ., ncol = 4)

multigeneplot

ggsave("gene_plots/EoD_FR_responsive_TFs-SI_data.tiff",
       plot = last_plot(),
       device = "tiff",
       scale = 1,
       width = 20,
       height = 12.4,
       units = "cm",
       dpi = 600,
       limitsize = TRUE)


# BBX6 - COL5, BBX17, BBX21, BBX23, BBX27, BBX28, BBX29, NF-YA2, NF-YA3, NF-YA5, NF-YA7, NF-YA8, NF-YA9,NF-YB6, NF-YB8, NF-YC12
### Plot sup fig TF genes (logFC)
gene_to_plot <- tp.de.data %>%
  filter(Gene_ID %in% gene.list2.de.data$Gene_ID)
gene_to_plot <- gene_to_plot %>% mutate(symbol = sapply(gene_to_plot$Gene_ID, function(x) toString(md$Symbol[which(md$Gene_ID == x)])))

### Multigene Plot: LogCPM in WL, EoD FR1 and EoD FR2 of core cell cycle genes (one plot per gene)
rm(tgc)
tgc <- summarySE(gene_to_plot, measurevar="logFC", groupvars=c("Gene_ID","symbol","treatment", "day"))
multigeneplot <- ggplot(data = tgc, aes(x = day, y = logFC, color = treatment, group = treatment)) +
  geom_line(size=1.5) +
  geom_point(size=2) +
  scale_color_manual(values=c("#A00000", "#f08080")) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(colour = "grey20", size = 12),
        text = element_text(size = 14))+
  xlab("D.A.S.") + ylab(label = expression("Log"[2]*"FC")) +
  scale_x_discrete(breaks = unique(gene_to_plot$day), labels = unique(gene_to_plot$day)) +
  ggtitle("Shade responsive TFs",
          subtitle = "Leaf blade development") +
  # scale_fill_discrete(name = "Symbol", labels = tgc$symbol) +
  facet_wrap(factor(symbol, levels  = c("BBX6 - COL5", "BBX17", "BBX21", "BBX23","BBX27", "BBX28",
                                        "BBX29", "NF-YA2", "NF-YA3", "NF-YA5", "NF-YA7", "NF-YA8",
                                        "NF-YA9", "NF-YB6", "NF-YB8", "NF-YC12")
  ) ~ ., ncol = 4)

multigeneplot


### Plot sup fig hormone genes
gene_to_plot <- gene.logCPM.data %>%
  filter(Gene_ID %in% gene.list3.de.data$Gene_ID)
gene_to_plot <- gene_to_plot %>% mutate(symbol = sapply(gene_to_plot$Gene_ID, function(x) toString(md$Symbol[which(md$Gene_ID == x)])))
for (i in unique(gene_to_plot$Gene_ID)) {
  gene_to_plot$logCPM[gene_to_plot$Gene_ID == i] <- scale(gene_to_plot$logCPM[gene_to_plot$Gene_ID == i], center = TRUE, scale = TRUE)
}

### Multigene Plot: LogCPM in WL, EoD FR1 and EoD FR2 of core cell cycle genes (one plot per gene)
rm(tgc)
tgc <- summarySE(gene_to_plot, measurevar="logCPM", groupvars=c("Gene_ID","symbol","treatment", "day"))
multigeneplot <- ggplot(data = tgc, aes(x = day, y = logCPM, color = treatment, group = treatment)) +
  geom_errorbar(aes(ymin=logCPM-se, ymax=logCPM+se), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=2) +
  scale_color_manual(values=c("#a9a9a9", "#A00000", "#f08080")) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(colour = "grey20", size = 12),
        text = element_text(size = 14))+
  xlab("D.A.S.") + ylab("Scaled normalised expression (logCPM)") +
  scale_x_discrete(breaks = unique(gene_to_plot$day), labels = unique(gene_to_plot$day)) +
  ggtitle("TF and light signalling genes",
          subtitle = "Leaf blade development - expression values") +
  scale_fill_discrete(name = "Symbol", labels = tgc$symbol) +
  facet_wrap(~symbol, ncol = 5)

multigeneplot

### Plot sup fig tRNAcyl genes
gene_to_plot <- gene.logCPM.data %>%
  filter(Gene_ID %in% gene.list4.de.data$Gene_ID)
gene_to_plot <- gene_to_plot %>% mutate(symbol = sapply(gene_to_plot$Gene_ID, function(x) toString(md$Symbol[which(md$Gene_ID == x)])))
for (i in unique(gene_to_plot$Gene_ID)) {
  gene_to_plot$logCPM[gene_to_plot$Gene_ID == i] <- scale(gene_to_plot$logCPM[gene_to_plot$Gene_ID == i], center = TRUE, scale = TRUE)
}

### Multigene Plot: LogCPM in WL, EoD FR1 and EoD FR2 of tRNAcyl genes (one plot per gene)
rm(tgc)
tgc <- summarySE(gene_to_plot, measurevar="logCPM", groupvars=c("Gene_ID","symbol","treatment", "day"))
multigeneplot <- ggplot(data = tgc, aes(x = day, y = logCPM, color = treatment, group = treatment)) +
  geom_errorbar(aes(ymin=logCPM-se, ymax=logCPM+se), width=.1) +
  geom_line(size=1) +
  geom_point(size=2) +
  scale_color_manual(values=c("#a9a9a9", "#A00000", "#f08080")) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(colour = "grey20", size = 12),
        text = element_text(size = 14))+
  xlab("D.A.S.") + ylab("Scaled normalised expression (logCPM)") +
  scale_x_discrete(breaks = unique(gene_to_plot$day), labels = unique(gene_to_plot$day)) +
  ggtitle("Aminoacyl tRNA biosynthesis genes",
          subtitle = "Leaf blade development - expression values") +
  scale_fill_discrete(name = "Symbol", labels = tgc$symbol) +
  facet_wrap(~symbol, ncol = 5)

multigeneplot

### Plot fig average auxin signalling hormone genes
gene_to_plot <- gene.logCPM.data %>%
  filter(Gene_ID %in% gene.list5.de.data$Gene_ID)
gene_to_plot <- gene_to_plot %>% mutate(symbol = sapply(gene_to_plot$Gene_ID, function(x) toString(md$Symbol[which(md$Gene_ID == x)])))
for (i in unique(gene_to_plot$Gene_ID)) {
  gene_to_plot$logCPM[gene_to_plot$Gene_ID == i] <- scale(gene_to_plot$logCPM[gene_to_plot$Gene_ID == i], center = TRUE, scale = TRUE)
}

### Multigene Plot: LogCPM in WL, EoD FR1 and EoD FR2 of core cell cycle genes (one plot per gene)
rm(tgc)
tgc <- summarySE(gene_to_plot, measurevar="logCPM", groupvars=c("treatment", "day"))
multigeneplot <- ggplot(data = tgc, aes(x = day, y = logCPM, color = treatment, group = treatment)) +
  geom_errorbar(aes(ymin=logCPM-se, ymax=logCPM+se), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=2) +
  scale_color_manual(values=c("#a9a9a9", "#A00000", "#f08080")) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(colour = "grey20", size = 12),
        text = element_text(size = 14))+
  xlab("D.A.S.") + ylab("Scaled normalised expression (logCPM)") +
  scale_x_discrete(breaks = unique(gene_to_plot$day), labels = unique(gene_to_plot$day)) +
  ggtitle("Auxin signalling genes",
          subtitle = "Leaf blade development - expression values")

multigeneplot

#####################################################################
## Reuse previous section for Cell cycle and associated processes ###
#####################################################################
# Load Core Cell Cycle DE data
cc.and.other.genes <- read_tsv("DE/selected_genes/cell_cycle_and_associated_processes_genes_for_manuscript.txt")

# Plot core cell cycle genes
gene_list <- unique(cc.and.other.genes$Gene_ID[which(cc.and.other.genes$Process == "Cell Cycle")])
gene_to_plot1 <- gene.logCPM.data %>%
  filter(Gene_ID %in% gene_list)
gene_to_plot1 <- gene_to_plot1 %>% mutate(symbol = sapply(gene_to_plot1$Gene_ID, function(x) toString(md$Symbol[which(md$Gene_ID == x)])))
gene_to_plot1 <- gene_to_plot1 %>% mutate(group = "Cell Cycle")

gene_list <- unique(cc.and.other.genes$Gene_ID[which(cc.and.other.genes$Process == "DNA replication")])
gene_to_plot2 <- gene.logCPM.data %>%
  filter(Gene_ID %in% gene_list)
gene_to_plot2 <- gene_to_plot2 %>% mutate(symbol = sapply(gene_to_plot2$Gene_ID, function(x) toString(md$Symbol[which(md$Gene_ID == x)])))
gene_to_plot2 <- gene_to_plot2 %>% mutate(group = "DNA replication")

gene_list <- unique(cc.and.other.genes$Gene_ID[which(cc.and.other.genes$Process == "DNA repair")])
gene_to_plot3 <- gene.logCPM.data %>%
  filter(Gene_ID %in% gene_list)
gene_to_plot3 <- gene_to_plot3 %>% mutate(symbol = sapply(gene_to_plot3$Gene_ID, function(x) toString(md$Symbol[which(md$Gene_ID == x)])))
gene_to_plot3 <- gene_to_plot3 %>% mutate(group = "DNA repair")

gene_list <- unique(cc.and.other.genes$Gene_ID[which(cc.and.other.genes$Process == "Cytokinesis")])
gene_to_plot4 <- gene.logCPM.data %>%
  filter(Gene_ID %in% gene_list)
gene_to_plot4 <- gene_to_plot4 %>% mutate(symbol = sapply(gene_to_plot4$Gene_ID, function(x) toString(md$Symbol[which(md$Gene_ID == x)])))
gene_to_plot4 <- gene_to_plot4 %>% mutate(group = "Cytokinesis")

gene_to_plot <- rbind(gene_to_plot1, gene_to_plot2, gene_to_plot3, gene_to_plot4)

for (i in unique(gene_to_plot$Gene_ID)) {
  gene_to_plot$logCPM[gene_to_plot$Gene_ID == i] <- scale(gene_to_plot$logCPM[gene_to_plot$Gene_ID == i], center = TRUE, scale = TRUE)
}


rm(tgc)
tgc <- summarySE(gene_to_plot, measurevar="logCPM", groupvars=c("treatment", "day", "group"))

multigeneplot <- ggplot(data = tgc, aes(x = day, y = logCPM, color = treatment, group = treatment)) +
  geom_errorbar(aes(ymin=logCPM-se, ymax=logCPM+se), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=2) +
  scale_color_manual(values=c("#a9a9a9", "#A00000", "#f08080")) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 45, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(colour = "grey20", size = 12),
        text = element_text(size = 14),
        strip.text = element_text(size = 14))+
  xlab("D.A.S.") + ylab("Scaled normalised expression (logCPM)") +
  scale_x_discrete(breaks = unique(gene_to_plot$day), labels = unique(gene_to_plot$day)) +
  ggtitle("Basic cellular processes",
          subtitle = "Leaf blade development - expression values")+
  facet_wrap(~group, ncol = 4)

multigeneplot

ggsave("gene_plots/EoD_FR_basic_cellular_processes.tiff",
       plot = last_plot(),
       device = "tiff",
       scale = 1,
       width = 28,
       height = 10,
       units = "cm",
       dpi = 600,
       limitsize = TRUE)

# Plot core cell cycle genes
gene_list <- unique(cc.and.other.genes$Gene_ID[which(cc.and.other.genes$Process == "Cell Cycle")])
gene_to_plot1 <- tp.de.data %>%
  filter(Gene_ID %in% gene_list)
gene_to_plot1 <- gene_to_plot1 %>% mutate(symbol = sapply(gene_to_plot1$Gene_ID, function(x) toString(md$Symbol[which(md$Gene_ID == x)])))
gene_to_plot1 <- gene_to_plot1 %>% mutate(group = "Cell Cycle")

gene_list <- unique(cc.and.other.genes$Gene_ID[which(cc.and.other.genes$Process == "DNA replication")])
gene_to_plot2 <- tp.de.data %>%
  filter(Gene_ID %in% gene_list)
gene_to_plot2 <- gene_to_plot2 %>% mutate(symbol = sapply(gene_to_plot2$Gene_ID, function(x) toString(md$Symbol[which(md$Gene_ID == x)])))
gene_to_plot2 <- gene_to_plot2 %>% mutate(group = "DNA replication")

gene_list <- unique(cc.and.other.genes$Gene_ID[which(cc.and.other.genes$Process == "DNA repair")])
gene_to_plot3 <- tp.de.data %>%
  filter(Gene_ID %in% gene_list)
gene_to_plot3 <- gene_to_plot3 %>% mutate(symbol = sapply(gene_to_plot3$Gene_ID, function(x) toString(md$Symbol[which(md$Gene_ID == x)])))
gene_to_plot3 <- gene_to_plot3 %>% mutate(group = "DNA repair")

gene_list <- unique(cc.and.other.genes$Gene_ID[which(cc.and.other.genes$Process == "Cytokinesis")])
gene_to_plot4 <- tp.de.data %>%
  filter(Gene_ID %in% gene_list)
gene_to_plot4 <- gene_to_plot4 %>% mutate(symbol = sapply(gene_to_plot4$Gene_ID, function(x) toString(md$Symbol[which(md$Gene_ID == x)])))
gene_to_plot4 <- gene_to_plot4 %>% mutate(group = "Cytokinesis")

gene_to_plot <- rbind(gene_to_plot1, gene_to_plot2, gene_to_plot3, gene_to_plot4)

for (i in unique(gene_to_plot$Gene_ID)) {
  gene_to_plot$logFC[gene_to_plot$Gene_ID == i] <- scale(gene_to_plot$logFC[gene_to_plot$Gene_ID == i], center = TRUE, scale = TRUE)
}

rm(tgc)
tgc <- summarySE(gene_to_plot, measurevar="logFC", groupvars=c("treatment", "day", "group"))

multigeneplot <- ggplot(data = tgc, aes(x = day, y = logFC, color = treatment, group = treatment)) +
  geom_line(size=1.5) +
  geom_point(size=2) +
  scale_color_manual(values=c("#A00000", "#f08080")) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(colour = "grey20", size = 12),
        text = element_text(size = 14),
        strip.text = element_text(size = 14))+
  xlab("D.A.S.") + ylab("Mean logFC") +
  scale_x_discrete(breaks = unique(gene_to_plot$day), labels = unique(gene_to_plot$day)) +
  ggtitle("Basic cellular processes",
          subtitle = "Leaf blade development - mean DGE values")+
  facet_wrap(~group, ncol = 2)

multigeneplot

#############################################
# Plot heatmaps of selected lists of genes  #
#############################################
# This list of genes was compiled based on manual curation from GO Terms and KEGG pathways analysis
# A total of 91 categories were chosen

#Set wd
setwd(basedir)

# get metadata
md <- read_tsv(file="araport11_metadata_17-11-2018_and_curator_summary_31-12-2018.txt")

# List of logFC of all misregulated genes
full.degs.d13 <- read.delim(file="DE/00-primordiaDE_full.tab", sep="\t", stringsAsFactors = FALSE, header = TRUE, dec= ".")
full.degs.d16 <- read.delim(file="DE/01-blade_d16_FR06_vs_WLDE_full.tab", sep="\t", stringsAsFactors = FALSE, header = TRUE, dec= ".")
full.degs.d20 <- read.delim(file="DE/02-blade_d20_FR06_vs_WLDE_full.tab", sep="\t", stringsAsFactors = FALSE, header = TRUE, dec= ".")
full.degs.d20.lt <- read.delim(file="DE/03-blade_d20_FR18_vs_WLDE_full.tab", sep="\t", stringsAsFactors = FALSE, header = TRUE, dec= ".")
degs.d13 <- full.degs.d13[which(abs(full.degs.d13$logFC) > 0.58 & full.degs.d13$PValue < p.value & full.degs.d13$fdr.gen < fdr.gen),] # 3046
degs.d16 <- full.degs.d16[which(abs(full.degs.d16$logFC) > 0.58 & full.degs.d16$PValue < p.value & full.degs.d16$fdr.gen < fdr.gen),] # 2069
degs.d20 <- full.degs.d20[which(abs(full.degs.d20$logFC) > 0.58 & full.degs.d20$PValue < p.value & full.degs.d20$fdr.gen < fdr.gen),] # 2529
degs.d20.lt <- full.degs.d20.lt[which(abs(full.degs.d20.lt$logFC) > 0.58 & full.degs.d20.lt$PValue < p.value & full.degs.d20.lt$fdr.gen < fdr.gen),] # 3011
all.degs <- unique(c(degs.d13$Gene_ID, degs.d16$Gene_ID, degs.d20$Gene_ID, degs.d20.lt$Gene_ID))

degs.logFC <- data.frame(Gene_ID = all.degs,
                         logFC.d13 = full.degs.d13$logFC[match(all.degs,full.degs.d13$Gene_ID)],
                         logFC.d16 = full.degs.d16$logFC[match(all.degs,full.degs.d16$Gene_ID)],
                         logFC.d20 = full.degs.d20$logFC[match(all.degs,full.degs.d20$Gene_ID)],
                         logFC.d20.lt = full.degs.d20.lt$logFC[match(all.degs,full.degs.d20.lt$Gene_ID)])

degs.logFC.et <- data.frame(Gene_ID = all.degs,
                            logFC.d13 = full.degs.d13$logFC[match(all.degs,full.degs.d13$Gene_ID)],
                            logFC.d16 = full.degs.d16$logFC[match(all.degs,full.degs.d16$Gene_ID)],
                            logFC.d20 =  full.degs.d20$logFC[match(all.degs,full.degs.d20$Gene_ID)])

# Get lists of genes
genes.table <- read.csv("DE/selected_genes/selected_genes.txt", sep = "\t", stringsAsFactors = FALSE)
gene.lists <- list()

for (i in seq_along(colnames(genes.table))){
  gene.lists[[i]] <- unique(genes.table[,i])
}
names(gene.lists) <- colnames(genes.table)

# Generate dataframes for proportions and filtered genes
proportions <- data.frame(Term = c(),
                          Misregulated = c(),
                          Total = c())
filtered.gene.lists <- list()
gene.list.df <- data.frame(Process = c(),
                           Gene_ID = c())

# Calculate proportions of genes of each category that are also DEGs
for (i in seq_along(names(gene.lists))){
  print(names(gene.lists)[i])
  genes_to_plot <- unique(unlist(gene.lists[[i]]))
  genes_to_plot <- genes_to_plot[which(genes_to_plot %in% c(degs.d13$Gene_ID, degs.d16$Gene_ID, degs.d20$Gene_ID, degs.d20.lt$Gene_ID))]
  filtered.gene.lists[[i]] <- genes_to_plot
  print(c(names(gene.lists)[i],
          length(genes_to_plot),
          length(unique(unlist(gene.lists[[i]])))
  )
  )
  
  proportions[i,1] <- names(gene.lists)[i]
  proportions[i,2] <- length(genes_to_plot)
  proportions[i,3] <- length(unique(unlist(gene.lists[[i]])))
  
  print(genes_to_plot)
  size <- length(rownames(gene.list.df))
  for (j in seq_along(genes_to_plot)){
    if (genes_to_plot[j] != ""){
      gene.list.df[(j + size),1] <- names(gene.lists)[i]
      gene.list.df[(j + size),2] <- genes_to_plot[j]
    }
  }
}  

colnames(gene.list.df) <- c("Term", "Gene_ID")
names(filtered.gene.lists) <- names(gene.lists)
colnames(proportions) <- c("Term", "Misregulated", "Total")
proportions <- cbind(proportions,
                     Ratio = proportions$Misregulated/proportions$Total
)

# Save data to files
write.table(proportions, file = "DE/selected_genes/proportions_of_selected_genes.txt", sep = "\t", col.names = TRUE, row.names = FALSE)
write.table(gene.list.df, file = "DE/selected_genes/compiled_filtered_selected_genes.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

# Reload saved filtered list
gene.list.df <- read.delim(file="DE/selected_genes/compiled_filtered_selected_genes.txt", sep="\t", stringsAsFactors = FALSE, header = TRUE, dec= ".")

# Add DGE data to the list
df <- data.frame(Gene_ID = gene.list.df$Gene_ID,
                 Term = gene.list.df$Term)

df <- cbind(df,
            logFC.d13 = full.degs.d13$logFC[match(df$Gene_ID, full.degs.d13$Gene_ID)],
            p.value.d13 = full.degs.d13$PValue[match(df$Gene_ID, full.degs.d13$Gene_ID)],
            FDR.d13 = full.degs.d13$fdr.gen[match(df$Gene_ID, full.degs.d13$Gene_ID)],
            
            logFC.d16 = full.degs.d16$logFC[match(df$Gene_ID, full.degs.d16$Gene_ID)],
            p.value.d16 = full.degs.d16$PValue[match(df$Gene_ID, full.degs.d16$Gene_ID)],
            FDR.d16 = full.degs.d16$fdr.gen[match(df$Gene_ID, full.degs.d16$Gene_ID)],
            
            logFC.d20 = full.degs.d20$logFC[match(df$Gene_ID, full.degs.d20$Gene_ID)],
            p.value.d20 = full.degs.d20$PValue[match(df$Gene_ID, full.degs.d20$Gene_ID)],
            FDR.d20 = full.degs.d20$fdr.gen[match(df$Gene_ID, full.degs.d20$Gene_ID)],
            
            logFC.d20.lt = full.degs.d20.lt$logFC[match(df$Gene_ID, full.degs.d20.lt$Gene_ID)],
            p.value.d20.lt = full.degs.d20.lt$PValue[match(df$Gene_ID, full.degs.d20.lt$Gene_ID)],
            FDR.d20.lt = full.degs.d20.lt$fdr.gen[match(df$Gene_ID, full.degs.d20.lt$Gene_ID)])

# Annotate the gene lists
df <- getAnnnot(df, md)
gene.list.df <- cbind(df[,15],
                      df[,1:14],
                      df[,16:length(colnames(df))])

colnames(gene.list.df)[1] <- "Term"
gene.list.df <- gene.list.df[order(gene.list.df$Term),]

# Save the table to a file
write.table(gene.list.df, file = "DE/selected_genes/compiled_filtered_selected_genes.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

# Load logCPM normalised dataset to use for the heatmaps
df2.rearranged <- read.delim(file="blade.genes.logCPM.normalised.values.rearranged.tab", sep="\t", stringsAsFactors = FALSE, header = TRUE, row.names = 1, dec= ".")
df2.rearranged <- t(scale(t(data.matrix(df2.rearranged))))
# get an index of rows from the full data that match the subset of genes
indexFD.rearranged <- which(rownames(df2.rearranged) %in% unique(c(degs.d13$Gene_ID, degs.d16$Gene_ID, degs.d20$Gene_ID, degs.d20.lt$Gene_ID)))  # 2765/6357 DEGs
df2.rearranged <- df2.rearranged[indexFD.rearranged,]

# Labels for the Subset of columns for the heatmap
labels.rearranged <- c(
  rep("WL - d13",2),   # Blade in WL
  rep("WL - d16",2),
  rep("WL - d20",2),
  rep(expression("EoD FR" ^ 1*" - d13"),2), # Blade in EoD FR06
  rep(expression("EoD FR" ^ 1*" - d16"),2),
  rep(expression("EoD FR" ^ 1*" - d20"),2)
)
rows.rearranged <- c(1:12)
labels.d20 <- c(
  rep("WL - d13",2),   # Blade in WL
  rep("WL - d16",2),
  rep("WL - d20",2),
  rep(expression("EoD FR"^1*" - d13"),2), # Blade in EoD FR06
  rep(expression("EoD FR"^1*" - d16"),2),
  rep(expression("EoD FR"^1*" - d20"),2),
  rep(expression("EoD FR"^2*" - d20"),2)  # Blade in EoD FR18
)
rows.d20 <- c(1:14)

# Graph each of the genes within the vector of gene lists
dir.create("heatmaps") 
dir.create("heatmaps/selected_genes/") 

for (i in seq_along(names(filtered.gene.lists))){
  genes_to_plot <- filtered.gene.lists[i]
  print(i)
  print(genes_to_plot)
  selected_dir <- paste0("heatmaps/selected_genes/")
  
  gene_subset <- unique(unlist(genes_to_plot)) # obtain all genes belonging to that category
  gene_subset <- gene_subset[which(gene_subset %in% rownames(df2.rearranged))] # obtain the subset that has been differentially regulated
  gene_subset <- gene_subset[which(gene_subset %in% c(degs.d13$Gene_ID, degs.d16$Gene_ID, degs.d20$Gene_ID, degs.d20.lt$Gene_ID))]
  
  if (length(gene_subset) >= 2){
    # Select subset of genes to be graphed
    df3 <- df2.rearranged[which(rownames(df2.rearranged) %in% gene_subset),1:14]
    # Remove rows with zero values
    # Go through each row and determine if a value is zero
    row_sub = apply(df3, 1, function(row) all(row = 0 ))
    # Subset rows with non-zero values
    df3 <- df3[!row_sub,]
    
    # Hierarchical clustering of rearranged logCPM data
    d <- dist(df3, method = "euclidean") # calculate euclidean distance
    clust <- hclust(d, method = "complete") # hierarchical cluster of the data
    dend <- as.dendrogram(clust) # as dendogram (for clustering rows, Rowv parameter of heatmap.2)
    
    # creates a own color palette from red to green
    my_palette <- colorRampPalette(c("blue", "#ffc000"))(n = 199)
    
    # (optional) defines the color breaks manually for a "skewed" color transition
    col_breaks = c(seq(min(df3),median(df3),length=100),  # for orange
                   seq(median(df3) + 0.1, max(df3),length=100))  # for blue
    
    png(paste0(selected_dir,
               "/",
               i,
               " - ",
               gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", names(filtered.gene.lists)[i]),
               ".png"),    # create png for the heat map        
        width = 6.5*1200,        
        if (length(rownames(df3)) > 5) {height = 5*1200
        } else {height = 3.5*1200},
        res = 1200,
        pointsize = 5)        # smaller font size
    heatmap.2(df3,
              main = paste0(gsub("_up","",gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", names(filtered.gene.lists)[i]))), # heat map title
              notecol="black",      # change font color of cell labels to black
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins =c(10,22),     # widens margins around plot
              col=my_palette,       # use on color palette defined earlier
              breaks=col_breaks,    # enable color transition at specified limits
              dendrogram="row",    # only draw a row dendrogram
              Colv="NA",            # turn off column clustering 
              Rowv = dend,
              #scale = "row",
              useRaster=TRUE,
              labRow = paste0(md$Symbol[match(rownames(df3),md$Gene_ID)], " (", rownames(df3), ") ", 
                              round_dec(degs.logFC$logFC.d13[match(rownames(df3),degs.logFC$Gene_ID)],1), " ",
                              round_dec(degs.logFC$logFC.d16[match(rownames(df3),degs.logFC$Gene_ID)],1), " ",
                              round_dec(degs.logFC$logFC.d20[match(rownames(df3),degs.logFC$Gene_ID)],1), " ",
                              round_dec(degs.logFC$logFC.d20.lt[match(rownames(df3),degs.logFC$Gene_ID)],1)) ,
              labCol = labels.d20,
              cexRow = 0.2 + 1/log10(length(rownames(df3))),
              cexCol = 0.2 + 1/log10(length(labels.rearranged)),
              ColSideColors = c(rep("gray",6),
                                rep("#A00000",6),
                                rep("#ff8080",2)),
              key.xlab = "Scaled logCPM")
    dev.off()
  }
  
  # Redo to graph logFC only
  gene_subset <- unique(unlist(genes_to_plot)) # obtain all genes belonging to that category
  gene_subset <- gene_subset[which(gene_subset %in% rownames(df2.rearranged))] # obtain the subset that has been differentially regulated
  gene_subset <- gene_subset[which(gene_subset %in% c(degs.d13$Gene_ID, degs.d16$Gene_ID, degs.d20$Gene_ID, degs.d20.lt$Gene_ID))] # subset to early treated plants
  
  if (length(gene_subset) >= 2){
    # Select subset of genes to be graphed
    df3 <- degs.logFC[which(degs.logFC$Gene_ID %in% gene_subset),1:5]
    rownames(df3) <- df3$Gene_ID
    df3 <- as.matrix(df3[,2:5])
    # Remove rows with zero values
    # Go through each row and determine if a value is zero
    row_sub = apply(df3, 1, function(row) all(row = 0 ))
    # Subset rows with non-zero values
    df3 <- df3[!row_sub,]
    
    # Hierarchical clustering of rearranged logCPM data
    d <- dist(df3, method = "euclidean") # calculate euclidean distance
    clust <- hclust(d, method = "complete") # hierarchical cluster of the data
    dend <- as.dendrogram(clust) # as dendogram (for clustering rows, Rowv parameter of heatmap.2)
    
    # creates a own color palette from blue to yellow
    my_palette <- c(colorRampPalette(c("blue", "#b3cde0"))(n = 99),
                    "#dddddd",
                    colorRampPalette(c("#FFF4B2", "#ffc000"))(n = 99))
    col_breaks = c(seq(-3,-0.58, length=99), # for blue
                   seq(-0.57, 0.57, length=1),  # for grey
                   seq(0.58, 3, length=100))  # for orange
    
    png(paste0(selected_dir,
               "/",
               i,
               " - ",
               gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", names(filtered.gene.lists)[i]),
               "_logFC",
               ".png"),    # create png for the heat map        
        width = 6.5*1200,        # 5 x 600 pixels
        if (length(rownames(df3)) > 5) {height = 5*1200
        } else {height = 3.5*1200},
        res = 1200,            # 600 pixels per inch
        #paper = "a4r",
        pointsize = 7)        # smaller font size
    heatmap.2(df3,
              main = paste0(gsub("_up","",gsub("([.|()\\^{}+$*?]|\\[|\\]|/)", "", names(filtered.gene.lists)[i]))), # heat map title
              notecol="black",      # change font color of cell labels to black
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins =c(10,22),     # widens margins around plot
              col=my_palette,       # use on color palette defined earlier
              breaks=col_breaks,    # enable color transition at specified limits
              dendrogram="row",    # only draw a row dendrogram
              Colv="NA",            # turn off column clustering 
              Rowv = dend,
              #scale = "row",
              useRaster=TRUE,
              labRow = paste0(md$Symbol[match(rownames(df3),md$Gene_ID)], " (", rownames(df3), ") ", 
                              round(degs.logFC$logFC.d13[match(rownames(df3),degs.logFC$Gene_ID)],1), " ",
                              round(degs.logFC$logFC.d16[match(rownames(df3),degs.logFC$Gene_ID)],1), " ",
                              round(degs.logFC$logFC.d20[match(rownames(df3),degs.logFC$Gene_ID)],1), " ",
                              round(degs.logFC$logFC.d20.lt[match(rownames(df3),degs.logFC$Gene_ID)],1)) ,
              labCol = c("day 13", "day 16", expression("day 20" ^ 06), expression("day 20" ^ 18)),
              cexRow = 0.2 + 1/log10(length(rownames(df3))),
              cexCol = 0.2 + 1/log10(length(labels.rearranged)),
              ColSideColors = c(rep("#A00000",3),
                                rep("#ff8080",1)),
              key.xlab = "logFC")
    dev.off()
  }
}


## Save Session Information
setwd(basedir)
sessionInfo()
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

#####################################################################################
#                               End of analysis                                     #
#####################################################################################