#-----------------#
# MOVICS pipeline #

# set working path and creat directory
tumor.path <- "G:/BLCA_MOVICS/movics_pipeline"; setwd(tumor.path) #create dir
res.path    <- file.path(tumor.path, "Results")
fig.path    <- file.path(tumor.path, "Figures")
data.path    <- file.path(tumor.path, "InputData")

if (!file.exists(tumor.path)) { dir.create(tumor.path) }
if (!file.exists(res.path)) { dir.create(res.path) }
if (!file.exists(fig.path)) { dir.create(fig.path) }
if (!file.exists(data.path)) { dir.create(data.path) }

# set colors
clust.col <- c("#DD492E","#40548A","#32A087","#EC7D21")

# load R package
library(MOVICS)
library(wateRmelon)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(magrittr)
library(readxl)
library(stringr)
library(forcats)
library(consensusMIBC)
library(preprocessCore)
library(ggpubr)
library(survival)
library(survminer)
library(survRM2)
library(ggpmisc)
library(clusterProfiler)
library(enrichplot)
library(estimate)
library(gplots)
library(ClassDiscovery)
library(ComplexHeatmap)
library(GSVA)
library(ggalluvial)

# customized fucntion
display.progress = function (index, totalN, breakN=20) {
  
  if ( index %% ceiling(totalN/breakN)  ==0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
  
}    

standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}

Ginfo <- read.table(file.path(data.path,"overlapTable27.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
Lids <- Ginfo[Ginfo$genetype %in% c("non_coding","3prime_overlapping_ncRNA","antisense_RNA","lincRNA","sense_intronic","sense_overlapping","macro_lncRNA","bidirectional_promoter_lncRNA"),"genename"]
Mids <- Ginfo[Ginfo$genetype == "protein_coding","genename"]

# load data
mrna.tpm <- read.table(file.path(data.path,"TCGA_BLCA_mRNA_TPM.txt"),sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
lncrna.tpm <- read.table(file.path(data.path,"TCGA_BLCA_lncRNA_TPM.txt"),sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
meth <- read.table(file.path(data.path,"TCGA_BLCA_Methpromoter.txt"),sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
mut <- read.table(file.path(data.path,"TCGA_BLCA_mut2.txt"),sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
cna <- read.table(file.path(data.path,"TCGA_BLCA_cna.txt"),sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
surv <- read.delim(file.path(data.path,"TCGA_BLCA_surv.txt"),sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
gc()

count <- read.table(file.path(data.path,"TCGA_BLCA_Count.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
tpm <- rbind.data.frame(mrna.tpm,lncrna.tpm)
segment <- read.table(file.path(data.path,"TCGA_BLCA_segment.txt"),row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
maf <- read_tsv(file.path(data.path,"TCGA_BLCA_MAF.txt"), comment = "#")
maf$Tumor_Sample_Barcode <- paste0(maf$Tumor_Sample_Barcode,"A")

# filter elites
mo.mrna <- getElites(dat = mrna.tpm,
                     method = "mad",
                     doLog2 = T,
                     lowpct = 0.1,
                     elite.num = 1500, # get top 1000 genes according mad values
                     scaleFlag = F,
                     centerFlag = F)

mo.lncrna <- getElites(dat = lncrna.tpm,
                       method = "mad",
                       doLog2 = T,
                       lowpct = 0.1,
                       elite.num = 1500, # get top 1000 genes according mad values
                       scaleFlag = F,
                       centerFlag = F)

mo.meth <- getElites(dat = meth,
                     method = "mad",
                     doLog2 = F,
                     elite.num = 1500, # get top 1000 genes according mad values
                     scaleFlag = F,
                     centerFlag = F)

mo.mut <- getElites(dat = mut,
                    method = "freq",
                    doLog2 = F,
                    elite.pct = 0.03, # get mutations that have mutated rate greater than 3%
                    scaleFlag = F,
                    centerFlag = F)

mo.cna <- getElites(dat = cna,
                    method = "mad",
                    doLog2 = F,
                    elite.num = 1500,
                    scaleFlag = F,
                    centerFlag = F)

# create input for MOVICS (very IMPORTANT!!!)
BLCA.tcga <- list(mRNA      = mo.mrna$elite.dat,
                  lncRNA    = mo.lncrna$elite.dat,
                  meth      = mo.meth$elite.dat,
                  cna       = mo.cna$elite.dat,
                  mut       = mo.mut$elite.dat,
                  count     = count,
                  tpm       = tpm,
                  maf       = maf,
                  segment   = segment,
                  surv.info = surv)

# print name of example data
names(BLCA.tcga)
# [1] "mRNA"      "lncRNA"    "meth"      "mut"       "cna"       "count"     "tpm"       "maf"       "segment"   "surv.info"

# extract mo.data
mo.data <- BLCA.tcga[1:5]

#-------------------------------------------------------#
# identify optimal clustering number (may take a while) #
optk.BLCA <- getClustNum(data        = mo.data,
                         is.binary   = c(F,F,F,F,T), # note: the last data is somatic mutation which is a binary matrix
                         try.N.clust = 2:6, # try cluster number from 2 to 6
                         fig.path    = fig.path,
                         fig.name    = "CLUSTER NUMBER OF TCGA-BLCA")

#-------------------------------------------------------------------#
# perform multi-omics integrative clustering with the 10 algorithms #
moic.res.list <- getMOIC(data        = mo.data,
                         methodslist = list("iClusterBayes","SNF", "PINSPlus", "NEMO", "COCA", "LRAcluster", "ConsensusClustering", "IntNMF", "CIMLR", "MoCluster"),
                         N.clust     = 4,
                         type        = c("gaussian", "gaussian", "gaussian", "gaussian", "binomial"))
save(moic.res.list, file = file.path(res.path,"moic.res.list.rda"))

iClusterBayes <- moic.res.list$iClusterBayes$clust.res
SNF <- moic.res.list$SNF$clust.res
PINSPlus <- moic.res.list$PINSPlus$clust.res
NEMO <- moic.res.list$NEMO$clust.res
COCA <- moic.res.list$COCA$clust.res
LRAcluster <- moic.res.list$LRAcluster$clust.res
ConsensusClustering <- moic.res.list$ConsensusClustering$clust.res
IntNMF <- moic.res.list$IntNMF$clust.res
CIMLR <- moic.res.list$CIMLR$clust.res
MoCluster <- moic.res.list$MoCluster$clust.res

moic.res <- cbind.data.frame(iClusterBayes = paste0("CS",iClusterBayes$clust),
                             SNF = paste0("CS",SNF$clust),
                             PINSPlus = paste0("CS",PINSPlus$clust),
                             NEMO = paste0("CS",NEMO$clust),
                             COCA = paste0("CS",COCA$clust),
                             LRAcluster = paste0("CS",LRAcluster$clust),
                             ConsensusClustering = paste0("CS",ConsensusClustering$clust),
                             IntNMF = paste0("CS",IntNMF$clust),
                             CIMLR = paste0("CS",CIMLR$clust),
                             MoCluster = paste0("CS",MoCluster$clust))
rownames(moic.res) <- rownames(iClusterBayes)

#---------------------------------------------------#
# get consensus clustering from different algorithm #
cmoic.BLCA <- getConsensusMOIC(moic.res.list = moic.res.list,
                               fig.name      = "CONSENSUS HEATMAP",
                               fig.path      = fig.path,
                               distance      = "euclidean",
                               linkage       = "ward.D2")

#------------------#
# compare survival #
surv.info <- BLCA.tcga$surv.info[,c("OS","OS.time")]
colnames(surv.info) <- c("fustat","futime")
surv.BLCA <- compSurv(moic.res         = cmoic.BLCA,
                      surv.info        = surv.info,
                      convt.time       = "m", # convert day unit to month
                      surv.median.line = "h", # draw horizontal line at median survival
                      xyrs.est         = c(5,10), # estimate 5 and 10-year survival
                      p.adjust.method  = "none",
                      clust.col        = clust.col,
                      fig.path         = fig.path,
                      fig.name         = "KAPLAN-MEIER CURVE OF CONSENSUSMOIC")

surv.info2 <- BLCA.tcga$surv.info[,c("PFI","PFI.time")]
colnames(surv.info2) <- c("fustat","futime")
surv.BLCA <- compSurv(moic.res         = cmoic.BLCA,
                      surv.info        = surv.info2,
                      convt.time       = "m", # convert day unit to month
                      surv.median.line = "h", # draw horizontal line at median survival
                      xyrs.est         = c(5,10), # estimate 5 and 10-year survival
                      p.adjust.method  = "none",
                      clust.col        = clust.col,
                      fig.path         = fig.path,
                      fig.name         = "KAPLAN-MEIER CURVE OF CONSENSUSMOIC PFS")
#---------------------------#
# compare clinical features #
clin.BLCA <- compClinvar(moic.res      = cmoic.BLCA,
                         var2comp      = annCol, # data.frame needs to summarize (must has row names of samples)
                         strata        = "Subtype", # stratifying variable (e.g., Subtype in this example)
                         factorVars    = colnames(annCol), # features that are considered categorical variables
                         nonnormalVars = NULL, # feature(s) that are considered using nonparametric test
                         exactVars     = colnames(annCol), # feature(s) that are considered using exact test
                         doWord        = TRUE, # generate .docx file in local path
                         res.path      = res.path,
                         includeNA     = FALSE,
                         tab.name      = "SUMMARIZATION OF CLINICAL FEATURES")

#---------------------------------#
# mutational frequency comparison #
mut.BLCA <- compMut(moic.res     = cmoic.BLCA,
                    mut.matrix   = BLCA.tcga$mut, # binary somatic mutation matrix
                    doWord       = TRUE, # generate table in .docx format
                    doPlot       = TRUE, # draw OncoPrint
                    freq.cutoff  = 0.03, # keep those genes that mutated in at least 3% of samples
                    p.cutoff     = 0.05, # keep those genes with nominal p value < 0.25 to draw OncoPrint
                    p.adj.cutoff = 0.25, # keep those genes with adjusted p value < 1 to draw OncoPrint
                    innerclust   = TRUE, # perform clustering within each subtype
                    annCol       = annCol, # same annotation for heatmap
                    annColors    = annColors, # same annotation color for heatmap
                    width        = 8, 
                    height       = 5,
                    fig.name     = "ONCOPRINT FOR SIGNIFICANT MUTATIONS",
                    tab.name     = "INDEPENDENT TEST BETWEEN SUBTYPE AND MUTATION",
                    res.path     = res.path,
                    fig.path     = fig.path)

#-------------------------------#
# compare total mutation burden #
tmb.BLCA <- compTMB(moic.res     = cmoic.BLCA,
                    maf          = BLCA.tcga$maf,
                    nonSyn       = NULL,
                    rmDup        = TRUE, # remove duplicated variants per sample
                    rmFLAGS      = FALSE, # keep FLAGS mutations
                    exome.size   = 38, # estimated exome size
                    test.method  = "nonparametric", # statistical testing method
                    fig.path     = fig.path,
                    fig.name     = "DISTRIBUTION OF TMB AND TITV")

#---------------------------------#
# compare fraction genome altered #
# change column names of segment data
segment2 <- segment[,setdiff(colnames(segment),"num.mark")]
colnames(segment2) <- c("sample","chrom","start","end","value")

# compare FGA, FGG, and FGL
fga.BLCA <- compFGA(moic.res     = cmoic.BLCA,
                    segment      = segment2,
                    iscopynumber = FALSE, # this is a segmented copy number file
                    cnathreshold = 0.3, # threshold to determine CNA gain or loss
                    test.method  = "nonparametric", # statistical testing method
                    clust.col    = clust.col,
                    width = 8,height = 3,
                    fig.path     = fig.path,
                    fig.name     = "BARPLOT OF FGA")

#--------------------------------------#
# run differential expression analysis #
# run DEA with limma
runDEA(dea.method = "limma",
       expr       = BLCA.tcga$tpm[intersect(rownames(BLCA.tcga$tpm),Mids),], # raw count data
       moic.res   = cmoic.BLCA,
       res.path   = res.path,
       prefix     = "TCGA-BLCA") # prefix of figure name

#----------------------------------------#
# run biomarker identification procedure #
# choose limma result to identify subtype-specific up-regulated biomarkers
marker.up <- runMarker(moic.res      = cmoic.BLCA,
                       dea.method    = "limma", # name of DEA method
                       prefix        = "TCGA-BLCA", # MUST be the same of argument in runDEA()
                       dat.path      = res.path, # path of DEA files
                       res.path      = res.path, # path to save marker files
                       p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                       dirct         = "up", # direction of dysregulation in expression
                       n.marker      = 30, # number of biomarkers for each subtype
                       doplot        = TRUE, # generate diagonal heatmap
                       norm.expr     = tpm, # use normalized expression as heatmap input
                       annCol        = annCol, # sample annotation in heatmap
                       annColors     = annColors, # colors for sample annotation
                       show_rownames = FALSE, # show no rownames (biomarker name)
                       fig.path      = fig.path,
                       width         = 12,
                       height        = 12,
                       fig.name      = "UPREGULATED BIOMARKER HEATMAP")
