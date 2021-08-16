#---------------------------------------------------#
# prepare input data for MOVICS (tumor sample only) #

# set working path 
workdir <- "G:/BLCA_MOVICS/data_preparing"; setwd(workdir)

# load R package
library(TCGAbiolinks)
library(data.table)
library(ChAMPdata)
library(tidyverse)
library(magrittr)
library(readxl)
library(stringr)
library(forcats)

# customized functions
countToFpkm <- function(counts, effLen){
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
# load annotation
Ginfo <- read.table("overlapTable27.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

#----------------------#
# Clinical information #
clinical <- GDCquery(project = "TCGA-BLCA", 
                     data.category = "Clinical", 
                     file.type = "xml")
GDCdownload(clinical,directory = "GDCdata")
cliquery <- GDCprepare_clinic(clinical, clinical.info = "patient",directory = "GDCdata")
colnames(cliquery)[1] <- "Tumor_Sample_Barcode"
write.table(cliquery,"TCGA_BLCA_Clinical.txt",  quote=F, row.names=F,sep = "\t")

#----------------#
# RNA expression #

# count data
expquery <- GDCquery(project = "TCGA-BLCA", 
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification",
                     workflow.type = "HTSeq - Counts"
)
GDCdownload(expquery,directory = "GDCdata")
expquery2 <- GDCprepare(expquery,directory = "GDCdata",summarizedExperiment = T)
save(expquery2,file = "blca.gdc.RData")
expMatrix <- TCGAanalyze_Preprocessing(expquery2)
colnames(expMatrix) <- substr(colnames(expMatrix), start = 1,stop = 16)
normsamples <- colnames(expMatrix)[which(substr(colnames(expMatrix),14,15) == "11")] # get normal samples
tumorsamples <- colnames(expMatrix)[which(substr(colnames(expMatrix),14,15) == "01")] # get tumor samples
expMatrix <- expMatrix[,tumorsamples]; expMatrix <- as.data.frame(expMatrix[rowSums(expMatrix) > 0,])
comgene <- intersect(rownames(expMatrix),rownames(Ginfo))
count <- as.data.frame(expMatrix)[comgene,]; Ginfo <- Ginfo[comgene,]

# extract mRNA and lncRNA
Lids <- Ginfo[Ginfo$genetype %in% c("non_coding","3prime_overlapping_ncRNA","antisense_RNA","lincRNA","sense_intronic","sense_overlapping","macro_lncRNA","bidirectional_promoter_lncRNA"),]
Mids <- Ginfo[Ginfo$genetype == "protein_coding",]

count <- count[c(rownames(Mids),rownames(Lids)),]

# FPKM
fpkm.mrna <- as.data.frame(apply(as.matrix(count[rownames(Mids),]), 2, countToFpkm, effLen = Ginfo[rownames(Mids),"unqlen"]))
fpkm.lncrna <- as.data.frame(apply(as.matrix(count[rownames(Lids),]), 2, countToFpkm, effLen = Ginfo[rownames(Lids),"unqlen"]))

# TPM
tpm.mrna <- as.data.frame(round(apply(fpkm.mrna,2,fpkmToTpm),2))
tpm.mrna$Gene <- Ginfo[rownames(tpm.mrna),"genename"]
tpm.mrna <- apply(tpm.mrna[,setdiff(colnames(tpm.mrna), "Gene")], 2, function(x) tapply(x, INDEX=factor(tpm.mrna$Gene), FUN=median, na.rm=TRUE))

tpm.lncrna <- as.data.frame(round(apply(fpkm.lncrna,2,fpkmToTpm),2))
tpm.lncrna$Gene <- Ginfo[rownames(tpm.lncrna),"genename"]
tpm.lncrna <- apply(tpm.lncrna[,setdiff(colnames(tpm.lncrna), "Gene")], 2, function(x) tapply(x, INDEX=factor(tpm.lncrna$Gene), FUN=median, na.rm=TRUE))

write.table(tpm.mrna, "TCGA_BLCA_mRNA_TPM.txt", quote=F, row.names=T,col.names = NA,sep = "\t")
write.table(tpm.lncrna, "TCGA_BLCA_lncRNA_TPM.txt", quote=F, row.names=T,col.names = NA,sep = "\t")

countsTable <- count
countsTable$Gene <- Ginfo[rownames(count),"genename"]
countsTable <- countsTable[!duplicated(countsTable$Gene),]
rownames(countsTable) <- countsTable$Gene
countsTable <- countsTable[,-ncol(countsTable)]
write.table(countsTable, "TCGA_BLCA_Count.txt", quote=F, row.names=T,col.names = NA,sep = "\t")

# normal sample
tmp <- TCGAanalyze_Preprocessing(expquery2)
colnames(tmp) <- substr(colnames(tmp), start = 1,stop = 16)
tmp <- tmp[,normsamples]; tmp <- as.data.frame(tmp[rownames(count),])
fpkm.mrna.norm <- as.data.frame(apply(as.matrix(tmp[rownames(Mids),]), 2, countToFpkm, effLen = Ginfo[rownames(Mids),"unqlen"]))
tpm.mrna.norm <- as.data.frame(round(apply(fpkm.mrna.norm,2,fpkmToTpm),2))
tpm.mrna.norm$Gene <- Ginfo[rownames(tpm.mrna.norm),"genename"]
tpm.mrna.norm <- apply(tpm.mrna.norm[,setdiff(colnames(tpm.mrna.norm), "Gene")], 2, function(x) tapply(x, INDEX=factor(tpm.mrna.norm$Gene), FUN=median, na.rm=TRUE))
write.table(tpm.mrna.norm,"TCGA_BLCAnorm_mRNA_TPM.txt", quote=F, row.names=T,col.names = NA,sep = "\t")

#-----------------#
# DNA methylation #

# please download this in XENA (GDC TCGA Bladder Cancer (BLCA))
# https://xenabrowser.net/datapages/?dataset=TCGA-BLCA.methylation450.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
meth <- fread("TCGA-BLCA.methylation450.tsv.gz",sep = "\t",check.names = F,stringsAsFactors = F,header = T)
meth <- as.data.frame(meth); rownames(meth) <- meth[,1]; meth <- as.data.frame(na.omit(meth[,-1]))
colnames(meth) <- substr(colnames(meth), start = 1,stop = 16)
meth <- meth[,which(substr(colnames(meth), 14, 15) == "01")]

# get promoter
data("probe.features")
promoter <- intersect(rownames(probe.features[which(probe.features$cgi == "island" & probe.features$feature %in% c("TSS1500","TSS200")),]),rownames(meth))
meth <- round(meth[promoter,],3)

# map to gene name to take median value if multiple mapping
meth$gene <- as.character(probe.features[promoter,"gene"])
meth <- apply(meth[,setdiff(colnames(meth), "gene")], 2, function(x) tapply(x, INDEX=factor(meth$gene), FUN=median, na.rm=TRUE)) # be patient because this will take a while
meth <- as.data.frame(round(meth,3))

#---------------------#
# Copy number segment #

# download from firehose http://firebrowse.org/?cohort=BLCA&download_dialog=true
segment <- read.table("BLCA.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
colnames(segment) <- c("sample","chromosome","start","end","num.mark","seg.mean")
segment$sample <- substr(segment$sample, start = 1,stop = 16)
segment <- segment[which(substr(segment$sample, 14, 15) == "01"),]
write.table(segment, "TCGA_BLCA_segment.txt", quote=F, row.names=F,col.names = T,sep = "\t")

# cnv <- read.delim("TCGA.BLCA.sampleMap_Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz",sep = "\t",header = T,check.names = F,stringsAsFactors = F)
library(iClusterPlus)
cna <- CNregions(segment,
                 epsilon = 0,
                 adaptive = FALSE,
                 rmCNV = FALSE,
                 cnv = FALSE,
                 frac.overlap = 0.5, 
                 rmSmallseg = TRUE, 
                 nProbes = 15)
cna <- as.data.frame(t(cna))

#-----#
# MAF #

# please download this from cBioPortal under Data Sets archive
# https://www.cbioportal.org/datasets
maf <- read_tsv("data_mutations_extended.txt", comment = "#")
flag <- read.table("Mutation Flags 100.txt",sep = "\t")$V1
label <- c("Tumor_Sample_Barcode",
           "Hugo_Symbol",
           "Chromosome",
           "Start_Position",
           "End_Position",
           "Variant_Classification",
           "Variant_Type",
           "Reference_Allele",
           "Tumor_Seq_Allele1",
           "Tumor_Seq_Allele2")
maf <- maf[-which(maf$Hugo_Symbol %in% flag),label]
maf$Hugo_Symbol <- toupper(maf$Hugo_Symbol) # transfer gene name to capital. (eg, C1orf198 to C1ORF198)
write.table(maf, "TCGA_BLCA_MAF.txt", quote=F, row.names=F,col.names = T,sep = "\t")

#-------------------------#
# binary somatic mutation #

# binary mutation
mut.binary <- matrix(0,nrow = length(unique(maf$Hugo_Symbol)),ncol = length(unique(maf$Tumor_Sample_Barcode)),dimnames = list(unique(maf$Hugo_Symbol),unique(maf$Tumor_Sample_Barcode)))
for (i in colnames(mut.binary)) {
  tmp <- maf[which(maf$Tumor_Sample_Barcode == i),]
  tmp <- tmp[which(tmp$Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation")),]
  for (j in tmp$Hugo_Symbol)
    mut.binary[j,i] <- 1
}
mut.binary <- as.data.frame(mut.binary)
mut.binary <- mut.binary[rowSums(mut.binary) > 0,]
mut.binary <- mut.binary[,which(substr(colnames(mut.binary), 14, 15) == "01")]
colnames(mut.binary) <- paste0(colnames(mut.binary),"A")

# load survival data
surv.info <- read.table("pancancerSurvivalData_XLu.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
surv.info <- as.data.frame(surv.info[which(surv.info$type == "BLCA"),c(3,6,25:36)])
surv.info <- surv.info[which(surv.info$OS.time > 0),] # get survival time greater than 0
rownames(surv.info) <- paste0(rownames(surv.info),"-01A")

clin.info <- read.delim("BLCA_ClinicalInfo_forDavid.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
clin.info <- clin.info[which(clin.info$SampleType == "PrimaryTumor"),]
rownames(clin.info) <- paste0(clin.info$bcr_patient_barcode,"-",sapply(strsplit(clin.info$sampleID_XSU,"-",fixed = T),"[",3))

# get common samples
comsam <- intersect(colnames(tpm.mrna),colnames(mut.binary))
comsam <- intersect(comsam, colnames(meth))
comsam <- intersect(comsam, rownames(surv.info))
comsam <- intersect(comsam,rownames(clin.info))
comsam <- intersect(comsam, colnames(cna))

surv.info <- cbind.data.frame(surv.info[comsam,],clin.info[comsam,])

# output file
write.table(tpm.mrna[,comsam], "TCGA_BLCA_mRNA_TPM.txt", quote=F, row.names=T,col.names = NA,sep = "\t")
write.table(tpm.lncrna[,comsam], "TCGA_BLCA_lncRNA_TPM.txt", quote=F, row.names=T,col.names = NA,sep = "\t")
write.table(meth[,comsam], "TCGA_BLCA_Methpromoter.txt", quote=F, row.names=T,col.names = NA,sep = "\t")
write.table(mut.binary[,comsam],"TCGA_BLCA_mut2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(cna[,comsam],"TCGA_BLCA_cna.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(surv.info[comsam,],"TCGA_BLCA_surv.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

# save RData
save.image("data_preparing.RData")
