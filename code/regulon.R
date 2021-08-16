# load R package
library(RTN)
library(snow)

# load expression and phenotype
tcgaBLCA <- read.table("tcga_tpm_forRegulon.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
pheno <- read.table("tcga_phenotype_forRegulon.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

# load TFs
tfs <- read.table("regulon.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
chr.tfs1 <- read.delim("Appendix A. AHAT mutations in cancer.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
chr.tfs2 <- read.delim("Appendix B.Histone deacetylases in cancer.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
chr.tfs3 <- read.delim("Appendix C.Selected histone methyltransferases in cancer.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
chr.tfs4 <- read.delim("Appendix D.Histone demethylases in cancer.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
chr.tfs <- c(unlist(strsplit(chr.tfs1$Gene,"/",fixed = T)),unlist(strsplit(chr.tfs1$`Common name`,"/",fixed = T)),
             unlist(strsplit(chr.tfs2$Name,"/",fixed = T)),
             unlist(strsplit(chr.tfs3$Name,"/",fixed = T)),unlist(strsplit(chr.tfs3$Synonyms,"/",fixed = T)),
             unlist(strsplit(chr.tfs4$Name,"/",fixed = T)),unlist(strsplit(chr.tfs4$Synonyms,"/",fixed = T)))
chr.tfs <- c(chr.tfs, unlist(strsplit(chr.tfs3$Synonyms,", ",fixed = T)), unlist(strsplit(chr.tfs4$Synonyms,", ",fixed = T)))
chr.tfs <- unlist(strsplit(chr.tfs,", ",fixed = T))
chr.tfs <- unlist(strsplit(chr.tfs,"/",fixed = T))
chr.tfs <- unique(chr.tfs)
chr.tfs <- c(chr.tfs,"KDM7A","JHDM1D")
chr.tfs <- setdiff(chr.tfs,c("JHDM1D (KDM7A)","MYST family","GNAT family","Orphan family"))

# get intersection of TFs
regulatoryElements <- intersect(colnames(tfs), rownames(tcgaBLCA))

# Run the TNI constructor
rtni_tcgaBLCA <- tni.constructor(expData = as.matrix(tcgaBLCA), 
                                 regulatoryElements = regulatoryElements)

# Compute the reference regulatory network by permutation and bootstrap analyses.
# Please set 'spec' according to your available hardware
options(cluster=snow::makeCluster(spec=12, "SOCK"))
rtni_tcgaBLCA <- tni.permutation(rtni_tcgaBLCA, pValueCutoff = 1e-5, nPermutations = 1000)
rtni_tcgaBLCA <- tni.bootstrap(rtni_tcgaBLCA, nBootstraps = 1000)
stopCluster(getOption("cluster"))

# Compute the DPI-filtered regulatory network
rtni_tcgaBLCA <- tni.dpi.filter(rtni_tcgaBLCA, eps = 0, sizeThreshold = TRUE, minRegulonSize = 15)

# Save the TNI object for subsequent analyses
save(rtni_tcgaBLCA, file="rtni_tcgaBLCA.RData")

# Compute regulon activity for individual samples
rtnigsea_tcgaBLCA <- tni.gsea2(rtni_tcgaBLCA, regulatoryElements = regulatoryElements)
MIBC_regact <- tni.get(rtnigsea_tcgaBLCA, what = "regulonActivity")
save(MIBC_regact,file = "MIBC_regact.RData")

# get intersection of TFs
regulatoryElements <- intersect(chr.tfs, rownames(tcgaBLCA))

# Run the TNI constructor
rtni_tcgaBLCA2 <- tni.constructor(expData = as.matrix(tcgaBLCA), 
                                  regulatoryElements = regulatoryElements)

# Compute the reference regulatory network by permutation and bootstrap analyses.
# Please set 'spec' according to your available hardware
options(cluster=snow::makeCluster(spec=12, "SOCK"))
rtni_tcgaBLCA2 <- tni.permutation(rtni_tcgaBLCA2, pValueCutoff = 1e-5, nPermutations = 1000)
rtni_tcgaBLCA2 <- tni.bootstrap(rtni_tcgaBLCA2, nBootstraps = 1000)
stopCluster(getOption("cluster"))

# Compute the DPI-filtered regulatory network
rtni_tcgaBLCA2 <- tni.dpi.filter(rtni_tcgaBLCA2, eps = 0, sizeThreshold = TRUE, minRegulonSize = 15)

# Save the TNI object for subsequent analyses
save(rtni_tcgaBLCA2, file="rtni_tcgaBLCA2.RData")

# Compute regulon activity for individual samples
rtnigsea_tcgaBLCA2 <- tni.gsea2(rtni_tcgaBLCA2, regulatoryElements = regulatoryElements)
chromatinRemodeling_regact <- tni.get(rtnigsea_tcgaBLCA2, what = "regulonActivity")
save(chromatinRemodeling_regact,file = "chromatinRemodeling_regact.RData")
