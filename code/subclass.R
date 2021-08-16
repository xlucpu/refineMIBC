# submap using SKCM dataset
generateInputFileForSubMap <- function(in_gct, gct_file, cls_file, sam_info, type_name = "type"){
  in_gct <- data.frame(GeneID=rownames(in_gct),
                       description="na",
                       in_gct, 
                       stringsAsFactors = F,
                       check.names = F)
  cat("#1.2\n", file = gct_file)
  cat(nrow(in_gct),"\t",ncol(in_gct)-2,"\n", file = gct_file, append = T)
  cat(paste(colnames(in_gct), collapse = "\t"),"\n", file = gct_file, append = T)
  for(i in 1:nrow(in_gct)) cat(paste(in_gct[i,], collapse = "\t"),"\n", file = gct_file, append = T)
  
  cat(nrow(sam_info),length(levels(factor(sam_info$rank))),1, "\n", file = cls_file )
  cat("#", paste0(levels(factor(sam_info[, type_name])), collapse = " " ), "\n", file = cls_file, sep = "", append = T)
  cat(as.numeric(factor(sam_info[, type_name])), file = cls_file, append = T)
}

GENELIST <- intersect(rownames(tpm),rownames(skcm.immunotherapy.logNC)) 
GENELIST <- intersect(GENELIST,rownames(combined.expr.illumina.combat))
GENELIST <- intersect(GENELIST,rownames(combined.expr.affy.combat))

skcm.immunotherapy.logNC <- read.table(file.path(data.path,"skcm.immunotherapy.47samples.log2CountsNorm.txt"),sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F) #原文提供的log2转化的标准化count值
rownames(skcm.immunotherapy.logNC) <- toupper(rownames(skcm.immunotherapy.logNC)) # 基因大写，因为我使用的数据是把基因名都大写的
skcm.immunotherapy.info <- read.table(file.path(data.path,"skcm.immunotherapy.47sampleInfo.txt"),sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
skcm.immunotherapy.info <- skcm.immunotherapy.info[order(skcm.immunotherapy.info$label),]
skcm.immunotherapy.info$rank <- rep(c(1,2,3,4),times=as.character(table(skcm.immunotherapy.info$label))) #1: CTLA4_noR 2: CTLA4_R 3:PD1_noR 4:PD1_R
sam_info <- skcm.immunotherapy.info
in_gct <- skcm.immunotherapy.logNC[GENELIST,rownames(skcm.immunotherapy.info)]
gct_file <- file.path(res.path,"skcm.immunotherapy.for.SubMap.gct")
cls_file <- file.path(res.path,"skcm.immunotherapy.for.SubMap.cls")
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

# tcga
sam_info <- data.frame("CMOIC"=blca.sinfo$CMOIC,row.names = rownames(blca.sinfo))
sam_info$rank <- as.numeric(gsub("CS","",sam_info$CMOIC))
sam_info <- sam_info[order(sam_info$rank),]
gct_file <- file.path(res.path,"tcga.cmoic.for.SubMap.gct")
cls_file <- file.path(res.path,"tcga.cmoic.for.SubMap.cls")
in_gct <- log2(tpm[GENELIST,rownames(sam_info)] + 1)
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

# illumina
sam_info <- data.frame("CMOIC"=combined.surv.illumina$NTP,row.names = rownames(combined.surv.illumina))
sam_info$rank <- as.numeric(gsub("CS","",sam_info$CMOIC))
sam_info <- sam_info[order(sam_info$rank),]
gct_file <- file.path(res.path,"illumina.cmoic.for.SubMap.gct")
cls_file <- file.path(res.path,"illumina.cmoic.for.SubMap.cls")
in_gct <- combined.expr.illumina.combat[GENELIST,rownames(sam_info)]
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

# affy
sam_info <- data.frame("CMOIC"=combined.surv.affy$NTP,row.names = rownames(combined.surv.affy))
sam_info$rank <- as.numeric(gsub("CS","",sam_info$CMOIC))
sam_info <- sam_info[order(sam_info$rank),]
gct_file <- file.path(res.path,"affy.cmoic.for.SubMap.gct")
cls_file <- file.path(res.path,"affy.cmoic.for.SubMap.cls")
in_gct <- combined.expr.affy.combat[GENELIST,rownames(sam_info)]
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

# submap using IMVIGOR210 dataset
source("twoclasslimma.R")
imm.backgroud <- read.table(file.path(data.path,"Immune.gene.background.txt"),sep = "\t",row.names = NULL,header = T,check.names = F)
imvigor210.expr <- read.table(file.path(data.path,"cds.tpm.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
imvigor210.surv <- read.table(file.path(data.path,"cds.pData.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

cds.blca <- imvigor210.surv[which(imvigor210.surv$`Best Confirmed Overall Response` %in% c("CR","PD")),]
saminfo <- data.frame(condition = as.factor(cds.blca$`Best Confirmed Overall Response`),
                      row.names = rownames(cds.blca),
                      stringsAsFactors = F)
saminfo <- saminfo[order(saminfo$condition),,drop = F]
twoclasslimma(subtype  = saminfo, # subtype information (must contain a column named 'condition')
              featmat  = imvigor210.expr[intersect(rownames(imvigor210.expr),imm.backgroud$`HUGO Name`),rownames(saminfo)], # expression file (fill detect data scale automatically)
              treatVar = "CR", # name of treatment group
              ctrlVar  = "PD", # name of control group
              prefix   = "imvigor210", # prefix of the DE file
              overwt   = T, # whether over write files that already existed
              sort.p   = TRUE, # if sorting the results by adjusted p value
              verbose  = TRUE, # if showing verbose result
              res.path = res.path) # path for result

blca.immunotherapy.logNC <- imvigor210.expr[intersect(rownames(imvigor210.expr),imm.backgroud$`HUGO Name`),rownames(saminfo)]
blca.immunotherapy.info <- read.table(file.path(data.path,"blca.imvrgor210.sampleInfo.txt"),sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
blca.immunotherapy.info <- blca.immunotherapy.info[order(blca.immunotherapy.info$label),]
blca.immunotherapy.info$rank <- rep(c(1,2),times=as.character(table(blca.immunotherapy.info$label))) #1:PDL1_noR 2:PDL1_R
sam_info <- blca.immunotherapy.info

GENELIST <- intersect(rownames(tpm),rownames(blca.immunotherapy.logNC)) 
GENELIST <- intersect(GENELIST,rownames(combined.expr.illumina.combat)) 
GENELIST <- intersect(GENELIST,rownames(combined.expr.affy.combat))

in_gct <- log2(blca.immunotherapy.logNC[GENELIST,rownames(sam_info)] + 1)
gct_file <- file.path(res.path,"IMvorblca.immunotherapy.for.SubMap.gct")
cls_file <- file.path(res.path,"IMvorblca.immunotherapy.for.SubMap.cls")
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

# tcga
sam_info <- data.frame("CMOIC"=blca.sinfo$CMOIC,row.names = rownames(blca.sinfo))
sam_info$rank <- as.numeric(gsub("CS","",sam_info$CMOIC))
sam_info <- sam_info[order(sam_info$rank),]
gct_file <- file.path(res.path,"tcga.cmoic.for.SubMap2.gct")
cls_file <- file.path(res.path,"tcga.cmoic.for.SubMap2.cls")
in_gct <- log2(tpm[GENELIST,rownames(sam_info)] + 1)
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

# illumina
sam_info <- data.frame("CMOIC"=combined.surv.illumina$NTP,row.names = rownames(combined.surv.illumina))
sam_info$rank <- as.numeric(gsub("CS","",sam_info$CMOIC))
sam_info <- sam_info[order(sam_info$rank),]
gct_file <- file.path(res.path,"illumina.cmoic.for.SubMap2.gct")
cls_file <- file.path(res.path,"illumina.cmoic.for.SubMap2.cls")
in_gct <- combined.expr.illumina.combat[GENELIST,rownames(sam_info)]
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

# affy
sam_info <- data.frame("CMOIC"=combined.surv.affy$NTP,row.names = rownames(combined.surv.affy))
sam_info$rank <- as.numeric(gsub("CS","",sam_info$CMOIC))
sam_info <- sam_info[order(sam_info$rank),]
gct_file <- file.path(res.path,"affy.cmoic.for.SubMap2.gct")
cls_file <- file.path(res.path,"affy.cmoic.for.SubMap2.cls")
in_gct <- combined.expr.affy.combat[GENELIST,rownames(sam_info)]
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")
