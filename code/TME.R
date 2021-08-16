#---------------------------------------------#
# immune microenvironment enrichment analysis #
library(GSVA)
library(ComplexHeatmap)
library(viridis)
library(gplots)
library(estimate)
immune.signature <- read.table(file.path(data.path,"CCR_Curated_Immune_Cell_Signature.txt"),sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = NULL)
cell.type <- unique(immune.signature$CellType)
immune.sig.ccr <- list()
for (i in cell.type) {
  immune.sig.ccr[[i]] <- immune.signature[which(immune.signature$CellType == i),"Symbol"]
}
imm.targets <- c("CD274","PDCD1","CD247","PDCD1LG2","CTLA4","TNFRSF9","TNFRSF4","TLR9")
immune.sig.ccr.order <- c("T.cells.CD8",
                          "T.cells.regulatory..Tregs.",
                          "T.cells.CD4.naive",
                          "T.cells.follicular.helper",
                          "B.cells.naive",
                          "B.cells.memory",
                          "T.cells.gamma.delta",
                          "Dendritic.cells.activated",
                          "Macrophages.M1",
                          "NK.cells.activated",
                          "Plasma.cells",
                          "T.cells.CD4.memory.resting",
                          "T.cells.CD4.memory.activated",
                          "Mast.cells.activated",
                          "NK.cells.resting",
                          "Macrophages.M0",
                          "Macrophages.M2",
                          "Eosinophils",
                          "Monocytes",
                          "Dendritic.cells.resting",
                          "Mast.cells.resting",
                          "Neutrophils",
                          "Endothelial cells",
                          "Fibroblasts")

# tcga
indata <- log2(tpm + 1)
write.table(indata,file = file.path(res.path,"TCGA_log2TPM_hugo.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
filterCommonGenes(input.f=file.path(res.path, "TCGA_log2TPM_hugo.txt") , output.f=file.path(res.path,"TCGA_log2TPM_hugo_ESTIMATE.txt"), id="GeneSymbol")
estimateScore(file.path(res.path,"TCGA_log2TPM_hugo_ESTIMATE.txt"), file.path(res.path,"TCGA_log2TPM_hugo_estimate_score.txt"), platform="affymetrix")
est.tcga <- read.table(file = file.path(res.path,"TCGA_log2TPM_hugo_estimate_score.txt"),header = T,row.names = NULL,check.names = F,stringsAsFactors = F,sep = "\t")
rownames(est.tcga) <- est.tcga[,2]; colnames(est.tcga) <- est.tcga[1,]; est.tcga <- est.tcga[-1,c(-1,-2)];
est.tcga <- sapply(est.tcga, as.numeric); rownames(est.tcga) <- c("StromalScore","ImmuneScore","ESTIMATEScore","TumorPurity"); est.tcga.backup = as.data.frame(est.tcga); colnames(est.tcga.backup) <- colnames(indata)
est.tcga <- annTrackScale(indata = est.tcga, halfwidth = 2, poolsd = F); est.tcga <- as.data.frame(t(est.tcga))
rownames(est.tcga) <- colnames(tpm)

tcga.immune.gsva <- gsva(as.matrix(log2(tpm + 1)),
                         immune.sig.ccr,
                         method = "gsva")
annCol.tcga <- annCol
annCol.tcga$CMOIC <- paste0("CS",cmoic.BLCA$clust.res$clust)
annCol.tcga$ImmuneScore <- as.numeric(est.tcga[rownames(annCol.tcga),"ImmuneScore"])
annCol.tcga$StromalScore <- as.numeric(est.tcga[rownames(annCol.tcga),"StromalScore"])
annCol.tcga <- annCol.tcga[order(annCol.tcga$CMOIC),]
annColors.tcga <- annColors
annColors.tcga[["CMOIC"]] <- c("CS1" = clust.col[1],
                               "CS2" = clust.col[2],
                               "CS3" = clust.col[3],
                               "CS4" = clust.col[4])
annColors.tcga[["ImmuneScore"]] <- inferno(64)
annColors.tcga[["StromalScore"]] <- viridis(64)
indata <- log2(tpm[intersect(rownames(tpm),imm.targets),] + 1)
hm1 <- pheatmap(standarize.fun(indata[imm.targets,rownames(annCol.tcga)],halfwidth = 2),
                border_color = NA,
                annotation_col = annCol.tcga[,c("CMOIC","StromalScore","ImmuneScore")],
                annotation_colors = annColors.tcga[c("CMOIC","StromalScore","ImmuneScore")],
                color = colorpanel(64,low=blue,mid = "black",high=gold),
                show_rownames = T,
                show_colnames = F,
                cellheight = 12,
                cellwidth = 0.6,
                name = "Immunotherapy targets",
                cluster_rows = F,
                cluster_cols = F)
hm1

hm2 <- pheatmap(standarize.fun(tcga.immune.gsva[immune.sig.ccr.order,rownames(annCol.tcga)],halfwidth = 1),
                border_color = NA,
                color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
                gaps_row = c(14,22),
                show_rownames = T,
                show_colnames = F,
                cellheight = 12,
                cellwidth = 0.6,
                name = "Tumor microenvironment",
                cluster_rows = F,
                cluster_cols = F)
hm2

annCol.tcga$MeTIL <- MeTIL
indata <- t(annCol.tcga[,"MeTIL",drop = F])
hm3 <- pheatmap(standarize.fun(indata,halfwidth = 1),
                border_color = NA,
                color = NMF:::ccRamp(c(cyan,"black","#F12AFE"),64),
                show_rownames = T,
                show_colnames = F,
                cellheight = 12,
                cellwidth = 0.6,
                name = "MeTIL markers",
                cluster_rows = F,
                cluster_cols = F)
hm3
pdf(file.path(fig.path,"immune checkpoint expression and TME enrichment in TCGA cohort.pdf"),width = 10,height = 10)
draw(hm1 %v% hm2 %v% hm3, heatmap_legend_side = "right", annotation_legend_side = "right")
invisible(dev.off())

# illumina
write.table(combined.expr.illumina.combat,file = file.path(res.path,"ILLUMINA_hugo.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
filterCommonGenes(input.f=file.path(res.path, "ILLUMINA_hugo.txt") , output.f=file.path(res.path,"ILLUMINA_hugo_ESTIMATE.txt"), id="GeneSymbol")
estimateScore(file.path(res.path,"ILLUMINA_hugo_ESTIMATE.txt"), file.path(res.path,"ILLUMINA_hugo_estimate_score.txt"), platform="affymetrix")
est.illumina <- read.table(file = file.path(res.path,"ILLUMINA_hugo_estimate_score.txt"),header = T,row.names = NULL,check.names = F,stringsAsFactors = F,sep = "\t")
rownames(est.illumina) <- est.illumina[,2]; colnames(est.illumina) <- est.illumina[1,]; est.illumina <- est.illumina[-1,c(-1,-2)];
est.illumina <- sapply(est.illumina, as.numeric); rownames(est.illumina) <- c("StromalScore","ImmuneScore","ESTIMATEScore","TumorPurity"); est.illumina.backup = as.data.frame(est.illumina); colnames(est.illumina.backup) <- colnames(combined.expr.illumina.combat)
est.illumina <- annTrackScale(indata = est.illumina, halfwidth = 2, poolsd = F); est.illumina <- as.data.frame(t(est.illumina))
rownames(est.illumina) <- colnames(combined.expr.illumina.combat)

cms.illumina <- getConsensusClass(x = combined.expr.illumina.combat, gene_id = "hgnc_symbol")
annCol.illumina <- data.frame(Consensus = cms.illumina$consensusClass,
                              CMOIC = paste0("CS",illumina.ntp.pred$clust.res$clust),
                              row.names = colnames(combined.expr.illumina.combat))
annCol.illumina <- cbind.data.frame(annCol.illumina,combined.surv.illumina[rownames(annCol.illumina),c("Age","Gender","T_Stage","fustat","Cohort")])
annCol.illumina <- annCol.illumina[order(annCol.illumina$CMOIC),]
colnames(annCol.illumina)[6] <- "OS"
annCol.illumina$ImmuneScore <- as.numeric(est.illumina[rownames(annCol.illumina),"ImmuneScore"])
annCol.illumina$StromalScore <- as.numeric(est.illumina[rownames(annCol.illumina),"StromalScore"])
annColors.illumina <- annColors.tcga
annColors.illumina[["T_Stage"]] <- c("T2" = alpha(darkred,0.3),
                                     "T3" = alpha(darkred,0.7),
                                     "T4" = alpha(darkred,1))
annColors.illumina[["Cohort"]] <- c("GSE13507" = alpha(cohort.col[2],0.7),
                                    "GSE32548" = alpha(cohort.col[4],0.7),
                                    "GSE32894" = alpha(cohort.col[6],0.7),
                                    "GSE48075" = alpha(cohort.col[8],0.7),
                                    "GSE48276" = alpha(cohort.col[10],0.7))
illumina.immune.gsva <- gsva(as.matrix(combined.expr.illumina.combat),
                             immune.sig.ccr,
                             method = "gsva")

indata <- combined.expr.illumina.combat[intersect(rownames(combined.expr.illumina.combat),imm.targets),]
hm1 <- pheatmap(standarize.fun(indata[imm.targets,rownames(annCol.illumina)],halfwidth = 2),
                border_color = NA,
                annotation_col = annCol.illumina[,c("CMOIC","StromalScore","ImmuneScore")],
                annotation_colors = annColors.illumina[c("CMOIC","StromalScore","ImmuneScore")],
                color = colorpanel(64,low=blue,mid = "black",high=gold),
                show_rownames = T,
                show_colnames = F,
                cellheight = 12,
                cellwidth = 1.2,
                name = "Immunotherapy targets",
                cluster_rows = F,
                cluster_cols = F)
hm1

hm2 <- pheatmap(standarize.fun(illumina.immune.gsva[immune.sig.ccr.order,rownames(annCol.illumina)],halfwidth = 1),
                border_color = NA,
                color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
                gaps_row = c(14,22),
                show_rownames = T,
                show_colnames = F,
                cellheight = 12,
                cellwidth = 1.2,
                name = "Tumor microenvironment",
                cluster_rows = F,
                cluster_cols = F)
hm2

pdf(file.path(fig.path,"immune checkpoint expression and TME enrichment in illumina cohort.pdf"),width = 10,height = 10)
draw(hm1 %v% hm2, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
invisible(dev.off())

# affy
write.table(combined.expr.affy.combat,file = file.path(res.path,"AFFY_hugo.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
filterCommonGenes(input.f=file.path(res.path, "AFFY_hugo.txt") , output.f=file.path(res.path,"AFFY_hugo_ESTIMATE.txt"), id="GeneSymbol")
estimateScore(file.path(res.path,"AFFY_hugo_ESTIMATE.txt"), file.path(res.path,"AFFY_hugo_estimate_score.txt"), platform="affymetrix")
est.affy <- read.delim(file = file.path(res.path,"AFFY_hugo_estimate_score.txt"),header = T,row.names = NULL,check.names = F,stringsAsFactors = F,sep = "\t")
rownames(est.affy) <- est.affy[,2]; colnames(est.affy) <- est.affy[2,]; est.affy <- est.affy[c(-1,-2),c(-1,-2)];
est.affy <- est.affy[1:2,]
est.affy <- sapply(est.affy, as.numeric); rownames(est.affy) <- c("StromalScore","ImmuneScore"); est.affy.backup = as.data.frame(est.affy); colnames(est.affy.backup) <- colnames(combined.expr.affy.combat)
est.affy <- annTrackScale(indata = est.affy, halfwidth = 2, poolsd = F); est.affy <- as.data.frame(t(est.affy))
rownames(est.affy) <- colnames(combined.expr.affy.combat)

cms.affy <- getConsensusClass(x = combined.expr.affy.combat, gene_id = "hgnc_symbol")
annCol.affy <- data.frame(Consensus = cms.affy$consensusClass,
                          CMOIC = paste0("CS",affy.ntp.pred$clust.res$clust),
                          row.names = colnames(combined.expr.affy.combat))
annCol.affy <- cbind.data.frame(annCol.affy,combined.surv.affy[rownames(annCol.affy),c("Age","Gender","T_Stage","fustat","Cohort")])
annCol.affy <- annCol.affy[order(annCol.affy$CMOIC),]
colnames(annCol.affy)[6] <- "OS"
annCol.affy$ImmuneScore <- as.numeric(est.affy[rownames(annCol.affy),"ImmuneScore"])
annCol.affy$StromalScore <- as.numeric(est.affy[rownames(annCol.affy),"StromalScore"])
annColors.affy <- annColors.tcga
annColors.affy[["T_Stage"]] <- c("T2" = alpha(darkred,0.3),
                                 "T3" = alpha(darkred,0.7),
                                 "T4" = alpha(darkred,1))
annColors.affy[["Cohort"]] <- c("E-MTAB-1803" = "black",
                                "GSE31684" = "grey90")
affy.immune.gsva <- gsva(as.matrix(combined.expr.affy.combat),
                         immune.sig.ccr,
                         method = "gsva")

indata <- combined.expr.affy.combat[intersect(rownames(combined.expr.affy.combat),imm.targets),]
hm1 <- pheatmap(standarize.fun(indata[imm.targets,rownames(annCol.affy)],halfwidth = 2),
                border_color = NA,
                annotation_col = annCol.affy[,c("CMOIC","StromalScore","ImmuneScore")],
                annotation_colors = annColors.affy[c("CMOIC","StromalScore","ImmuneScore")],
                color = colorpanel(64,low=blue,mid = "black",high=gold),
                show_rownames = T,
                show_colnames = F,
                cellheight = 12,
                cellwidth = 2.5,
                name = "Immunotherapy targets",
                cluster_rows = F,
                cluster_cols = F)
hm1

hm2 <- pheatmap(standarize.fun(affy.immune.gsva[immune.sig.ccr.order,rownames(annCol.affy)],halfwidth = 1),
                border_color = NA,
                color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
                gaps_row = c(14,22),
                show_rownames = T,
                show_colnames = F,
                cellheight = 12,
                cellwidth = 2.5,
                name = "Tumor microenvironment",
                cluster_rows = F,
                cluster_cols = F)
hm2

pdf(file.path(fig.path,"immune checkpoint expression and TME enrichment in affy cohort.pdf"),width = 10,height = 10)
draw(hm1 %v% hm2, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
invisible(dev.off())

# imvor
indata <- log2(imvigor210.expr + 1)
write.table(indata,file = file.path(res.path,"imvigor210_log2TPM_hugo.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
filterCommonGenes(input.f=file.path(res.path, "imvigor210_log2TPM_hugo.txt") , output.f=file.path(res.path,"imvigor210_log2TPM_hugo_ESTIMATE.txt"), id="GeneSymbol")
estimateScore(file.path(res.path,"imvigor210_log2TPM_hugo_ESTIMATE.txt"), file.path(res.path,"imvigor210_log2TPM_hugo_estimate_score.txt"), platform="affymetrix")
est.imvigor210 <- read.table(file = file.path(res.path,"imvigor210_log2TPM_hugo_estimate_score.txt"),header = T,row.names = NULL,check.names = F,stringsAsFactors = F,sep = "\t")
rownames(est.imvigor210) <- est.imvigor210[,2]; colnames(est.imvigor210) <- est.imvigor210[1,]; est.imvigor210 <- est.imvigor210[-1,c(-1,-2)];
est.imvigor210 <- sapply(est.imvigor210, as.numeric); rownames(est.imvigor210) <- c("StromalScore","ImmuneScore","ESTIMATEScore","TumorPurity"); est.imvigor210.backup = as.data.frame(est.imvigor210); colnames(est.imvigor210.backup) <- colnames(indata)
est.imvigor210 <- annTrackScale(indata = est.imvigor210, halfwidth = 2, poolsd = F); est.imvigor210 <- as.data.frame(t(est.imvigor210))
rownames(est.imvigor210) <- colnames(imvigor210.expr)

cms.imvigor210 <- getConsensusClass(x = log2(imvigor210.expr + 1), gene_id = "hgnc_symbol")
annCol.imvigor210 <- data.frame(Response = imvigor210.surv$`Best Confirmed Overall Response`,
                                Response2 = imvigor210.surv$binaryResponse,
                                OS = imvigor210.surv$censOS,
                                IC = imvigor210.surv$`IC Level`,
                                TC = imvigor210.surv$`TC Level`,
                                Immunephenotype = imvigor210.surv$`Immune phenotype`,
                                Consensus = cms.imvigor210$consensusClass,
                                Lund = imvigor210.surv$Lund2,
                                CMOIC = paste0("CS",imvigor210.ntp.pred$clust.res$clust),
                                row.names = colnames(imvigor210.expr))
annCol.imvigor210 <- annCol.imvigor210[order(annCol.imvigor210$CMOIC),]
annCol.imvigor210[is.na(annCol.imvigor210)] <- "Missing"
annCol.imvigor210$ImmuneScore <- as.numeric(est.imvigor210[rownames(annCol.imvigor210),"ImmuneScore"])
annCol.imvigor210$StromalScore <- as.numeric(est.imvigor210[rownames(annCol.imvigor210),"StromalScore"])
annCol.imvigor210$Response <- factor(annCol.imvigor210$Response,levels = c("CR","PR","SD","PD"))

annColors.imvigor210 <- annColors.tcga
annColors.imvigor210[["CMOIC"]] <-  c("CS1" = clust.col[1],"CS2" = clust.col[2],"CS3" = clust.col[3],"CS4" = clust.col[4])
annColors.imvigor210[["Lund"]] <-  c("Genomically unstable" = mycol[2],"Infiltrated" = mycol[3],"Basal/SCC-like" = mycol[1],"UroA" = mycol[4],"UroB" = mycol[5])
annColors.imvigor210[["TCGA"]] <-  c("Basal_squamous" = mycol[1], "Luminal" = mycol[3],"Luminal_infiltrated" = mycol[2],"Luminal_papillary" = mycol[4],"Neuronal" = mycol[5])
annColors.imvigor210[["Consensus"]] <- c("Ba/Sq" = mycol[1], "LumNS" = mycol[3], "LumP" = mycol[4], "LumU" = mycol[5], "NE-like" = lightred, "Stroma-rich" = mycol[6])
annColors.imvigor210[["OS"]] <- c("0" = "grey90","1" = "black")
annColors.imvigor210[["IC"]] <- c("IC0" = clust.col[3],"IC1" = clust.col[4], "IC2+" = clust.col[1], Missing = "grey90")
annColors.imvigor210[["TC"]] <- c("TC0" = clust.col[3],"TC1" = clust.col[4], "TC2+" = clust.col[1], Missing = "grey90")
annColors.imvigor210[["Immunephenotype"]] <- c("desert" = clust.col[3],"excluded" = clust.col[4], "inflamed" = clust.col[1], Missing = "grey90")
annColors.imvigor210[["Response"]] <- c("CR" = sun,"PR" = lightred, "SD" = blue,"PD" = darkblue)
annColors.imvigor210[["Response2"]] <- c("CR/PR" = jco[1],"SD/PD" = jco[2], Missing = "grey90")

imvigor210.immune.gsva <- gsva(as.matrix(log2(imvigor210.expr + 1)),
                               immune.sig.ccr,
                               method = "gsva")
imvor.order <- rownames(annCol.imvigor210[order(annCol.imvigor210$CMOIC,annCol.imvigor210$Response),])
indata <- log2(imvigor210.expr + 1)[intersect(rownames(log2(imvigor210.expr + 1)),imm.targets),]
hm1 <- pheatmap(standarize.fun(indata[imm.targets,imvor.order],halfwidth = 2),
                border_color = NA,
                annotation_col = annCol.imvigor210[imvor.order,c("CMOIC","Response","Consensus","IC","TC","Immunephenotype","StromalScore","ImmuneScore")],
                annotation_colors = annColors.imvigor210[c("CMOIC","Response","Consensus","IC","TC","Immunephenotype","StromalScore","ImmuneScore")],
                color = colorpanel(64,low=blue,mid = "black",high=gold),
                show_rownames = T,
                show_colnames = F,
                cellheight = 12,
                cellwidth = 0.8,
                name = "Immunotherapy targets",
                cluster_rows = F,
                cluster_cols = F)
hm1

hm2 <- pheatmap(standarize.fun(imvigor210.immune.gsva[immune.sig.ccr.order,imvor.order],halfwidth = 1),
                border_color = NA,
                color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
                gaps_row = c(14,22),
                show_rownames = T,
                show_colnames = F,
                cellheight = 12,
                cellwidth = 0.8,
                name = "Tumor microenvironment",
                cluster_rows = F,
                cluster_cols = F)
hm2

pdf(file.path(fig.path,"immune checkpoint expression and TME enrichment in imvigor210 cohort.pdf"),width = 10,height = 10)
draw(hm1 %v% hm2, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
invisible(dev.off())
