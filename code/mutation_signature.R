#------------------------------#
# Mutation signature of APOBEC #
library(maftools)
library(deconstructSigs)
maf2 <- read.maf(maf[which(maf$Tumor_Sample_Barcode %in% rownames(annCol.tcga)),])
maf.tnm = trinucleotideMatrix(maf = maf2, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
APOBEC <- as.data.frame(maf.tnm$APOBEC_scores)
APOBEC$CMOIC <- annCol.tcga[APOBEC$Tumor_Sample_Barcode,"CMOIC"] 

tmp <- maf[which(maf$Tumor_Sample_Barcode %in% rownames(annCol.tcga)),]
tmp$Chromosome <- paste0("chr",tmp$Chromosome)
unique(tmp$Variant_Classification)
# [1] "Nonsense_Mutation" "Missense_Mutation" "Splice_Site"       "Frame_Shift_Del"  
# [5] "In_Frame_Ins"      "Frame_Shift_Ins"   "In_Frame_Del" 
sigs.input.sbs <- mut.to.sigs.input(mut.ref = as.data.frame(tmp), 
                                    sample.id = "Tumor_Sample_Barcode", 
                                    chr = "Chromosome", 
                                    pos = "Start_Position", 
                                    ref = "Reference_Allele", 
                                    alt = "Tumor_Seq_Allele2",
                                    bsg = BSgenome.Hsapiens.UCSC.hg19,
                                    sig.type = "SBS")
write.csv(sigs.input.sbs,file.path(res.path,"mutation.sig.input.bydeconstructSigs.csv"),row.names = T,quote = F)

cut.off <- 0
mut.wt.cosmic2013 <- mut.wt.cosmic2019 <- data.frame()
sigs.out.cosmic2013.list <- sigs.out.cosmic2019.list <- list()
for (sample in rownames(sigs.input.sbs)) {
  tmp <- whichSignatures(tumor.ref = sigs.input.sbs, 
                         signatures.ref = signatures.cosmic, 
                         sample.id = sample, 
                         contexts.needed = TRUE,
                         tri.counts.method = 'exome2genome',
                         signature.cutoff = cut.off)
  
  sigs.out.cosmic2013.list[[sample]] <- tmp
  tmp <- data.frame(c(tmp$weights,unknown=tmp$unknown),row.names = sample)
  mut.wt.cosmic2013 <- rbind.data.frame(mut.wt.cosmic2013,tmp)
  
  tmp <- whichSignatures(tumor.ref = sigs.input.sbs, 
                         signatures.ref = signatures.exome.cosmic.v3.may2019, 
                         sample.id = sample, 
                         contexts.needed = TRUE,
                         tri.counts.method = 'exome2genome',
                         signature.cutoff = cut.off)
  
  sigs.out.cosmic2019.list[[sample]] <- tmp
  tmp <- data.frame(c(tmp$weights,unknown=tmp$unknown),row.names = sample)
  mut.wt.cosmic2019 <- rbind.data.frame(mut.wt.cosmic2019,tmp)
}
write.csv(mut.wt.cosmic2013,file.path(res.path,"mutation.snp.signature.weightMatrix.cosmic2013.csv"),row.names = T,quote = F)
write.csv(mut.wt.cosmic2019,file.path(res.path,"mutation.snp.signature.weightMatrix.cosmic2019.csv"),row.names = T,quote = F)

mut.wt.cosmic2013.trim <- mut.wt.cosmic2013[,1:30]
mut.wt.cosmic2013.trim <- mut.wt.cosmic2013.trim[,colSums(mut.wt.cosmic2013.trim) > 0]
mut.wt.cosmic2013.trim.backup <- mut.wt.cosmic2013.trim

pheatmap(t(mut.wt.cosmic2013.trim[rownames(annCol.tcga),]),
         border_color = NA,
         annotation_col = annCol.tcga,
         annotation_colors = annColors.tcga,
         color = NMF:::ccRamp(x = c("#EAF0FA","#6081C3","#3454A7"),n = 64),
         cluster_cols = F,
         cluster_rows = T,
         show_rownames = T,
         show_colnames = F,
         cellwidth = 0.8,
         cellheight = 10)
rm(mut.wt.cosmic2013.trim)
rm(mut.wt.cosmic2013.trim.backup)
rm(maf2)
rm(APOBEC)
rm(maf.tnm)
gc()

mut.wt.cosmic2013 <- read.csv(file.path(res.path,"mutation.snp.signature.weightMatrix.cosmic2013.csv"),row.names = 1,check.names = F,stringsAsFactors = F,header = T)
mut.wt.cosmic2013 <- mut.wt.cosmic2013[,c("Signature.1","Signature.2","Signature.5","Signature.13")]
mut.wt.cosmic2013$APOBEC <- mut.wt.cosmic2013$Signature.2 + mut.wt.cosmic2013$Signature.13
mut.wt.cosmic2013$CMOIC <- annCol.tcga[rownames(mut.wt.cosmic2013),"CMOIC"]
mut.wt.cosmic2013 <- mut.wt.cosmic2013[order(mut.wt.cosmic2013$CMOIC,-mut.wt.cosmic2013$APOBEC,decreasing = F),]
