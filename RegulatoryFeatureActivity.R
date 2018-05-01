library(GenomicRanges)
load("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/SE/RData/all_cells_SE_type_and_PSI.RData")
load("/mnt/data5/BGI/UCB/tangchao/DSU/cell_by_cell/SE/RData/all_cells_SE_GTF.RData")
motiffeatures <- read.table(gzfile("/mnt/data5/BGI/UCB/tangchao/data/RegulatoryFeatureActivity/homo_sapiens.GRCh38.motiffeatures.20161111.gff.gz"), stringsAsFactors=F)


SE_loci <- SE_GTF[SE_GTF$AStype == "exon-skipping_exactly" & SE_GTF$OverType == "equal", c("loci", "V1", "V4", "V5", "V7")]

gr <- GRanges(
  seqnames = SE_loci$V1,
  ranges = IRanges(start = SE_loci$V4, end = SE_loci$V5, names = SE_loci$loci),
  strand = SE_loci$V7
)
gr

motiffeatures <- motiffeatures[motiffeatures$V1 %in% c(1:22,"X","Y"),]
motiffeatures$loci <- paste(motiffeatures$V1,":",motiffeatures$V4,"-",motiffeatures$V5,sep = "")



gri <- GRanges(seqnames = motiffeatures$V1[i], 
               ranges = IRanges(start = motiffeatures$V4[i], end = motiffeatures$V5[i]), 
               strand = motiffeatures$V7[i], score = motiffeatures$V6[i])

countOverlaps(gri, gr)

library(parallel)
counts <- mcmapply(function(i){
  countOverlaps(GRanges(seqnames = motiffeatures$V1[i], 
                        ranges = IRanges(start = motiffeatures$V4[i], end = motiffeatures$V5[i]), 
                        strand = motiffeatures$V7[i], score = motiffeatures$V6[i]), gr)
}, 1:nrow(motiffeatures), mc.cores = 30)

sum(counts>0)
[1] 223
length(counts)
[1] 558299

which(counts==3)[1]
[1] 40054
i = which(counts==3)[1]
gri <- GRanges(seqnames = motiffeatures$V1[i], 
               ranges = IRanges(start = motiffeatures$V4[i], end = motiffeatures$V5[i]), 
               strand = motiffeatures$V7[i], score = motiffeatures$V6[i])

countOverlaps(gri, gr)
findOverlapPairs(gri, gr)

which(counts>0)

test <- list()
for(j in 1:sum(counts>0)){
	i <- which(counts>0)[j]
	gri <- GRanges(seqnames = motiffeatures$V1[i], 
               ranges = IRanges(start = motiffeatures$V4[i], end = motiffeatures$V5[i]), 
               strand = motiffeatures$V7[i], score = motiffeatures$V6[i])
	test[j] <- findOverlapPairs(gri, gr)
}

findOverlaps(gri, gr)
subjectHits(findOverlaps(gri, gr))


test2 <- list()
for(j in 1:sum(counts>0)){
	i <- which(counts>0)[j]
	gri <- GRanges(seqnames = motiffeatures$V1[i], 
               ranges = IRanges(start = motiffeatures$V4[i], end = motiffeatures$V5[i]), 
               strand = motiffeatures$V7[i], score = motiffeatures$V6[i])
	test2[[j]] <- subjectHits(findOverlaps(gri, gr))
}

unlist(test2)

rep(which(counts>0),counts[counts>0])

search_Hits <- data.frame(Hits_re = rep(which(counts>0),counts[counts>0]), Hits_se = unlist(test2))

motiffeatures[rep(which(counts>0),counts[counts>0]),]

SE_GTF[SE_GTF$AStype == "exon-skipping_exactly" & SE_GTF$OverType == "equal", ][unlist(test2),]

tang <- cbind(motiffeatures[rep(which(counts>0),counts[counts>0]),],SE_GTF[SE_GTF$AStype == "exon-skipping_exactly" & SE_GTF$OverType == "equal", ][unlist(test2),])

colnames(tang) <- c("seqname_TF", "source_TF", "feature_TF", "start_TF", "end_TF", "score_TF", "strand_TF", "frame_TF", "attribute_TF", "loci_TF",
"loci", "left", "right", "chr", "dist_left", "dist_right", "dist_insert", "ASrange", "AStype", "OverlapCount", "ASlen", "OverType",
"seqname_GTF", "source_GTF", "feature_GTF", "start_GTF", "end_GTF", "score_GTF", "strand_GTF", "frame_GTF", "attribute_GTF")


tang$binding_matrix <- do.call(rbind,strsplit(tang$attribute_TF, split="[=;]"))[,2]
tang$motif_feature_type <- do.call(rbind,strsplit(tang$attribute_TF, split="[=;]"))[,4]
tang <- tang[,c("seqname_TF", "source_TF", "feature_TF", "start_TF", "end_TF", "score_TF", "strand_TF", "frame_TF", "attribute_TF", "binding_matrix", "motif_feature_type", "loci_TF",
"loci", "left", "right", "chr", "dist_left", "dist_right", "dist_insert", "ASrange", "AStype", "OverlapCount", "ASlen", "OverType",
"seqname_GTF", "source_GTF", "feature_GTF", "start_GTF", "end_GTF", "score_GTF", "strand_GTF", "frame_GTF", "attribute_GTF")]

test <- as.data.frame(do.call(rbind, strsplit(tang$attribute_GTF, split = "[; ]")))
valueFind <- function(x, type = "gene_id"){
  for(i in 1:length(x)){
    if(sum(x == type) == 0){
      return(NA)
    }else{
      return(x[which(x == type)[1]+1])
    }
  }
}
gtf_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),
                        gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),             
                        transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
                        transcript_name = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_name")})),
                        gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})),
                        exon_number = as.vector(apply(test,1,function(x){valueFind(x,type = "exon_number")})), stringsAsFactors=F)

tang <- cbind(tang, gtf_infor)

save(tang, file = "/mnt/data5/BGI/UCB/tangchao/Regulatory/SE/SE_exon_Overlap_homo_sapiens.GRCh38.motiffeatures.20161111.RData")
write.table(tang, file = "/mnt/data5/BGI/UCB/tangchao/Regulatory/SE_exon_Overlap_homo_sapiens.GRCh38.motiffeatures.20161111.txt", row.names=F, sep = "\t", quote = F)









