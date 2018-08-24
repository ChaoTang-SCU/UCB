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






#### sashimi plot of alternative exons, which has TF binding sites.
# /mnt/data5/BGI/UCB/tangchao/Regulatory/SE/Figure
length(unique(tang$ASrange))
[1] 135

#cells <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, sep = "\t", stringsAsFactors = F)

#### sashimi2
sashimi <- function(junction, cell, left_expand = 1000, right_expand = 1000, tsv = "/mnt/data5/BGI/UCB/tangchao/data/BAM_tsv/", outDir = "~/", min = 10, ann_height = 4){
  if(dir.exists(outDir) == FALSE) dir.create(outDir)
  Coordinates <- data.frame(do.call(rbind,strsplit(junction, split = "[:-]")),stringsAsFactors=F)
  Coordinates$X2 <- abs(as.numeric(Coordinates$X2) - left_expand)
  Coordinates$X3 <- abs(as.numeric(Coordinates$X3) + right_expand)
  window <- paste(Coordinates$X1, ":", Coordinates$X2, "-", Coordinates$X3, sep = "") 
  if(length(cell) == 1){
    Dir_out <- paste(outDir, cell, "_", junction, sep = "")
    work <- paste("system('bash /mnt/data5/BGI/UCB/tangchao/novel_SS/sashimi/sashimi.sh ", window, cell, Dir_out, min, ann_height, "')")
    eval(parse(text = work))
  }else{
    Dir_out <- paste(outDir, junction, sep = "")
    tmpdir <- tempdir()
    work1 <- paste("system('", "cat ", paste(tsv,cell,".tsv",sep = "", collapse = " "), " > ", tmpdir,"/tmp.tsv", "')", sep = "")
    work2 <- paste("system('bash /mnt/data5/BGI/UCB/tangchao/novel_SS/sashimi/sashimi2.sh ", window, paste(tmpdir,"/tmp.tsv", sep = ""), Dir_out, min, ann_height, "')")
    work3 <- paste("system('", "rm",paste(tmpdir,"/tmp.tsv", "')", sep = ""))
    eval(parse(text = work1))
    eval(parse(text = work2))
    eval(parse(text = work3))
  }
}



test <- all_cells_SE_type_and_PSI[all_cells_SE_type_and_PSI$loci %in% tang$loci & !duplicated(all_cells_SE_type_and_PSI$ASrange), ]

apply(test[,12:ncol(test)],1,function(x) {!is.na(x) & names(x) %in% cells$V1}) -> test2

cell <- vector()
for(i in 1:nrow(test)){
  tmp <- test[i, 12:ncol(test)]
  if(sum(!is.na(tmp))==1){
    cell[i] <- names(tmp[which(!is.na(tmp))])
  }else{
    tmp <- tmp[, !is.na(tmp)]
    cell[i] <- names(which.min(abs(tmp-.5)))
  }
}


for(i in 1:length(cell)){
  sashimi(junction = substr(test$ASrange,4,100)[i], cell = cell[i], 
          outDir = "/mnt/data5/BGI/UCB/tangchao/Regulatory/SE/Figure/", 
          ann_height = 10, min = 2,
          left_expand = round(test$dist_left*1.2)[i],
          right_expand = round(test$dist_right*1.2)[i])
}


cell <- list()
for(i in 1:nrow(test)){
  tmp <- test[i, 12:ncol(test)]
  if(sum(!is.na(tmp))==1){
    cell[[i]] <- names(tmp[which(!is.na(tmp))])
  }else{
    tmp <- tmp[, !is.na(tmp)]
    cell[[i]] <- c(names(which.min(abs(tmp-.5))),names(which.max(abs(tmp-.5))))
  }
}

for(i in 1:length(cell)){
  sashimi(junction = substr(test$ASrange,4,100)[i], cell = cell[[i]], 
          outDir = "/mnt/data5/BGI/UCB/tangchao/Regulatory/SE/Figure2/", 
          ann_height = 10, min = 1,
          left_expand = round(test$dist_left*1.2)[i],
          right_expand = round(test$dist_right*1.2)[i])
}


for(i in 127){
  sashimi(junction = substr(test$ASrange,4,100)[i], cell = cell[[i]], 
          outDir = "/mnt/data5/BGI/UCB/tangchao/Regulatory/SE/Figure2/", 
          ann_height = 10, min = 1,
          left_expand = round(test$dist_left*4)[i],
          right_expand = round(test$dist_right*1.2)[i])
}

grep("103526066",test$ASrange)
[1] 26

for(i in 26){
  sashimi(junction = substr(test$ASrange,4,100)[i], cell = cell[[i]], 
          outDir = "/mnt/data5/BGI/UCB/tangchao/Regulatory/SE/Figure2/", 
          ann_height = 10, min = 1,
          left_expand = round(test$dist_left)[i],
          right_expand = round(test$dist_right*4)[i])
}


grep("187041916",test$ASrange)
[1] 130

for(i in 130){
  sashimi(junction = substr(test$ASrange,4,100)[i], cell = cell[[i]], 
          outDir = "/mnt/data5/BGI/UCB/tangchao/Regulatory/SE/Figure2/", 
          ann_height = 10, min = 1,
          left_expand = round(test$dist_left*8)[i],
          right_expand = round(test$dist_right*2)[i])
}

grep("48094076", test$ASrange)
[1] 34
for(i in 34){
  sashimi(junction = substr(test$ASrange,4,100)[i], cell = cell[[i]], 
          outDir = "/mnt/data5/BGI/UCB/tangchao/Regulatory/SE/Figure2/", 
          ann_height = 10, min = 1,
          left_expand = round(test$dist_left*5)[i],
          right_expand = round(test$dist_right*1.2)[i])
}

for(i in 38){
  sashimi(junction = substr(test$ASrange,4,100)[i], cell = cell[[i]], 
          outDir = "/mnt/data5/BGI/UCB/tangchao/Regulatory/SE/Figure2/", 
          ann_height = 10, min = 1,
          left_expand = round(test$dist_left)[i],
          right_expand = round(test$dist_right*10)[i])
}

loci_tang <- paste(do.call(rbind,strsplit(test$loci, "[:-]"))[,1], ":", do.call(rbind,strsplit(test$loci, "[:-]"))[,2], "-", do.call(rbind,strsplit(test$loci, "[:-]"))[,5], sep = "")


for(i in 84){
  sashimi(junction = "loci_tang[i]", cell = cell[[i]], 
          outDir = "/mnt/data5/BGI/UCB/tangchao/Regulatory/SE/Figure2/", 
          ann_height = 10, min = 1,
          left_expand = round(test$dist_left*1)[i],
          right_expand = round(test$dist_right*1)[i])
}





