#### Single cell DSU new vertion

#### 0. Basic settings ==========================================================================================================================

require(parallel)
library(dplyr)
library(plyr)
library(pryr)
require(parallel)
library(ggplot2)
library(cowplot)
library(GenomicRanges)
library(reshape2)
library(vioplot)

depar <- par()

## path of gtf
pfgtf <- file.path("/mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.87.gtf")

### load SJs merge file =========================================================================================================================

load("/mnt/data5/BGI/UCB/tangchao/data/SJ/SJ_merged_raw_te(all_more_than_10).RData")

Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

te_table <- data.frame(as.data.frame(te[,-1]), row.names = te[[1]])
rm("te");gc()

colnames(te_table) <- do.call(rbind,strsplit(colnames(te_table),"_SJ"))[,1]
dim(te_table)
# [1] 461646    3574
te_table <- te_table[, row.names(Cell_type)]

require(data.table)  # v1.6.6
require(gdata)
f_dowle3 = function(DT) {
  # either of the following for loops
  
  # by name :
  for (j in names(DT))
    set(DT,which(is.na(DT[[j]])),j,0)
  
  # or by number (slightly faster than by name) :
  for (j in seq_len(ncol(DT)))
    set(DT,which(is.na(DT[[j]])),j,0)
}

f_dowle3(te_table)

## We only care about chr1-22 and X Y
chr_tu <- do.call(rbind, strsplit(rownames(te_table), split = ":"))[,1] %in% c(1:22,"X","Y")
sum(chr_tu)
# [1] 460327

SJ_tu <- te_table[chr_tu,]
dim(SJ_tu)
# [1] 460327   3039
rm("te_table");gc()



#### 2. Identify alternative splicing events ====================================================================================================


#junc.names=do.call(rbind,(strsplit(sub(x=rownames(junc),
#                                       pattern="^([[:alnum:]]+):([[:digit:]]+)-([[:digit:]]+)$",
#                                       replace="\\1;\\2;\\3"),split=";")))
junc.names=do.call(rbind,strsplit(rownames(SJ_tu),split="[:-]"))
#alnum-- Letters and Numbers；digit -- Numbers
colnames(junc.names) <- c("chr","start","end")
rownames(junc.names) <- rownames(SJ_tu)
junc.names=data.frame(junc.names,stringsAsFactors=F)
junc.names$start=as.integer(junc.names$start)
junc.names$end=as.integer(junc.names$end)

junc.names=junc.names[order(junc.names$chr,junc.names$start,junc.names$end),]
junc.names$names=rownames(junc.names)


junc.as=do.call(c,mclapply(unique(junc.names$chr),function(chr) {
  junc.chr=junc.names[junc.names$chr==chr,]
  same.start=dlply(junc.chr,c("start"),function(x) x)
  same.end=dlply(junc.chr,c("end"),function(x) {if(nrow(x)>1) {return(x)} else {return(NULL)}})
  print(chr)
  return(c(same.start,same.end[!sapply(same.end,is.null)]))
},mc.cores=10))


juncs <- junc.as

junc.as <- junc.as[sapply(junc.as,nrow)>1]
names(junc.as) <- sapply(junc.as,function(x) paste(x$names,collapse="_"))

## junction just has two alternative starts/ends:
junc.as.2=junc.as[sapply(junc.as,nrow)==2]

table(sapply(junc.as.2,nrow))



#### 3.calculate PSI of all junction 2 ==========================================================================================================



#psi_sj2=mclapply(junc.as.2,function(sjs){
#  sjs.names=sjs$names[order(as.numeric(sjs$end)-as.numeric(sjs$start))]# find the splicing_in isoform to calculate the psi
#  tab=as.matrix(t(SJ_tu[sjs.names,]))
#  distc=max(c(diff(as.numeric(sjs$start)),diff(as.numeric(sjs$end))))
#  rs<-apply(tab, 1, function(a) sum(a))
#  tsrs<-tab/rs
#  return(tsrs[,1])
#},mc.cores=10)#

#psi.co<-t(do.call(cbind,psi_sj2))



#### 4.find SE ==================================================================================================================================



junc.as2.df=do.call(rbind,junc.as.2)
dim(junc.as2.df)

junc.as2.df[order(junc.as2.df$chr,junc.as2.df$start,junc.as2.df$end),][1:20,]

junc.as.2.ot <- junc.as2.df[order(junc.as2.df$chr,junc.as2.df$start,junc.as2.df$end),]
junc.as.2.ot[1:20,]
dim(junc.as.2.ot)


temp <- dlply(data.frame(chr = unique(junc.as.2.ot$chr)),"chr")


for(i in 1:length(unique(junc.as.2.ot$chr))){
  ## Calculated separately by chr
  junc.chr=junc.as.2.ot[junc.as.2.ot$chr==unique(junc.as.2.ot$chr)[i],]
  if(length(names(table(junc.chr$start)[table(junc.chr$start) == 3])) == 0 ){
    print(paste("The chr",unique(junc.as.2.ot$chr)[i],"have no SE", sep = " "))
  }else{
    ## left:
    ## because the longest SJ has been counted twice in SE, so the same starts/ends will appear three times
    se_start_left <- names(table(junc.chr$start)[table(junc.chr$start) == 3])
    if(length(se_start_left)==0) next()
    se_left_all <- junc.chr[junc.chr$start %in% se_start_left, ]
    
    ## because we have ordered the data table by chr, start and end, 
    ## so the order of same start and end in SE will appear in a certain order: 1,2,1,2
    se_start_left_yn <- as.vector(rep(NA,length(unique(se_left_all$start))))
    for(j in 1:length(unique(se_left_all$start))){
      se_left_all_j <- se_left_all[se_left_all$start == unique(se_left_all$start)[j],]
      se_start_left_yn[j] <- do.call(rbind, strsplit(rownames(se_left_all_j), split = "\\."))[1,2] == "1" &
        do.call(rbind, strsplit(rownames(se_left_all_j), split = "\\."))[2,2] == "2" &
        do.call(rbind, strsplit(rownames(se_left_all_j), split = "\\."))[3,2] == "1"
    }
    if(sum(se_start_left_yn)==0) next()
    se_left_all <- se_left_all[rep(se_start_left_yn,rep(3,length(se_start_left_yn))),]
    
    se_left <- se_left_all[seq(2,nrow(se_left_all),3),]
    se_left <- se_left[do.call(rbind,strsplit(rownames(se_left), split = "\\."))[,2] == 2, ]
    rownames(se_left) <- do.call(rbind,strsplit(rownames(se_left), split = "\\."))[,1]
    ## right:
    se_end_right <- names(table(junc.chr$end)[table(junc.chr$end) == 3])
    if(length(se_end_right)==0) next()
    se_right_all <- junc.chr[junc.chr$end %in% se_end_right, ]
    
    se_end_right_yn <- as.vector(rep(NA,length(unique(se_right_all$end))))
    for(j in 1:length(unique(se_right_all$end))){
      se_right_all_j <- se_right_all[se_right_all$end == unique(se_right_all$end)[j],]
      se_end_right_yn[j] <- do.call(rbind, strsplit(rownames(se_right_all_j), split = "\\."))[1,2] == "2" &
        do.call(rbind, strsplit(rownames(se_right_all_j), split = "\\."))[2,2] == "1" &
        do.call(rbind, strsplit(rownames(se_right_all_j), split = "\\."))[3,2] == "2"
    }
    se_right_all <- se_right_all[rep(se_end_right_yn,rep(3,length(se_end_right_yn))),]
    
    se_right <- se_right_all[seq(2,nrow(se_right_all),3),]
    
    se_right <- se_right[do.call(rbind,strsplit(rownames(se_right), split = "\\."))[,2] == 1, ]
    
    rownames(se_right) <- do.call(rbind,strsplit(rownames(se_right), split = "\\."))[,1]
    ## merge:
    se_left_tu <- se_left[se_left$names %in% se_right$names,]
    se_right_tu <- se_right[se_right$names %in% se_left$names,]
    
    se_raw_name <- paste(row.names(se_left_tu),row.names(se_right_tu),sep = "_")
    
    se_names <- do.call(rbind, strsplit(se_raw_name,split = "[-:_\\.]"))[,c(1,2,3,11,12)]
    if(length(se_names) == 5){
      se_name <- paste(se_names[1], paste(se_names[2],se_names[3],se_names[4],se_names[5], sep = "-"), sep = ":")
    }else{
      se_name <- paste(se_names[,1], paste(se_names[,2],se_names[,3],se_names[,4],se_names[,5], sep = "-"), sep = ":")
    }
    se <- data.frame(loci = se_name, left = rownames(se_left_tu), right = rownames(se_right_tu), row.names = se_left_tu$names)
    temp[[i]] <- se
  }
}

temp <- temp[lapply(temp,ncol)>1]

SE <- do.call(rbind, temp)
dim(SE)
# [1] 7142    3

rownames(SE) <- do.call(rbind, strsplit(rownames(SE), split = "\\."))[,2]

loci <- do.call(rbind, strsplit(as.vector(SE$loci), split = "[-:]"))

SE$chr <- loci[,1]

SE$dist_left <- as.numeric(loci[,5]) - as.numeric(loci[,3])

SE$dist_right <- as.numeric(loci[,4]) - as.numeric(loci[,2])

SE$dist_insert <- as.numeric(loci[,4]) - as.numeric(loci[,3])

SE <- SE[SE$dist_insert > 0,]
dim(SE)
# [1] 6848   7



#### 5.SE type ==================================================================================================================================



gtf <- read.table(pfgtf, head = F, sep = "\t", stringsAsFactors = F)
exon <- subset(gtf,gtf$V3=="exon")
exon <- exon[!duplicated(exon[,c(1,4,5)]),]
exon <- subset(exon,exon$V1 %in% c(1:22,"X","Y"))
exon$gene_id <- substr(exon$V9,9,23)
row.names(exon) <- 1:nrow(exon)


exon_range <- GRanges(seqnames = exon$V1, 
               		  ranges = IRanges(start = exon$V4, end = exon$V5), 
               		  strand = exon$V7, gene_id = exon$gene_id)


SE$AE_up <- as.numeric(do.call(rbind,strsplit(as.character(SE$loci), split = "[:-]"))[,3])+1
SE$AE_dn <- as.numeric(do.call(rbind,strsplit(as.character(SE$loci), split = "[:-]"))[,4])-1


counts_equal <- mcmapply(function(i){
  countOverlaps(GRanges(seqnames = SE$chr[i], ranges = IRanges(start = SE$AE_up[i], end = SE$AE_dn[i])), 
                exon_range, type="equal")
}, 1:nrow(SE), mc.cores = 20)

counts_any <- mcmapply(function(i){
  countOverlaps(GRanges(seqnames = SE$chr[i], ranges = IRanges(start = SE$AE_up[i], end = SE$AE_dn[i])), 
                exon_range, type="any")
}, 1:nrow(SE), mc.cores = 20)

head(data.frame(equal = counts_equal, any = counts_any))


SE$counts_equal <- counts_equal
SE$counts_any <- counts_any

SE$SE_type <- "Other"
SE[SE$counts_equal > 0, "SE_type"] <- "exon-skipping_exactly"
SE[SE$counts_any == 0, "SE_type"] <- "novel_exon-skipping"

table(SE$SE_type)
# exon-skipping_exactly   novel_exon-skipping                 Other
#                  4170                  1592                  1086



#### 6. Host gene of SE =========================================================================================================================



SE$ASrange <- paste(SE$chr, ":", SE$AE_up, "-", SE$AE_dn, sep = "")
exon_gene_id <- data.frame(ID = paste(exon$V1, ":", exon$V4, "-", exon$V5, sep = ""), gene_id = exon$gene_id, stringsAsFactors = F)

SE <- merge(x = SE, y = exon_gene_id, by.x = "ASrange", by.y = "ID", all.x = T)
length(table(na.omit(SE)[,"gene_id"]))
#[1] 3174
table(table(na.omit(SE)[,"gene_id"]))
#    1    2    3    4    5    6
# 2407  587  145   25    6    4

## merge PSI and SE information
#for(i in colnames(psi.co)){
#  SE[,paste(i, ".L_psi", sep = "")] <- psi.co[as.character(SE$right),i]
#  SE[,paste(i, ".R_psi", sep = "")] <- psi.co[as.character(SE$left),i]
#  SE[,i] <- rowMeans(cbind(psi.co[as.character(SE$left),i], psi.co[as.character(SE$right),i]), na.rm = F)
#}# as.vector() is necessary

dim(SE)
#[1] 6848   14

save(SE, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/SE/SE_type.RData")



#### 7. find A3SS ===============================================================================================================================



head(junc.as.2)
head(junc.as.2.ot)

junc.as.2.ot_uniq <- junc.as.2.ot[!duplicated(junc.as.2.ot),]
## 去除重复是因为要将所有的 AS SJ 都用起来，主要目的是搜集那些在 SE 的时候因为是出现 3 次但其实并不是 SE 的那些 SJ, 
## 它们之前出现三次，我们认为有可能是 SE ， 但是之后因为 AS insert 等原因判断出并不是 SE，但也出现了三次，
## 其中还有一次是重复的，所以把重复去掉。

sum(table(junc.as.2.ot_uniq$start)==2)

## 把真正的 SE 去掉：
head(junc.as.2.ot_uniq)
junc.as.2.ot_uniq$event <- do.call(rbind, strsplit(rownames(junc.as.2.ot_uniq), split = "\\."))[,1]

a3_raw <- junc.as.2.ot_uniq[junc.as.2.ot_uniq$start %in% names(table(junc.as.2.ot_uniq$start)[table(junc.as.2.ot_uniq$start)==2]), ]

sum(a3_raw$event %in% SE$left)
sum(a3_raw$event %in% SE$left == FALSE)

a3_raw <- a3_raw[a3_raw$event %in% SE$left == FALSE, ]
## 把真正的 SE 去掉


a3_raw <- dlply(a3_raw, c("chr","start"))
sum(unlist(lapply(a3_raw,function(x) length(unique(x$chr)) != 1)))
# if != 0, there are some problem!!!

names(a3_raw) <- sapply(a3_raw,function(x) paste(x$names,collapse="_"))

names(a3_raw)
as.vector(sapply(a3_raw,function(x) diff(as.numeric(x$end))))

cbind(names(a3_raw),as.vector(sapply(a3_raw,function(x) diff(x$end))))

sapply(a3_raw, function(x) mean(as.numeric(x$start)))
sapply(a3_raw, function(x) min(as.numeric(x$end)))
sapply(a3_raw, function(x) max(as.numeric(x$end)))
sapply(a3_raw, function(x) unique(x$chr))



a3 <- data.frame(row.names = names(a3_raw), name = names(a3_raw),
                 distan = as.character(sapply(a3_raw,function(x) diff(as.numeric(x$end)))),
                 chr = as.character(sapply(a3_raw, function(x) unique(x$chr))),
                 start = sapply(a3_raw, function(x) mean(as.numeric(x$start))),
                 end1 = sapply(a3_raw, function(x) min(as.numeric(x$end))),
                 end2 = sapply(a3_raw, function(x) max(as.numeric(x$end))),stringsAsFactors=F)



dim(a3)
head(a3)

a3SS <- a3[rownames(a3) %in% names(junc.as.2),]

counts_equal <- mcmapply(function(i){
  sum(countOverlaps(exon_range, GRanges(seqnames = a3SS$chr[i], ranges = IRanges(start = a3SS$end1[i]+1, end = a3SS$end2[i]+1)), 
                type="equal"))
}, 1:nrow(a3SS), mc.cores = 20)

counts_start <- mcmapply(function(i){
  sum(countOverlaps(exon_range, GRanges(seqnames = a3SS$chr[i], ranges = IRanges(start = a3SS$end1[i]+1, end = a3SS$end2[i]+1)), 
                type="start"))
}, 1:nrow(a3SS), mc.cores = 20)

counts_within <- mcmapply(function(i){
  sum(countOverlaps(exon_range, GRanges(seqnames = a3SS$chr[i], ranges = IRanges(start = a3SS$end1[i]+1, end = a3SS$end2[i]+1)), 
                type="within"))
}, 1:nrow(a3SS), mc.cores = 20)


length(counts_start)
#[1] 29612
sum(counts_start>0&counts_within==0&counts_equal==0)
#[1] 4083
sum(counts_start==1&counts_within==0&counts_equal==0)
#[1] 2745

a3SS$counts_start <- counts_start
a3SS$counts_within <- counts_within

A3SS <- a3SS[a3SS$counts_start == 1 & a3SS$counts_within == 0, ]
dim(A3SS)
#[1] 2745    8



#### 8. Host gene of A3SS =======================================================================================================================



queryHits_start <- mcmapply(function(i){
  queryHits(findOverlaps(exon_range, GRanges(seqnames = A3SS$chr[i], ranges = IRanges(start = A3SS$end1[i]+1, end = A3SS$end2[i]+1)), type="start"))
}, 1:nrow(A3SS), mc.cores = 20)

A3SS <- cbind(A3SS, data.frame(exon_range[queryHits_start,], stringsAsFactors=F))

save(A3SS, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/A3SS/A3SS_type.RData")



#### 9. find A5SS ===============================================================================================================================



junc.as.2.ot_uniq <- junc.as.2.ot[!duplicated(junc.as.2.ot, fromLast = TRUE),]

junc.as.2.ot_uniq$event <- do.call(rbind, strsplit(rownames(junc.as.2.ot_uniq), split = "\\."))[,1]

a5_raw <- junc.as.2.ot_uniq[junc.as.2.ot_uniq$end %in% names(table(junc.as.2.ot_uniq$end)[table(junc.as.2.ot_uniq$end)==2]), ]

a5_raw <- a5_raw[order(a5_raw$chr, a5_raw$end, a5_raw$start),]

a5_raw <- a5_raw[a5_raw$event %in% SE$right == FALSE, ]

a5_raw <- dlply(a5_raw, c("chr","end"))
sum(unlist(lapply(a5_raw,function(x) length(unique(x$chr)) != 1)))
# if != 0, there are some problem!!!


names(a5_raw) <- sapply(a5_raw,function(x) paste(x$names,collapse="_"))

a5 <- data.frame(name = names(a5_raw), row.names = names(a5_raw),
                 distan = as.integer(sapply(a5_raw,function(x) diff(as.numeric(x$start)))), 
                 chr = as.character(sapply(a5_raw, function(x) unique(x$chr))),
                 start1 = sapply(a5_raw, function(x) min(as.numeric(x$start))),
                 start2 = sapply(a5_raw, function(x) max(as.numeric(x$start))),
                 end = sapply(a5_raw, function(x) mean(as.numeric(x$end))), stringsAsFactors = F)


a5SS <- a5[rownames(a5) %in% names(junc.as.2),]

#identical(rownames(psi.co[row.names(A5SS),]),rownames(A5SS))
#identical(row.names(psi.co)[which(row.names(psi.co) %in% row.names(A5SS))],rownames(A5SS))


identical(names(junc.as.2)[which(names(junc.as.2) %in% row.names(a5SS))],rownames(a5SS))


counts_equal <- mcmapply(function(i){
  sum(countOverlaps(exon_range, GRanges(seqnames = a5SS$chr[i], ranges = IRanges(start = a5SS$start1[i]-1, end = a5SS$start2[i]-1)), 
                type="equal"))
}, 1:nrow(a5SS), mc.cores = 20)

counts_end <- mcmapply(function(i){
  sum(countOverlaps(exon_range, GRanges(seqnames = a5SS$chr[i], ranges = IRanges(start = a5SS$start1[i]-1, end = a5SS$start2[i]-1)), 
                type="end"))
}, 1:nrow(a5SS), mc.cores = 20)

counts_within <- mcmapply(function(i){
  sum(countOverlaps(exon_range, GRanges(seqnames = a5SS$chr[i], ranges = IRanges(start = a5SS$start1[i]-1, end = a5SS$start2[i]-1)), 
                type="within"))
}, 1:nrow(a5SS), mc.cores = 20)


length(counts_end)
#[1] 29419
sum(counts_end>0&counts_within==0&counts_equal==0)
#[1] 4059
sum(counts_end==1&counts_within==0&counts_equal==0)
#[1] 2744


a5SS$counts_end <- counts_end
a5SS$counts_within <- counts_within

A5SS <- a5SS[a5SS$counts_end == 1 & a5SS$counts_within == 0, ]
dim(A5SS)
#[1] 2744    8



#### 10. Host gene of A5SS ======================================================================================================================



queryHits_end <- mcmapply(function(i){
  queryHits(findOverlaps(exon_range, GRanges(seqnames = A5SS$chr[i], ranges = IRanges(start = A5SS$start1[i]-1, end = A5SS$start2[i]-1)), type="end"))
}, 1:nrow(A5SS), mc.cores = 20)

A5SS <- cbind(A5SS, data.frame(exon_range[queryHits_end,], stringsAsFactors=F))

dim(A5SS)
#[1] 2744   14

save(A5SS, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/A5SS/A5SS_type.RData")



#### 11. Find MXE ===============================================================================================================================



sum(table(junc.as.2.ot$start)==2)
sum(table(junc.as.2.ot$end)==2)

start_pos <- names(table(junc.as.2.ot$start)[table(junc.as.2.ot$start)==2])
end_pos <- names(table(junc.as.2.ot$end)[table(junc.as.2.ot$end)==2])

junc.as.2.ot[junc.as.2.ot$start %in% start_pos | junc.as.2.ot$end %in% end_pos, ]

## start1 < end2 < start2
## end1 < start1 < end2
## start < end1 < start1 < end2 < start2 < end

head(a3[order(a3$chr, a3$start, a3$end1, a3$end2),])
head(a5[order(a5$chr, a5$start1, a5$start2, a5$end),])


merge.data.frame(x = a5, y = a3, by = "name", all = T)

merge(x = data.table(a3, key = "name"),y = data.table(a5, key = "name"),all = T, by = "name")

a3$name %in% a5$name

a3_for_mxe <- data.frame(a3[order(a3$chr, a3$start, a3$end1, a3$end2),])
a5_for_mxe <- data.frame(a5[order(a5$chr, a5$start1, a5$start2, a5$end),])

## 思路是将 A3 与 A5 每两个都进行两两组合，然后根据位置（start < end1 < start1 < end2 < start2 < end）判断是否为 MXE

library(dplyr)
library(tidyr)

a3_for_mxe_split <- dlply(a3_for_mxe, .variables = "chr")
a5_for_mxe_split <- dlply(a5_for_mxe, .variables = "chr")

mxe_raw <- list()
for(i in 1:length(a3_for_mxe_split)){
	print(paste(i, "of", length(a3_for_mxe_split)))
	if(names(a5_for_mxe_split)[i] == names(a3_for_mxe_split)[i]){
		x = a3_for_mxe_split[[i]]
		y = a5_for_mxe_split[[i]]
		a5_for_mxe_rep <- bind_rows(replicate(nrow(x), y, simplify = FALSE))
		a3_for_mxe_rep <- bind_rows(replicate(nrow(y), x, simplify = FALSE))
		a3_for_mxe_rep <- data.frame(a3_for_mxe_rep[order(a3_for_mxe_rep$chr, a3_for_mxe_rep$start, a3_for_mxe_rep$end1, a3_for_mxe_rep$end2),])

		mxe_raw[[i]] <- cbind(a3_for_mxe_rep, a5_for_mxe_rep)
	}
}

mxe_raw <- lapply(mxe_raw, function(x){
	x[x[,"start1"] < x[,"end2"] & x[,"start1"] > x[,"end1"] & x[,"end2"] < x[,"start2"] & x[,"end2"] > x[,"start1"], ]
	})

MXE <- do.call(rbind, mxe_raw)

## we just neeed the A3 and A5 from the same chr
#MXE <- MXE[apply(MXE[,colnames(MXE) == "chr"], 1, function(x) sum(duplicated(as.character(x))) > 0),]




#### 12. validate MXE using GTF =================================================================================================================




counts_equal1 <- mcmapply(function(i){
  sum(countOverlaps(exon_range, GRanges(seqnames = MXE$chr[i], ranges = IRanges(start = MXE$end1[i]+1, end = MXE$start1[i]-1)), 
                type="equal"))
}, 1:nrow(MXE), mc.cores = 20)

counts_equal2 <- mcmapply(function(i){
  sum(countOverlaps(exon_range, GRanges(seqnames = MXE$chr[i], ranges = IRanges(start = MXE$end2[i]+1, end = MXE$start2[i]-1)), 
                type="equal"))
}, 1:nrow(MXE), mc.cores = 20)


counts_any1 <- mcmapply(function(i){
  sum(countOverlaps(exon_range, GRanges(seqnames = MXE$chr[i], ranges = IRanges(start = MXE$end1[i]+1, end = MXE$start1[i]-1)), 
                type="any"))
}, 1:nrow(MXE), mc.cores = 20)

counts_any2 <- mcmapply(function(i){
  sum(countOverlaps(exon_range, GRanges(seqnames = MXE$chr[i], ranges = IRanges(start = MXE$end2[i]+1, end = MXE$start2[i]-1)), 
                type="any"))
}, 1:nrow(MXE), mc.cores = 20)


MXE$counts_equal1 <- counts_equal1
MXE$counts_equal2 <- counts_equal2
MXE$counts_any1 <- counts_any1
MXE$counts_any2 <- counts_any2


queryHits_equal1 <- mcmapply(function(i){
  queryHits(findOverlaps(exon_range, GRanges(seqnames = MXE$chr[i], ranges = IRanges(start = MXE$end1[i]+1, end = MXE$start1[i]-1)), type="equal"))
}, 1:nrow(MXE), mc.cores = 20)

queryHits_equal2 <- mcmapply(function(i){
  queryHits(findOverlaps(exon_range, GRanges(seqnames = MXE$chr[i], ranges = IRanges(start = MXE$end2[i]+1, end = MXE$start2[i]-1)), type="equal"))
}, 1:nrow(MXE), mc.cores = 20)

data.frame(exon_range[unlist(queryHits_equal1),], stringsAsFactors=F)

MXE[MXE$counts_equal1 == 1, "gene_id1"] <- data.frame(exon_range[unlist(queryHits_equal1),], stringsAsFactors=F)$gene_id
MXE[MXE$counts_equal2 == 1, "gene_id2"] <- data.frame(exon_range[unlist(queryHits_equal2),], stringsAsFactors=F)$gene_id

sum(MXE$counts_equal1 == 1 & MXE$counts_equal2 == 1 | MXE$counts_equal1 == 1 & MXE$counts_any2 == 0 | MXE$counts_any1 == 0 & MXE$counts_equal2 == 1 | MXE$counts_any1 == 0 & MXE$counts_any2 == 0)
# [1] 907

MXE <- MXE[MXE$counts_equal1 == 1 & MXE$counts_equal2 == 1 | MXE$counts_equal1 == 1 & MXE$counts_any2 == 0 | MXE$counts_any1 == 0 & MXE$counts_equal2 == 1 | MXE$counts_any1 == 0 & MXE$counts_any2 == 0, ]
dim(MXE)
# [1] 907   18

sum(MXE[MXE$counts_equal1 == 1 & MXE$counts_equal2 == 1, "gene_id1"] != MXE[MXE$counts_equal1 == 1 & MXE$counts_equal2 == 1, "gene_id2"])
# [1] 9

sum(paste(MXE$chr, ":", MXE$start1, "-", MXE$end2, sep = "") %in% row.names(SJ_tu))
#[1] 771
sum(!paste(MXE$chr, ":", MXE$start1, "-", MXE$end2, sep = "") %in% row.names(SJ_tu))
#[1] 136

MXE <- MXE[!paste(MXE$chr, ":", MXE$start1, "-", MXE$end2, sep = "") %in% row.names(SJ_tu), ]

save(MXE, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/MXE/MXE_type.RData")



#### 13. AFE ====================================================================================================================================




load(file = "/mnt/data5/BGI/UCB/tangchao/DSU2/Homo_sapiens.GRCh38.87.gtf_AFE_index.RData")

length(junc.as.2)
#[1] 72769
length(AFE_index)
#[1] 12799

afe <- mclapply(junc.as.2, function(y){
	sum(unlist(lapply(AFE_index, function(x) sum(!y$names %in% x$AFSJ)==0)))
	}, mc.cores = 10)

table(unlist(afe))
#     0     1     2
# 71339  1412    18

afe.as <- junc.as.2[afe==1]
afe_gene <- mclapply(afe.as, function(y){
	unique(AFE_index[[which(lapply(AFE_index, function(x) sum(!y$names %in% x$AFSJ)==0)==1)]]$gene_id)
}, mc.cores = 10)
afe_gene_name <- mclapply(afe.as, function(y){
	unique(AFE_index[[which(lapply(AFE_index, function(x) sum(!y$names %in% x$AFSJ)==0)==1)]]$gene_name)
}, mc.cores = 10)
afe_gene_biotype <- mclapply(afe.as, function(y){
	unique(AFE_index[[which(lapply(AFE_index, function(x) sum(!y$names %in% x$AFSJ)==0)==1)]]$gene_biotype)
}, mc.cores = 10)
afe_gene_strand <- mclapply(afe.as, function(y){
	unique(AFE_index[[which(lapply(AFE_index, function(x) sum(!y$names %in% x$AFSJ)==0)==1)]]$strand)
}, mc.cores = 10)


AFE <- data.frame(loci = names(afe.as), 
				  gene_id = as.character(unlist(afe_gene)),
				  gene_name = as.character(unlist(afe_gene_name)),
				  strand = as.character(unlist(afe_gene_strand)),
				  biotype = as.character(unlist(afe_gene_biotype)), stringsAsFactors = F)

save(AFE, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/AFE/AFE_type.RData")



#### 14. ALE ====================================================================================================================================



load(file = "/mnt/data5/BGI/UCB/tangchao/DSU2/Homo_sapiens.GRCh38.87.gtf_ALE_index.RData")

ale <- mclapply(junc.as.2, function(y){
	sum(unlist(lapply(ALE_index, function(x) sum(!y$names %in% x$ALSJ)==0)))
	}, mc.cores = 10)

table(unlist(ale))
#     0     1     2
# 71892   870     7

ale.as <- junc.as.2[ale==1]
ale_gene <- mclapply(ale.as, function(y){
	unique(ALE_index[[which(lapply(ALE_index, function(x) sum(!y$names %in% x$ALSJ)==0)==1)]]$gene_id)
}, mc.cores = 10)

ale_gene_name <- mclapply(ale.as, function(y){
	unique(ALE_index[[which(lapply(ALE_index, function(x) sum(!y$names %in% x$ALSJ)==0)==1)]]$gene_name)
}, mc.cores = 10)

ale_gene_biotype <- mclapply(ale.as, function(y){
	unique(ALE_index[[which(lapply(ALE_index, function(x) sum(!y$names %in% x$ALSJ)==0)==1)]]$gene_biotype)
}, mc.cores = 10)

ale_gene_strand <- mclapply(ale.as, function(y){
	unique(ALE_index[[which(lapply(ALE_index, function(x) sum(!y$names %in% x$ALSJ)==0)==1)]]$strand)
}, mc.cores = 10)


ALE <- data.frame(loci = names(ale.as), 
				  gene_id = as.character(unlist(ale_gene)),
				  gene_name = as.character(unlist(ale_gene_name)),
				  strand = as.character(unlist(ale_gene_strand)),
				  biotype = as.character(unlist(ale_gene_biotype)), stringsAsFactors = F)

save(ALE, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/ALE/ALE_type.RData")



#### 15. PSI of SE ==============================================================================================================================



load("/mnt/data5/BGI/UCB/tangchao/data/SJ/SJ_merged_raw_te(all_more_than_10).RData")
te_table <- data.frame(as.data.frame(te[,-1]), row.names = te[[1]])
colnames(te_table) <- do.call(rbind,strsplit(colnames(te_table),"_SJ"))[,1]
te_table <- te_table[, row.names(Cell_type)]
dim(te_table)
# [1] 461646   3039
f_dowle3(te_table)

SE_loci <- data.frame(do.call(rbind,strsplit(as.character(SE$loci), "[:-]")), stringsAsFactors=F)
head(SE_loci)
SE_loci$SJ12 <- paste(SE_loci[,1], ":", SE_loci[,2], "-", SE_loci[,3], sep = "")
SE_loci$SJ23 <- paste(SE_loci[,1], ":", SE_loci[,4], "-", SE_loci[,5], sep = "")
SE_loci$SJ13 <- paste(SE_loci[,1], ":", SE_loci[,2], "-", SE_loci[,5], sep = "")

sum(!SE_loci$SJ12 %in% row.names(te_table))
sum(!SE_loci$SJ23 %in% row.names(te_table))
sum(!SE_loci$SJ13 %in% row.names(te_table))

sj_list <- apply(SE_loci, 1, function(x) te_table[x[6:8],])
psi_list <- lapply(sj_list, function(x) colSums(x[1:2,])/colSums(x))
table(unlist(lapply(psi_list, function(x) sum(!is.na(x)))))

psi <- do.call(rbind, psi_list)

SE_psi <- cbind(SE, psi)
save(SE_psi, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/SE/SE_psi.RData")



#### 16. PSI of A3SS ============================================================================================================================



A3SS_loci <- data.frame(do.call(rbind,strsplit(as.character(A3SS$name), "_")), stringsAsFactors=F)
head(A3SS_loci)

sum(!A3SS_loci$X1 %in% row.names(te_table))
sum(!A3SS_loci$X2 %in% row.names(te_table))

sj_list <- apply(A3SS_loci, 1, function(x) te_table[x,])
psi_list <- lapply(sj_list, function(x) x[1,]/colSums(x))
head(table(unlist(lapply(psi_list, function(x) sum(!is.na(x))))))

psi <- do.call(rbind, psi_list)

A3SS_psi <- cbind(A3SS, psi)
save(A3SS_psi, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/A3SS/A3SS_psi.RData")



#### 17. PSI of A5SS ============================================================================================================================



A5SS_loci <- data.frame(do.call(rbind,strsplit(as.character(A5SS$name), "_")), stringsAsFactors=F)
head(A5SS_loci)

sum(!A5SS_loci$X1 %in% row.names(te_table))
sum(!A5SS_loci$X2 %in% row.names(te_table))

sj_list <- apply(A5SS_loci, 1, function(x) te_table[x,])
psi_list <- lapply(sj_list, function(x) x[2,]/colSums(x))
head(table(unlist(lapply(psi_list, function(x) sum(!is.na(x))))))

psi <- do.call(rbind, psi_list)

A5SS_psi <- cbind(A5SS, psi)
save(A5SS_psi, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/A5SS/A5SS_psi.RData")



#### 18. PSI of AFE =============================================================================================================================



AFE_loci <- data.frame(do.call(rbind,strsplit(as.character(AFE$loci), "_")), stringsAsFactors=F)
head(AFE_loci)

sum(!AFE_loci$X1 %in% row.names(te_table))
sum(!AFE_loci$X2 %in% row.names(te_table))

sj_list <- apply(AFE_loci, 1, function(x) te_table[x,])
psi_list <- list()
for(i in 1:length(sj_list)){
	print(paste(i, "of", length(sj_list)))
	if(AFE[i,"strand"]=="-"){
		psi_list[[i]] <- sj_list[[i]][1,]/colSums(sj_list[[i]])
	}else{
		psi_list[[i]] <- sj_list[[i]][2,]/colSums(sj_list[[i]])
	}
}

head(table(unlist(lapply(psi_list, function(x) sum(!is.na(x))))))

psi <- do.call(rbind, psi_list)

AFE_psi <- cbind(AFE, psi)
save(AFE_psi, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/AFE/AFE_psi.RData")



#### 19. PSI of ALE =============================================================================================================================



ALE_loci <- data.frame(do.call(rbind,strsplit(as.character(ALE$loci), "_")), stringsAsFactors=F)
head(ALE_loci)

sum(!ALE_loci$X1 %in% row.names(te_table))
sum(!ALE_loci$X2 %in% row.names(te_table))

sj_list <- apply(ALE_loci, 1, function(x) te_table[x,])
psi_list <- list()
for(i in 1:length(sj_list)){
	print(paste(i, "of", length(sj_list)))
	if(AFE[i,"strand"]=="-"){
		psi_list[[i]] <- sj_list[[i]][2,]/colSums(sj_list[[i]])
	}else{
		psi_list[[i]] <- sj_list[[i]][1,]/colSums(sj_list[[i]])
	}
}

head(table(unlist(lapply(psi_list, function(x) sum(!is.na(x))))))

psi <- do.call(rbind, psi_list)

ALE_psi <- cbind(ALE, psi)
save(ALE_psi, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/ALE/ALE_psi.RData")



#### 20. PSI of MXE =============================================================================================================================



mxe_sj <- data.frame(SJ12 = do.call(rbind,strsplit(MXE$name, "_"))[,1], SJ13 = do.call(rbind,strsplit(MXE$name, "_"))[,2],
		   			 SJ24 = do.call(rbind,strsplit(MXE$name.1, "_"))[,1], SJ34 = do.call(rbind,strsplit(MXE$name.1, "_"))[,2], stringsAsFactors = F)

colSums(SJ_tu[as.character(mxe_sj[1,c(2,4)]),])/colSums(SJ_tu[as.character(mxe_sj[1,]),])

sj_tab <- list()
for(i in 1:nrow(MXE)){
	sj_tab[[i]] <- SJ_tu[as.character(mxe_sj[i,]),]
}

for(i in 1:length(sj_tab)){
	tmp <- sj_tab[[i]]
	tmp[,colSums(tmp!=0)==1] <- 0
	sj_tab[[i]] <- tmp
}

for(i in 1:length(sj_tab)){
	tmp <- sj_tab[[i]]
	tmp[,colSums(tmp[1:2,]!=0)==0] <- 0
	sj_tab[[i]] <- tmp
}

mxe_psi <- list()
for(i in 1:length(sj_tab)){
	mxe_psi[[i]] <- colSums(sj_tab[[i]][c(2,4),])/colSums(sj_tab[[i]])
}

MXE_psi <- cbind(MXE,do.call(rbind, mxe_psi))
dim(MXE_psi)
#[1] 136 3057
row.names(MXE_psi) <- 1:nrow(MXE_psi)
save(MXE_psi, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/MXE/MXE_psi.RData")













