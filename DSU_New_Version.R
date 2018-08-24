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
library(data.table)
depar <- par()

## path of gtf
pfgtf <- file.path("/mnt/data1/reference/ensembl/human/Homo_sapiens.GRCh38.87.gtf")

### load SJs merge file =========================================================================================================================

load("/mnt/data5/BGI/UCB/tangchao/data/SJ/SJ_merged_raw_te(all_more_than_10).RData")

Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)

names(te) <- do.call(rbind,strsplit(names(te),"_SJ"))[,1]
te <- te[, c("name", row.names(Cell_type)), with=F]


#require(data.table)  # v1.6.6
#require(gdata)
#f_dowle3 = function(DT) {
#  # either of the following for loops
#  
#  # by name :
#  for (j in names(DT))
#    set(DT,which(is.na(DT[[j]])),j,0)
#  
#  # or by number (slightly faster than by name) :
#  for (j in seq_len(ncol(DT)))
#    set(DT,which(is.na(DT[[j]])),j,0)
#}#

#f_dowle3(te_table)

## We only care about chr1-22 and X Y
chr_tu <- do.call(rbind, strsplit(te[[1]], split = ":"))[,1] %in% c(1:22,"X","Y")
sum(chr_tu)
# [1] 460327

te <- te[chr_tu,]
dim(te)
# [1] 460327   3040



#### 2. Identify alternative splicing events ====================================================================================================



#junc.names=do.call(rbind,(strsplit(sub(x=rownames(junc),
#                                       pattern="^([[:alnum:]]+):([[:digit:]]+)-([[:digit:]]+)$",
#                                       replace="\\1;\\2;\\3"),split=";")))
junc.names=do.call(rbind,strsplit(te[[1]],split="[:-]"))
#alnum-- Letters and Numbers；digit -- Numbers
colnames(junc.names) <- c("chr","start","end")
rownames(junc.names) <- te[[1]]
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

## junction just has two alternative starts/ends:
junc.as2 <- junc.as[sapply(junc.as,nrow)==2]
names(junc.as2) <- sapply(junc.as2,function(x) paste(x$names,collapse="_"))

table(sapply(junc.as.2,nrow))

te[!is.na(te[[2]]), name]

junc.as2 <- do.call(list, mcmapply(function(x){
  print(x)
  junc.names=do.call(rbind,strsplit(te[!is.na(te[[x]]), name],split="[:-]"))
  colnames(junc.names) <- c("chr","start","end")
  rownames(junc.names) <- te[!is.na(te[[x]]), name]
  junc.names=data.frame(junc.names,stringsAsFactors=F)
  junc.names$start=as.integer(junc.names$start)
  junc.names$end=as.integer(junc.names$end)
  junc.names=junc.names[order(junc.names$chr,junc.names$start,junc.names$end),]
  junc.names$names=rownames(junc.names)

  junc.as=do.call(c,mclapply(unique(junc.names$chr),function(chr) {
    junc.chr=junc.names[junc.names$chr==chr,]
    same.start=dlply(junc.chr,c("start"),function(x) x)
    same.end=dlply(junc.chr,c("end"),function(x) {if(nrow(x)>1) {return(x)} else {return(NULL)}})
    #print(chr)
    return(c(same.start,same.end[!sapply(same.end,is.null)]))
    },mc.cores=2))
  junc.as2 <- junc.as[sapply(junc.as,nrow)==2]
  names(junc.as2) <- sapply(junc.as2,function(x) paste(x$names,collapse="_"))
  return(junc.as2)
}, 1:ncol(te), mc.cores=10))


save(junc.as2, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/junc.as2.RData")



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



SE_list <- mclapply(junc.as2, function(tang){
  tang2 <- data.frame(unlist(lapply(tang, function(x) {if(diff(x$start)==0){return("SameStart")}else{return("SameEnd")}})), stringsAsFactors=F)
  colnames(tang2) <- "type"
  tang2$as2 <- names(tang)
  tang2$chr <- as.character(unlist(lapply(tang, function(x) {unique(x$chr)})))
  tang3 <- dlply(tang2, .variables = "chr")
  if(sum(unlist(lapply(tang3, nrow))>1)==0){
    tang3_se <- NULL
  }else{
    tang3_se <- mcmapply(function(x){
      if(length(unique(x$type))==1){
        return(NULL)
      }else{
        samestart <- subset.data.frame(x, x$type=="SameStart")
        samesend <- subset.data.frame(x, x$type=="SameEnd")
        samestart_rep <- bind_rows(replicate(nrow(samesend), samestart, simplify = FALSE))
        sameend_rep <- bind_rows(replicate(nrow(samestart), samesend, simplify = FALSE))
        sameend_rep <- data.frame(sameend_rep[order(sameend_rep$as2),], stringsAsFactors=F)
        se_rep <- cbind(samestart_rep, sameend_rep)
        se_rep_loci <- data.frame(do.call(rbind, strsplit(paste(se_rep[,2], se_rep[,5], sep = "_"), split="[:_-]")), stringsAsFactors=F)
        if(sum(se_rep_loci$X5 == se_rep_loci$X8 & se_rep_loci$X6 == se_rep_loci$X9 & as.numeric(se_rep_loci$X3) < as.numeric(se_rep_loci$X11))==0){
          return(NULL)
        }else{
          loci_tu <- se_rep_loci[se_rep_loci$X5 == se_rep_loci$X8 & se_rep_loci$X6 == se_rep_loci$X9 & as.numeric(se_rep_loci$X3) < as.numeric(se_rep_loci$X11),]
          se_loci <- paste(loci_tu$X1,":",loci_tu$X2, "-", loci_tu$X3, "-", loci_tu$X11, "-", loci_tu$X12, sep="")
          return(se_loci)
        }   
      }
      
    }, tang3, mc.cores = 1)
    tang3_se <- data.frame(loci = as.character(unlist(tang3_se)), stringsAsFactors=F)
    if(nrow(tang3_se)==0){
      tang3_se <- NULL
    }else{
      tmp <- do.call(rbind, strsplit(tang3_se$loci, "[:-]"))
      tang3_se$sj12 <- paste(tmp[,1], ":", tmp[,2], "-", tmp[,3], sep = "")
      tang3_se$sj23 <- paste(tmp[,1], ":", tmp[,4], "-", tmp[,5], sep = "")
      tang3_se$sj13 <- paste(tmp[,1], ":", tmp[,2], "-", tmp[,5], sep = "")
    }
  }
  return(tang3_se)
}, mc.cores = 10)

SE <- do.call(rbind, SE_list)
head(SE)
dim(SE[!duplicated(SE),])
# [1] 17670     4
SE <- SE[!duplicated(SE),]
sum(SE$sj12 %in% te[[1]])
#[1] 17670
sum(SE$sj13 %in% te[[1]])
#[1] 17670
sum(SE$sj23 %in% te[[1]])
#[1] 17670





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

SE$chr <- as.character(do.call(rbind,strsplit(as.character(SE$loci), split = "[:-]"))[,1])
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
#                 10372                  3398                  3900


#### 6. Host gene of SE =========================================================================================================================



SE$ASrange <- paste(SE$chr, ":", SE$AE_up, "-", SE$AE_dn, sep = "")
exon_gene_id <- data.frame(ID = paste(exon$V1, ":", exon$V4, "-", exon$V5, sep = ""), gene_id = exon$gene_id, stringsAsFactors = F)

SE <- merge(x = SE, y = exon_gene_id, by.x = "ASrange", by.y = "ID", all.x = T)
length(table(na.omit(SE)[,"gene_id"]))
#[1] 5896
table(table(na.omit(SE)[,"gene_id"]))
#    1    2    3    4    5    6    7    8    9   11   12   24
# 3364 1457  620  236  123   56   15   16    4    2    2    1
## merge PSI and SE information
#for(i in colnames(psi.co)){
#  SE[,paste(i, ".L_psi", sep = "")] <- psi.co[as.character(SE$right),i]
#  SE[,paste(i, ".R_psi", sep = "")] <- psi.co[as.character(SE$left),i]
#  SE[,i] <- rowMeans(cbind(psi.co[as.character(SE$left),i], psi.co[as.character(SE$right),i]), na.rm = F)
#}# as.vector() is necessary

dim(SE)
#[1] 17670   12

save(SE, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/SE/SE_type_new.RData")



#### 7. find A3SS ===============================================================================================================================



all_same_start <- data.frame(unlist(lapply(junc.as2, function(x) names(x)[unlist(lapply(x, function(x) diff(x$start)==0))] )), stringsAsFactors=F)
colnames(all_same_start) <- "as2"
all_same_start$sj12 <- do.call(rbind, strsplit(all_same_start$as2, split="_"))[,1]
all_same_start$sj13 <- do.call(rbind, strsplit(all_same_start$as2, split="_"))[,2]
all_same_start$chr <- do.call(rbind, strsplit(all_same_start$as2, split="[:_-]"))[,1]

all_same_start$AE_up <- as.numeric(do.call(rbind, strsplit(all_same_start$as2, split="[:_-]"))[,3])
all_same_start$AE_dn <- as.numeric(do.call(rbind, strsplit(all_same_start$as2, split="[:_-]"))[,6])

sum(all_same_start$sj12 %in% c(SE$sj12, SE$sj13, SE$sj23))
#[1] 96812
dim(all_same_start)
#[1] 152018      6
sum(all_same_start$sj13 %in% c(SE$sj12, SE$sj13, SE$sj23))
#[1] 92302
sum(all_same_start$sj13 %in% c(SE$sj12, SE$sj13, SE$sj23) | all_same_start$sj12 %in% c(SE$sj12, SE$sj13, SE$sj23))
#[1] 106257

a3SS <- all_same_start[!(all_same_start$sj13 %in% c(SE$sj12, SE$sj13, SE$sj23) | all_same_start$sj12 %in% c(SE$sj12, SE$sj13, SE$sj23)), ]
a3SS <- unique(a3SS)
dim(a3SS)
head(a3SS)

#counts_equal <- vector()
#for(i in 1:nrow(a3SS)){
#  print(paste(i, "of", nrow(a3SS)))
#  counts_equal[i] <- sum(countOverlaps(exon_range, GRanges(seqnames = a3SS$chr[i], ranges = IRanges(start = a3SS$AE_up[i]+1, end = a3SS$AE_dn[i]+1)), type="equal"))
#}

counts_equal <- mcmapply(function(i){
  return(sum(countOverlaps(exon_range, GRanges(seqnames = a3SS$chr[i], ranges = IRanges(start = a3SS$AE_up[i]+1, end = a3SS$AE_dn[i]+1)), 
                type="equal")))
}, 1:nrow(a3SS), mc.cores = 20)

counts_start <- mcmapply(function(i){
  return(sum(countOverlaps(exon_range, GRanges(seqnames = a3SS$chr[i], ranges = IRanges(start = a3SS$AE_up[i]+1, end = a3SS$AE_dn[i]+1)), type="start")))
}, 1:nrow(a3SS), mc.cores = 10)

counts_within <- mcmapply(function(i){
  return(sum(countOverlaps(exon_range, GRanges(seqnames = a3SS$chr[i], ranges = IRanges(start = a3SS$AE_up[i]+1, end = a3SS$AE_dn[i]+1)), type="within")))
}, 1:nrow(a3SS), mc.cores = 10)


length(counts_start)
#[1] 27832
sum(counts_start>0&counts_within==0&counts_equal==0)
#[1] 4507
sum(counts_start==1&counts_within==0&counts_equal==0)
#[1] 3050

a3SS$counts_start <- counts_start
a3SS$counts_within <- counts_within

A3SS <- a3SS[a3SS$counts_start == 1 & a3SS$counts_within == 0, ]
dim(A3SS)
#[1] 3050    8



#### 8. Host gene of A3SS =======================================================================================================================



queryHits_start <- mcmapply(function(i){
  return(queryHits(findOverlaps(exon_range, GRanges(seqnames = A3SS$chr[i], ranges = IRanges(start = A3SS$AE_up[i]+1, end = A3SS$AE_dn[i]+1)), type="start")))
}, 1:nrow(A3SS), mc.cores = 20)

A3SS <- cbind(A3SS, data.frame(exon_range[queryHits_start,], stringsAsFactors=F))

save(A3SS, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/A3SS/A3SS_type_new.RData")



#### 9. find A5SS ===============================================================================================================================



all_same_end <- data.frame(unlist(lapply(junc.as2, function(x) names(x)[unlist(lapply(x, function(x) diff(x$end)==0))] )), stringsAsFactors=F)
colnames(all_same_end) <- "as2"
all_same_end$sj13 <- do.call(rbind, strsplit(all_same_end$as2, split="_"))[,1]
all_same_end$sj23 <- do.call(rbind, strsplit(all_same_end$as2, split="_"))[,2]
all_same_end$chr <- do.call(rbind, strsplit(all_same_end$as2, split="[:_-]"))[,1]

all_same_end$AE_up <- as.numeric(do.call(rbind, strsplit(all_same_end$as2, split="[:_-]"))[,2])
all_same_end$AE_dn <- as.numeric(do.call(rbind, strsplit(all_same_end$as2, split="[:_-]"))[,5])

sum(all_same_end$sj13 %in% c(SE$sj12, SE$sj13, SE$sj23))
#[1] 95294
dim(all_same_end)
#[1] 166009      6
sum(all_same_end$sj13 %in% c(SE$sj12, SE$sj13, SE$sj23))
#[1] 95294
sum(all_same_end$sj13 %in% c(SE$sj12, SE$sj13, SE$sj23) | all_same_end$sj23 %in% c(SE$sj12, SE$sj13, SE$sj23))
#[1] 114601

a5SS <- all_same_end[!(all_same_end$sj13 %in% c(SE$sj12, SE$sj13, SE$sj23) | all_same_end$sj23 %in% c(SE$sj12, SE$sj13, SE$sj23)), ]
dim(a5SS)
# [1] 51408     6
a5SS <- unique(a5SS)
dim(a5SS)
#[1] 27776     6
head(a5SS)


#identical(rownames(psi.co[row.names(A5SS),]),rownames(A5SS))
#identical(row.names(psi.co)[which(row.names(psi.co) %in% row.names(A5SS))],rownames(A5SS))

counts_equal <- mcmapply(function(i){
  return(sum(countOverlaps(exon_range, GRanges(seqnames = a5SS$chr[i], ranges = IRanges(start = a5SS$AE_up[i]-1, end = a5SS$AE_dn[i]-1)), 
                type="equal")))
}, 1:nrow(a5SS), mc.cores = 20)

counts_end <- mcmapply(function(i){
  return(sum(countOverlaps(exon_range, GRanges(seqnames = a5SS$chr[i], ranges = IRanges(start = a5SS$AE_up[i]-1, end = a5SS$AE_dn[i]-1)), 
                type="end")))
}, 1:nrow(a5SS), mc.cores = 20)

counts_within <- mcmapply(function(i){
  return(sum(countOverlaps(exon_range, GRanges(seqnames = a5SS$chr[i], ranges = IRanges(start = a5SS$AE_up[i]-1, end = a5SS$AE_dn[i]-1)), 
                type="within")))
}, 1:nrow(a5SS), mc.cores = 20)


length(counts_end)
#[1] 27776
sum(counts_end>0&counts_within==0&counts_equal==0)
#[1] 4451
sum(counts_end==1&counts_within==0&counts_equal==0)
#[1] 3036


a5SS$counts_end <- counts_end
a5SS$counts_within <- counts_within

A5SS <- a5SS[a5SS$counts_end == 1 & a5SS$counts_within == 0, ]
dim(A5SS)
#[1] 3036    8



#### 10. Host gene of A5SS ======================================================================================================================



queryHits_end <- mcmapply(function(i){
  return(queryHits(findOverlaps(exon_range, GRanges(seqnames = A5SS$chr[i], ranges = IRanges(start = A5SS$AE_up[i]-1, end = A5SS$AE_dn[i]-1)), type="end")))
}, 1:nrow(A5SS), mc.cores = 20)

A5SS <- cbind(A5SS, data.frame(exon_range[queryHits_end,], stringsAsFactors=F))

dim(A5SS)
#[1] 3036   14

save(A5SS, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/A5SS/A5SS_type_new.RData")



#### 11. Find MXE ===============================================================================================================================



a3_for_mxe <- a3SS[!a3SS$as2 %in% A3SS$as2, 1:4]
a3_for_mxe$start <- as.numeric(do.call(rbind, strsplit(a3_for_mxe$as2, "[:_-]"))[,2])
a3_for_mxe$end1 <- as.numeric(do.call(rbind, strsplit(a3_for_mxe$as2, "[:_-]"))[,3])
a3_for_mxe$end2 <- as.numeric(do.call(rbind, strsplit(a3_for_mxe$as2, "[:_-]"))[,6])
a3_for_mxe <- data.frame(a3_for_mxe[order(a3_for_mxe$chr, a3_for_mxe$start, a3_for_mxe$end1, a3_for_mxe$end2),])


a5_for_mxe <- a5SS[!a5SS$as2 %in% A5SS$as2, 1:4]
a5_for_mxe$start1 <- as.numeric(do.call(rbind, strsplit(a5_for_mxe$as2, "[:_-]"))[,2])
a5_for_mxe$start2 <- as.numeric(do.call(rbind, strsplit(a5_for_mxe$as2, "[:_-]"))[,5])
a5_for_mxe$end <- as.numeric(do.call(rbind, strsplit(a5_for_mxe$as2, "[:_-]"))[,6])
a5_for_mxe <- data.frame(a5_for_mxe[order(a5_for_mxe$chr, a5_for_mxe$start1, a5_for_mxe$start2, a5_for_mxe$end),])

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
dim(MXE)
# [1] 6368   14
## we just neeed the A3 and A5 from the same chr
#MXE <- MXE[apply(MXE[,colnames(MXE) == "chr"], 1, function(x) sum(duplicated(as.character(x))) > 0),]




#### 12. validate MXE using GTF =================================================================================================================




counts_equal1 <- mcmapply(function(i){
  return(sum(countOverlaps(exon_range, GRanges(seqnames = MXE$chr[i], ranges = IRanges(start = MXE$end1[i]+1, end = MXE$start1[i]-1)), 
                type="equal")))
}, 1:nrow(MXE), mc.cores = 20)

counts_equal2 <- mcmapply(function(i){
  return(sum(countOverlaps(exon_range, GRanges(seqnames = MXE$chr[i], ranges = IRanges(start = MXE$end2[i]+1, end = MXE$start2[i]-1)), 
                type="equal")))
}, 1:nrow(MXE), mc.cores = 20)


counts_any1 <- mcmapply(function(i){
  return(sum(countOverlaps(exon_range, GRanges(seqnames = MXE$chr[i], ranges = IRanges(start = MXE$end1[i]+1, end = MXE$start1[i]-1)), 
                type="any")))
}, 1:nrow(MXE), mc.cores = 20)

counts_any2 <- mcmapply(function(i){
  return(sum(countOverlaps(exon_range, GRanges(seqnames = MXE$chr[i], ranges = IRanges(start = MXE$end2[i]+1, end = MXE$start2[i]-1)), 
                type="any")))
}, 1:nrow(MXE), mc.cores = 20)


MXE$counts_equal1 <- counts_equal1
MXE$counts_equal2 <- counts_equal2
MXE$counts_any1 <- counts_any1
MXE$counts_any2 <- counts_any2


queryHits_equal1 <- mcmapply(function(i){
  return(queryHits(findOverlaps(exon_range, GRanges(seqnames = MXE$chr[i], ranges = IRanges(start = MXE$end1[i]+1, end = MXE$start1[i]-1)), type="equal")))
}, 1:nrow(MXE), mc.cores = 20)

queryHits_equal2 <- mcmapply(function(i){
  return(queryHits(findOverlaps(exon_range, GRanges(seqnames = MXE$chr[i], ranges = IRanges(start = MXE$end2[i]+1, end = MXE$start2[i]-1)), type="equal")))
}, 1:nrow(MXE), mc.cores = 20)

data.frame(exon_range[unlist(queryHits_equal1),], stringsAsFactors=F)

MXE[MXE$counts_equal1 == 1, "gene_id1"] <- data.frame(exon_range[unlist(queryHits_equal1),], stringsAsFactors=F)$gene_id
MXE[MXE$counts_equal2 == 1, "gene_id2"] <- data.frame(exon_range[unlist(queryHits_equal2),], stringsAsFactors=F)$gene_id

sum(MXE$counts_equal1 == 1 & MXE$counts_equal2 == 1 | MXE$counts_equal1 == 1 & MXE$counts_any2 == 0 | MXE$counts_any1 == 0 & MXE$counts_equal2 == 1 | MXE$counts_any1 == 0 & MXE$counts_any2 == 0)
# [1] 808

MXE <- MXE[MXE$counts_equal1 == 1 & MXE$counts_equal2 == 1 | MXE$counts_equal1 == 1 & MXE$counts_any2 == 0 | MXE$counts_any1 == 0 & MXE$counts_equal2 == 1 | MXE$counts_any1 == 0 & MXE$counts_any2 == 0, ]
dim(MXE)
# [1] 808   20

sum(MXE[MXE$counts_equal1 == 1 & MXE$counts_equal2 == 1, "gene_id1"] != MXE[MXE$counts_equal1 == 1 & MXE$counts_equal2 == 1, "gene_id2"])
# [1] 5

sum(paste(MXE$chr, ":", MXE$start1, "-", MXE$end2, sep = "") %in% te[[1]])
#[1] 202
sum(!paste(MXE$chr, ":", MXE$start1, "-", MXE$end2, sep = "") %in% te[[1]])
#[1] 606

MXE <- MXE[!paste(MXE$chr, ":", MXE$start1, "-", MXE$end2, sep = "") %in% te[[1]], ]
dim(MXE)
#[1] 606  20
save(MXE, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/MXE/MXE_type_new.RData")



#### 13. AFE ====================================================================================================================================



load(file = "/mnt/data5/BGI/UCB/tangchao/DSU2/Homo_sapiens.GRCh38.87.gtf_AFE_index.RData")

length(AFE_index)
#[1] 12799
names(AFE_index) <- lapply(AFE_index, function(x) paste(x$AFSJ[order(x$AFSJ)], collapse="_", sep = ""))

all_as2 <- rbind(all_same_start[, c(1,4:6)], all_same_end[, c(1,4:6)])
all_as2 <- unique(all_as2)
dim(all_as2)
#[1] 100776      4

all_as2$sj1 <- do.call(rbind, strsplit(all_as2$as2, "_"))[,1]
all_as2$sj2 <- do.call(rbind, strsplit(all_as2$as2, "_"))[,2]

chao1 <- mcmapply(function(i){
  print(paste(i, "of", nrow(all_as2)))
  sj1_in <- grep(all_as2[i, "sj1"], names(AFE_index))
  sj2_in <- grep(all_as2[i, "sj2"], names(AFE_index))
  if(length(sj1_in) != 1 | length(sj2_in) != 1){
    res <- NA
  }else{
    if(sj1_in != sj2_in){
      res <- NA
      }else{
        res <- sj1_in
      }
  }
  return(res)
  }, 1:nrow(all_as2), mc.cores = 30)

sum(!is.na(chao1))
# [1] 2873
afe_as2 <- all_as2[which(!is.na(chao1)), ]
dim(afe_as2)
#[1] 2873     6

afe_as2$AFE_index <- names(AFE_index[chao1[!is.na(chao1)]])

afe_gene <- mcmapply(function(y){
  return(unique(y$gene_id))
}, AFE_index[chao1[!is.na(chao1)]], mc.cores = 10)

afe_gene_name <- mcmapply(function(y){
  return(unique(y$gene_name))
}, AFE_index[chao1[!is.na(chao1)]], mc.cores = 10)

afe_gene_biotype <- mcmapply(function(y){
  return(unique(y$gene_biotype))
}, AFE_index[chao1[!is.na(chao1)]], mc.cores = 10)

afe_strand <- mcmapply(function(y){
  return(unique(y$strand))
}, AFE_index[chao1[!is.na(chao1)]], mc.cores = 10)


AFE <- data.frame(afe_as2, 
				  gene_id = as.character(afe_gene),
				  gene_name = as.character(afe_gene_name),
				  strand = as.character(afe_strand),
				  biotype = as.character(afe_gene_biotype), stringsAsFactors = F)
dim(AFE)
#[1] 2873   11

save(AFE, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/AFE/AFE_type_new.RData")



#### 14. ALE ====================================================================================================================================



load(file = "/mnt/data5/BGI/UCB/tangchao/DSU2/Homo_sapiens.GRCh38.87.gtf_ALE_index.RData")
names(ALE_index) <- lapply(ALE_index, function(x) paste(x$ALSJ[order(x$ALSJ)], collapse="_", sep = ""))

chao2 <- mcmapply(function(i){
  print(paste(i, "of", nrow(all_as2)))
  sj1_in <- grep(all_as2[i, "sj1"], names(ALE_index))
  sj2_in <- grep(all_as2[i, "sj2"], names(ALE_index))
  if(length(sj1_in) != 1 | length(sj2_in) != 1){
    res <- NA
  }else{
    if(sj1_in != sj2_in){
      res <- NA
      }else{
        res <- sj1_in
      }
  }
  return(res)
  }, 1:nrow(all_as2), mc.cores = 30)

sum(!is.na(chao2))
# [1] 1390
ale_as2 <- all_as2[which(!is.na(chao2)), ]
dim(ale_as2)
#[1] 1390     6

ale_as2$ALE_index <- names(ALE_index[chao2[!is.na(chao2)]])

ale_gene <- mcmapply(function(y){
  return(unique(y$gene_id))
}, ALE_index[chao2[!is.na(chao2)]], mc.cores = 10)

ale_gene_name <- mcmapply(function(y){
  return(unique(y$gene_name))
}, ALE_index[chao2[!is.na(chao2)]], mc.cores = 10)

ale_gene_biotype <- mcmapply(function(y){
  return(unique(y$gene_biotype))
}, ALE_index[chao2[!is.na(chao2)]], mc.cores = 10)

ale_strand <- mcmapply(function(y){
  return(unique(y$strand))
}, ALE_index[chao2[!is.na(chao2)]], mc.cores = 10)



ALE <- data.frame(ale_as2, 
				          gene_id = as.character(ale_gene),
				          gene_name = as.character(ale_gene_name),
				          strand = as.character(ale_strand),
				          biotype = as.character(ale_gene_biotype), stringsAsFactors = F)

save(ALE, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/ALE/ALE_type_new.RData")



#### 15. PSI of SE ==============================================================================================================================



load("/mnt/data5/BGI/UCB/tangchao/DSU/RData/psi_list_left_right_table.RData")
sum(!SE$sj12 %in% row.names(psi_sj_same_start_table))
#[1] 0
sum(!SE$sj23 %in% row.names(psi_sj_same_end_table))
#[1] 0

# test <- vector()
# for (i in 1:nrow(PSI_SE2)) {
#   if(is.na(rowSums(PSI_SE2[i,]))){
#     test[i] <- NaN
#   }else{
#     if(sum(PSI_SE2[i,] == 0) == 1){
#       test[i] <- max(PSI_SE2[i,])
#     }else{
#       test[i] <- rowMeans(PSI_SE2[i,])
#     }
#   }
# }

# test2 <- list()
# for(j in 1:nrow(SE)){
#   print(paste(j, "of", nrow(SE)))
#   PSI_SEi <- data.frame(psi12 = as.numeric(psi_sj_same_start_table[SE$sj12[j],]), 
#                         psi23 = as.numeric(psi_sj_same_end_table[SE$sj23[j],]), 
#                         row.names = colnames(psi_sj_same_start_table))
#   test <- vector()
#   for (i in 1:nrow(PSI_SEi)) {
#     if(is.na(rowSums(PSI_SEi[i,]))){
#       test[i] <- NaN
#     }else{
#       if(sum(PSI_SEi[i,] == 0) == 1){
#         test[i] <- max(PSI_SEi[i,])
#       }else{
#         test[i] <- rowMeans(PSI_SEi[i,])
#       }
#     }
#   }
#   test2[[j]] <- test
# }


psi_mat <- mcmapply(function(j){
  print(paste(j, "of", nrow(SE)))
  PSI_SEi <- data.frame(psi12 = as.numeric(psi_sj_same_start_table[SE$sj12[j],]), 
                        psi23 = as.numeric(psi_sj_same_end_table[SE$sj23[j],]), 
                        row.names = colnames(psi_sj_same_start_table))
  test <- vector()
  for (i in 1:nrow(PSI_SEi)) {
    if(is.na(rowSums(PSI_SEi[i,]))){
      test[i] <- NaN
    }else{
      if(sum(PSI_SEi[i,] == 0) == 1){
        test[i] <- max(PSI_SEi[i,])
      }else{
        test[i] <- rowMeans(PSI_SEi[i,])
      }
    }
  }
  return(test)
}, 1:nrow(SE), mc.cores = 10)

dim(psi_mat)
#[1]  3574 17670
row.names(psi_mat) <- colnames(psi_sj_same_start_table)
SE_psi <- cbind(SE, t(psi_mat))

save(SE_psi, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/SE/SE_psi_new.RData")



#### 16. PSI of A3SS ============================================================================================================================



A3SS_psi <- cbind(A3SS, psi_sj_same_start_table[A3SS$sj12,])
dim(A3SS_psi)
#[1] 3050 3588
save(A3SS_psi, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/A3SS/A3SS_psi_new.RData")



#### 17. PSI of A5SS ============================================================================================================================



A5SS_psi <- cbind(A5SS, psi_sj_same_end_table[A5SS$sj23,])
dim(A5SS_psi)
#[1] 3036 3588
save(A5SS_psi, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/A5SS/A5SS_psi_new.RData")



#### 18. PSI of AFE =============================================================================================================================



afe_sj <- list()
for(i in 1:nrow(AFE)){
  if(AFE[i, "strand"] == "-"){
    afe_sj[[i]] <- psi_sj_same_start_table[AFE[i, "sj1"], ]
  }else{
    afe_sj[[i]] <- psi_sj_same_end_table[AFE[i, "sj2"], ]
  }
}

psi <- do.call(rbind, afe_sj)
AFE_psi <- cbind(AFE, psi)
save(AFE_psi, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/AFE/AFE_psi_new.RData")



#### 19. PSI of ALE =============================================================================================================================



ale_sj <- list()
for(i in 1:nrow(ALE)){
  if(ALE[i, "strand"] == "-"){
    ale_sj[[i]] <- psi_sj_same_end_table[ALE[i, "sj2"], ]
  }else{
    ale_sj[[i]] <- psi_sj_same_start_table[ALE[i, "sj1"], ]
  }
}

psi <- do.call(rbind, ale_sj)
ALE_psi <- cbind(ALE, psi)
save(ALE_psi, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/ALE/ALE_psi_new.RData")



#### 20. PSI of MXE =============================================================================================================================



mxe_psi <- list()
for(i in 1:nrow(MXE)){
  mxe_psi[[i]] <- rbind(psi_sj_same_start_table[MXE[i, "sj12"], ], psi_sj_same_end_table[MXE[i, "sj13.1"], ])
}

MXE_psi <- lapply(mxe_psi, function(x){
  tmp <- vector()
  for(i in 1:ncol(x)){
    tmp[i] <- mean(x[,i])
    if(sum(is.na(x[,i])) > 0){ tmp[i] <- NaN }else{
          if(sum(x[,i] == 0) == 1){ tmp[i] <- max(x[,i]) }
          if(sum(x[,i] == 1) == 1){ tmp[i] <- min(x[,i]) }
    }
  }
  return(tmp)
})
MXE_psi <- do.call(rbind, MXE_psi)
colnames(MXE_psi) <- colnames(psi_sj_same_start_table)

MXE_psi <- cbind(MXE, MXE_psi)
MXE_psi <- MXE_psi[!duplicated(paste(MXE_psi$as2,MXE_psi$as2.1,sep="_")),]
save(MXE_psi, file = "/mnt/data5/BGI/UCB/tangchao/DSU2/MXE/MXE_psi_new.RData")



table(SE$SE_type)
# exon-skipping_exactly   novel_exon-skipping                 Other
#                 10372                  3398                  3900
dim(A3SS)
#[1] 3050   14
dim(A5SS)
#[1] 3036   14
dim(AFE)
#[1] 2873   11
dim(ALE)
#[1] 1390   11
table(rowSums(is.na(MXE[,19:20])))
#  0   1   2
#451  92  63

# 10372 + 3398 + 3900 + 3050 + 3036 + 2873 + 1390 + 116
# [1] 28625











