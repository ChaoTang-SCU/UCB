load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/RData/sj11.RData")

gtf <- read.table("/mnt/data5/BGI/UCB/tangchao/Homo_sapiens.GRCh38.87.gtf", header = FALSE, sep = "\t", stringsAsFactors = F)

all_gtf <- gtf[, c(1,4,5,7)]
colnames(all_gtf) <- c("chr", "start", "end", "strand")
all_gtf$gene_id <- substr(gtf$V9,9,23)

sum(!duplicated(all_gtf))
[1] 1013331
dim(all_gtf)
[1] 2575494       5

all_gtf <- all_gtf[!duplicated(all_gtf), ]
row.names(all_gtf) <- 1:nrow(all_gtf)

all_gtf[all_gtf$strand=="+",]$strand <- 1
all_gtf[all_gtf$strand=="-",]$strand <- 2


library(parallel)
SJ_start_exon <- mclapply(which(sj$annotation == 0 | sj$annotation == 1), function(i){
  if(sum(all_gtf$chr == sj[i, "chr"] & all_gtf$strand == sj[i, "strand"] & all_gtf$end == sj[i, "start"]-1) > 0){
    paste(unique(all_gtf[which(all_gtf$chr == sj[i, "chr"] & all_gtf$strand == sj[i, "strand"] & all_gtf$end == sj[i, "start"]-1),]$gene_id),collapse="|")
  }
}, mc.cores = 40)

SJ_end_exon <- mclapply(which(sj$annotation == 0 | sj$annotation == 1), function(i){
  if(sum(all_gtf$chr == sj[i, "chr"] & all_gtf$strand == sj[i, "strand"] & all_gtf$start == sj[i, "end"]+1) > 0){
    paste(unique(all_gtf[which(all_gtf$chr == sj[i, "chr"] & all_gtf$strand == sj[i, "strand"] & all_gtf$start == sj[i, "end"]+1),]$gene_id),collapse="|")
  }
}, mc.cores = 40)


SJ_start_exon_vector <- unlist(lapply(SJ_start_exon,function(x) if(is.null(x)) {x <- NA} else {x <- x}))
SJ_end_exon_vector <- unlist(lapply(SJ_end_exon,function(x) if(is.null(x)) {x <- NA} else {x <- x}))

sum(is.na(SJ_start_exon_vector))
[1] 633016

sum(is.na(SJ_end_exon_vector))
[1] 633094


sj$SJ_start_exon <- SJ_start_exon_vector
sj$SJ_end_exon <- SJ_end_exon_vector

table(nchar(SJ_start_exon_vector))
table(nchar(SJ_end_exon_vector))


gene_gtf <- gtf[gtf$V3 == "gene", c(1,4,5,7)]
colnames(gene_gtf) <- c("chr", "start", "end", "strand")
gene_gtf$gene_id <- substr(gtf[gtf$V3 == "gene", ]$V9,9,23)

row.names(gene_gtf) <- 1:nrow(gene_gtf)



SJ_start_exon_novel <- mclapply(1:nrow(sj), function(i){
	if(sum(gene_gtf$chr == sj[i, "chr"] & gene_gtf$start < sj[i, "start"] & gene_gtf$end > sj[i, "start"]) > 0){
		paste(unique(gene_gtf[which(gene_gtf$chr == sj[i, "chr"] & gene_gtf$start < sj[i, "start"] & gene_gtf$end > sj[i, "start"]),]$gene_id),collapse="|")
		
	}
}, mc.cores = 40)


SJ_end_exon_novel <- mclapply(1:nrow(sj), function(i){
	if(sum(gene_gtf$chr == sj[i, "chr"] & gene_gtf$start < sj[i, "end"] & gene_gtf$end > sj[i, "end"]) > 0){
		paste(unique(gene_gtf[which(gene_gtf$chr == sj[i, "chr"] & gene_gtf$start < sj[i, "end"] & gene_gtf$end > sj[i, "end"]),]$gene_id),collapse="|")
		
	}
}, mc.cores = 40)


SJ_start_exon_novel_vector <- unlist(lapply(SJ_start_exon_novel,function(x) if(is.null(x)) {x <- NA} else {x <- x}))
SJ_end_exon_novel_vector <- unlist(lapply(SJ_end_exon_novel,function(x) if(is.null(x)) {x <- NA} else {x <- x}))

sum(is.na(SJ_start_exon_novel_vector))
[1] 136097
sum(is.na(SJ_end_exon_novel_vector))
[1] 135276

sj$SJ_start_exon_novel <- SJ_start_exon_novel_vector
sj$SJ_end_exon_novel <- SJ_end_exon_novel_vector


length(SJ_start_exon_vector[which(sj$annotation==1)])
[1] 176529
sum(is.na(SJ_start_exon_vector[which(sj$annotation==1)]))
[1] 753
dim(sj[sj$annotation == 1 & is.na(sj$SJ_start_exon),c(1:7,42:45)])
[1] 753  11
head(sj[sj$annotation == 1 & is.na(sj$SJ_start_exon),c(1:7,42:45)])

i = 2600
sum(all_gtf$chr == sj[i, "chr"] & all_gtf$strand == sj[i, "strand"] & all_gtf$end == sj[i, "start"]-1) > 0
## Because this annotated SJ's strand are not 1 or 2 (are 0).

for(i in which(sj$annotation == 1 & is.na(sj$SJ_start_exon))){
	if(sum(all_gtf$chr == sj[i, "chr"] & all_gtf$end == sj[i, "start"]-1) > 0){
    sj[i, "SJ_start_exon"] <- paste(unique(all_gtf[which(all_gtf$chr == sj[i, "chr"] & all_gtf$end == sj[i, "start"]-1),]$gene_id),collapse="|")
  		}
  	if(sum(all_gtf$chr == sj[i, "chr"] & all_gtf$start == sj[i, "end"]+1) > 0){
    sj[i, "SJ_end_exon"] <- paste(unique(all_gtf[which(all_gtf$chr == sj[i, "chr"] & all_gtf$start == sj[i, "end"]+1),]$gene_id),collapse="|")
  		}
}

dim(sj[sj$annotation == 1 & is.na(sj$SJ_start_exon),c(1:7,42:45)])
[1]  0 11
save(sj, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj12.RData")




sj[,c(1:7,42:45)] -> test
dim(test)
[1] 870046     11
table(test$annotation)
     0      1
693517 176529
table(test$strand)
     0      1      2
151776 353137 365133
sum(is.na(test$SJ_start_exon) & is.na(test$SJ_end_exon))
[1] 594054
sum(is.na(test$SJ_start_exon) & !is.na(test$SJ_end_exon))
[1] 38209
sum(!is.na(test$SJ_start_exon) & is.na(test$SJ_end_exon))
[1] 38287
sum(test$annotation == 0 & !is.na(test$SJ_start_exon) & !is.na(test$SJ_end_exon))
[1] 22967

sum(test$annotation == 0 & is.na(test$SJ_start_exon) & is.na(test$SJ_end_exon))
[1] 594054

Known junction       176529
Unannotated junction 594054
Unannotated left     38209
Unannotated right    38287
Unannotated pattern  22967


sum(!is.na(test$SJ_start_exon) & !is.na(test$SJ_end_exon))
[1] 199496
test2 <- test[!is.na(test$SJ_start_exon) & !is.na(test$SJ_end_exon),]
teseyn <- apply(test2,1,function(x){
  !(grepl(x[8],x[9]) | grepl(x[9],x[8]))
  })

sum(teseyn)
[1] 1699
head(test2[teseyn,])
table(test2[teseyn,]$annotation)
   0
1699
trans_splicing 1699


load(file = "/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Cell_sep_PSI.RData")
sj_Bcell <- row.names(Bcell_Sep[Bcell_Sep$adj_pval < 0.01 & Bcell_Sep$Cell_p > .1 & Bcell_Sep$Cell_p > 2*Bcell_Sep$Other_p, ])
sj_CD4 <- row.names(CD4_Sep[CD4_Sep$adj_pval < 0.01 & CD4_Sep$Cell_p > .1 & CD4_Sep$Cell_p > 2*CD4_Sep$Other_p, ])
sj_CD8 <- row.names(CD8_Sep[CD8_Sep$adj_pval < 0.01 & CD8_Sep$Cell_p > .1 & CD8_Sep$Cell_p > 2*CD8_Sep$Other_p, ])
sj_MK <- row.names(MK_Sep[MK_Sep$adj_pval < 0.01 & MK_Sep$Cell_p > .1 & MK_Sep$Cell_p > 2*MK_Sep$Other_p, ])
sj_Mono <- row.names(Mono_Sep[Mono_Sep$adj_pval < 0.01 & Mono_Sep$Cell_p > .1 & Mono_Sep$Cell_p > 2*Mono_Sep$Other_p, ])
sj_NK <- row.names(NK_Sep[NK_Sep$adj_pval < 0.01 & NK_Sep$Cell_p > .1 & NK_Sep$Cell_p > 2*NK_Sep$Other_p, ])

library(gdata)
marker <- read.xls('/mnt/data7/huangfei/single-cell/UCB/PCA/Markers.xlsx',sheet = 2,header=T)
read.xls('/mnt/data7/huangfei/single-cell/UCB/PCA/Markers.xlsx',sheet = 1,header=T)
head(sj[sj$sj %in% sj_Bcell,42:45])

gtf <- read.table("/mnt/data5/BGI/UCB/tangchao/Homo_sapiens.GRCh38.87.gtf", header = FALSE, sep = "\t", stringsAsFactors = F)
gtf <- gtf[gtf$V3=="gene",]
test <- as.data.frame(do.call(rbind, strsplit(as.character(gtf$V9), split = "[; ]")))

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
           gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})), stringsAsFactors=F)

gene_gtf <- cbind(gtf[,1:8],gtf_infor)

sj[,c(1:7,42:45)] -> test
dim(test)

dim(test[test$sj %in% sj_Bcell,])
[1] 1003   11
sj_Bcell_info <- test[test$sj %in% sj_Bcell,]
table(sj_Bcell_info$annotation)
  0   1
 45 958
table(nchar(sj_Bcell_info$SJ_start_exon))
 15  31
941  34
table(nchar(sj_Bcell_info$SJ_end_exon))
 15  31
938  39
dim(sj_Bcell_info[nchar(sj_Bcell_info$SJ_end_exon)>15 | nchar(sj_Bcell_info$SJ_start_exon)>15,])
[1] 76 11

apply(sj_Bcell_info,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])})
sum(apply(sj_Bcell_info,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])}),na.rm=T)
[1] 958
dim(sj_Bcell_info[which(apply(sj_Bcell_info,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])})),])
[1] 958  11
head(sj_Bcell_info[which(apply(sj_Bcell_info,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])})),1:9],20)

dim(na.omit(sj_Bcell_info[,8:9]))
[1] 970   2

library(reshape2)
sj_Bcell_gene <- unique(unlist(strsplit(c(sj_Bcell_info[which(apply(sj_Bcell_info,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])})),8],sj_Bcell_info[which(apply(sj_Bcell_info,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])})),9]), split="\\|")))

read.xls('/mnt/data7/huangfei/single-cell/UCB/PCA/Markers.xlsx',sheet = 1,header=T)
dim(marker[marker$cluster==1,])
[1] 63  7

Bcell_gene <- gene_gtf[gene_gtf$gene_name %in% marker[marker$cluster==1,]$gene,]$gene_id

length(Bcell_gene)
[1] 63
length(sj_Bcell_gene)
[1] 216
length(intersect(Bcell_gene,sj_Bcell_gene))
[1] 55


sj_CD4_info <- test[test$sj %in% sj_CD4,]
table(sj_CD4_info$annotation)
 0  1
 3 81
sum(apply(sj_CD4_info,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])}),na.rm=T)
[1] 81
sj_CD4_gene <- unique(unlist(strsplit(c(sj_CD4_info[which(apply(sj_CD4_info,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])})),8],sj_CD4_info[which(apply(sj_CD4_info,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])})),9]), split="\\|")))
dim(marker[marker$cluster==0,])
[1] 2 7
CD4_gene <- gene_gtf[gene_gtf$gene_name %in% marker[marker$cluster==0,]$gene,]$gene_id

length(CD4_gene)
[1] 2
length(sj_CD4_gene)
[1] 22
length(intersect(CD4_gene,sj_CD4_gene))
[1] 2



sj_CD8_info <- test[test$sj %in% sj_CD8,]
table(sj_CD8_info$annotation)
  0   1
  2 140
sum(apply(sj_CD8_info,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])}),na.rm=T)
[1] 140
sj_CD8_gene <- unique(unlist(strsplit(c(sj_CD8_info[which(apply(sj_CD8_info,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])})),8],sj_CD8_info[which(apply(sj_CD8_info,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])})),9]), split="\\|")))
dim(marker[marker$cluster==2,])
[1] 3 7
CD8_gene <- gene_gtf[gene_gtf$gene_name %in% marker[marker$cluster==2,]$gene,]$gene_id

length(CD8_gene)
[1] 3
length(sj_CD8_gene)
[1] 44
length(intersect(CD8_gene,sj_CD8_gene))
[1] 3



sj_NK_info <- test[test$sj %in% sj_NK,]
table(sj_NK_info$annotation)
  0   1
 26 921
sum(apply(sj_NK_info,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])}),na.rm=T)
[1] 923
sj_NK_gene <- unique(unlist(strsplit(c(sj_NK_info[which(apply(sj_NK_info,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])})),8],sj_NK_info[which(apply(sj_NK_info,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])})),9]), split="\\|")))
dim(marker[marker$cluster==3,])
[1] 42 7
NK_gene <- gene_gtf[gene_gtf$gene_name %in% marker[marker$cluster==3,]$gene,]$gene_id

length(NK_gene)
[1] 42
length(sj_NK_gene)
[1] 256
length(intersect(NK_gene,sj_NK_gene))
[1] 40



sj_Mono_info <- test[test$sj %in% sj_Mono,]
table(sj_Mono_info$annotation)
   0    1
  98 1549
sum(apply(sj_Mono_info,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])}),na.rm=T)
[1] 1550
sj_Mono_gene <- unique(unlist(strsplit(c(sj_Mono_info[which(apply(sj_Mono_info,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])})),8],sj_Mono_info[which(apply(sj_Mono_info,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])})),9]), split="\\|")))
dim(marker[marker$cluster==4,])
[1] 116 7
Mono_gene <- gene_gtf[gene_gtf$gene_name %in% marker[marker$cluster==4,]$gene,]$gene_id

length(Mono_gene)
[1] 116
length(sj_Mono_gene)
[1] 512
length(intersect(Mono_gene,sj_Mono_gene))
[1] 98



sj_MK_info <- test[test$sj %in% sj_MK,]
table(sj_MK_info$annotation)
   0    1
  54 2640
sum(apply(sj_MK_info,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])}),na.rm=T)
[1] 2642
sj_MK_gene <- unique(unlist(strsplit(c(sj_MK_info[which(apply(sj_MK_info,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])})),8],sj_MK_info[which(apply(sj_MK_info,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])})),9]), split="\\|")))
dim(marker[marker$cluster==5,])
[1] 72 7
MK_gene <- gene_gtf[gene_gtf$gene_name %in% marker[marker$cluster==5,]$gene,]$gene_id

length(MK_gene)
[1] 71
length(sj_MK_gene)
[1] 599
length(intersect(MK_gene,sj_MK_gene))
[1] 63

save(Bcell_gene,sj_Bcell_gene,CD4_gene,sj_CD4_gene,CD8_gene,sj_CD8_gene,NK_gene,sj_NK_gene,Mono_gene,sj_Mono_gene,MK_gene,sj_MK_gene,
  file = "/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/Cell_sepecific_gene_and_sj_host_gene.RData")




### Host gene of cell type specific unannotated SJ
load(file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/sj12.RData")
load("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/cell_type_specific_unannotated_SJ.RData")
sj[,c(1:7,42:45)] -> test
dim(test)
[1] 870046     11
dim(novel_sj_maxreads)
[1] 210   3
sum(test$sj %in% novel_sj_maxreads$SJ)
[1] 191
unannotated <- test[test$sj %in% novel_sj_maxreads$SJ,]

unique(unlist(strsplit(c(unannotated[which(apply(unannotated,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])})),8],unannotated[which(apply(unannotated,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])})),9]), split="\\|")))

cell_sep_unannotated_SJ_gene <- unique(unlist(strsplit(c(unannotated[,8],unannotated[,9]), split="\\|")))[-1]

save(cell_sep_unannotated_SJ_gene, file = "/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/cell_sep_unannotated_SJ_gene.RData")




### Host gene of all unannotated SJ
annotated <- test[sj$annotation == 1, ]
unique(unlist(strsplit(c(annotated[which(apply(annotated,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])})),8],annotated[which(apply(annotated,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])})),9]), split="\\|")))
annotated_SJ_gene <- unique(unlist(strsplit(c(annotated[which(apply(annotated,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])})),8],annotated[which(apply(annotated,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])})),9]), split="\\|")))


### Host gene of all unannotated SJ
unannotated <- test[sj$annotation == 0 & sj$novel == "N",]

unique(unlist(strsplit(c(unannotated[which(apply(unannotated,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])})),8],unannotated[which(apply(unannotated,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])})),9]), split="\\|")))

unannotated_SJ_gene <- unique(unlist(strsplit(c(unannotated[which(apply(unannotated,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])})),8],unannotated[which(apply(unannotated,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])})),9]), split="\\|")))


### Host gene of all novel SJ
novel <- test[sj$annotation == 0 & sj$novel == "Y",]

unique(unlist(strsplit(c(novel[which(apply(novel,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])})),8],novel[which(apply(novel,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])})),9]), split="\\|")))

novel_SJ_gene <- unique(unlist(strsplit(c(novel[which(apply(novel,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])})),8],novel[which(apply(novel,1,function(x){grepl(x[8],x[9]) | grepl(x[9],x[8])})),9]), split="\\|")))

save(annotated_SJ_gene, unannotated_SJ_gene, novel_SJ_gene, file = "/mnt/data5/BGI/UCB/tangchao/Cell_spe_SJ/all_annotated_unannotated_and_novel_SJ_gene.RData")


library(org.Hs.eg.db)##### data base target to human-geneIDs
library(AnnotationDbi)
library(DOSE)
library(GO.db)
library(topGO)
library(GSEABase)
library(clusterProfiler)

background <- as.character(read.table("/mnt/data5/BGI/UCB/tangchao/AS_stat/Background_for_GO.txt", header = F, stringsAsFactors = F)[,1])

all_Novel_BP_ego <- enrichGO(gene = novel_SJ_gene,
                            universe = background,
                            OrgDb  = org.Hs.eg.db,
                            keyType= 'ENSEMBL',
                            ont  = "BP", 
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05)

interested_gene <- c(novel_SJ_gene[novel_SJ_gene %in% all_Novel_BP_ego@geneSets$`GO:0050852`],novel_SJ_gene[novel_SJ_gene %in% all_Novel_BP_ego@geneSets$`GO:0050851`])

dim(novel[novel$SJ_start_exon %in% interested_gene & novel$SJ_end_exon %in% interested_gene,])
[1] 112  11

interested_sj <- novel[novel$SJ_start_exon %in% interested_gene & novel$SJ_end_exon %in% interested_gene,]$sj

load("/mnt/data5/BGI/UCB/tangchao/data/SJ/SJ_merged_raw_te(no_NA).RData")

te_sub <- te_table[interested_sj,]

tabn <- vector()
for(i in 1:nrow(te_sub)){
  test <- data.frame(R = as.character(te_sub[i,]), I = substr(colnames(te_sub),1,4), stringsAsFactors = F)
  test <- test[test$R > 0, ]
       tabn[i] <- length(unique(test$I))
}
table(tabn)
tabn
 1  2  3
97 12  3

dim(te_sub2[tabn >= 2,])
#[1]  15 3574

#### sashimi plot ====================================================================================================================

novel_sj_reads <- rowSums(te_sub[tabn >= 2,])
novel_sj_maxreads <- as.data.frame(apply(te_sub[tabn >= 2,], 1, function(x) colnames(te_sub)[which.max(x)]))
colnames(novel_sj_maxreads) <- "Cell"
novel_sj_maxreads$SJ <- row.names(te_sub[tabn >= 2,])
save(novel_sj_maxreads, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/interested_novel_SJ_and_max_cell.RData")

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


for(i in 1:nrow(novel_sj_maxreads)){
  print(i)
  
    sashimi(junction = novel_sj_maxreads[i,]$SJ, cell = novel_sj_maxreads[i,]$Cell, outDir = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Figure/interested_novel_SJ_sashimi/", ann_height = 6)
  
}


novel_sj_reads <- rowSums(te_sub[tabn == 1,])
novel_sj_maxreads <- as.data.frame(apply(te_sub[tabn == 1,], 1, function(x) colnames(te_sub)[which.max(x)]))
colnames(novel_sj_maxreads) <- "Cell"
novel_sj_maxreads$SJ <- row.names(te_sub[tabn == 1,])
save(novel_sj_maxreads, file = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/interested_novel_SJ_and_max_cell_only_in_one_people.RData")


for(i in 1:nrow(novel_sj_maxreads)){
  print(i)
  sashimi(junction = novel_sj_maxreads[i,]$SJ, cell = novel_sj_maxreads[i,]$Cell, outDir = "/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Figure/interested_novel_SJ_sashimi/onepeople/", ann_height = 6)  
}




#### gene biotype
sort(table(gene_gtf[gene_gtf$gene_id %in% annotated_SJ_gene, ]$gene_biotype), decreasing = T)/length(annotated_SJ_gene)
sort(table(gene_gtf[gene_gtf$gene_id %in% unannotated_SJ_gene, ]$gene_biotype), decreasing = T)/length(unannotated_SJ_gene)
sort(table(gene_gtf[gene_gtf$gene_id %in% novel_SJ_gene, ]$gene_biotype), decreasing = T)/length(novel_SJ_gene)

bio_novel <- as.data.frame(t(t(sort(table(gene_gtf[gene_gtf$gene_id %in% novel_SJ_gene, ]$gene_biotype), decreasing = T)/length(novel_SJ_gene))))
bio_annotated <- as.data.frame(t(t(sort(table(gene_gtf[gene_gtf$gene_id %in% annotated_SJ_gene, ]$gene_biotype), decreasing = T)/length(annotated_SJ_gene))))
bio_unannotated <- as.data.frame(t(t(sort(table(gene_gtf[gene_gtf$gene_id %in% unannotated_SJ_gene, ]$gene_biotype), decreasing = T)/length(unannotated_SJ_gene))))

bio_gtf <- as.data.frame(t(t(sort(table(gene_gtf$gene_biotype), decreasing = T)/nrow(gene_gtf))))

bio <- merge(merge(merge(bio_annotated, bio_unannotated, by = "Var1"),bio_novel,by = "Var1"),bio_gtf,by = "Var1")
bio <- bio[,c(1,3,5,7,9)]
colnames(bio) <- c("TYPE","Annotated","Unannotated","Novel","GRCh38.87")
bio <- bio[order(bio$Annotated, decreasing = T),]

bio <- bio[-1,]
pdf("/mnt/data5/BGI/UCB/tangchao/novel_SS/All_info_table/Figure/SJ_biotype.pdf")
par(mar = c(16,4,2,1))
plot(x = 1:nrow(bio), y = bio$GRCh38.87, type = "l", xaxt="n", xlab = "", ylab = "")
lines(x = 1:nrow(bio), y = bio$Annotated, col = 2)
lines(x = 1:nrow(bio), y = bio$Unannotated, col = 3)
lines(x = 1:nrow(bio), y = bio$Novel, col = 4)
axis(side = 1,at = 1:13, labels = bio$TYPE, tick = TRUE, las=2)
legend("topright", col = 1:4, legend = c("GRCh38.87","Annotated","Unannotated","Novel"), bty = "n",
       lty = 1, cex = .5)

dev.off()









