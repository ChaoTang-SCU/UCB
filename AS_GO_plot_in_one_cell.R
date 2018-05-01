library(org.Hs.eg.db)##### data base target to human-geneIDs
library(AnnotationDbi)
library(DOSE)
library(GO.db)
library(topGO)
library(GSEABase)
library(clusterProfiler)

background <- as.character(read.table("/Users/tangchao/CloudStation/tangchao/project/UCB/AS_GO/Background_for_GO.txt", header = F, stringsAsFactors = F)[,1])


#### SE ===================================================================================================
load("/Users/tangchao/CloudStation/tangchao/project/UCB/all_cells_SE_type_and_PSI.RData")
load("/Users/tangchao/CloudStation/tangchao/project/UCB/all_cells_SE_GTF.RData")
rowSums(!is.na(all_cells_SE_type_and_PSI[,12:ncol(all_cells_SE_type_and_PSI)])) -> se_psi_rowsum
table(all_cells_SE_type_and_PSI[se_psi_rowsum>1,]$AStype)

all_cells_SE_type_and_PSI[se_psi_rowsum == 1 & all_cells_SE_type_and_PSI$AStype == "exon-skipping_exactly",]$loci -> se_loci

test <- as.data.frame(do.call(rbind, strsplit(SE_GTF[SE_GTF$loci %in% se_loci, ]$V9, split = "[; ]")))

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
                        transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
                        gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
                        transcript_name = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_name")})),
                        gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})),
                        transcript_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_biotype")})), stringsAsFactors=F)

dim(gtf_infor)
gtf_infor <- gtf_infor[!duplicated(gtf_infor),]
dim(gtf_infor)
# [1] 6591    6
length(unique(gtf_infor$gene_id))
# [1] 2864
length(unique(gtf_infor$transcript_id))
# [1] 6591

table(gtf_infor$gene_biotype)

gtf_infor -> SE_gene_info
unique(as.character(gtf_infor$gene_id)) -> SE_gene_list

#### GO

SE_BP_ego <- enrichGO(gene = SE_gene_list,
                      universe = background,
                      OrgDb  = org.Hs.eg.db,
                      keytype= 'ENSEMBL',
                      ont  = "BP", 
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/AS_GO/one cell/SE_GO.pdf",height = 10,width = 15)
print(dotplot(SE_BP_ego,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(SE_BP_ego,firstSigNodes = 15, useInfo = "all"))
enrichMap(SE_BP_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


#### A3SS ===========================================================================================================================

load("/Users/tangchao/CloudStation/tangchao/project/UCB/Parsed_A3SS_PSI_GTF.RData")
dim(A3SS_type_and_PSI)
# [1] 5874 3584
rowSums(!is.na(A3SS_type_and_PSI[,11:ncol(A3SS_type_and_PSI)])) -> a3ss_psi_rowsum
sum(a3ss_psi_rowsum>1)
# [1] 3179

length(a3ss_psi_rowsum)
sum(a3ss_psi_rowsum==1)
A3SS_type_and_PSI$name[a3ss_psi_rowsum==1] -> a3ss_loci

test <- as.data.frame(do.call(rbind, strsplit(na.omit(A3SS_GTF[A3SS_GTF$name %in% a3ss_loci,]$V9), split = "[; ]")))

gtf_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),             
                        transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
                        gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
                        transcript_name = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_name")})),
                        gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})),
                        transcript_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_biotype")})), stringsAsFactors=F)

gtf_infor <- gtf_infor[!duplicated(gtf_infor),]
dim(gtf_infor)
# [1] 4044    6
length(unique(gtf_infor$gene_id))
# [1] 2183
length(unique(gtf_infor$transcript_id))
# [1] 4044

unique(gtf_infor$gene_id) -> A3SS_gene_list

#### GO

A3SS_BP_ego <- enrichGO(gene = A3SS_gene_list,
                        universe = background,
                        OrgDb  = org.Hs.eg.db,
                        keytype= 'ENSEMBL',
                        ont  = "BP", 
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/AS_GO/one cell/A3SS_GO.pdf",height = 10,width = 15)
print(dotplot(A3SS_BP_ego,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(A3SS_BP_ego,firstSigNodes = 15, useInfo = "all"))
enrichMap(A3SS_BP_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


#### A5SS ===========================================================================================================================

load("/Users/tangchao/CloudStation/tangchao/project/UCB/Parsed_A5SS_PSI_GTF.RData")
dim(A5SS_type_and_PSI)
# [1] 3908 3584
rowSums(!is.na(A5SS_type_and_PSI[,11:ncol(A5SS_type_and_PSI)])) -> a5ss_psi_rowsum
sum(a5ss_psi_rowsum>1)
# [1] 2256

length(a5ss_psi_rowsum)
sum(a5ss_psi_rowsum==1)
A5SS_type_and_PSI$name[a5ss_psi_rowsum==1] -> a5ss_loci

test <- as.data.frame(do.call(rbind, strsplit(na.omit(A5SS_GTF[A5SS_GTF$name %in% a5ss_loci,]$V9), split = "[; ]")))

gtf_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),             
                        transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
                        gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
                        transcript_name = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_name")})),
                        gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})),
                        transcript_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_biotype")})), stringsAsFactors=F)

gtf_infor <- gtf_infor[!duplicated(gtf_infor),]
dim(gtf_infor)
# [1] 2588    6
length(unique(gtf_infor$gene_id))
# [1] 1476
length(unique(gtf_infor$transcript_id))
# [1] 2588

unique(gtf_infor$gene_id) -> A5SS_gene_list

#### GO

A5SS_BP_ego <- enrichGO(gene = A5SS_gene_list,
                        universe = background,
                        OrgDb  = org.Hs.eg.db,
                        keytype= 'ENSEMBL',
                        ont  = "BP", 
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/AS_GO/one cell/A5SS_GO.pdf",height = 10,width = 15)
print(dotplot(A5SS_BP_ego,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(A5SS_BP_ego,firstSigNodes = 5, useInfo = "all"))
enrichMap(A5SS_BP_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


#### AFE ===========================================================================================================================

load("/Users/tangchao/CloudStation/tangchao/project/UCB/all_cells_AFE_PSI_GTF.RData")
dim(all_cells_AFE_PSI)
# [1] 2289 3584
rowSums(!is.na(all_cells_AFE_PSI[,11:ncol(all_cells_AFE_PSI)])) -> afe_psi_rowsum
sum(afe_psi_rowsum==1)
# [1] 928

all_cells_AFE_PSI$name[afe_psi_rowsum==1] -> afe_loci
length(afe_loci)
# [1] 928


test <- as.data.frame(do.call(rbind, strsplit(na.omit(AFE_GTF[AFE_GTF$name %in% afe_loci,]$V9), split = "[; ]")))

gtf_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),             
                        transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
                        gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
                        transcript_name = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_name")})),
                        gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})),
                        transcript_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_biotype")})), stringsAsFactors=F)

gtf_infor <- gtf_infor[!duplicated(gtf_infor),]
dim(gtf_infor)
# [1] 1128    6
length(unique(gtf_infor$gene_id))
# [1] 871
length(unique(gtf_infor$transcript_id))
# [1] 1128

unique(gtf_infor$gene_id) -> AFE_gene_list

#### GO

AFE_BP_ego <- enrichGO(gene = AFE_gene_list,
                       universe = background,
                       OrgDb  = org.Hs.eg.db,
                       keytype= 'ENSEMBL',
                       ont  = "BP", 
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/AS_GO/one cell/AFE_GO.pdf",height = 10,width = 15)
print(dotplot(AFE_BP_ego,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(AFE_BP_ego,firstSigNodes = 15, useInfo = "all"))
enrichMap(AFE_BP_ego,n = 30, vertex.label.font = 1,cex = .1, font.size = .1)
dev.off()


 #### ALE ===========================================================================================================================

load("/Users/tangchao/CloudStation/tangchao/project/UCB/all_cells_ALE_PSI_GTF.RData")
dim(all_cells_ALE_PSI)
# [1] 1211 3584
rowSums(!is.na(all_cells_ALE_PSI[,11:ncol(all_cells_ALE_PSI)])) -> ale_psi_rowsum
sum(ale_psi_rowsum==1)
# [1] 587

all_cells_ALE_PSI$name[ale_psi_rowsum==1] -> ale_loci
length(ale_loci)
# [1] 587


test <- as.data.frame(do.call(rbind, strsplit(na.omit(ALE_GTF[ALE_GTF$name %in% ale_loci,]$V9), split = "[; ]")))

gtf_infor <- data.frame(gene_id = as.vector(apply(test,1,valueFind)),             
                        transcript_id = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_id")})),
                        gene_name = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_name")})),
                        transcript_name = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_name")})),
                        gene_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "gene_biotype")})),
                        transcript_biotype = as.vector(apply(test,1,function(x){valueFind(x,type = "transcript_biotype")})), stringsAsFactors=F)

gtf_infor <- gtf_infor[!duplicated(gtf_infor),]
dim(gtf_infor)
# [1] 846    6
length(unique(gtf_infor$gene_id))
# [1] 545
length(unique(gtf_infor$transcript_id))
# [1] 846

unique(gtf_infor$gene_id) -> ALE_gene_list

#### GO

ALE_BP_ego <- enrichGO(gene = ALE_gene_list,
                       universe = background,
                       OrgDb  = org.Hs.eg.db,
                       keytype= 'ENSEMBL',
                       ont  = "BP", 
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/AS_GO/one cell/ALE_GO.pdf", height = 10,width = 15)
print(dotplot(ALE_BP_ego,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(ALE_BP_ego,firstSigNodes = 15, useInfo = "all"))
enrichMap(ALE_BP_ego,n = 30, vertex.label.font = 1, cex = .1, font.size = .1)
dev.off()


#### circro plot ===========================================================================================================

library(circlize)
library(RColorBrewer)
SE_BP_ego@result -> SE
A3SS_BP_ego@result -> A3SS
A5SS_BP_ego@result -> A5SS
AFE_BP_ego@result -> AFE
ALE_BP_ego@result -> ALE

lapply(list(SE,A3SS,A5SS,AFE,ALE), dim)

lapply(list(SE_gene_list,A3SS_gene_list,A5SS_gene_list,AFE_gene_list,ALE_gene_list), length)
Reduce(union,list(SE_gene_list,A3SS_gene_list,A5SS_gene_list,AFE_gene_list,ALE_gene_list))

length(unique(c(SE$ID, A3SS$ID, A5SS$ID, AFE$ID, ALE$ID)))
mat <- matrix(0, nrow = 5, ncol = length(unique(c(SE$ID, A3SS$ID, A5SS$ID, AFE$ID, ALE$ID))))
row.names(mat) <- c("SE","A3SS","A5SS","AFE","ALE")
colnames(mat) <- unique(c(SE$ID, A3SS$ID, A5SS$ID, AFE$ID, ALE$ID))

for(i in 1:nrow(SE)){
  mat[1, SE[i,"ID"]] <- SE[i,"Count"]
}

for(i in 1:nrow(A3SS)){
  mat[2, A3SS[i,"ID"]] <- A3SS[i,"Count"]
}

for(i in 1:nrow(A5SS)){
  mat[3, A5SS[i,"ID"]] <- A5SS[i,"Count"]
}

for(i in 1:nrow(AFE)){
  mat[4, AFE[i,"ID"]] <- AFE[i,"Count"]
}

for(i in 1:nrow(ALE)){
  mat[5, ALE[i,"ID"]] <- ALE[i,"Count"]
}



circos.par("canvas.xlim" = c(-2, 2), "canvas.ylim" = c(-2, 2), start.degree = 270, gap.after = c(rep(2, nrow(mat)-1), 10, rep(0, ncol(mat)-1), 10))# 块儿与块儿之间的空格
circos.par(start.degree = 270, gap.after = c(rep(2, nrow(mat)-1), 10, rep(0, ncol(mat)-1), 10))# 块儿与块儿之间的空格
chordDiagram(mat, annotationTrack = "grid", transparency = 0.5, 
             grid.col = c("#FF0000","#00FF00","#0000FF","#FFFF00","#00FFFF",colorRampPalette(brewer.pal(7, "Paired"))(ncol(mat))),
             row.col = c("#FF0000","#00FF00","#0000FF","#FFFF00","#00FFFF"),
             preAllocateTracks = list(track.height = 0.1))
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  if(sector.name %in% get.all.sector.index()[1:5]){
    circos.text(mean(xlim), 0.3*nchar(sector.name), sector.name, niceFacing = T, facing = "clockwise")
  }
}, bg.border = NA)
circos.clear()


library(pheatmap)
pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/AS_GO/one cell/AS_GO_Count_Heatmap.pdf")
pheatmap(mat, show_colnames = F)
dev.off()

pheatmap(mat, silent = T) -> cluster_tab
mat_sort <- mat[,colnames(mat)[cluster_tab$tree_col$order]]

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/AS_GO/one cell/AS_GO_Count_circular.pdf")
circos.par(start.degree = 270, gap.after = c(rep(2, nrow(mat)-1), 10, rep(0, ncol(mat)-1), 10))# 块儿与块儿之间的空格
chordDiagram(mat_sort, annotationTrack = "grid", transparency = 0.5, 
             grid.col = c("#FF0000","#00FF00","#0000FF","#FFFF00","#00FFFF",colorRampPalette(brewer.pal(7, "Paired"))(ncol(mat))),
             row.col = c("#FF0000","#00FF00","#0000FF","#FFFF00","#00FFFF"),
             preAllocateTracks = list(track.height = 0.1))
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  if(sector.name %in% get.all.sector.index()[1:5]){
    circos.text(mean(xlim), 0.3*nchar(sector.name), sector.name, niceFacing = T, facing = "clockwise")
  }
}, bg.border = NA)
circos.clear()
dev.off()


#### 

library(ggsci)
library(scales)
rev(colorRampPalette(pal_material("deep-purple")(10))(max(mat_sort)))
colorRampPalette(pal_material("deep-purple")(10))(max(mat_sort))
show_col(colorRampPalette(pal_material("deep-purple")(10))(max(mat_sort)))

colorRampPalette(c("#FFFFFF","#FF0000"))(100)

show_col(rev(colorRampPalette(c("#FFFFFF","#FF0000"))(max(mat_sort))))

rev(colorRampPalette(pal_material("deep-purple")(10))(max(mat_sort)))

rev(colorRampPalette(c("#FFFFFF","#FF0000"))(max(mat_sort)))

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/AS_GO/one cell/AS_GO_circos.pdf")
circos.par(start.degree = 270, gap.after = c(rep(2, nrow(mat)-1), 10, rep(0, ncol(mat)-1), 10))# 块儿与块儿之间的空格
chordDiagram(mat_sort, annotationTrack = "grid", transparency = 0.5, 
             grid.col = c("#FF0000","#00FF00","#0000FF","#FFFF00","#00FFFF",rep("#FFFFFF", ncol(mat_sort))),
             row.col = c("#FF0000","#00FF00","#0000FF","#FFFF00","#00FFFF"),
             preAllocateTracks = list(track.height = 0.1))
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  if(sector.name %in% get.all.sector.index()[1:5]){
    circos.text(mean(xlim), 0.3*nchar(sector.name), sector.name, niceFacing = T, facing = "clockwise")
  }
}, bg.border = NA)

circos.track(track.index = 2, panel.fun = function(x, y) {
  sector.name = get.cell.meta.data("sector.index")
  if(sector.name %in% colnames(mat_sort)) {
    circos.rect(get.cell.meta.data("xlim")[1], 
                get.cell.meta.data("ylim")[1], 
                get.cell.meta.data("xlim")[2], 
                get.cell.meta.data("ylim")[2],  
                col = rev(colorRampPalette(c("#FFFFFF","#FF0000"))(max(mat_sort)))[as.numeric(mat_sort[1,sector.name])],
                border = NA)
  }
}, bg.border = NA)

circos.track(track.index = 2, panel.fun = function(x, y) {
  sector.name = get.cell.meta.data("sector.index")
  if(sector.name %in% colnames(mat_sort)) {
    circos.rect(get.cell.meta.data("xlim")[1], 
                get.cell.meta.data("ylim")[2], 
                get.cell.meta.data("xlim")[2], 
                get.cell.meta.data("ylim")[2]*2,  
                col = rev(colorRampPalette(c("#FFFFFF","#00FF00"))(max(mat_sort)))[as.numeric(mat_sort[2,sector.name])],
                border = NA)
  }
}, bg.border = NA)

circos.track(track.index = 2, panel.fun = function(x, y) {
  sector.name = get.cell.meta.data("sector.index")
  if(sector.name %in% colnames(mat_sort)) {
    circos.rect(get.cell.meta.data("xlim")[1], 
                get.cell.meta.data("ylim")[2]*2, 
                get.cell.meta.data("xlim")[2], 
                get.cell.meta.data("ylim")[2]*3,  
                col = rev(colorRampPalette(c("#FFFFFF","#0000FF"))(max(mat_sort)))[as.numeric(mat_sort[3,sector.name])],
                border = NA)
  }
}, bg.border = NA)

circos.track(track.index = 2, panel.fun = function(x, y) {
  sector.name = get.cell.meta.data("sector.index")
  if(sector.name %in% colnames(mat_sort)) {
    circos.rect(get.cell.meta.data("xlim")[1], 
                get.cell.meta.data("ylim")[2]*3, 
                get.cell.meta.data("xlim")[2], 
                get.cell.meta.data("ylim")[2]*4,  
                col = rev(colorRampPalette(c("#FFFFFF","#FFFF00"))(max(mat_sort)))[as.numeric(mat_sort[4,sector.name])],
                border = NA)
  }
}, bg.border = NA)

circos.track(track.index = 2, panel.fun = function(x, y) {
  sector.name = get.cell.meta.data("sector.index")
  if(sector.name %in% colnames(mat_sort)) {
    circos.rect(get.cell.meta.data("xlim")[1], 
                get.cell.meta.data("ylim")[2]*4, 
                get.cell.meta.data("xlim")[2], 
                get.cell.meta.data("ylim")[2]*5,  
                col = rev(colorRampPalette(c("#FFFFFF","#00FFFF"))(max(mat_sort)))[as.numeric(mat_sort[5,sector.name])],
                border = NA)
  }
}, bg.border = NA)

circos.clear()
dev.off()



#### heatmap of adj.pvalue ============================================================================================================

matp <- matrix(0, nrow = 5, ncol = length(unique(c(SE$ID, A3SS$ID, A5SS$ID, AFE$ID, ALE$ID))))
row.names(matp) <- c("SE","A3SS","A5SS","AFE","ALE")
colnames(matp) <- unique(c(SE$ID, A3SS$ID, A5SS$ID, AFE$ID, ALE$ID))

for(i in 1:nrow(SE)){
  matp[1, SE[i,"ID"]] <- SE[i,"p.adjust"]
}

for(i in 1:nrow(A3SS)){
  matp[2, A3SS[i,"ID"]] <- A3SS[i,"p.adjust"]
}

for(i in 1:nrow(A5SS)){
  matp[3, A5SS[i,"ID"]] <- A5SS[i,"p.adjust"]
}

for(i in 1:nrow(AFE)){
  matp[4, AFE[i,"ID"]] <- AFE[i,"p.adjust"]
}

for(i in 1:nrow(ALE)){
  matp[5, ALE[i,"ID"]] <- ALE[i,"p.adjust"]
}


-log10(matp) -> logmatp
logmatp[is.infinite(logmatp)] <- -1
pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/AS_GO/one cell/AS_GO_p.adjust_Heatmap.pdf")
pheatmap(logmatp, show_colnames = F)
dev.off()

pheatmap(logmatp, show_colnames = F, cutree_cols = 2)

test2 <- pheatmap(logmatp, silent = T)
cutree(tree = test2$tree_col, 2)
table(cutree(tree = test2$tree_col, 2))

cluster1_GO <- names(cutree(tree = test2$tree_col, 2)[cutree(tree = test2$tree_col, 2)==1])
sapply(cluster1_GO, function(x) Term(GOTERM[[x]]))

cluster2_GO <- names(cutree(tree = test2$tree_col, 2)[cutree(tree = test2$tree_col, 2)==2])
sapply(cluster2_GO, function(x) Term(GOTERM[[x]]))


Term(GOTERM[["GO:0006413"]])
GOBPANCESTOR[["GO:0006413"]]

tang <- do.call(rbind,lapply(as.list(cluster1_GO), function(x) {select(GO.db, keys = GOBPANCESTOR[[x]], columns = "TERM")}))

table(tang$TERM)[order(table(tang$TERM))]
tail(table(tang$TERM)[order(table(tang$TERM))],11)

tang2 <- lapply(as.list(cluster2_GO), function(x) {select(GO.db, keys = GOBPANCESTOR[[x]], columns = "TERM")})
tang2 <- do.call(rbind,tang2)
table(tang2$TERM)[order(table(tang2$TERM))]
tail(table(tang2$TERM)[order(table(tang2$TERM))],11)


names(tail(table(tang$TERM)[order(table(tang$TERM))],11)) -> n1
names(tail(table(tang2$TERM)[order(table(tang2$TERM))],11)) -> n2

n1[!n1 %in% Reduce(intersect, list(n1,n2))]
n2[!n2 %in% Reduce(intersect, list(n1,n2))]



#### GO slim at high level =============================================================================================

myGOslim <- function(x){
  requireNamespace("topGO") || stop("package topGO is required")
  groupGOTerms <- rvcheck::get_fun_from_pkg("topGO", "groupGOTerms")
  annFUN.gene2GO <- rvcheck::get_fun_from_pkg("topGO", "annFUN.gene2GO")
  if (!class(x) %in% c("gseaResult", "enrichResult")) {
    stop("x should be output of gseGO or enrichGO...")
  }
  gene2GO <- inverseList(x@geneSets)
  if (is(x, "gseaResult")) {
    ont <- x@setType
    allgenes <- x@geneList
    core_genes <- unique(unlist(geneInCategory(x)))
    allgenes[!names(allgenes) %in% core_genes] <- -1
    allgenes[core_genes] <- 1
  }
  else {
    ont <- x@ontology
    universe <- x@universe
    allgenes <- numeric(length(universe))
    names(allgenes) <- universe
    allgenes[x@gene] <- 1
  }
  selector <- function(scores) return(scores == 1)
  if (!ont %in% c("BP", "MF", "CC")) {
    stop("ontology should be one of 'BP', 'MF' or 'CC'...")
  }
  pvalue <- x@result$p.adjust
  names(pvalue) <- x@result$ID
  groupGOTerms()
  GOdata <- new("topGOdata", description = "clusterProfiler enrichment results", 
                ontology = ont, 
                allGenes = allgenes, 
                geneSel = selector, 
                annot = annFUN.gene2GO, 
                gene2GO = gene2GO)
  result <- buildLevels(graph(GOdata))
  return(result)
}

SE_DAG <- myGOslim(SE_BP_ego)
A3SS_DAG <- myGOslim(A3SS_BP_ego)
A5SS_DAG <- myGOslim(A5SS_BP_ego)
AFE_DAG <- myGOslim(AFE_BP_ego)
ALE_DAG <- myGOslim(ALE_BP_ego)

getNoOfLevels(SE_DAG)
getNoOfLevels(A3SS_DAG)
getNoOfLevels(A5SS_DAG)
getNoOfLevels(AFE_DAG)
getNoOfLevels(ALE_DAG)

SE_DAG$level2nodes$`5`
A3SS_DAG$level2nodes$`5`
A5SS_DAG$level2nodes$`5`
AFE_DAG$level2nodes$`5`
ALE_DAG$level2nodes$`5`

lapply(FUN = length, list(SE_DAG$level2nodes$`2`,
                          A3SS_DAG$level2nodes$`2`,
                          A5SS_DAG$level2nodes$`2`,
                          AFE_DAG$level2nodes$`2`,
                          ALE_DAG$level2nodes$`2`))

length(Reduce(intersect, list(SE_DAG$level2nodes$`2`,
                              A3SS_DAG$level2nodes$`2`,
                              A5SS_DAG$level2nodes$`2`,
                              AFE_DAG$level2nodes$`2`,
                              ALE_DAG$level2nodes$`2`)))


select(GO.db, keys = SE_DAG$level2nodes$`2`, columns = "TERM")


SE$Description[1]
Term(GOTERM[[SE$ID[1]]])

select(GO.db, keys = GOBPANCESTOR[[SE$ID[1]]][GOBPANCESTOR[[SE$ID[1]]] %in% SE_DAG$level2nodes$`3`], columns = "TERM")
select(GO.db, keys = GOBPANCESTOR[[SE$ID[2]]][GOBPANCESTOR[[SE$ID[2]]] %in% SE_DAG$level2nodes$`3`], columns = "TERM")
select(GO.db, keys = GOBPANCESTOR[[SE$ID[12]]][GOBPANCESTOR[[SE$ID[12]]] %in% SE_DAG$level2nodes$`3`], columns = "TERM")




#### enrichment at level 3

SE[SE$ID %in% SE_DAG$level2nodes$`3`,]$Description
#[1] "multi-organism metabolic process"    "antigen processing and presentation"
#[3] "methylation"                        
A3SS[A3SS$ID %in% A3SS_DAG$level2nodes$`3`,]$Description
#[1] "multi-organism metabolic process"    "antigen processing and presentation"
#[3] "maintenance of cell number"          "protein folding"                    
A5SS[A5SS$ID %in% A5SS_DAG$level2nodes$`3`,]$Description
#[1] "multi-organism metabolic process"    "methylation"                        
#[3] "protein folding"                     "antigen processing and presentation"
AFE[AFE$ID %in% AFE_DAG$level2nodes$`3`,]$Description
# [1] "production of molecular mediator of immune response"
ALE[ALE$ID %in% ALE_DAG$level2nodes$`4`,]$Description
# [1] "ribonucleoprotein complex biogenesis"



SE_L3 <-SE[SE$ID %in% SE_DAG$level2nodes$`3`,]
A3SS_L3 <- A3SS[A3SS$ID %in% A3SS_DAG$level2nodes$`3`,]
A5SS_L3 <- A5SS[A5SS$ID %in% A5SS_DAG$level2nodes$`3`,]
AFE_L3 <- AFE[AFE$ID %in% AFE_DAG$level2nodes$`3`,]
ALE_L3 <- ALE[ALE$ID %in% ALE_DAG$level2nodes$`4`,]


length(unique(c(SE_L3$ID, A3SS_L3$ID, A5SS_L3$ID, AFE_L3$ID, ALE_L3$ID)))
mat_l3 <- matrix(0, nrow = 5, ncol = length(unique(c(SE_L3$ID, A3SS_L3$ID, A5SS_L3$ID, AFE_L3$ID, ALE_L3$ID))))
row.names(mat_l3) <- c("SE","A3SS","A5SS","AFE","ALE")
colnames(mat_l3) <- unique(c(SE_L3$Description, A3SS_L3$Description, A5SS_L3$Description, AFE_L3$Description, ALE_L3$Description))

for(i in 1:nrow(SE_L3)){
  mat_l3[1, SE_L3[i,"Description"]] <- SE_L3[i,"Count"]
}

for(i in 1:nrow(A3SS_L3)){
  mat_l3[2, A3SS_L3[i,"Description"]] <- A3SS_L3[i,"Count"]
}

for(i in 1:nrow(A5SS_L3)){
  mat_l3[3, A5SS_L3[i,"Description"]] <- A5SS_L3[i,"Count"]
}

for(i in 1:nrow(AFE_L3)){
  mat_l3[4, AFE_L3[i,"Description"]] <- AFE_L3[i,"Count"]
}

for(i in 1:nrow(ALE_L3)){
  mat_l3[5, ALE_L3[i,"Description"]] <- ALE_L3[i,"Count"]
}

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/AS_GO/one cell/AS_GO_Count_circos_at_levels3.pdf")
#circos.par("canvas.xlim" = c(-2, 2), "canvas.ylim" = c(-2, 2), start.degree = 230, gap.after = c(rep(2, nrow(mat_l3)-1), 10, rep(2, ncol(mat_l3)-1), 10))#
circos.par(start.degree = 240, gap.after = c(rep(2, nrow(mat_l3)-1), 60, rep(2, ncol(mat_l3)-1), 60))# 块儿与块儿之间的空格
chordDiagram(mat_l3, annotationTrack = "grid", transparency = 0.5, 
             grid.col = c("#FF0000","#00FF00","#0000FF","#FFFF00","#00FFFF",colorRampPalette(brewer.pal(7, "Paired"))(ncol(mat_l3))),
             row.col = c("#FF0000","#00FF00","#0000FF","#FFFF00","#00FFFF"),
             preAllocateTracks = list(track.height = 0.1))
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  if(sector.name %in% get.all.sector.index()[1:5]){
    circos.text(mean(xlim), 0.3*nchar(sector.name), sector.name, niceFacing = T, facing = "clockwise")
  }
}, bg.border = NA)

circos.clear()

library(ComplexHeatmap)
lgd_points = Legend(at = colnames(mat_l3), type = "points", 
                    legend_gp = gpar(col = colorRampPalette(brewer.pal(7, "Paired"))(ncol(mat_l3))), 
                    background = colorRampPalette(brewer.pal(7, "Paired"))(ncol(mat_l3)),
                    title_position = "topleft", 
                    title = "GO Term", 
                    nrow = 7)

lgd_list_vertical = packLegend(lgd_points, direction = "vertical")

pushViewport(viewport(x = unit(105, "mm"), y = unit(20, "mm"), 
                      width = grobWidth(lgd_list_vertical), 
                      height = grobHeight(lgd_list_vertical)))
grid.draw(lgd_list_vertical)
upViewport()
dev.off()









