library(org.Hs.eg.db)##### data base target to human-geneIDs
library(AnnotationDbi)
library(DOSE)
library(GO.db)
library(topGO)
library(GSEABase)
library(clusterProfiler)


load("/Users/tangchao/CloudStation/tangchao/project/UCB/Cell_sepecific_gene_and_sj_host_gene.RData")
background <- as.character(read.table("/Users/tangchao/CloudStation/tangchao/project/UCB/AS_GO/Background_for_GO.txt", header = F, stringsAsFactors = F)[,1])

length(Bcell_gene);length(sj_Bcell_gene);length(intersect(sj_Bcell_gene,Bcell_gene))
length(CD4_gene);length(sj_CD4_gene);length(intersect(sj_CD4_gene,CD4_gene))
length(CD8_gene);length(sj_CD8_gene);length(intersect(sj_CD8_gene,CD8_gene))
length(NK_gene);length(sj_NK_gene);length(intersect(sj_NK_gene,NK_gene))
length(Mono_gene);length(sj_Mono_gene);length(intersect(sj_Mono_gene,Mono_gene))
length(MK_gene);length(sj_MK_gene);length(intersect(sj_MK_gene,MK_gene))

library(VennDiagram)
library(VennDiagram)
venn_plot <- venn.diagram(list(`Gene Expression` = Bcell_gene, `Intron-centric PSI` = sj_Bcell_gene), 
                          fill = rainbow(2), filename = NULL, main.cex = 2, 
                          main = "Venn plot of B cell specific gene and SJ host gene", 
                          height = 6, width = 20)

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/Cell_type_specific/Venn plot of B cell specific gene and SJ host gene.pdf", 
    width = 19, height = 12)
grid.draw(venn_plot)
grid.draw(venn_plot)
dev.off()


venn_plot <- venn.diagram(list(`Gene Expression` = CD4_gene, `Intron-centric PSI` = sj_CD4_gene), 
                          fill = rainbow(2), filename = NULL, main.cex = 2, 
                          main = "Venn plot of CD4+ T cell specific gene and SJ host gene", 
                          height = 6, width = 20)
pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/Cell_type_specific/Venn plot of CD4 T cell specific gene and SJ host gene.pdf", 
    width = 19, height = 12)
grid.draw(venn_plot)
grid.draw(venn_plot)
dev.off()


venn_plot <- venn.diagram(list(`Gene Expression` = CD8_gene, `Intron-centric PSI` = sj_CD8_gene), 
                          fill = rainbow(2), filename = NULL, main.cex = 2, 
                          main = "Venn plot of CD8+ T cell specific gene and SJ host gene", 
                          height = 6, width = 20)
pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/Cell_type_specific/Venn plot of CD8 T cell specific gene and SJ host gene.pdf", 
    width = 19, height = 12)
grid.draw(venn_plot)
grid.draw(venn_plot)
dev.off()


venn_plot <- venn.diagram(list(`Gene Expression` = NK_gene, `Intron-centric PSI` = sj_NK_gene), 
                          fill = rainbow(2), filename = NULL, main.cex = 2, 
                          main = "Venn plot of NK cell specific gene and SJ host gene", 
                          height = 6, width = 20)
pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/Cell_type_specific/Venn plot of NK cell specific gene and SJ host gene.pdf", 
    width = 19, height = 12)
grid.draw(venn_plot)
grid.draw(venn_plot)
dev.off()


venn_plot <- venn.diagram(list(`Gene Expression` = Mono_gene, `Intron-centric PSI` = sj_Mono_gene), 
                          fill = rainbow(2), filename = NULL, main.cex = 2, 
                          main = "Venn plot of Monocytes specific gene and SJ host gene", 
                          height = 6, width = 20)
pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/Cell_type_specific/Venn plot of Monocytes specific gene and SJ host gene.pdf", 
    width = 19, height = 12)
grid.draw(venn_plot)
grid.draw(venn_plot)
dev.off()


venn_plot <- venn.diagram(list(`Gene Expression` = MK_gene, `Intron-centric PSI` = sj_MK_gene), 
                          fill = rainbow(2), filename = NULL, main.cex = 2, 
                          main = "Venn plot of MK cell specific gene and SJ host gene", 
                          height = 6, width = 20)
pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/Cell_type_specific/Venn plot of MK cell specific gene and SJ host gene.pdf", 
    width = 19, height = 12)
grid.draw(venn_plot)
grid.draw(venn_plot)
dev.off()


#### GO ==================================================================================================================
#### Bcell
Bcell_BP_ego <- enrichGO(gene = sj_Bcell_gene[!sj_Bcell_gene %in% Bcell_gene],
                         universe = background,
                         OrgDb  = org.Hs.eg.db,
                         keytype= 'ENSEMBL',
                         ont  = "BP", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)

#### CD4
CD4_BP_ego <- enrichGO(gene = sj_CD4_gene[!sj_CD4_gene %in% CD4_gene],
                         universe = background,
                         OrgDb  = org.Hs.eg.db,
                         keytype= 'ENSEMBL',
                         ont  = "BP", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)

#### CD8
CD8_BP_ego <- enrichGO(gene = sj_CD8_gene[!sj_CD8_gene %in% CD8_gene],
                       universe = background,
                       OrgDb  = org.Hs.eg.db,
                       keytype= 'ENSEMBL',
                       ont  = "BP", 
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)

#### NK
NK_BP_ego <- enrichGO(gene = sj_NK_gene[!sj_NK_gene %in% NK_gene],
                       universe = background,
                       OrgDb  = org.Hs.eg.db,
                       keytype= 'ENSEMBL',
                       ont  = "BP", 
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)

#### Mono
Mono_BP_ego <- enrichGO(gene = sj_Mono_gene[!sj_Mono_gene %in% Mono_gene],
                      universe = background,
                      OrgDb  = org.Hs.eg.db,
                      keytype= 'ENSEMBL',
                      ont  = "BP", 
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

#### MK
MK_BP_ego <- enrichGO(gene = sj_MK_gene[!sj_MK_gene %in% MK_gene],
                      universe = background,
                      OrgDb  = org.Hs.eg.db,
                      keytype= 'ENSEMBL',
                      ont  = "BP", 
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/Cell_type_specific/Bcell_GO.pdf",height = 10,width = 15)
print(dotplot(Bcell_BP_ego,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(Bcell_BP_ego,firstSigNodes = 15, useInfo = "all"))
enrichMap(Bcell_BP_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/Cell_type_specific/CD4_GO.pdf",height = 10,width = 15)
print(dotplot(CD4_BP_ego,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(CD4_BP_ego,firstSigNodes = 15, useInfo = "all"))
enrichMap(CD4_BP_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/Cell_type_specific/CD8_GO.pdf",height = 10,width = 15)
print(dotplot(CD8_BP_ego,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(CD8_BP_ego,firstSigNodes = 15, useInfo = "all"))
enrichMap(CD8_BP_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/Cell_type_specific/NK_GO.pdf",height = 10,width = 15)
print(dotplot(NK_BP_ego,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(NK_BP_ego,firstSigNodes = 15, useInfo = "all"))
enrichMap(NK_BP_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/Cell_type_specific/Mono_GO.pdf",height = 10,width = 15)
print(dotplot(Mono_BP_ego,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(Mono_BP_ego,firstSigNodes = 15, useInfo = "all"))
enrichMap(Mono_BP_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/Cell_type_specific/MK_GO.pdf",height = 10,width = 15)
print(dotplot(MK_BP_ego,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(MK_BP_ego,firstSigNodes = 15, useInfo = "all"))
enrichMap(MK_BP_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


#### All gene GO ==========================================================================================================
#### Bcell
Bcell_BP_ego <- enrichGO(gene = sj_Bcell_gene,
                         universe = background,
                         OrgDb  = org.Hs.eg.db,
                         keytype= 'ENSEMBL',
                         ont  = "BP", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)

#### CD4
CD4_BP_ego <- enrichGO(gene = sj_CD4_gene,
                       universe = background,
                       OrgDb  = org.Hs.eg.db,
                       keytype= 'ENSEMBL',
                       ont  = "BP", 
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)

#### CD8
CD8_BP_ego <- enrichGO(gene = sj_CD8_gene,
                       universe = background,
                       OrgDb  = org.Hs.eg.db,
                       keytype= 'ENSEMBL',
                       ont  = "BP", 
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)

#### NK
NK_BP_ego <- enrichGO(gene = sj_NK_gene,
                      universe = background,
                      OrgDb  = org.Hs.eg.db,
                      keytype= 'ENSEMBL',
                      ont  = "BP", 
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

#### Mono
Mono_BP_ego <- enrichGO(gene = sj_Mono_gene,
                        universe = background,
                        OrgDb  = org.Hs.eg.db,
                        keytype= 'ENSEMBL',
                        ont  = "BP", 
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

#### MK
MK_BP_ego <- enrichGO(gene = sj_MK_gene,
                      universe = background,
                      OrgDb  = org.Hs.eg.db,
                      keytype= 'ENSEMBL',
                      ont  = "BP", 
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/Cell_type_specific/Bcell_all_host_gene_GO.pdf",height = 10,width = 15)
print(dotplot(Bcell_BP_ego,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(Bcell_BP_ego,firstSigNodes = 15, useInfo = "all"))
enrichMap(Bcell_BP_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/Cell_type_specific/CD4_all_host_gene_GO.pdf",height = 10,width = 15)
print(dotplot(CD4_BP_ego,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(CD4_BP_ego,firstSigNodes = 15, useInfo = "all"))
enrichMap(CD4_BP_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/Cell_type_specific/CD8_all_host_gene_GO.pdf",height = 10,width = 15)
print(dotplot(CD8_BP_ego,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(CD8_BP_ego,firstSigNodes = 15, useInfo = "all"))
enrichMap(CD8_BP_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/Cell_type_specific/NK_all_host_gene_GO.pdf",height = 10,width = 15)
print(dotplot(NK_BP_ego,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(NK_BP_ego,firstSigNodes = 15, useInfo = "all"))
enrichMap(NK_BP_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/Cell_type_specific/Mono_all_host_gene_GO.pdf",height = 10,width = 15)
print(dotplot(Mono_BP_ego,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(Mono_BP_ego,firstSigNodes = 15, useInfo = "all"))
enrichMap(Mono_BP_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/Cell_type_specific/MK_all_host_gene_GO.pdf",height = 10,width = 15)
print(dotplot(MK_BP_ego,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(MK_BP_ego,firstSigNodes = 15, useInfo = "all"))
enrichMap(MK_BP_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()







#### cell_sep_unannotated_SJ_gene ====
load("/Users/tangchao/CloudStation/tangchao/project/UCB/cell_sep_unannotated_SJ_gene.RData")
Unan_BP_ego <- enrichGO(gene = cell_sep_unannotated_SJ_gene,
                      universe = background,
                      OrgDb  = org.Hs.eg.db,
                      keytype= 'ENSEMBL',
                      ont  = "BP", 
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/Cell_type_specific/Unannotated_all_SJ_host_gene_GO.pdf",height = 10,width = 15)
print(dotplot(Unan_BP_ego,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(Unan_BP_ego,firstSigNodes = 15, useInfo = "all"))
enrichMap(Unan_BP_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()




#### All unannotated & novel SJ host gene GO ====
load("/Users/tangchao/CloudStation/tangchao/project/UCB/all_unannotated_and_novel_SJ_gene.RData")
all_Unan_BP_ego <- enrichGO(gene = unannotated_SJ_gene,
                        universe = background,
                        OrgDb  = org.Hs.eg.db,
                        keytype= 'ENSEMBL',
                        ont  = "BP", 
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)

all_Novel_BP_ego <- enrichGO(gene = novel_SJ_gene,
                            universe = background,
                            OrgDb  = org.Hs.eg.db,
                            keytype= 'ENSEMBL',
                            ont  = "BP", 
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05)

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/NSS4/all_Unannotated_SJ_host_gene_GO.pdf",height = 10,width = 15)
print(dotplot(all_Unan_BP_ego,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(all_Unan_BP_ego,firstSigNodes = 15, useInfo = "all"))
enrichMap(all_Unan_BP_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/Users/tangchao/CloudStation/tangchao/project/UCB/NSS4/all_Novel_SJ_host_gene_GO.pdf",height = 10,width = 15)
print(dotplot(all_Novel_BP_ego,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(all_Novel_BP_ego,firstSigNodes = 15, useInfo = "all"))
enrichMap(all_Novel_BP_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


View(all_Novel_BP_ego@result)
novel_go_result <- all_Novel_BP_ego@result
novel_go_result[novel_go_result$ID %in% c("GO:0050852","GO:0050851"),]
as.character(novel_go_result[novel_go_result$ID %in% c("GO:0050852","GO:0050851"),]$geneID)
columns(GO.db)

novel_SJ_gene[novel_SJ_gene %in% all_Novel_BP_ego@geneSets$`GO:0050852`]
novel_SJ_gene[novel_SJ_gene %in% all_Novel_BP_ego@geneSets$`GO:0050851`]
