#### AS GO
library(org.Hs.eg.db)##### data base target to human-geneIDs
library(AnnotationDbi)
library(DOSE)
library(GO.db)
library(topGO)
library(GSEABase)
library(clusterProfiler)

#background <- as.character(read.table("/Users/tangchao/CloudStation/tangchao/project/UCB/AS_GO/Background_for_GO.txt", header = F, stringsAsFactors = F)[,1])
background <- as.character(read.table("/mnt/data5/BGI/UCB/tangchao/DSU2/Document/Background_for_GO.txt", header = F, stringsAsFactors = F)[,1])
Cell_type <- read.table("/mnt/data5/BGI/UCB/ExpMat_NewID/Cell_type.txt", header = F, row.names = 1, sep = "\t", stringsAsFactors=F)
colnames(Cell_type) <- "CellType"
Cell_type <- data.frame(CellType = Cell_type[order(Cell_type),], row.names = row.names(Cell_type)[order(Cell_type)])
Cell_type$Individual <- substr(row.names(Cell_type), 1, 4)



#### SE =========================================================================================================================================



load("/mnt/data5/BGI/UCB/tangchao/DSU2/SE/SE_psi_new.RData")
PSI <- SE_psi[, row.names(Cell_type)]
row.names(PSI) <- SE_psi$loci
rs <- rowSums(!is.na(PSI))
summary(rs)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    0.0    17.0    76.0   256.3   273.8  3039.0
sum(rs==0)
#[1] 53
sum(rs==1)
#[1] 513
PSI <- PSI[rs!=0, ]

rs <- rowSums(!is.na(PSI) & PSI != 0 & PSI != 1)
# summary(rs)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   0.00    0.00    1.00   12.62    3.00 2600.00
sum(rs==0)
# [1] 7802
sum(rs==1)
# [1] 3404
table(SE_psi[SE_psi$loci %in% names(rs)[rs>=1], "SE_type"])
# exon-skipping_exactly   novel_exon-skipping                 Other
#                  5272                  1859                  2684
table(SE_psi[SE_psi$loci %in% names(rs)[rs==0], "SE_type"])
# exon-skipping_exactly   novel_exon-skipping                 Other
#                  5088                  1529                  1185

all_SE_gene <- unique(SE_psi[SE_psi$loci %in% row.names(PSI), "gene_id"])
length(all_SE_gene)
#[1] 5890

SE_gene_0 <- unique(SE_psi[SE_psi$loci %in% names(rs)[rs==0], "gene_id"])
SE_gene_1 <- unique(SE_psi[SE_psi$loci %in% names(rs)[rs>=1], "gene_id"])
length(SE_gene_0)
#[1] 3686
length(SE_gene_1)
#[1] 3045

#### GO

SE_BP_ego <- enrichGO(gene          = all_SE_gene,
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
SE_MF_ego <- enrichGO(gene          = all_SE_gene,
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
SE_CC_ego <- enrichGO(gene          = all_SE_gene,
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/SE/all_SE_gene_GO.pdf",height = 10,width = 15)
print(dotplot(SE_BP_ego,showCategory = 10,font.size = 16,title = "BP"))
print(dotplot(SE_MF_ego,showCategory = 10,font.size = 16,title = "MF"))
print(dotplot(SE_CC_ego,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(SE_BP_ego,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(SE_MF_ego,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(SE_CC_ego,firstSigNodes = 15, useInfo = "all"))
enrichMap(SE_BP_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(SE_MF_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(SE_CC_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()



SE_BP_ego_0 <- enrichGO(gene          = SE_gene_0,
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "BP", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
SE_MF_ego_0 <- enrichGO(gene          = SE_gene_0,
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "MF", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
SE_CC_ego_0 <- enrichGO(gene          = SE_gene_0,
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "CC", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/SE/Pseudo_SE_gene_GO.pdf",height = 10,width = 15)
print(dotplot(SE_BP_ego_0,showCategory = 10,font.size = 16,title = "BP"))
print(dotplot(SE_MF_ego_0,showCategory = 10,font.size = 16,title = "MF"))
print(dotplot(SE_CC_ego_0,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(SE_BP_ego_0,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(SE_MF_ego_0,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(SE_CC_ego_0,firstSigNodes = 15, useInfo = "all"))
enrichMap(SE_BP_ego_0,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(SE_MF_ego_0,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(SE_CC_ego_0,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


SE_BP_ego_1 <- enrichGO(gene          = SE_gene_1,
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "BP", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
SE_MF_ego_1 <- enrichGO(gene          = SE_gene_1,
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "MF", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
SE_CC_ego_1 <- enrichGO(gene          = SE_gene_1,
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "CC", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/SE/Real_SE_gene_GO.pdf",height = 10,width = 15)
print(dotplot(SE_BP_ego_1,showCategory = 10,font.size = 16,title = "BP"))
print(dotplot(SE_MF_ego_1,showCategory = 10,font.size = 16,title = "MF"))
print(dotplot(SE_CC_ego_1,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(SE_BP_ego_1,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(SE_MF_ego_1,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(SE_CC_ego_1,firstSigNodes = 15, useInfo = "all"))
enrichMap(SE_BP_ego_1,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(SE_MF_ego_1,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(SE_CC_ego_1,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


background_enid <- as.character(na.omit(mapIds(org.Hs.eg.db, keys = background, column = "ENTREZID", keytype = "ENSEMBL")))
all_SE_gene_enid <- as.character(na.omit(mapIds(org.Hs.eg.db, keys = all_SE_gene, column = "ENTREZID", keytype = "ENSEMBL")))
all_SE_gene_enid <- all_SE_gene_enid[!is.na(all_SE_gene_enid)]

keytypes(org.Hs.eg.db)
background_enid <- bitr(background, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
all_SE_gene_enid <- bitr(all_SE_gene, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID

SE_gene_0_enid <- bitr(SE_gene_0, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
SE_gene_1_enid <- bitr(SE_gene_1, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID


SE_KEGG <- enrichKEGG(gene         = all_SE_gene_enid,
                 	  organism     = 'hsa',
                 	  universe	   = background_enid,
                 	  pvalueCutoff = 0.05)
head(SE_KEGG)

SE_0_KEGG <- enrichKEGG(gene         = SE_gene_0_enid,
                 	  	organism     = 'hsa',
                 	  	universe	 = background_enid,
                 	  	pvalueCutoff = 0.05)

SE_1_KEGG <- enrichKEGG(gene         = SE_gene_1_enid,
                 	  	organism     = 'hsa',
                 	  	universe	 = background_enid,
                 	  	pvalueCutoff = 0.05)



SE_mkk <- enrichMKEGG(gene     = all_SE_gene_enid,
                   	  organism = 'hsa',
                   	  universe = background_enid)
SE_0_mkk <- enrichMKEGG(gene     = SE_gene_0_enid,
                   	  	organism = 'hsa',
                   	  	universe = background_enid)
SE_1_mkk <- enrichMKEGG(gene     = SE_gene_1_enid,
                   	  	organism = 'hsa',
                   	  	universe = background_enid)
head(SE_mkk)



SE_do <- DOSE::enrichDO(gene    	  = all_SE_gene_enid, 
				  		ont           = "DO", 
				  		pvalueCutoff  = 0.05, 
				  		pAdjustMethod = "BH",
				  		#universe 	  = background_enid,
				  		minGSSize     = 5, 
				  		maxGSSize     = 500, 
				  		qvalueCutoff  = 0.05,
				  		readable      = FALSE)
SE_0_do <- DOSE::enrichDO(gene    	    = SE_gene_0_enid, 
				  		  ont           = "DO", 
				  		  pvalueCutoff  = 0.05, 
				  		  pAdjustMethod = "BH",
				  		  #universe 	 	= background_enid,
				  		  minGSSize     = 5, 
				  		  maxGSSize     = 500, 
				  		  qvalueCutoff  = 0.05,
				  		  readable      = FALSE)
SE_1_do <- DOSE::enrichDO(gene    	    = SE_gene_1_enid, 
				  		  ont           = "DO", 
				  		  pvalueCutoff  = 0.05, 
				  		  pAdjustMethod = "BH",
				  		  #universe 	 	= background_enid,
				  		  minGSSize     = 5, 
				  		  maxGSSize     = 500, 
				  		  qvalueCutoff  = 0.05,
				  		  readable      = FALSE)
head(SE_do)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/SE/all_SE_gene_DO.pdf",height = 10,width = 15)
print(dotplot(SE_do,showCategory = 10,font.size = 16, title = "DO"))
enrichMap(SE_do,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/SE/Pseudo_SE_gene_DO.pdf",height = 10,width = 15)
print(dotplot(SE_0_do,showCategory = 10,font.size = 16, title = "DO"))
enrichMap(SE_0_do,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/SE/Real_SE_gene_DO.pdf",height = 10,width = 15)
print(dotplot(SE_1_do,showCategory = 10,font.size = 16, title = "DO"))
enrichMap(SE_1_do,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


dgn <- enrichDGN(gene         = all_SE_gene_enid,
				 #universe     = background_enid,
				 qvalueCutoff = 0.05,
				 readable     = TRUE)
dgn_0 <- enrichDGN(gene         = SE_gene_0_enid,
				   #universe     = background_enid,
				   qvalueCutoff = 0.05,
				   readable     = TRUE)
dgn_1 <- enrichDGN(gene         = SE_gene_1_enid,
				   #universe     = background_enid,
				   qvalueCutoff = 0.05,
				   readable     = TRUE)
head(dgn)
pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/SE/all_SE_gene_DGN.pdf",height = 10,width = 15)
print(dotplot(dgn,showCategory = 10,font.size = 16, title = "DO"))
enrichMap(dgn,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/SE/Pseudo_SE_gene_DGN.pdf",height = 10,width = 15)
print(dotplot(dgn_0,showCategory = 10,font.size = 16, title = "DO"))
enrichMap(dgn_0,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/SE/Real_SE_gene_DGN.pdf",height = 10,width = 15)
print(dotplot(dgn_1,showCategory = 10,font.size = 16, title = "DO"))
enrichMap(dgn_1,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()





library(ReactomePA)
SE_Path <- enrichPathway(gene = all_SE_gene_enid, 
						 organism = "human", 
						 pvalueCutoff = 0.05,
						 pAdjustMethod = "BH", 
						 qvalueCutoff = 0.2, 
						 universe = background_enid, 
						 minGSSize = 10,
						 maxGSSize = 500, 
						 readable = FALSE)

SE_0_Path <- enrichPathway(gene = SE_gene_0_enid, 
						   organism = "human", 
						   pvalueCutoff = 0.05,
						   pAdjustMethod = "BH", 
						   qvalueCutoff = 0.2, 
						   universe = background_enid, 
						   minGSSize = 10,
						   maxGSSize = 500, 
						   readable = FALSE)

SE_1_Path <- enrichPathway(gene = SE_gene_1_enid, 
						   organism = "human", 
						   pvalueCutoff = 0.05,
						   pAdjustMethod = "BH", 
						   qvalueCutoff = 0.2, 
						   universe = background_enid, 
						   minGSSize = 10,
						   maxGSSize = 500, 
						   readable = FALSE)



SE_BP_ego_1
SE_1_KEGG
dgn_1
print(dotplot(SE_BP_ego_1,showCategory = 10,font.size = 16,title = "BP"))
print(dotplot(SE_1_KEGG,showCategory = 10,font.size = 16,title = "BP"))
print(dotplot(dgn_1,showCategory = 10,font.size = 16,title = "BP"))




pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/SE/all_SE_gene_KEGG.pdf",height = 10,width = 15)
print(dotplot(SE_KEGG,showCategory = 10,font.size = 16, title = "BP"))
#suppressMessages(plotGOgraph(SE_KEGG,firstSigNodes = 15, useInfo = "all"))
enrichMap(SE_KEGG,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/SE/Pseudo_SE_gene_KEGG.pdf",height = 10,width = 15)
print(dotplot(SE_0_KEGG,showCategory = 10,font.size = 16, title = "BP"))
enrichMap(SE_0_KEGG,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/SE/Real_SE_gene_KEGG.pdf",height = 10,width = 15)
print(dotplot(SE_1_KEGG,showCategory = 10,font.size = 16, title = "BP"))
enrichMap(SE_1_KEGG,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/SE/Real_SE_gene_DGN.pdf",height = 10,width = 15)
print(dotplot(dgn_1,showCategory = 10,font.size = 16, title = "BP"))
#suppressMessages(plotGOgraph(SE_KEGG,firstSigNodes = 15, useInfo = "all"))
enrichMap(dgn_1,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/SE/all_SE_gene_Pathway.pdf",height = 10,width = 15)
print(dotplot(SE_Path,showCategory = 10,font.size = 16, title = "BP"))
#suppressMessages(plotGOgraph(SE_Path,firstSigNodes = 15, useInfo = "all"))
enrichMap(SE_Path,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/SE/Pseudo_SE_gene_Pathway.pdf",height = 10,width = 15)
print(dotplot(SE_0_Path,showCategory = 10,font.size = 16, title = "BP"))
enrichMap(SE_0_Path,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/SE/Real_SE_gene_Pathway.pdf",height = 10,width = 15)
print(dotplot(SE_1_Path,showCategory = 10,font.size = 16, title = "BP"))
enrichMap(SE_1_Path,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()



CompareGO_BP <- compareCluster(geneCluster   = list(All = all_SE_gene, Pseudo = SE_gene_0, Real = SE_gene_1), 
                               universe      = background,
                               fun 		  ="enrichGO", 
                               OrgDb  		  = org.Hs.eg.db, 
                               keyType 	  = 'ENSEMBL', 
                               ont  		  = "BP", 
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)
CompareGO_MF <- compareCluster(geneCluster   = list(All = all_SE_gene, Pseudo = SE_gene_0, Real = SE_gene_1), 
                               universe      = background,
                               fun 		  ="enrichGO", 
                               OrgDb  		  = org.Hs.eg.db, 
                               keyType 	  = 'ENSEMBL', 
                               ont  		  = "MF", 
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)
CompareGO_CC <- compareCluster(geneCluster   = list(All = all_SE_gene, Pseudo = SE_gene_0, Real = SE_gene_1), 
                               universe      = background,
                               fun 		  ="enrichGO", 
                               OrgDb  		  = org.Hs.eg.db, 
                               keyType 	  = 'ENSEMBL', 
                               ont  		  = "CC", 
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/SE/SE_CompareGO.pdf", height = 8, width = 20)
print(dotplot(CompareGO_BP,showCategory = 10,font.size = 16,title = "SE GO BP"))
print(dotplot(CompareGO_MF,showCategory = 10,font.size = 16,title = "SE GO MF"))
print(dotplot(CompareGO_CC,showCategory = 10,font.size = 16,title = "SE GO CC"))
dev.off()




CompareKEGG <- compareCluster(geneCluster  = list(All = all_SE_gene_enid, Pseudo = SE_gene_0_enid, Real = SE_gene_1_enid),
							  universe 	   = background_enid,
							  fun          = "enrichKEGG",
							  organism	   = "hsa", 
							  pvalueCutoff = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/SE/SE_CompareKEGG.pdf", height = 8, width = 20)
print(dotplot(CompareKEGG,showCategory = 10,font.size = 16,title = "SE KEGG"))
dev.off()



CompareDO <- compareCluster(geneCluster   = list(All = all_SE_gene_enid, Pseudo = SE_gene_0_enid, Real = SE_gene_1_enid),
							fun 		  = "enrichDO",
				  	  		ont           = "DO", 
				  	  		pvalueCutoff  = 0.05, 
				  	  		pAdjustMethod = "BH",
				  	  		#universe 	  = background_enid,
				  	  		minGSSize     = 5, 
				  	  		maxGSSize     = 500, 
				  	  		qvalueCutoff  = 0.05,
				  	  		readable      = FALSE)
pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/SE/SE_CompareDO.pdf", height = 8, width = 20)
print(dotplot(CompareDO,showCategory = 10,font.size = 16,title = "SE KEGG"))
dev.off()


library(ReactomePA)
ComparePathway <- compareCluster(geneCluster   = list(All = all_SE_gene_enid, Pseudo = SE_gene_0_enid, Real = SE_gene_1_enid),
								 fun="enrichPathway",
								 organism = "human",
								 pvalueCutoff = 0.05, 
								 pAdjustMethod = "BH",
								 qvalueCutoff = 0.2,
								 universe = background_enid,
								 minGSSize = 10,
								 maxGSSize = 500,
								 readable = FALSE)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/SE/SE_ComparePathway.pdf", height = 8, width = 20)
print(dotplot(ComparePathway,showCategory = 10,font.size = 16,title = "SE Pathway"))
dev.off()






#### A3SS =======================================================================================================================================



load("/mnt/data5/BGI/UCB/tangchao/DSU2/A3SS/A3SS_psi_new.RData")

PSI <- A3SS_psi[, row.names(Cell_type)]
row.names(PSI) <- A3SS_psi$as2
rs <- rowSums(!is.na(PSI))
summary(rs)

PSI <- PSI[rs!=0, ]

rs <- rowSums(!is.na(PSI) & PSI != 0 & PSI != 1)
summary(rs)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.000   0.000   0.000   2.791   1.000 575.000

sum(rs==0)
#[1] 2056
sum(rs==1)
#[1] 470
sum(rs>=1)
#[1] 987

all_A3SS_gene <- unique(A3SS_psi[A3SS_psi$as2 %in% row.names(PSI), "gene_id"])
length(all_A3SS_gene)
#[1] 2412
A3SS_gene_0 <- unique(A3SS_psi[A3SS_psi$as2 %in% names(rs)[rs==0], "gene_id"])
A3SS_gene_1 <- unique(A3SS_psi[A3SS_psi$as2 %in% names(rs)[rs>=1], "gene_id"])
length(A3SS_gene_0)
#[1] 1756
length(A3SS_gene_1)
#[1] 838


#### GO

A3SS_BP_ego <- enrichGO(gene         = all_A3SS_gene,
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
A3SS_MF_ego <- enrichGO(gene          = all_A3SS_gene,
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
A3SS_CC_ego <- enrichGO(gene          = all_A3SS_gene,
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A3SS/all_A3SS_gene_GO.pdf",height = 10,width = 15)
print(dotplot(A3SS_BP_ego,showCategory = 10,font.size = 16,title = "BP"))
print(dotplot(A3SS_MF_ego,showCategory = 10,font.size = 16,title = "MF"))
print(dotplot(A3SS_CC_ego,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(A3SS_BP_ego,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(A3SS_MF_ego,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(A3SS_CC_ego,firstSigNodes = 15, useInfo = "all"))
enrichMap(A3SS_BP_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(A3SS_MF_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(A3SS_CC_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()



A3SS_BP_ego_0 <- enrichGO(gene          = A3SS_gene_0,
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "BP", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
A3SS_MF_ego_0 <- enrichGO(gene          = A3SS_gene_0,
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "MF", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
A3SS_CC_ego_0 <- enrichGO(gene          = A3SS_gene_0,
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "CC", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A3SS/Pseudo_A3SS_gene_GO.pdf",height = 10,width = 15)
print(dotplot(A3SS_BP_ego_0,showCategory = 10,font.size = 16,title = "BP"))
print(dotplot(A3SS_MF_ego_0,showCategory = 10,font.size = 16,title = "MF"))
print(dotplot(A3SS_CC_ego_0,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(A3SS_BP_ego_0,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(A3SS_MF_ego_0,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(A3SS_CC_ego_0,firstSigNodes = 15, useInfo = "all"))
enrichMap(A3SS_BP_ego_0,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(A3SS_MF_ego_0,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(A3SS_CC_ego_0,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


A3SS_BP_ego_1 <- enrichGO(gene          = A3SS_gene_1,
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "BP", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
A3SS_MF_ego_1 <- enrichGO(gene          = A3SS_gene_1,
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "MF", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
A3SS_CC_ego_1 <- enrichGO(gene          = A3SS_gene_1,
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "CC", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A3SS/Real_A3SS_gene_GO.pdf",height = 10,width = 15)
print(dotplot(A3SS_BP_ego_1,showCategory = 10,font.size = 16,title = "BP"))
print(dotplot(A3SS_MF_ego_1,showCategory = 10,font.size = 16,title = "MF"))
print(dotplot(A3SS_CC_ego_1,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(A3SS_BP_ego_1,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(A3SS_MF_ego_1,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(A3SS_CC_ego_1,firstSigNodes = 15, useInfo = "all"))
enrichMap(A3SS_BP_ego_1,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(A3SS_MF_ego_1,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(A3SS_CC_ego_1,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()



background_enid <- as.character(na.omit(mapIds(org.Hs.eg.db, keys = background, column = "ENTREZID", keytype = "ENSEMBL")))
all_A3SS_gene_enid <- as.character(na.omit(mapIds(org.Hs.eg.db, keys = all_A3SS_gene, column = "ENTREZID", keytype = "ENSEMBL")))
all_A3SS_gene_enid <- all_A3SS_gene_enid[!is.na(all_A3SS_gene_enid)]

keytypes(org.Hs.eg.db)
background_enid <- bitr(background, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
all_A3SS_gene_enid <- bitr(all_A3SS_gene, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID

A3SS_gene_0_enid <- bitr(A3SS_gene_0, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
A3SS_gene_1_enid <- bitr(A3SS_gene_1, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID


A3SS_KEGG <- enrichKEGG(gene         = all_A3SS_gene_enid,
                      organism     = 'hsa',
                      universe	   = background_enid,
                      pvalueCutoff = 0.05)
head(A3SS_KEGG)

A3SS_0_KEGG <- enrichKEGG(gene         = A3SS_gene_0_enid,
                        organism     = 'hsa',
                        universe	 = background_enid,
                        pvalueCutoff = 0.05)

A3SS_1_KEGG <- enrichKEGG(gene         = A3SS_gene_1_enid,
                        organism     = 'hsa',
                        universe	 = background_enid,
                        pvalueCutoff = 0.05)



A3SS_mkk <- enrichMKEGG(gene     = all_A3SS_gene_enid,
                      organism = 'hsa',
                      universe = background_enid)
A3SS_0_mkk <- enrichMKEGG(gene     = A3SS_gene_0_enid,
                        organism = 'hsa',
                        universe = background_enid)
A3SS_1_mkk <- enrichMKEGG(gene     = A3SS_gene_1_enid,
                        organism = 'hsa',
                        universe = background_enid)
head(A3SS_mkk)



A3SS_do <- DOSE::enrichDO(gene    	  = all_A3SS_gene_enid, 
                        ont           = "DO", 
                        pvalueCutoff  = 0.05, 
                        pAdjustMethod = "BH",
                        #universe 	  = background_enid,
                        minGSSize     = 5, 
                        maxGSSize     = 500, 
                        qvalueCutoff  = 0.05,
                        readable      = FALSE)
A3SS_0_do <- DOSE::enrichDO(gene    	    = A3SS_gene_0_enid, 
                          ont           = "DO", 
                          pvalueCutoff  = 0.05, 
                          pAdjustMethod = "BH",
                          #universe 	 	= background_enid,
                          minGSSize     = 5, 
                          maxGSSize     = 500, 
                          qvalueCutoff  = 0.05,
                          readable      = FALSE)
A3SS_1_do <- DOSE::enrichDO(gene    	    = A3SS_gene_1_enid, 
                          ont           = "DO", 
                          pvalueCutoff  = 0.05, 
                          pAdjustMethod = "BH",
                          #universe 	 	= background_enid,
                          minGSSize     = 5, 
                          maxGSSize     = 500, 
                          qvalueCutoff  = 0.05,
                          readable      = FALSE)
head(A3SS_do)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A3SS/all_A3SS_gene_DO.pdf",height = 10,width = 15)
print(dotplot(A3SS_do,showCategory = 10,font.size = 16, title = "DO"))
enrichMap(A3SS_do,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A3SS/Pseudo_A3SS_gene_DO.pdf",height = 10,width = 15)
print(dotplot(A3SS_0_do,showCategory = 10,font.size = 16, title = "DO"))
enrichMap(A3SS_0_do,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A3SS/Real_A3SS_gene_DO.pdf",height = 10,width = 15)
print(dotplot(A3SS_1_do,showCategory = 10,font.size = 16, title = "DO"))
enrichMap(A3SS_1_do,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


dgn <- enrichDGN(gene         = all_A3SS_gene_enid,
                 #universe     = background_enid,
                 qvalueCutoff = 0.05,
                 readable     = TRUE)
dgn_0 <- enrichDGN(gene         = A3SS_gene_0_enid,
                   #universe     = background_enid,
                   qvalueCutoff = 0.05,
                   readable     = TRUE)
dgn_1 <- enrichDGN(gene         = A3SS_gene_1_enid,
                   #universe     = background_enid,
                   qvalueCutoff = 0.05,
                   readable     = TRUE)
head(dgn)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A3SS/all_A3SS_gene_DGN.pdf",height = 10,width = 15)
print(dotplot(dgn,showCategory = 10,font.size = 16, title = "DO"))
enrichMap(dgn,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A3SS/Pseudo_A3SS_gene_DGN.pdf",height = 10,width = 15)
print(dotplot(dgn_0,showCategory = 10,font.size = 16, title = "DO"))
enrichMap(dgn_0,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A3SS/Real_A3SS_gene_DGN.pdf",height = 10,width = 15)
print(dotplot(dgn_1,showCategory = 10,font.size = 16, title = "DO"))
enrichMap(dgn_1,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()





library(ReactomePA)
A3SS_Path <- enrichPathway(gene = all_A3SS_gene_enid, 
                         organism = "human", 
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH", 
                         qvalueCutoff = 0.2, 
                         universe = background_enid, 
                         minGSSize = 10,
                         maxGSSize = 500, 
                         readable = FALSE)

A3SS_0_Path <- enrichPathway(gene = A3SS_gene_0_enid, 
                           organism = "human", 
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH", 
                           qvalueCutoff = 0.2, 
                           universe = background_enid, 
                           minGSSize = 10,
                           maxGSSize = 500, 
                           readable = FALSE)

A3SS_1_Path <- enrichPathway(gene = A3SS_gene_1_enid, 
                           organism = "human", 
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH", 
                           qvalueCutoff = 0.2, 
                           universe = background_enid, 
                           minGSSize = 10,
                           maxGSSize = 500, 
                           readable = FALSE)



A3SS_BP_ego_1
A3SS_1_KEGG
dgn_1
print(dotplot(A3SS_BP_ego_1,showCategory = 10,font.size = 16,title = "BP"))
print(dotplot(A3SS_1_KEGG,showCategory = 10,font.size = 16,title = "BP"))
print(dotplot(dgn_1,showCategory = 10,font.size = 16,title = "BP"))




pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A3SS/all_A3SS_gene_KEGG.pdf",height = 10,width = 15)
print(dotplot(A3SS_KEGG,showCategory = 10,font.size = 16, title = "BP"))
#suppressMessages(plotGOgraph(A3SS_KEGG,firstSigNodes = 15, useInfo = "all"))
enrichMap(A3SS_KEGG,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A3SS/Pseudo_A3SS_gene_KEGG.pdf",height = 10,width = 15)
print(dotplot(A3SS_0_KEGG,showCategory = 10,font.size = 16, title = "BP"))
enrichMap(A3SS_0_KEGG,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A3SS/Real_A3SS_gene_KEGG.pdf",height = 10,width = 15)
print(dotplot(A3SS_1_KEGG,showCategory = 10,font.size = 16, title = "BP"))
enrichMap(A3SS_1_KEGG,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A3SS/Real_A3SS_gene_DGN.pdf",height = 10,width = 15)
print(dotplot(dgn_1,showCategory = 10,font.size = 16, title = "BP"))
#suppressMessages(plotGOgraph(A3SS_KEGG,firstSigNodes = 15, useInfo = "all"))
enrichMap(dgn_1,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A3SS/all_A3SS_gene_Pathway.pdf",height = 10,width = 15)
print(dotplot(A3SS_Path,showCategory = 10,font.size = 16, title = "BP"))
#suppressMessages(plotGOgraph(A3SS_Path,firstSigNodes = 15, useInfo = "all"))
enrichMap(A3SS_Path,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A3SS/Pseudo_A3SS_gene_Pathway.pdf",height = 10,width = 15)
print(dotplot(A3SS_0_Path,showCategory = 10,font.size = 16, title = "BP"))
enrichMap(A3SS_0_Path,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A3SS/Real_A3SS_gene_Pathway.pdf",height = 10,width = 15)
print(dotplot(A3SS_1_Path,showCategory = 10,font.size = 16, title = "BP"))
enrichMap(A3SS_1_Path,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()



CompareGO_BP <- compareCluster(geneCluster   = list(All = all_A3SS_gene, Pseudo = A3SS_gene_0, Real = A3SS_gene_1), 
                               universe      = background,
                               fun 		  ="enrichGO", 
                               OrgDb  		  = org.Hs.eg.db, 
                               keyType 	  = 'ENSEMBL', 
                               ont  		  = "BP", 
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)
CompareGO_MF <- compareCluster(geneCluster   = list(All = all_A3SS_gene, Pseudo = A3SS_gene_0, Real = A3SS_gene_1), 
                               universe      = background,
                               fun 		  ="enrichGO", 
                               OrgDb  		  = org.Hs.eg.db, 
                               keyType 	  = 'ENSEMBL', 
                               ont  		  = "MF", 
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)
CompareGO_CC <- compareCluster(geneCluster   = list(All = all_A3SS_gene, Pseudo = A3SS_gene_0, Real = A3SS_gene_1), 
                               universe      = background,
                               fun 		  ="enrichGO", 
                               OrgDb  		  = org.Hs.eg.db, 
                               keyType 	  = 'ENSEMBL', 
                               ont  		  = "CC", 
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A3SS/A3SS_CompareGO.pdf", height = 8, width = 20)
print(dotplot(CompareGO_BP,showCategory = 10,font.size = 16,title = "A3SS GO BP"))
print(dotplot(CompareGO_MF,showCategory = 10,font.size = 16,title = "A3SS GO MF"))
print(dotplot(CompareGO_CC,showCategory = 10,font.size = 16,title = "A3SS GO CC"))
dev.off()





CompareKEGG <- compareCluster(geneCluster  = list(All = all_A3SS_gene_enid, Pseudo = A3SS_gene_0_enid, Real = A3SS_gene_1_enid),
                              universe 	   = background_enid,
                              fun          = "enrichKEGG",
                              organism	   = "hsa", 
                              pvalueCutoff = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A3SS/A3SS_CompareKEGG.pdf", height = 8, width = 20)
print(dotplot(CompareKEGG,showCategory = 10,font.size = 16,title = "A3SS KEGG"))
dev.off()



CompareDO <- compareCluster(geneCluster   = list(All = all_A3SS_gene_enid, Pseudo = A3SS_gene_0_enid, Real = A3SS_gene_1_enid),
                            fun 		  = "enrichDO",
                            ont           = "DO", 
                            pvalueCutoff  = 0.05, 
                            pAdjustMethod = "BH",
                            universe 	  = background_enid,
                            minGSSize     = 5, 
                            maxGSSize     = 500, 
                            qvalueCutoff  = 0.05,
                            readable      = FALSE)


library(ReactomePA)
ComparePathway <- compareCluster(geneCluster   = list(All = all_A3SS_gene_enid, Pseudo = A3SS_gene_0_enid, Real = A3SS_gene_1_enid),
                                 fun="enrichPathway",
                                 organism = "human",
                                 pvalueCutoff = 0.05, 
                                 pAdjustMethod = "BH",
                                 qvalueCutoff = 0.2,
                                 universe = background_enid,
                                 minGSSize = 10,
                                 maxGSSize = 500,
                                 readable = FALSE)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A3SS/A3SS_ComparePathway.pdf", height = 8, width = 20)
print(dotplot(ComparePathway,showCategory = 10,font.size = 16,title = "A3SS Pathway"))
dev.off()




















#### A5SS =======================================================================================================================================



load("/mnt/data5/BGI/UCB/tangchao/DSU2/A5SS/A5SS_psi_new.RData")

PSI <- A5SS_psi[, row.names(Cell_type)]
row.names(PSI) <- A5SS_psi$as2
rs <- rowSums(!is.na(PSI))
summary(rs)

PSI <- PSI[rs!=0, ]

rs <- rowSums(!is.na(PSI) & PSI != 0 & PSI != 1)
summary(rs)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.000   0.000   0.000   3.632   1.000 960.000

sum(rs==0)
#[1] 2107
sum(rs==1)
#[1] 470
sum(rs>=1)
#[1] 923

all_A5SS_gene <- unique(A5SS_psi[A5SS_psi$as2 %in% row.names(PSI), "gene_id"])
length(all_A5SS_gene)
#[1] 2416
A5SS_gene_0 <- unique(A5SS_psi[A5SS_psi$as2 %in% names(rs)[rs==0], "gene_id"])
A5SS_gene_1 <- unique(A5SS_psi[A5SS_psi$as2 %in% names(rs)[rs>=1], "gene_id"])
length(A5SS_gene_0)
#[1] 1790
length(A5SS_gene_1)
#[1] 809


#### GO

A5SS_BP_ego <- enrichGO(gene          = all_A5SS_gene,
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
A5SS_MF_ego <- enrichGO(gene          = all_A5SS_gene,
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
A5SS_CC_ego <- enrichGO(gene          = all_A5SS_gene,
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A5SS/all_A5SS_gene_GO.pdf",height = 10,width = 15)
print(dotplot(A5SS_BP_ego,showCategory = 10,font.size = 16,title = "BP"))
print(dotplot(A5SS_MF_ego,showCategory = 10,font.size = 16,title = "MF"))
print(dotplot(A5SS_CC_ego,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(A5SS_BP_ego,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(A5SS_MF_ego,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(A5SS_CC_ego,firstSigNodes = 15, useInfo = "all"))
enrichMap(A5SS_BP_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(A5SS_MF_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(A5SS_CC_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()



A5SS_BP_ego_0 <- enrichGO(gene          = A5SS_gene_0,
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "BP", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
A5SS_MF_ego_0 <- enrichGO(gene          = A5SS_gene_0,
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "MF", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
A5SS_CC_ego_0 <- enrichGO(gene          = A5SS_gene_0,
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "CC", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A5SS/Pseudo_A5SS_gene_GO.pdf",height = 10,width = 15)
print(dotplot(A5SS_BP_ego_0,showCategory = 10,font.size = 16,title = "BP"))
print(dotplot(A5SS_MF_ego_0,showCategory = 10,font.size = 16,title = "MF"))
print(dotplot(A5SS_CC_ego_0,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(A5SS_BP_ego_0,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(A5SS_MF_ego_0,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(A5SS_CC_ego_0,firstSigNodes = 15, useInfo = "all"))
enrichMap(A5SS_BP_ego_0,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(A5SS_MF_ego_0,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(A5SS_CC_ego_0,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


A5SS_BP_ego_1 <- enrichGO(gene          = A5SS_gene_1,
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "BP", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
A5SS_MF_ego_1 <- enrichGO(gene          = A5SS_gene_1,
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "MF", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
A5SS_CC_ego_1 <- enrichGO(gene          = A5SS_gene_1,
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "CC", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A5SS/Real_A5SS_gene_GO.pdf",height = 10,width = 15)
print(dotplot(A5SS_BP_ego_1,showCategory = 10,font.size = 16,title = "BP"))
print(dotplot(A5SS_MF_ego_1,showCategory = 10,font.size = 16,title = "MF"))
print(dotplot(A5SS_CC_ego_1,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(A5SS_BP_ego_1,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(A5SS_MF_ego_1,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(A5SS_CC_ego_1,firstSigNodes = 15, useInfo = "all"))
enrichMap(A5SS_BP_ego_1,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(A5SS_MF_ego_1,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(A5SS_CC_ego_1,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()



background_enid <- as.character(na.omit(mapIds(org.Hs.eg.db, keys = background, column = "ENTREZID", keytype = "ENSEMBL")))
all_A5SS_gene_enid <- as.character(na.omit(mapIds(org.Hs.eg.db, keys = all_A5SS_gene, column = "ENTREZID", keytype = "ENSEMBL")))
all_A5SS_gene_enid <- all_A5SS_gene_enid[!is.na(all_A5SS_gene_enid)]

keytypes(org.Hs.eg.db)
background_enid <- bitr(background, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
all_A5SS_gene_enid <- bitr(all_A5SS_gene, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID

A5SS_gene_0_enid <- bitr(A5SS_gene_0, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
A5SS_gene_1_enid <- bitr(A5SS_gene_1, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID


A5SS_KEGG <- enrichKEGG(gene         = all_A5SS_gene_enid,
                        organism     = 'hsa',
                        universe	 = background_enid,
                        pvalueCutoff = 0.05)
head(A5SS_KEGG)

A5SS_0_KEGG <- enrichKEGG(gene         = A5SS_gene_0_enid,
                          organism     = 'hsa',
                          universe	   = background_enid,
                          pvalueCutoff = 0.05)

A5SS_1_KEGG <- enrichKEGG(gene         = A5SS_gene_1_enid,
                          organism     = 'hsa',
                          universe	   = background_enid,
                          pvalueCutoff = 0.05)



A5SS_mkk <- enrichMKEGG(gene     = all_A5SS_gene_enid,
                        organism = 'hsa',
                        universe = background_enid)
A5SS_0_mkk <- enrichMKEGG(gene     = A5SS_gene_0_enid,
                          organism = 'hsa',
                          universe = background_enid)
A5SS_1_mkk <- enrichMKEGG(gene     = A5SS_gene_1_enid,
                          organism = 'hsa',
                          universe = background_enid)
head(A5SS_mkk)



A5SS_do <- DOSE::enrichDO(gene    	  = all_A5SS_gene_enid, 
                          ont           = "DO", 
                          pvalueCutoff  = 0.05, 
                          pAdjustMethod = "BH",
                          #universe 	  = background_enid,
                          minGSSize     = 5, 
                          maxGSSize     = 500, 
                          qvalueCutoff  = 0.05,
                          readable      = FALSE)
A5SS_0_do <- DOSE::enrichDO(gene    	    = A5SS_gene_0_enid, 
                            ont           = "DO", 
                            pvalueCutoff  = 0.05, 
                            pAdjustMethod = "BH",
                            #universe 	 	= background_enid,
                            minGSSize     = 5, 
                            maxGSSize     = 500, 
                            qvalueCutoff  = 0.05,
                            readable      = FALSE)
A5SS_1_do <- DOSE::enrichDO(gene    	    = A5SS_gene_1_enid, 
                            ont           = "DO", 
                            pvalueCutoff  = 0.05, 
                            pAdjustMethod = "BH",
                            #universe 	 	= background_enid,
                            minGSSize     = 5, 
                            maxGSSize     = 500, 
                            qvalueCutoff  = 0.05,
                            readable      = FALSE)
head(A5SS_do)
dgn <- enrichDGN(gene         = all_A5SS_gene_enid,
                 #universe     = background_enid,
                 qvalueCutoff = 0.05,
                 readable     = TRUE)
dgn_0 <- enrichDGN(gene         = A5SS_gene_0_enid,
                   #universe     = background_enid,
                   qvalueCutoff = 0.05,
                   readable     = TRUE)
dgn_1 <- enrichDGN(gene         = A5SS_gene_1_enid,
                   #universe     = background_enid,
                   qvalueCutoff = 0.05,
                   readable     = TRUE)

head(dgn)

library(ReactomePA)
A5SS_Path <- enrichPathway(gene          = all_A5SS_gene_enid, 
                           organism      = "human", 
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH", 
                           qvalueCutoff  = 0.2, 
                           universe      = background_enid, 
                           minGSSize     = 10,
                           maxGSSize     = 500, 
                           readable      = FALSE)

A5SS_0_Path <- enrichPathway(gene          = A5SS_gene_0_enid, 
                             organism      = "human", 
                             pvalueCutoff  = 0.05,
                             pAdjustMethod = "BH", 
                             qvalueCutoff  = 0.2, 
                             universe      = background_enid, 
                             minGSSize     = 10,
                             maxGSSize     = 500, 
                             readable      = FALSE)

A5SS_1_Path <- enrichPathway(gene          = A5SS_gene_1_enid, 
                             organism      = "human", 
                             pvalueCutoff  = 0.05,
                             pAdjustMethod = "BH", 
                             qvalueCutoff  = 0.2, 
                             universe      = background_enid, 
                             minGSSize     = 10,
                             maxGSSize     = 500, 
                             readable      = FALSE)



A5SS_BP_ego_1
A5SS_1_KEGG
dgn_1
print(dotplot(A5SS_BP_ego_1,showCategory = 10,font.size = 16,title = "BP"))
print(dotplot(A5SS_1_KEGG,showCategory = 10,font.size = 16,title = "BP"))
print(dotplot(dgn_1,showCategory = 10,font.size = 16,title = "BP"))




pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A5SS/all_A5SS_gene_KEGG.pdf",height = 10,width = 15)
print(dotplot(A5SS_KEGG,showCategory = 10,font.size = 16, title = "BP"))
#suppressMessages(plotGOgraph(A5SS_KEGG,firstSigNodes = 15, useInfo = "all"))
enrichMap(A5SS_KEGG,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A5SS/Pseudo_A5SS_gene_KEGG.pdf",height = 10,width = 15)
print(dotplot(A5SS_0_KEGG,showCategory = 10,font.size = 16, title = "BP"))
enrichMap(A5SS_0_KEGG,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A5SS/Real_A5SS_gene_KEGG.pdf",height = 10,width = 15)
print(dotplot(A5SS_1_KEGG,showCategory = 10,font.size = 16, title = "BP"))
enrichMap(A5SS_1_KEGG,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A5SS/Real_A5SS_gene_DGN.pdf",height = 10,width = 15)
print(dotplot(dgn_1,showCategory = 10,font.size = 16, title = "BP"))
#suppressMessages(plotGOgraph(A5SS_KEGG,firstSigNodes = 15, useInfo = "all"))
enrichMap(dgn_1,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A5SS/all_A5SS_gene_Pathway.pdf",height = 10,width = 15)
print(dotplot(A5SS_Path,showCategory = 10,font.size = 16, title = "BP"))
#suppressMessages(plotGOgraph(A5SS_Path,firstSigNodes = 15, useInfo = "all"))
enrichMap(A5SS_Path,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A5SS/Pseudo_A5SS_gene_Pathway.pdf",height = 10,width = 15)
print(dotplot(A5SS_0_Path,showCategory = 10,font.size = 16, title = "BP"))
enrichMap(A5SS_0_Path,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A5SS/Real_A5SS_gene_Pathway.pdf",height = 10,width = 15)
print(dotplot(A5SS_1_Path,showCategory = 10,font.size = 16, title = "BP"))
enrichMap(A5SS_1_Path,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()



CompareGO_BP <- compareCluster(geneCluster   = list(All = all_A5SS_gene, Pseudo = A5SS_gene_0, Real = A5SS_gene_1), 
                               universe      = background,
                               fun 		  ="enrichGO", 
                               OrgDb  		  = org.Hs.eg.db, 
                               keyType 	  = 'ENSEMBL', 
                               ont  		  = "BP", 
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)
CompareGO_MF <- compareCluster(geneCluster   = list(All = all_A5SS_gene, Pseudo = A5SS_gene_0, Real = A5SS_gene_1), 
                               universe      = background,
                               fun 		  ="enrichGO", 
                               OrgDb  		  = org.Hs.eg.db, 
                               keyType 	  = 'ENSEMBL', 
                               ont  		  = "MF", 
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)
CompareGO_CC <- compareCluster(geneCluster   = list(All = all_A5SS_gene, Pseudo = A5SS_gene_0, Real = A5SS_gene_1), 
                               universe      = background,
                               fun 		  ="enrichGO", 
                               OrgDb  		  = org.Hs.eg.db, 
                               keyType 	  = 'ENSEMBL', 
                               ont  		  = "CC", 
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A5SS/A5SS_CompareGO.pdf", height = 8, width = 20)
print(dotplot(CompareGO_BP,showCategory = 10,font.size = 16,title = "A5SS GO BP"))
print(dotplot(CompareGO_MF,showCategory = 10,font.size = 16,title = "A5SS GO MF"))
print(dotplot(CompareGO_CC,showCategory = 10,font.size = 16,title = "A5SS GO CC"))
dev.off()




CompareKEGG <- compareCluster(geneCluster  = list(All = all_A5SS_gene_enid, Pseudo = A5SS_gene_0_enid, Real = A5SS_gene_1_enid),
                              universe 	   = background_enid,
                              fun          = "enrichKEGG",
                              organism	   = "hsa", 
                              pvalueCutoff = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A5SS/A5SS_CompareKEGG.pdf", height = 8, width = 20)
print(dotplot(CompareKEGG,showCategory = 10,font.size = 16,title = "A5SS KEGG"))
dev.off()



CompareDO <- compareCluster(geneCluster   = list(All = all_A5SS_gene_enid, Pseudo = A5SS_gene_0_enid, Real = A5SS_gene_1_enid),
                            fun 		  = "enrichDO",
                            ont           = "DO", 
                            pvalueCutoff  = 0.05, 
                            pAdjustMethod = "BH",
                            universe 	  = background_enid,
                            minGSSize     = 5, 
                            maxGSSize     = 500, 
                            qvalueCutoff  = 0.05,
                            readable      = FALSE)


library(ReactomePA)
ComparePathway <- compareCluster(geneCluster   = list(All = all_A5SS_gene_enid, Pseudo = A5SS_gene_0_enid, Real = A5SS_gene_1_enid),
                                 fun="enrichPathway",
                                 organism = "human",
                                 pvalueCutoff = 0.05, 
                                 pAdjustMethod = "BH",
                                 qvalueCutoff = 0.2,
                                 universe = background_enid,
                                 minGSSize = 10,
                                 maxGSSize = 500,
                                 readable = FALSE)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/A5SS/A5SS_ComparePathway.pdf", height = 8, width = 20)
print(dotplot(ComparePathway,showCategory = 10,font.size = 16,title = "A5SS Pathway"))
dev.off()





























#### AFE =======================================================================================================================================



load("/mnt/data5/BGI/UCB/tangchao/DSU2/AFE/AFE_psi_new.RData")

PSI <- AFE_psi[, row.names(Cell_type)]
row.names(PSI) <- AFE_psi$as2
rs <- rowSums(!is.na(PSI))
summary(rs)

PSI <- PSI[rs!=0, ]

rs <- rowSums(!is.na(PSI) & PSI != 0 & PSI != 1)
summary(rs)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  	   0       0       1      16       5    2600

sum(rs==0)
#[1] 1093
sum(rs==1)
#[1] 537
sum(rs>=1)
#[1] 1780

all_AFE_gene <- unique(AFE_psi[AFE_psi$as2 %in% row.names(PSI), "gene_id"])
length(all_AFE_gene)
#[1] 2251
all_AFE_gene <- all_AFE_gene[mapply(nchar,all_AFE_gene)==15]
length(all_AFE_gene)
#[1] 2250

AFE_gene_0 <- unique(AFE_psi[AFE_psi$as2 %in% names(rs)[rs==0], "gene_id"])
AFE_gene_1 <- unique(AFE_psi[AFE_psi$as2 %in% names(rs)[rs>=1], "gene_id"])
length(AFE_gene_0)
#[1] 1055
length(AFE_gene_1)
#[1] 1258
AFE_gene_1 <- AFE_gene_1[mapply(nchar,AFE_gene_1)==15]
length(AFE_gene_1)
#[1] 1257

#### GO

AFE_BP_ego <- enrichGO(gene          = all_AFE_gene,
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
AFE_MF_ego <- enrichGO(gene          = all_AFE_gene,
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
AFE_CC_ego <- enrichGO(gene          = all_AFE_gene,
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/AFE/all_AFE_gene_GO.pdf",height = 10,width = 15)
print(dotplot(AFE_BP_ego,showCategory = 10,font.size = 16,title = "BP"))
print(dotplot(AFE_MF_ego,showCategory = 10,font.size = 16,title = "MF"))
print(dotplot(AFE_CC_ego,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(AFE_BP_ego,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(AFE_MF_ego,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(AFE_CC_ego,firstSigNodes = 15, useInfo = "all"))
enrichMap(AFE_BP_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(AFE_MF_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(AFE_CC_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()



AFE_BP_ego_0 <- enrichGO(gene          = AFE_gene_0,
                          universe      = background,
                          OrgDb         = org.Hs.eg.db,
                          keyType       = 'ENSEMBL',
                          ont           = "BP", 
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05)
AFE_MF_ego_0 <- enrichGO(gene          = AFE_gene_0,
                          universe      = background,
                          OrgDb         = org.Hs.eg.db,
                          keyType       = 'ENSEMBL',
                          ont           = "MF", 
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05)
AFE_CC_ego_0 <- enrichGO(gene          = AFE_gene_0,
                          universe      = background,
                          OrgDb         = org.Hs.eg.db,
                          keyType       = 'ENSEMBL',
                          ont           = "CC", 
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/AFE/Pseudo_AFE_gene_GO.pdf",height = 10,width = 15)
print(dotplot(AFE_BP_ego_0,showCategory = 10,font.size = 16,title = "BP"))
print(dotplot(AFE_MF_ego_0,showCategory = 10,font.size = 16,title = "MF"))
print(dotplot(AFE_CC_ego_0,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(AFE_BP_ego_0,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(AFE_MF_ego_0,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(AFE_CC_ego_0,firstSigNodes = 15, useInfo = "all"))
enrichMap(AFE_BP_ego_0,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(AFE_MF_ego_0,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(AFE_CC_ego_0,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


AFE_BP_ego_1 <- enrichGO(gene          = AFE_gene_1,
                          universe      = background,
                          OrgDb         = org.Hs.eg.db,
                          keyType       = 'ENSEMBL',
                          ont           = "BP", 
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05)
AFE_MF_ego_1 <- enrichGO(gene          = AFE_gene_1,
                          universe      = background,
                          OrgDb         = org.Hs.eg.db,
                          keyType       = 'ENSEMBL',
                          ont           = "MF", 
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05)
AFE_CC_ego_1 <- enrichGO(gene          = AFE_gene_1,
                          universe      = background,
                          OrgDb         = org.Hs.eg.db,
                          keyType       = 'ENSEMBL',
                          ont           = "CC", 
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/AFE/Real_AFE_gene_GO.pdf",height = 10,width = 15)
print(dotplot(AFE_BP_ego_1,showCategory = 10,font.size = 16,title = "BP"))
print(dotplot(AFE_MF_ego_1,showCategory = 10,font.size = 16,title = "MF"))
print(dotplot(AFE_CC_ego_1,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(AFE_BP_ego_1,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(AFE_MF_ego_1,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(AFE_CC_ego_1,firstSigNodes = 15, useInfo = "all"))
enrichMap(AFE_BP_ego_1,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(AFE_MF_ego_1,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(AFE_CC_ego_1,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()



background_enid <- as.character(na.omit(mapIds(org.Hs.eg.db, keys = background, column = "ENTREZID", keytype = "ENSEMBL")))
all_AFE_gene_enid <- as.character(na.omit(mapIds(org.Hs.eg.db, keys = all_AFE_gene, column = "ENTREZID", keytype = "ENSEMBL")))
all_AFE_gene_enid <- all_AFE_gene_enid[!is.na(all_AFE_gene_enid)]

keytypes(org.Hs.eg.db)
background_enid <- bitr(background, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
all_AFE_gene_enid <- bitr(all_AFE_gene, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID

AFE_gene_0_enid <- bitr(AFE_gene_0, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
AFE_gene_1_enid <- bitr(AFE_gene_1, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID


AFE_KEGG <- enrichKEGG(gene         = all_AFE_gene_enid,
                        organism     = 'hsa',
                        universe	   = background_enid,
                        pvalueCutoff = 0.05)
head(AFE_KEGG)

AFE_0_KEGG <- enrichKEGG(gene         = AFE_gene_0_enid,
                          organism     = 'hsa',
                          universe	 = background_enid,
                          pvalueCutoff = 0.05)

AFE_1_KEGG <- enrichKEGG(gene         = AFE_gene_1_enid,
                          organism     = 'hsa',
                          universe	 = background_enid,
                          pvalueCutoff = 0.05)



AFE_mkk <- enrichMKEGG(gene     = all_AFE_gene_enid,
                        organism = 'hsa',
                        universe = background_enid)
AFE_0_mkk <- enrichMKEGG(gene     = AFE_gene_0_enid,
                          organism = 'hsa',
                          universe = background_enid)
AFE_1_mkk <- enrichMKEGG(gene     = AFE_gene_1_enid,
                          organism = 'hsa',
                          universe = background_enid)
head(AFE_mkk)



AFE_do <- DOSE::enrichDO(gene    	  = all_AFE_gene_enid, 
                          ont           = "DO", 
                          pvalueCutoff  = 0.05, 
                          pAdjustMethod = "BH",
                          #universe 	  = background_enid,
                          minGSSize     = 5, 
                          maxGSSize     = 500, 
                          qvalueCutoff  = 0.05,
                          readable      = FALSE)
AFE_0_do <- DOSE::enrichDO(gene    	    = AFE_gene_0_enid, 
                            ont           = "DO", 
                            pvalueCutoff  = 0.05, 
                            pAdjustMethod = "BH",
                            #universe 	 	= background_enid,
                            minGSSize     = 5, 
                            maxGSSize     = 500, 
                            qvalueCutoff  = 0.05,
                            readable      = FALSE)
AFE_1_do <- DOSE::enrichDO(gene    	    = AFE_gene_1_enid, 
                            ont           = "DO", 
                            pvalueCutoff  = 0.05, 
                            pAdjustMethod = "BH",
                            #universe 	 	= background_enid,
                            minGSSize     = 5, 
                            maxGSSize     = 500, 
                            qvalueCutoff  = 0.05,
                            readable      = FALSE)
head(AFE_do)
dgn <- enrichDGN(gene         = all_AFE_gene_enid,
                 #universe     = background_enid,
                 qvalueCutoff = 0.05,
                 readable     = TRUE)
dgn_0 <- enrichDGN(gene         = AFE_gene_0_enid,
                   #universe     = background_enid,
                   qvalueCutoff = 0.05,
                   readable     = TRUE)
dgn_1 <- enrichDGN(gene         = AFE_gene_1_enid,
                   #universe     = background_enid,
                   qvalueCutoff = 0.05,
                   readable     = TRUE)
head(dgn)
pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/AFE/all_AFE_gene_DGN.pdf",height = 10,width = 15)
print(dotplot(dgn,showCategory = 10,font.size = 16, title = "BP"))
enrichMap(dgn,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()
pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/AFE/Pseudo_AFE_gene_DGN.pdf",height = 10,width = 15)
print(dotplot(dgn_0,showCategory = 10,font.size = 16, title = "BP"))
enrichMap(dgn_0,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()
pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/AFE/Real_AFE_gene_DGN.pdf",height = 10,width = 15)
print(dotplot(dgn_1,showCategory = 10,font.size = 16, title = "BP"))
enrichMap(dgn_1,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


library(ReactomePA)
AFE_Path <- enrichPathway(gene = all_AFE_gene_enid, 
                           organism = "human", 
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH", 
                           qvalueCutoff = 0.2, 
                           universe = background_enid, 
                           minGSSize = 10,
                           maxGSSize = 500, 
                           readable = FALSE)

AFE_0_Path <- enrichPathway(gene = AFE_gene_0_enid, 
                             organism = "human", 
                             pvalueCutoff = 0.05,
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.2, 
                             universe = background_enid, 
                             minGSSize = 10,
                             maxGSSize = 500, 
                             readable = FALSE)

AFE_1_Path <- enrichPathway(gene = AFE_gene_1_enid, 
                             organism = "human", 
                             pvalueCutoff = 0.05,
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.2, 
                             universe = background_enid, 
                             minGSSize = 10,
                             maxGSSize = 500, 
                             readable = FALSE)



AFE_BP_ego_1
AFE_1_KEGG
dgn_1
print(dotplot(AFE_BP_ego_1,showCategory = 10,font.size = 16,title = "BP"))
print(dotplot(AFE_1_KEGG,showCategory = 10,font.size = 16,title = "BP"))
print(dotplot(dgn_1,showCategory = 10,font.size = 16,title = "BP"))




pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/AFE/all_AFE_gene_KEGG.pdf",height = 10,width = 15)
print(dotplot(AFE_KEGG,showCategory = 10,font.size = 16, title = "BP"))
#suppressMessages(plotGOgraph(AFE_KEGG,firstSigNodes = 15, useInfo = "all"))
enrichMap(AFE_KEGG,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/AFE/Pseudo_AFE_gene_KEGG.pdf",height = 10,width = 15)
print(dotplot(AFE_0_KEGG,showCategory = 10,font.size = 16, title = "BP"))
enrichMap(AFE_0_KEGG,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/AFE/Real_AFE_gene_KEGG.pdf",height = 10,width = 15)
print(dotplot(AFE_1_KEGG,showCategory = 10,font.size = 16, title = "BP"))
enrichMap(AFE_1_KEGG,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/AFE/all_AFE_gene_Pathway.pdf",height = 10,width = 15)
print(dotplot(AFE_Path,showCategory = 10,font.size = 16, title = "BP"))
#suppressMessages(plotGOgraph(AFE_Path,firstSigNodes = 15, useInfo = "all"))
enrichMap(AFE_Path,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/AFE/Pseudo_AFE_gene_Pathway.pdf",height = 10,width = 15)
print(dotplot(AFE_0_Path,showCategory = 10,font.size = 16, title = "BP"))
enrichMap(AFE_0_Path,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/AFE/Real_AFE_gene_Pathway.pdf",height = 10,width = 15)
print(dotplot(AFE_1_Path,showCategory = 10,font.size = 16, title = "BP"))
enrichMap(AFE_1_Path,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()



CompareGO_BP <- compareCluster(geneCluster   = list(All = all_AFE_gene, Pseudo = AFE_gene_0, Real = AFE_gene_1), 
                            universe      = background,
                            fun 		  ="enrichGO", 
                            OrgDb  		  = org.Hs.eg.db, 
                            keyType 	  = 'ENSEMBL', 
                            ont  		  = "BP", 
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05)
CompareGO_MF <- compareCluster(geneCluster   = list(All = all_AFE_gene, Pseudo = AFE_gene_0, Real = AFE_gene_1), 
                            universe      = background,
                            fun 		  ="enrichGO", 
                            OrgDb  		  = org.Hs.eg.db, 
                            keyType 	  = 'ENSEMBL', 
                            ont  		  = "MF", 
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05)
CompareGO_CC <- compareCluster(geneCluster   = list(All = all_AFE_gene, Pseudo = AFE_gene_0, Real = AFE_gene_1), 
                            universe      = background,
                            fun 		  ="enrichGO", 
                            OrgDb  		  = org.Hs.eg.db, 
                            keyType 	  = 'ENSEMBL', 
                            ont  		  = "CC", 
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/AFE/AFE_CompareGO.pdf", height = 8, width = 20)
print(dotplot(CompareGO_BP,showCategory = 10,font.size = 16,title = "AFE GO BP"))
print(dotplot(CompareGO_MF,showCategory = 10,font.size = 16,title = "AFE GO MF"))
print(dotplot(CompareGO_CC,showCategory = 10,font.size = 16,title = "AFE GO CC"))
dev.off()





CompareKEGG <- compareCluster(geneCluster  = list(All = all_AFE_gene_enid, Pseudo = AFE_gene_0_enid, Real = AFE_gene_1_enid),
                              universe 	   = background_enid,
                              fun          = "enrichKEGG",
                              organism	   = "hsa", 
                              pvalueCutoff = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/AFE/AFE_CompareKEGG.pdf", height = 8, width = 20)
print(dotplot(CompareKEGG,showCategory = 10,font.size = 16,title = "AFE KEGG"))
dev.off()



CompareDO <- compareCluster(geneCluster   = list(All = all_AFE_gene_enid, Pseudo = AFE_gene_0_enid, Real = AFE_gene_1_enid),
                            fun 		  = "enrichDO",
                            ont           = "DO", 
                            pvalueCutoff  = 0.05, 
                            pAdjustMethod = "BH",
                            #universe 	  = background_enid,
                            minGSSize     = 5, 
                            maxGSSize     = 500, 
                            qvalueCutoff  = 0.05,
                            readable      = FALSE)


library(ReactomePA)
ComparePathway <- compareCluster(geneCluster   = list(All = all_AFE_gene_enid, Pseudo = AFE_gene_0_enid, Real = AFE_gene_1_enid),
                                 fun="enrichPathway",
                                 organism = "human",
                                 pvalueCutoff = 0.05, 
                                 pAdjustMethod = "BH",
                                 qvalueCutoff = 0.2,
                                 universe = background_enid,
                                 minGSSize = 10,
                                 maxGSSize = 500,
                                 readable = FALSE)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/AFE/AFE_ComparePathway.pdf", height = 8, width = 20)
print(dotplot(ComparePathway,showCategory = 10,font.size = 16,title = "AFE Pathway"))
dev.off()

















#### ALE =======================================================================================================================================



load("/mnt/data5/BGI/UCB/tangchao/DSU2/ALE/ALE_psi_new.RData")

PSI <- ALE_psi[, row.names(Cell_type)]
row.names(PSI) <- ALE_psi$as2
rs <- rowSums(!is.na(PSI))
summary(rs)

PSI <- PSI[rs!=0, ]

rs <- rowSums(!is.na(PSI) & PSI != 0 & PSI != 1)
summary(rs)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.00    0.00    1.00   11.97    3.00 1646.00

sum(rs==0)
#[1] 644
sum(rs==1)
#[1] 245
sum(rs>=1)
#[1] 743

all_ALE_gene <- unique(ALE_psi[ALE_psi$as2 %in% row.names(PSI), "gene_id"])
length(all_ALE_gene)
#[1] 1228
all_ALE_gene <- all_ALE_gene[mapply(nchar,all_ALE_gene)==15]
length(all_ALE_gene)
#[1] 1228

ALE_gene_0 <- unique(ALE_psi[ALE_psi$as2 %in% names(rs)[rs==0], "gene_id"])
ALE_gene_1 <- unique(ALE_psi[ALE_psi$as2 %in% names(rs)[rs>=1], "gene_id"])
length(ALE_gene_0)
#[1] 622
length(ALE_gene_1)
#[1] 627
ALE_gene_1 <- ALE_gene_1[mapply(nchar,ALE_gene_1)==15]
length(ALE_gene_1)
#[1] 627

#### GO

ALE_BP_ego <- enrichGO(gene          = all_ALE_gene,
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
ALE_MF_ego <- enrichGO(gene          = all_ALE_gene,
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
ALE_CC_ego <- enrichGO(gene          = all_ALE_gene,
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/ALE/all_ALE_gene_GO.pdf",height = 10,width = 15)
print(dotplot(ALE_BP_ego,showCategory = 10,font.size = 16,title = "BP"))
print(dotplot(ALE_MF_ego,showCategory = 10,font.size = 16,title = "MF"))
print(dotplot(ALE_CC_ego,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(ALE_BP_ego,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(ALE_MF_ego,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(ALE_CC_ego,firstSigNodes = 15, useInfo = "all"))
enrichMap(ALE_BP_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(ALE_MF_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(ALE_CC_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()



ALE_BP_ego_0 <- enrichGO(gene          = ALE_gene_0,
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "BP", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
ALE_MF_ego_0 <- enrichGO(gene          = ALE_gene_0,
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "MF", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
ALE_CC_ego_0 <- enrichGO(gene          = ALE_gene_0,
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "CC", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/ALE/Pseudo_ALE_gene_GO.pdf",height = 10,width = 15)
print(dotplot(ALE_BP_ego_0,showCategory = 10,font.size = 16,title = "BP"))
print(dotplot(ALE_MF_ego_0,showCategory = 10,font.size = 16,title = "MF"))
print(dotplot(ALE_CC_ego_0,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(ALE_BP_ego_0,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(ALE_MF_ego_0,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(ALE_CC_ego_0,firstSigNodes = 15, useInfo = "all"))
enrichMap(ALE_BP_ego_0,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(ALE_MF_ego_0,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(ALE_CC_ego_0,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


ALE_BP_ego_1 <- enrichGO(gene          = ALE_gene_1,
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "BP", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
ALE_MF_ego_1 <- enrichGO(gene          = ALE_gene_1,
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "MF", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
ALE_CC_ego_1 <- enrichGO(gene          = ALE_gene_1,
                         universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "CC", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/ALE/Real_ALE_gene_GO.pdf",height = 10,width = 15)
print(dotplot(ALE_BP_ego_1,showCategory = 10,font.size = 16,title = "BP"))
print(dotplot(ALE_MF_ego_1,showCategory = 10,font.size = 16,title = "MF"))
print(dotplot(ALE_CC_ego_1,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(ALE_BP_ego_1,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(ALE_MF_ego_1,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(ALE_CC_ego_1,firstSigNodes = 15, useInfo = "all"))
enrichMap(ALE_BP_ego_1,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(ALE_MF_ego_1,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(ALE_CC_ego_1,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()



background_enid <- as.character(na.omit(mapIds(org.Hs.eg.db, keys = background, column = "ENTREZID", keytype = "ENSEMBL")))
all_ALE_gene_enid <- as.character(na.omit(mapIds(org.Hs.eg.db, keys = all_ALE_gene, column = "ENTREZID", keytype = "ENSEMBL")))
all_ALE_gene_enid <- all_ALE_gene_enid[!is.na(all_ALE_gene_enid)]

keytypes(org.Hs.eg.db)
background_enid <- bitr(background, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
all_ALE_gene_enid <- bitr(all_ALE_gene, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID

ALE_gene_0_enid <- bitr(ALE_gene_0, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
ALE_gene_1_enid <- bitr(ALE_gene_1, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID


ALE_KEGG <- enrichKEGG(gene         = all_ALE_gene_enid,
                       organism     = 'hsa',
                       universe	   = background_enid,
                       pvalueCutoff = 0.05)
head(ALE_KEGG)

ALE_0_KEGG <- enrichKEGG(gene         = ALE_gene_0_enid,
                         organism     = 'hsa',
                         universe	 = background_enid,
                         pvalueCutoff = 0.05)

ALE_1_KEGG <- enrichKEGG(gene         = ALE_gene_1_enid,
                         organism     = 'hsa',
                         universe	 = background_enid,
                         pvalueCutoff = 0.05)



ALE_mkk <- enrichMKEGG(gene     = all_ALE_gene_enid,
                       organism = 'hsa',
                       universe = background_enid)
ALE_0_mkk <- enrichMKEGG(gene     = ALE_gene_0_enid,
                         organism = 'hsa',
                         universe = background_enid)
ALE_1_mkk <- enrichMKEGG(gene     = ALE_gene_1_enid,
                         organism = 'hsa',
                         universe = background_enid)
head(ALE_mkk)



ALE_do <- DOSE::enrichDO(gene    	  = all_ALE_gene_enid, 
                         ont           = "DO", 
                         pvalueCutoff  = 0.05, 
                         pAdjustMethod = "BH",
                         #universe 	  = background_enid,
                         minGSSize     = 5, 
                         maxGSSize     = 500, 
                         qvalueCutoff  = 0.05,
                         readable      = FALSE)
ALE_0_do <- DOSE::enrichDO(gene    	    = ALE_gene_0_enid, 
                           ont           = "DO", 
                           pvalueCutoff  = 0.05, 
                           pAdjustMethod = "BH",
                           #universe 	 	= background_enid,
                           minGSSize     = 5, 
                           maxGSSize     = 500, 
                           qvalueCutoff  = 0.05,
                           readable      = FALSE)
ALE_1_do <- DOSE::enrichDO(gene    	    = ALE_gene_1_enid, 
                           ont           = "DO", 
                           pvalueCutoff  = 0.05, 
                           pAdjustMethod = "BH",
                           #universe 	 	= background_enid,
                           minGSSize     = 5, 
                           maxGSSize     = 500, 
                           qvalueCutoff  = 0.05,
                           readable      = FALSE)
head(ALE_do)
dgn <- enrichDGN(gene         = all_ALE_gene_enid,
                 #universe     = background_enid,
                 qvalueCutoff = 0.05,
                 readable     = TRUE)
dgn_0 <- enrichDGN(gene         = ALE_gene_0_enid,
                   #universe     = background_enid,
                   qvalueCutoff = 0.05,
                   readable     = TRUE)
dgn_1 <- enrichDGN(gene         = ALE_gene_1_enid,
                   #universe     = background_enid,
                   qvalueCutoff = 0.05,
                   readable     = TRUE)
head(dgn)

library(ReactomePA)
ALE_Path <- enrichPathway(gene = all_ALE_gene_enid, 
                          organism = "human", 
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH", 
                          qvalueCutoff = 0.2, 
                          universe = background_enid, 
                          minGSSize = 10,
                          maxGSSize = 500, 
                          readable = FALSE)

ALE_0_Path <- enrichPathway(gene = ALE_gene_0_enid, 
                            organism = "human", 
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH", 
                            qvalueCutoff = 0.2, 
                            universe = background_enid, 
                            minGSSize = 10,
                            maxGSSize = 500, 
                            readable = FALSE)

ALE_1_Path <- enrichPathway(gene = ALE_gene_1_enid, 
                            organism = "human", 
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH", 
                            qvalueCutoff = 0.2, 
                            universe = background_enid, 
                            minGSSize = 10,
                            maxGSSize = 500, 
                            readable = FALSE)



ALE_BP_ego_1
ALE_1_KEGG
dgn_1
print(dotplot(ALE_BP_ego_1,showCategory = 10,font.size = 16,title = "BP"))
print(dotplot(ALE_1_KEGG,showCategory = 10,font.size = 16,title = "BP"))
print(dotplot(dgn_1,showCategory = 10,font.size = 16,title = "BP"))




pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/ALE/all_ALE_gene_KEGG.pdf",height = 10,width = 15)
print(dotplot(ALE_KEGG,showCategory = 10,font.size = 16, title = "BP"))
#suppressMessages(plotGOgraph(ALE_KEGG,firstSigNodes = 15, useInfo = "all"))
enrichMap(ALE_KEGG,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/ALE/Pseudo_ALE_gene_KEGG.pdf",height = 10,width = 15)
print(dotplot(ALE_0_KEGG,showCategory = 10,font.size = 16, title = "BP"))
enrichMap(ALE_0_KEGG,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/ALE/Real_ALE_gene_KEGG.pdf",height = 10,width = 15)
print(dotplot(ALE_1_KEGG,showCategory = 10,font.size = 16, title = "BP"))
enrichMap(ALE_1_KEGG,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()


pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/ALE/Real_ALE_gene_DGN.pdf",height = 10,width = 15)
print(dotplot(dgn_1,showCategory = 10,font.size = 16, title = "BP"))
#suppressMessages(plotGOgraph(ALE_KEGG,firstSigNodes = 15, useInfo = "all"))
enrichMap(dgn_1,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/ALE/all_ALE_gene_Pathway.pdf",height = 10,width = 15)
print(dotplot(ALE_Path,showCategory = 10,font.size = 16, title = "BP"))
#suppressMessages(plotGOgraph(ALE_Path,firstSigNodes = 15, useInfo = "all"))
enrichMap(ALE_Path,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/ALE/Pseudo_ALE_gene_Pathway.pdf",height = 10,width = 15)
print(dotplot(ALE_0_Path,showCategory = 10,font.size = 16, title = "BP"))
enrichMap(ALE_0_Path,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/ALE/Real_ALE_gene_Pathway.pdf",height = 10,width = 15)
print(dotplot(ALE_1_Path,showCategory = 10,font.size = 16, title = "BP"))
enrichMap(ALE_1_Path,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()



CompareGO_BP <- compareCluster(geneCluster   = list(All = all_ALE_gene, Pseudo = ALE_gene_0, Real = ALE_gene_1), 
                               universe      = background,
                               fun 		  ="enrichGO", 
                               OrgDb  		  = org.Hs.eg.db, 
                               keyType 	  = 'ENSEMBL', 
                               ont  		  = "BP", 
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)
CompareGO_MF <- compareCluster(geneCluster   = list(All = all_ALE_gene, Pseudo = ALE_gene_0, Real = ALE_gene_1), 
                               universe      = background,
                               fun 		  ="enrichGO", 
                               OrgDb  		  = org.Hs.eg.db, 
                               keyType 	  = 'ENSEMBL', 
                               ont  		  = "MF", 
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)
CompareGO_CC <- compareCluster(geneCluster   = list(All = all_ALE_gene, Pseudo = ALE_gene_0, Real = ALE_gene_1), 
                               universe      = background,
                               fun 		  ="enrichGO", 
                               OrgDb  		  = org.Hs.eg.db, 
                               keyType 	  = 'ENSEMBL', 
                               ont  		  = "CC", 
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05,
                               qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/ALE/ALE_CompareGO.pdf", height = 8, width = 20)
print(dotplot(CompareGO_BP,showCategory = 10,font.size = 16,title = "ALE GO BP"))
print(dotplot(CompareGO_MF,showCategory = 10,font.size = 16,title = "ALE GO MF"))
print(dotplot(CompareGO_CC,showCategory = 10,font.size = 16,title = "ALE GO CC"))
dev.off()





CompareKEGG <- compareCluster(geneCluster  = list(All = all_ALE_gene_enid, Pseudo = ALE_gene_0_enid, Real = ALE_gene_1_enid),
                              universe 	   = background_enid,
                              fun          = "enrichKEGG",
                              organism	   = "hsa", 
                              pvalueCutoff = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/ALE/ALE_CompareKEGG.pdf", height = 8, width = 20)
print(dotplot(CompareKEGG,showCategory = 10,font.size = 16,title = "ALE KEGG"))
dev.off()



CompareDO <- compareCluster(geneCluster   = list(All = all_ALE_gene_enid, Pseudo = ALE_gene_0_enid, Real = ALE_gene_1_enid),
                            fun 		  = "enrichDO",
                            ont           = "DO", 
                            pvalueCutoff  = 0.05, 
                            pAdjustMethod = "BH",
                            #universe 	  = background_enid,
                            minGSSize     = 5, 
                            maxGSSize     = 500, 
                            qvalueCutoff  = 0.05,
                            readable      = FALSE)


library(ReactomePA)
ComparePathway <- compareCluster(geneCluster   = list(All = all_ALE_gene_enid, Pseudo = ALE_gene_0_enid, Real = ALE_gene_1_enid),
                                 fun="enrichPathway",
                                 organism = "human",
                                 pvalueCutoff = 0.05, 
                                 pAdjustMethod = "BH",
                                 qvalueCutoff = 0.2,
                                 universe = background_enid,
                                 minGSSize = 10,
                                 maxGSSize = 500,
                                 readable = FALSE)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/ALE/ALE_ComparePathway.pdf", height = 8, width = 20)
print(dotplot(ComparePathway,showCategory = 10,font.size = 16,title = "ALE Pathway"))
dev.off()













#### MXE ========================================================================================================================================



load("/mnt/data5/BGI/UCB/tangchao/DSU2/MXE/MXE_psi_new.RData")

PSI <- MXE_psi[, row.names(Cell_type)]
rs <- rowSums(!is.na(PSI))
summary(rs)

sum(rs>=1)
[1] 75
c(MXE_psi[rs>=1, "gene_id1"], MXE_psi[rs>=1, "gene_id2"])
all_MXE_gene <- c(MXE_psi[rs>=1, "gene_id1"], MXE_psi[rs>=1, "gene_id2"])[!is.na(c(MXE_psi[rs>=1, "gene_id1"], MXE_psi[rs>=1, "gene_id2"]))]


#### GO

MXE_BP_ego <- enrichGO(gene          = all_MXE_gene,
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
MXE_MF_ego <- enrichGO(gene          = all_MXE_gene,
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
MXE_CC_ego <- enrichGO(gene          = all_MXE_gene,
                       universe      = background,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/MXE/all_MXE_gene_GO.pdf",height = 10,width = 15)
print(dotplot(MXE_BP_ego,showCategory = 10,font.size = 16,title = "BP"))
print(dotplot(MXE_MF_ego,showCategory = 10,font.size = 16,title = "MF"))
print(dotplot(MXE_CC_ego,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(MXE_BP_ego,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(MXE_MF_ego,firstSigNodes = 15, useInfo = "all"))
suppressMessages(plotGOgraph(MXE_CC_ego,firstSigNodes = 15, useInfo = "all"))
enrichMap(MXE_BP_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(MXE_MF_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
enrichMap(MXE_CC_ego,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
dev.off()



all_MXE_gene_enid <- bitr(all_MXE_gene, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID

MXE_KEGG <- enrichKEGG(gene         = all_MXE_gene_enid,
                       organism     = 'hsa',
                       universe	   = background_enid,
                       pvalueCutoff = 0.05)

MXE_mkk <- enrichMKEGG(gene     = all_MXE_gene_enid,
                       organism = 'hsa',
                       universe = background_enid)

MXE_do <- DOSE::enrichDO(gene    	  = all_MXE_gene_enid, 
                         ont           = "DO", 
                         pvalueCutoff  = 0.05, 
                         pAdjustMethod = "BH",
                         universe 	  = background_enid,
                         minGSSize     = 5, 
                         maxGSSize     = 500, 
                         qvalueCutoff  = 0.05,
                         readable      = FALSE)

pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/MXE/all_MXE_gene_DO.pdf",height = 10,width = 15)
print(dotplot(MXE_do,showCategory = 10,font.size = 16,title = "DO"))
dev.off()


dgn <- enrichDGN(gene         = all_MXE_gene_enid,
                 universe     = background_enid,
                 qvalueCutoff = 0.05,
                 readable     = TRUE)


MXE_Path <- enrichPathway(gene = all_MXE_gene_enid, 
                          organism = "human", 
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH", 
                          qvalueCutoff = 0.2, 
                          universe = background_enid, 
                          minGSSize = 10,
                          maxGSSize = 500, 
                          readable = FALSE)






BG_BP_ego_1 <- enrichGO(gene          = background,
                         #universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "BP", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
BG_MF_ego_1 <- enrichGO(gene          = background,
                         #universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "MF", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
BG_CC_ego_1 <- enrichGO(gene          = background,
                         #universe      = background,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "CC", 
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)


BG_KEGG <- enrichKEGG(gene         = background_enid,
                       organism     = 'hsa',
                       #universe	   = background_enid,
                       pvalueCutoff = 0.05)

BG_do <- DOSE::enrichDO(gene    	  = background_enid, 
                         ont           = "DO", 
                         pvalueCutoff  = 0.05, 
                         pAdjustMethod = "BH",
                         #universe 	  = background_enid,
                         minGSSize     = 5, 
                         maxGSSize     = 500, 
                         qvalueCutoff  = 0.05,
                         readable      = FALSE)

BG_dgn <- enrichDGN(gene         = background_enid,
                   #universe     = background_enid,
                   qvalueCutoff = 0.05,
                   readable     = TRUE)

library(ReactomePA)
BG_Path <- enrichPathway(gene = background_enid, 
                          organism = "human", 
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH", 
                          qvalueCutoff = 0.2, 
                          #universe = background_enid, 
                          minGSSize = 10,
                          maxGSSize = 500, 
                          readable = FALSE)


pdf("/mnt/data5/BGI/UCB/tangchao/DSU2_GO/Background_gene_GO.pdf",height = 10,width = 15)
print(dotplot(BG_BP_ego_1,showCategory = 10,font.size = 16,title = "BP"))
suppressMessages(plotGOgraph(BG_BP_ego_1,firstSigNodes = 15, useInfo = "all"))
enrichMap(BG_BP_ego_1,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
print(dotplot(BG_MF_ego_1,showCategory = 10,font.size = 16,title = "MF"))
suppressMessages(plotGOgraph(BG_MF_ego_1,firstSigNodes = 15, useInfo = "all"))
enrichMap(BG_MF_ego_1,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)
print(dotplot(BG_CC_ego_1,showCategory = 10,font.size = 16,title = "CC"))
suppressMessages(plotGOgraph(BG_CC_ego_1,firstSigNodes = 15, useInfo = "all"))
enrichMap(BG_CC_ego_1,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(BG_KEGG,showCategory = 10,font.size = 16, title = "KEGG"))
enrichMap(BG_KEGG,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(BG_do,showCategory = 10,font.size = 16, title = "DO"))
enrichMap(BG_do,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(BG_dgn,showCategory = 10,font.size = 16, title = "DGN"))
enrichMap(BG_dgn,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

print(dotplot(BG_Path,showCategory = 10,font.size = 16, title = "PathWay"))
enrichMap(BG_Path,n = 30, vertex.label.font = 1,cex = 1, font.size = 1)

dev.off()







